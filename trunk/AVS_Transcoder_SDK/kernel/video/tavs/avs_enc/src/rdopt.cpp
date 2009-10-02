/*
*****************************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2003, Advanced Audio Video Coding Standard, Part II
*
* DISCLAIMER OF WARRANTY
*
* The contents of this file are subject to the Mozilla Public License
* Version 1.1 (the "License"); you may not use this file except in
* compliance with the License. You may obtain a copy of the License at
* http://www.mozilla.org/MPL/
*
* Software distributed under the License is distributed on an "AS IS"
* basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
* License for the specific language governing rights and limitations under
* the License.
*                     
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
* The AVS Working Group doesn't represent or warrant that the programs
* furnished here under are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy for standardization procedure is available at 
* AVS Web site http://www.avs.org.cn. Patent Licensing is outside
* of AVS Working Group.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

/*
*************************************************************************************
* File name: 
* Function: 
*
*************************************************************************************
*/
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h>

#include "global.h"
#include "const_data.h"
#include "transcoding_type.h"
#include "sse_header.h"
#ifdef FastME
#include "fast_me.h"
#endif


void c_avs_enc::clear_rdopt ()
{
  free_mem_DCcoeff (cofDC);
  free_mem_ACcoeff (cofAC);
  free_mem_ACcoeff (cofAC8x8);
  free_mem_ACcoeff (cofAC4x4intern);
  free_mem4Dint(chromacofAC4x4, 2, 4 );
  // structure for saving the coding state
  delete_coding_state (cs_mb);
  delete_coding_state (cs_b8);
  delete_coding_state (cs_cm);

}


void c_avs_enc::init_rdopt ()
{
  get_mem_DCcoeff (&cofDC);
  get_mem_ACcoeff (&cofAC);
  get_mem_ACcoeff (&cofAC8x8);
  get_mem_ACcoeff (&cofAC4x4intern);
  cofAC4x4 = cofAC4x4intern[0];
  get_mem4Dint(&(chromacofAC4x4),2,4,2,17);

  // structure for saving the coding state
  cs_mb  = create_coding_state ();
  cs_b8  = create_coding_state ();
  cs_cm  = create_coding_state ();
}

/*
*************************************************************************
* Function:Get RD cost for AVS intra block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
double c_avs_enc::RDCost_for_AVSIntraBlocks(int_32_t *nonzero,
                                 int_32_t b8,
                                 int_32_t ipmode,
                                 double lambda,
                                 double  min_rdcost,
                                 int_32_t mostProbableMode)
{
  int_32_t x,y,rate,tmp_cbp,tmp_cbp_blk;
  int_32_t distortion;
  int_32_t block_x;
  int_32_t block_y;
  int_32_t pic_pix_x;
  int_32_t pic_pix_y;
  int_32_t pic_block_x;
  int_32_t pic_block_y;
  int_32_t even_block;
  int_16_t tmp_block_88_inv[8][8];
  Macroblock  *currMB;
  SyntaxElement *currSE;
  block_x    =8*(b8%2);
  block_y    =8*(b8/2);
  pic_pix_x  =img->pix_x+block_x;
  pic_pix_y  =img->pix_y+block_y;
  pic_block_x=pic_pix_x/4;
  pic_block_y=pic_pix_y/4;
  even_block =0;
  currMB   =&img->mb_data[img->current_mb_nr];
  currSE   =&img->MB_SyntaxElements[currMB->currSEnr];
  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====

  for(y=0;y<8;y++)
    for(x=0;x<8;x++)
      tmp_block_88_inv[y][x]=img->m7[y][x];        //the subblock to be processed is in the top left corner of img->m7[][].

  tmp_cbp=tmp_cbp_blk=0;
  avs_dct_sse(tmp_block_88_inv);
  scanquant_B8(img->qp, 4, b8, tmp_block_88_inv, 0, &tmp_cbp, &tmp_cbp_blk);
  *nonzero=(tmp_cbp!=0);
  //===== get distortion (SSD) of 4x4 block =====
  distortion =0;
  for (y=0; y<8; y++)
    for (x=0; x<8; x++)
      distortion += img->quad [imgY_org[pic_pix_y+y][pic_pix_x+x] - imgY[pic_pix_y+y][pic_pix_x+x]];

  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
  currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  currSE->context = 4*b8;
  currSE->type    = SE_INTRAPREDMODE;
  //--- encode and update rate ---
  writeSyntaxElement_Intra4x4PredictionMode(currSE, currBitStream);

  rate = currSE->len;
  currSE++;
  currMB->currSEnr++;

  //===== RATE for LUMINANCE COEFFICIENTS =====

  x=currMB->cbp;
  currMB->cbp=tmp_cbp;//writeLumaCoeffAVS_B8 needs a correct Macroblock->cbp .
  rate+=writeLumaCoeffAVS_B8(b8,1);
  currMB->cbp=x;
  //calc RD and return it.
  return (double)distortion+lambda*(double)rate;
}

/*
*************************************************************************
* Function:4x4 Intra mode decision for an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::Mode_Decision_for_AVSIntraMacroblock (double lambda,  int_32_t* total_cost)
{
  int_32_t  cbp=0, b8;
  int_32_t  min_cost;

  int_32_t ipmode,best_ipmode,x,y;
  int_32_t c_nz, nonzero;
  __declspec(align(16)) byte rec4x4[8][8];
  double rdcost=0;
  int_32_t block_x;
  int_32_t block_y;
  int_32_t pic_pix_x;
  int_32_t pic_pix_y;
  int_32_t pic_block_x;
  int_32_t pic_block_y;
  int_32_t upMode;
  int_32_t leftMode;
  int_32_t mostProbableMode;
  Macroblock *currMB;
  SyntaxElement *currSE;
  double min_rdcost ;
  int_32_t MBRowSize = img->width / MB_BLOCK_SIZE;
  byte  *d1, *d2, *d3, *d4;
  int_16_t *d5, *d6;
  int_32_t loc_cost = (int_32_t)floor(6.0 * lambda + 0.4999);//,best_AVS_cost;

  int_32_t rate,tmp_cbp,tmp_cbp_blk;
  int_32_t distortion=0;
  __declspec(align(16)) int_16_t tmp_block_88_inv[8][8];
  currMB =  img->mb_data+img->current_mb_nr;
  currSE = &img->MB_SyntaxElements[currMB->currSEnr];

  for (min_cost=0, b8=0; b8<4; b8++)
  {
    block_x    = (b8&1) << 3;
    block_y    = (b8>>1) << 3;
    pic_pix_x  =img->pix_x+block_x;
    pic_pix_y  =img->pix_y+block_y;
    pic_block_x=pic_pix_x>>3;
    pic_block_y=pic_pix_y>>3;
    min_rdcost =1e30;
    min_cost = (1<<20);

    upMode           = img->ipredmode[pic_block_x+1][pic_block_y  ];
    leftMode         = img->ipredmode[pic_block_x  ][pic_block_y+1];
    mostProbableMode = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

    //===== INTRA PREDICTION FOR 4x4 BLOCK =====
    intrapred_luma_AVS(pic_pix_x,pic_pix_y);

    //===== LOOP OVER ALL INTRA PREDICTION MODES =====
    for (ipmode=0;ipmode<NO_INTRA_PMODE;ipmode++)
    {
      if(img->available_intra_mode[ipmode] == 1)
      {
        d1 = &imgY_org[pic_pix_y][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+1][pic_pix_x];
        d3 = &img->mprr[ipmode][0][0];
        d4 = &img->mprr[ipmode][1][0];
        d5 = &img->mpr[block_y][block_x];
        d6 = &img->mpr[block_y+1][block_x];

        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm1, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm3, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm0, xmm7
          punpcklbw  xmm1, xmm7
          punpcklbw  xmm2, xmm7
          punpcklbw  xmm3, xmm7
          psubw       xmm0, xmm2
          psubw       xmm1, xmm3
          movdqa      tmp_block_88_inv,    xmm0
          movdqa      tmp_block_88_inv+16, xmm1
          mov      esi,  dword ptr [d5]
          movdqa      xmmword ptr [esi],   xmm2
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi],   xmm3
        }

        d1 = &imgY_org[pic_pix_y+2][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+3][pic_pix_x];
        d3 = &img->mprr[ipmode][2][0];
        d4 = &img->mprr[ipmode][3][0];
        d5 = &img->mpr[block_y+2][block_x];
        d6 = &img->mpr[block_y+3][block_x];

        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm1, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm3, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm0, xmm7
          punpcklbw  xmm1, xmm7
          punpcklbw  xmm2, xmm7
          punpcklbw  xmm3, xmm7
          psubw       xmm0, xmm2
          psubw       xmm1, xmm3
          movdqa      tmp_block_88_inv+32, xmm0
          movdqa      tmp_block_88_inv+48, xmm1
          mov      esi,  dword ptr [d5]
          movdqa      xmmword ptr [esi],   xmm2
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi],   xmm3
        }

        d1 = &imgY_org[pic_pix_y+4][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+5][pic_pix_x];
        d3 = &img->mprr[ipmode][4][0];
        d4 = &img->mprr[ipmode][5][0];
        d5 = &img->mpr[block_y+4][block_x];
        d6 = &img->mpr[block_y+5][block_x];

        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm1, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm3, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm0, xmm7
          punpcklbw  xmm1, xmm7
          punpcklbw  xmm2, xmm7
          punpcklbw  xmm3, xmm7
          psubw       xmm0, xmm2
          psubw       xmm1, xmm3
          movdqa      tmp_block_88_inv+64, xmm0
          movdqa      tmp_block_88_inv+80, xmm1
          mov      esi,  dword ptr [d5]
          movdqa      xmmword ptr [esi],   xmm2
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi],   xmm3
        }

        d1 = &imgY_org[pic_pix_y+6][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+7][pic_pix_x];
        d3 = &img->mprr[ipmode][6][0];
        d4 = &img->mprr[ipmode][7][0];
        d5 = &img->mpr[block_y+6][block_x];
        d6 = &img->mpr[block_y+7][block_x];

        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm1, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm3, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm0, xmm7
          punpcklbw  xmm1, xmm7
          punpcklbw  xmm2, xmm7
          punpcklbw  xmm3, xmm7
          psubw       xmm0, xmm2
          psubw       xmm1, xmm3
          movdqa      tmp_block_88_inv+96,  xmm0
          movdqa      tmp_block_88_inv+112, xmm1
          mov      esi,  dword ptr [d5]
          movdqa      xmmword ptr [esi],    xmm2
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi],    xmm3
        }

        store_coding_state (cs_cm);
        //展开RDCost_for_AVSIntraBlocks()
        //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
        tmp_cbp=tmp_cbp_blk=0;
        avs_dct_sse(tmp_block_88_inv);
        scanquant_B8(img->qp, 4, b8, tmp_block_88_inv, 0, &tmp_cbp, &tmp_cbp_blk);
        c_nz=(tmp_cbp!=0);
        //===== get distortion (SSD) of 4x4 block =====
        d1 = &imgY_org[pic_pix_y][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+1][pic_pix_x];
        d3 = &imgY[pic_pix_y][pic_pix_x];
        d4 = &imgY[pic_pix_y+1][pic_pix_x];
        __asm
        {
          mov      esi,  dword ptr [d1]    //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm4, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in imgY[pic_pix_y][pic_pix_x]
          movdqu    xmm5, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm6, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm0, xmm7
          punpcklbw  xmm4, xmm7
          punpcklbw  xmm5, xmm7
          punpcklbw  xmm6, xmm7
          psubw       xmm0, xmm5
          psubw       xmm4, xmm6
          pmullw      xmm0, xmm0
          pmullw      xmm4, xmm4
          paddw       xmm0, xmm4
        }

        d1 = &imgY_org[pic_pix_y+2][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+3][pic_pix_x];
        d3 = &imgY[pic_pix_y+2][pic_pix_x];
        d4 = &imgY[pic_pix_y+3][pic_pix_x];
        __asm
        {
          mov      esi,  dword ptr [d1]    //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm1, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm4, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in imgY[pic_pix_y][pic_pix_x]
          movdqu    xmm5, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm6, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm1, xmm7
          punpcklbw  xmm4, xmm7
          punpcklbw  xmm5, xmm7
          punpcklbw  xmm6, xmm7
          psubw       xmm1, xmm5
          psubw       xmm4, xmm6
          pmullw      xmm1, xmm1
          pmullw      xmm4, xmm4
          paddw       xmm1, xmm4
        }

        d1 = &imgY_org[pic_pix_y+4][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+5][pic_pix_x];
        d3 = &imgY[pic_pix_y+4][pic_pix_x];
        d4 = &imgY[pic_pix_y+5][pic_pix_x];
        __asm
        {
          mov      esi,  dword ptr [d1]    //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm4, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in imgY[pic_pix_y][pic_pix_x]
          movdqu    xmm5, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm6, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm2, xmm7
          punpcklbw  xmm4, xmm7
          punpcklbw  xmm5, xmm7
          punpcklbw  xmm6, xmm7
          psubw       xmm2, xmm5
          psubw       xmm4, xmm6
          pmullw      xmm2, xmm2
          pmullw      xmm4, xmm4
          paddw       xmm2, xmm4
        }

        d1 = &imgY_org[pic_pix_y+6][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+7][pic_pix_x];
        d3 = &imgY[pic_pix_y+6][pic_pix_x];
        d4 = &imgY[pic_pix_y+7][pic_pix_x];
        __asm
        {
          mov      esi,  dword ptr [d1]    //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm3, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm4, xmmword ptr [esi]

          mov      esi,  dword ptr [d3]  //read in imgY[pic_pix_y][pic_pix_x]
          movdqu    xmm5, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm6, xmmword ptr [esi]

          pxor    xmm7, xmm7        //byte -> int_16_t
          punpcklbw  xmm3, xmm7
          punpcklbw  xmm4, xmm7
          punpcklbw  xmm5, xmm7
          punpcklbw  xmm6, xmm7
          psubw       xmm3, xmm5
          psubw       xmm4, xmm6
          pmullw      xmm3, xmm3
          pmullw      xmm4, xmm4
          paddw       xmm3, xmm4
          paddw       xmm0, xmm1
          paddw       xmm0, xmm2
          paddw       xmm0, xmm3
          psadbw      xmm0, xmm7
          movdqa    xmm7, xmm0
          punpckhqdq  xmm0, xmm0
          paddw       xmm0, xmm7
          pextrw      eax,  xmm0, 0
          mov      distortion, eax
        }

        //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
        currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

        //--- set position and type ---
        currSE->context = 4*b8;
        currSE->type    = SE_INTRAPREDMODE;
        //--- encode and update rate ---
        writeSyntaxElement_Intra4x4PredictionMode(currSE, currBitStream);

        rate = currSE->len;
        currSE++;
        currMB->currSEnr++;

        //===== RATE for LUMINANCE COEFFICIENTS =====

        x=currMB->cbp;
        currMB->cbp=tmp_cbp;//writeLumaCoeffAVS_B8 needs a correct Macroblock->cbp .
        rate+=writeLumaCoeffAVS_B8(b8,1);
        currMB->cbp=x;
        //calc RD and return it.
        //return (double)distortion+lambda*(double)rate;

        // get and check rate-distortion cost
        //展开RDCost_for_AVSIntraBlocks()
        //if ((rdcost = RDCost_for_AVSIntraBlocks(&c_nz,b8,ipmode,lambda,min_rdcost, mostProbableMode)) < min_rdcost)
        if ((rdcost = (double)distortion+lambda*(double)rate) < min_rdcost)
        {
          for(y=0;y<4;y++)
          {
            for(x=0;x<2;x++)
            {
              d5 = &img->cofAC[b8][y][x][0];
              d6 = &cofAC4x4[y][x][0];
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              d5 += 8;
              d6 += 8;
              __asm
              {
                mov      esi,  dword ptr [d5]
                movdqa      xmm0, xmmword ptr [esi]
                mov      esi,  dword ptr [d6]
                movdqa      xmmword ptr [esi], xmm0
              }
              cofAC4x4[y][x][64]=img->cofAC[b8][y][x][64];
            }
          }
          //--- set reconstruction ---
          d1 = &imgY[pic_pix_y  ][pic_pix_x];
          d2 = &imgY[pic_pix_y+1][pic_pix_x];
          d3 = &imgY[pic_pix_y+2][pic_pix_x];
          d4 = &imgY[pic_pix_y+3][pic_pix_x];
          __asm
          {
            mov      esi,  dword ptr [d1]    //read in imgY[pic_pix_y][pic_pix_x]
            movdqu    xmm0, xmmword ptr [esi]
            mov      esi,  dword ptr [d2]
            movdqu    xmm1, xmmword ptr [esi]
            mov      esi,  dword ptr [d3]    //read in imgY[pic_pix_y][pic_pix_x]
            movdqu    xmm2, xmmword ptr [esi]
            mov      esi,  dword ptr [d4]
            movdqu    xmm3, xmmword ptr [esi]

            punpcklqdq  xmm0, xmm1
            punpcklqdq  xmm2, xmm3
            movdqa      rec4x4, xmm0
            movdqa      rec4x4+16, xmm2
          }

          d1 = &imgY[pic_pix_y+4][pic_pix_x];
          d2 = &imgY[pic_pix_y+5][pic_pix_x];
          d3 = &imgY[pic_pix_y+6][pic_pix_x];
          d4 = &imgY[pic_pix_y+7][pic_pix_x];
          __asm
          {
            mov      esi,  dword ptr [d1]    //read in imgY[pic_pix_y][pic_pix_x]
            movdqu    xmm0, xmmword ptr [esi]
            mov      esi,  dword ptr [d2]
            movdqu    xmm1, xmmword ptr [esi]
            mov      esi,  dword ptr [d3]    //read in imgY[pic_pix_y][pic_pix_x]
            movdqu    xmm2, xmmword ptr [esi]
            mov      esi,  dword ptr [d4]
            movdqu    xmm3, xmmword ptr [esi]

            punpcklqdq  xmm0, xmm1
            punpcklqdq  xmm2, xmm3
            movdqa      rec4x4+32, xmm0
            movdqa      rec4x4+48, xmm2
          }

          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          min_rdcost  = rdcost;
          min_cost=(int_32_t)min_rdcost;
          best_ipmode = ipmode;
        }
        reset_coding_state (cs_cm);
      }
    }

    //===== set intra mode prediction =====

    img->ipredmode[pic_block_x+1][pic_block_y+1] = best_ipmode;
    currMB->intra_pred_modes[b8] = mostProbableMode == best_ipmode ? -1 : best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1;

    for(y=0;y<4;y++)
    {
      for(x=0;x<2;x++)
      {
        d6 = &img->cofAC[b8][y][x][0];
        d5 = &cofAC4x4[y][x][0];
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        d5 += 8;
        d6 += 8;
        __asm
        {
          mov      esi,  dword ptr [d5]
          movdqa      xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d6]
          movdqa      xmmword ptr [esi], xmm0
        }
        img->cofAC[b8][y][x][64] = cofAC4x4[y][x][64];
      }
    }

  /*  for(k=0;k<4;k++)
      for(j=0;j<2;j++)
        for(i=0;i<65;i++)
          img->cofAC[b8][k][j][i] = cofAC4x4[k][j][i];
  */
    //--- set reconstruction ---
    d1 = &imgY[pic_pix_y][pic_pix_x];
    d2 = &imgY[pic_pix_y+1][pic_pix_x];
    d3 = &imgY[pic_pix_y+2][pic_pix_x];
    d4 = &imgY[pic_pix_y+3][pic_pix_x];
    __asm
    {
      movq        mm0, rec4x4
      movq        mm1, rec4x4+8
      movq        mm2, rec4x4+16
      movq        mm3, rec4x4+24

      mov      esi,  dword ptr [d1]
      movq    xmmword ptr [esi], mm0
      mov      esi,  dword ptr [d2]
      movq    xmmword ptr [esi], mm1
      mov      esi,  dword ptr [d3]
      movq    xmmword ptr [esi], mm2
      mov      esi,  dword ptr [d4]
      movq    xmmword ptr [esi], mm3
    }

    d1 = &imgY[pic_pix_y+4][pic_pix_x];
    d2 = &imgY[pic_pix_y+5][pic_pix_x];
    d3 = &imgY[pic_pix_y+6][pic_pix_x];
    d4 = &imgY[pic_pix_y+7][pic_pix_x];
    __asm
    {
      movq        mm0, rec4x4+32
      movq        mm1, rec4x4+40
      movq        mm2, rec4x4+48
      movq        mm3, rec4x4+56

      mov      esi,  dword ptr [d1]
      movq    xmmword ptr [esi], mm0
      mov      esi,  dword ptr [d2]
      movq    xmmword ptr [esi], mm1
      mov      esi,  dword ptr [d3]
      movq    xmmword ptr [esi], mm2
      mov      esi,  dword ptr [d4]
      movq    xmmword ptr [esi], mm3
      emms
    }

  for(y=0; y<8; y++)
    {
      for(x=0; x<8; x++)
      {
        img->mpr[block_y+y][block_x+x] = img->mprr[best_ipmode][y][x];
      }
    }
    min_cost += loc_cost;

    if (nonzero)
    {
      cbp |= (1<<b8);
    }
    *total_cost += min_cost;
  }
  return cbp;
}
  // xzhao }

/*
*************************************************************************
* Function:R-D Cost for an 8x8 Partition
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

double c_avs_enc:: RDCost_for_8x8blocks (int_32_t*    cnt_nonz,   // --> number of nonzero coefficients
                             int_32_t*    cbp_blk,    // --> cbp blk
                             double  lambda,     // <-- Lagrange multiplier
                             int_32_t     block,      // <-- 8x8 block number
                             int_32_t     mode,       // <-- partitioning mode
                             int_32_t     pdir,       // <-- prediction direction
                             int_32_t     ref,        // <-- reference frame
                             int_32_t     bwd_ref)   // <-- abp type
{
  int_32_t  i, j;
  int_32_t  rate=0, distortion=0;
  int_32_t  dummy, mrate;
  int_32_t  fw_mode, bw_mode;
  int_32_t  cbp     = 0;
  int_32_t  pax     = 8*(block%2);
  int_32_t  pay     = 8*(block/2);
  int_32_t  i0      = pax/8;
  int_32_t  j0      = pay/8;
  int_32_t  bframe  = (img->type==B_IMG);
  int_32_t  direct  = (bframe && mode==0);
  int_32_t  b8value = B8Mode2Value (mode, pdir);

  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t  **frefarr = refFrArr;  // For MB level field/frame
  int_32_t  block_y = img->block_y;
  int_32_t  block_x = img->block_x;
  int_32_t pix_y = img->pix_y;
  int_32_t pix_x = img->pix_x;
  byte **imgY_original  =  imgY_org;

  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  if (direct)
  {
    *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, 0, 0,
      max(0,frefarr[block_y+j0][block_x+i0]), 0);
  }
  else
  {
    fw_mode   = (pdir==0||pdir==2 ? mode : 0);
    bw_mode   = (pdir==1||pdir==2 ? mode : 0);
    *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, fw_mode, bw_mode, ref, bwd_ref);
  }

  for (j=0; j<8; j++)
  {
    for (i=pax; i<pax+8; i++)
    {
      distortion += img->quad [imgY_original[pix_y+pay+j][pix_x+i] - imgY[img->pix_y+pay+j][img->pix_x+i]];
    }
  }

  //=====   GET RATE
  //----- block 8x8 mode -----
  ue_linfo (b8value, dummy, &mrate, &dummy);
  rate += mrate;

  //----- motion information -----
  if (!direct)
  {
    if (pdir==0 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+1, j0+1, ref, 0, 1, mode);
    }
    if (pdir==1 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+1, j0+1, bwd_ref, 0, 0, mode);
    }
  }
  //----- luminance coefficients -----
  if (*cnt_nonz)
  {
    currMB->cbp = cbp;
    rate += writeLumaCoeff8x8 (block, 0);
  }
  return (double)distortion + lambda * (double)rate;
}

/*
*************************************************************************
* Function:Sets modes and reference frames for an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/


void c_avs_enc:: SetModesAndRefframeForBlocks (int_32_t mode)
{
  int_32_t i,j,k;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int_32_t  bframe  = (img->type==B_IMG);
  int_32_t  **fwrefarr = fw_refFrArr;
  int_32_t  **bwrefarr = bw_refFrArr;
  int_32_t  **frefarr = refFrArr;
  int_32_t  block_y = img->block_y;
  int_32_t  block_x = img->block_x;
  //--- macroblock type ---
  currMB->mb_type = mode;
  //--- block 8x8 mode and prediction direction ---
  switch (mode)
    {
    case 0:
      for(i=0;i<4;i++)
        {
        currMB->b8mode[i] = 0;
        currMB->b8pdir[i] = (bframe ? 2 : 0);
        }
      break;
    case 1:
    case 2:
    case 3:
      for(i=0;i<4;i++)
        {
        currMB->b8mode[i] = mode;
        currMB->b8pdir[i] = best8x8pdir[mode][i];
        }
      break;
    case P8x8:
      for(i=0;i<4;i++)
        {
        currMB->b8mode[i]   = best8x8mode[i];
        currMB->b8pdir[i]   = best8x8pdir[mode][i];
        }
      break;
    case I4MB:
      for(i=0;i<4;i++)
        {
        currMB->b8mode[i] = IBLOCK;
        currMB->b8pdir[i] = -1;
        }
      break;
    default:
      printf ("Unsupported mode in SetModesAndRefframeForBlocks!\n");
      exit (1);
    }

  //--- reference frame arrays ---
  if (mode==0 || mode==I4MB)
    {
    if (bframe)
      {
      for (j=0;j<2;j++)
        {
        for (i=0;i<2;i++)
          {
          k = 2*j+i;
          if(!mode)
            {
            if(!img->picture_structure)
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[mode][k][0];
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[mode][k][1];
              }
            else
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = 0;
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = 0;
              }
            }
          else
            {
            fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = -1;
            bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = -1;
            }
          }
        }
      }
    else
      {
      for (j=0;j<2;j++)
        {
        for (i=0;i<2;i++)
          {
          frefarr [(block_y>>1)+j][(block_x>>1)+i] = (mode==0?0:-1);
          }
        }
      }
    }
  else
    {
    if (bframe)
      {
      for (j=0;j<2;j++)
        {
        for (i=0;i<2;i++)
          {
          k = 2*j+i;
          if((mode==P8x8) && (best8x8mode[k]==0))
            {
            if(!img->picture_structure)
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[0][k][0];
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[0][k][1];
              }
            else
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = 0;
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = 0;
              }
            }
          else
            {
            if(IS_FW&&IS_BW)
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[mode][k][0];
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = best8x8symref[mode][k][1];
              }
            else
              {
              fwrefarr[(block_y>>1)+j][(block_x>>1)+i] = (IS_FW ? best8x8ref[mode][k] : -1);
              bwrefarr[(block_y>>1)+j][(block_x>>1)+i] = (IS_BW ? best8x8bwref[mode][k] : -1);
              }
            }
          }
        }
      }
    else
      {
      for (j=0;j<2;j++)
        {
        for (i=0;i<2;i++)
          {
          k = 2*j+i;
          frefarr [(block_y>>1)+j][(block_x>>1)+i] = (IS_FW ? best8x8ref[mode][k] : -1);
          }
        }
      }
    }
  }
/*
*************************************************************************
* Function:Sets Coefficients and reconstruction for an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: SetCoeffAndReconstruction8x8 (Macroblock* currMB)
{
  int_32_t block, k, j, i;
  int_32_t block_x = img->block_x;
  int_32_t block_y = img->block_y;
  int_32_t **ipredmodes = img->ipredmode;
  int_32_t  stage_block8x8_pos=0;
  int_32_t  incr_y=1,off_y=0;

  //--- restore coefficients ---
  for (block=0; block<6; block++)
    for (k=0; k<4; k++)
      for (j=0; j<2; j++)
        for (i=0; i<65; i++)
          img->cofAC[block][k][j][i] = cofAC8x8[block][k][j][i];

  if (cnt_nonz_8x8<=5)
  {
    currMB->cbp     = 0;
    currMB->cbp_blk = 0;

    for (j=0; j<16; j++)
      for (i=0; i<16; i++)
        imgY[img->pix_y+j][img->pix_x+i] = (byte)mpr8x8[j][i];

  }
  else
  {
    currMB->cbp     = cbp8x8;
    currMB->cbp_blk = cbp_blk8x8;

    for(k=0;k<2;k++)
    {
      for (j=0; j<8; j++)
        for (i=0; i<16; i++)
          imgY[img->pix_y+k*8+j][img->pix_x+i] = rec_mbY8x8[k*8+j][i];
    }
  }

  //===== restore intra prediction modes for 8x8+ macroblock mode =====
  for (j=img->block8_y+1; j<img->block8_y+3; j++)
    {
    for (i=img->block8_x+1; i<img->block8_x+3; i++)
      {
      k = (j-img->block8_y-1)*2 + i -img->block8_x -1;
      ipredmodes[i][j] = b8_ipredmode[k];
      currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
      }
    }
}

void c_avs_enc:: SetMotionVectorsMB (Macroblock* currMB, int_32_t bframe)
{
  int_32_t i, j, k, mode8, pdir8, ref, by, bx, bxr, dref;
  int_32_t ***mv        = tmp_mv;
  int_32_t *****all_mv  = img->all_mv;
  int_32_t *****all_bmv = img->all_bmv;
  int_32_t ***fwMV      = tmp_fwMV;
  int_32_t ***bwMV      = tmp_bwMV;
  int_32_t *****imgmv   = img->mv;
  int_32_t *****p_fwMV  = img->p_fwMV;
  int_32_t *****p_bwMV  = img->p_bwMV;
  int_32_t **refar      = refFrArr;
  int_32_t **frefar     = fw_refFrArr;
  int_32_t block_y      = img->block_y>>1;
  int_32_t block_x      = img->block_x>>1;
  int_32_t **brefar     = bw_refFrArr;
  int_32_t  bw_ref;

  for (j=0;j<2;j++)
  {
    for (i=0; i<2; i++)
    {
      k=2*j+i;
      mode8 = currMB->b8mode[k];
      pdir8 = currMB->b8pdir[k];
      if(pdir8 == 2 && mode8 != 0)
      {
        all_mv = img->all_omv;
        p_fwMV = img->omv;
      }
      else
      {
        all_mv = img->all_mv;
        imgmv = img->mv;
      }

      by     = block_y + j;
      bxr    = block_x + i;
      bx     = block_x + i + 4;
      ref    = (bframe?frefar:refar)[by][bxr];
      bw_ref = (bframe?brefar:refar)[by][bxr];

      if (!bframe)
      {
        if (mode8!=IBLOCK && ref != -1)
        {
          mv[0][by][bx] = all_mv[i][j][ref][mode8][0];
          mv[1][by][bx] = all_mv[i][j][ref][mode8][1];
        }
        else
        {
          mv[0][by][bx] = 0;
          mv[1][by][bx] = 0;
        }
      }
      else
      {
        if (pdir8==-1) // intra
        {
          fwMV [0][by][bx] = 0;
          fwMV [1][by][bx] = 0;
          bwMV [0][by][bx] = 0;
          bwMV [1][by][bx] = 0;
          dfMV [0][by][bx] = 0;
          dfMV [1][by][bx] = 0;
          dbMV [0][by][bx] = 0;
          dbMV [1][by][bx] = 0;
        }
        else if (pdir8==0) // forward
        {
          fwMV [0][by][bx] = all_mv [i][j][ ref][mode8][0];
          fwMV [1][by][bx] = all_mv [i][j][ ref][mode8][1];
          bwMV [0][by][bx] = 0;
          bwMV [1][by][bx] = 0;
          dfMV [0][by][bx] = 0;
          dfMV [1][by][bx] = 0;
          dbMV [0][by][bx] = 0;
          dbMV [1][by][bx] = 0;
        }
        else if (pdir8==1) // backward
        {
          fwMV [0][by][bx] = 0;
          fwMV [1][by][bx] = 0;
          {
            bwMV [0][by][bx] = all_bmv[i][j][bw_ref][mode8][0];
            bwMV [1][by][bx] = all_bmv[i][j][bw_ref][mode8][1];
          }
          dfMV[0][by][bx] = 0;
          dfMV[1][by][bx] = 0;
          dbMV[0][by][bx] = 0;
          dbMV[1][by][bx] = 0;
        }
        else if (mode8!=0) // bidirect
        {
          fwMV[0][by][bx] = all_mv [i][j][ ref][mode8][0];
          fwMV[1][by][bx] = all_mv [i][j][ ref][mode8][1];
          {
            int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,refframe,delta_PB;
            refframe = ref;
            delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
            //rm52j
            delta_P = (delta_P + 512)%512;
            if(img->picture_structure)
              TRp = (refframe+1)*delta_P;  //the lates backward reference
            else
            {
              TRp = delta_P;//refframe == 0 ? delta_P-1 : delta_P+1;
            }
            //rm52c
            //delta_PB = 2*(img->tr - img->imgtr_last_P_frm);
            //rm52j begin
            delta_PB = 2*(picture_distance - img->imgtr_last_P_frm);    // Tsinghua 200701
            TRp = (TRp + 512)%512;
            delta_PB = (delta_PB + 512)%512;  // Added by Xiaozhen ZHENG, 2007.05.05
            //      end
            if(!img->picture_structure)
            {
              if(img->current_mb_nr_fld < img->total_number_mb) //top field
                DistanceIndexFw =  refframe == 0 ? delta_PB-1:delta_PB;
              else
                DistanceIndexFw =  refframe == 0 ? delta_PB:delta_PB+1;
            }
            else
              DistanceIndexFw = delta_PB;
            DistanceIndexBw    = (TRp - DistanceIndexFw+512)%512; // Added by Zhijie Yang, 20070419, Broadcom
            bwMV [0][by][bx] = - ((all_mv[i][j][ref][mode8][0]*DistanceIndexBw*(512/DistanceIndexFw)+256)>>9);
            bwMV [1][by][bx] = - ((all_mv[i][j][ref][mode8][1]*DistanceIndexBw*(512/DistanceIndexFw)+256)>>9);
          }

          dfMV[0][by][bx] = 0;
          dfMV[1][by][bx] = 0;
          dbMV[0][by][bx] = 0;
          dbMV[1][by][bx] = 0;
        }
        else // direct
        {
          if(!img->picture_structure)
          {
            dref = max(0,fw_refFrArr[by][bxr]);
            fwMV [0][by][bx] = dfMV[0][by][bx] = all_mv [i][j][dref][0][0];
            fwMV [1][by][bx] = dfMV[1][by][bx] = all_mv [i][j][dref][0][1];
            dref = max(0,bw_refFrArr[by][bxr]);
            bwMV [0][by][bx] = dbMV[0][by][bx] = all_bmv[i][j][dref][0][0];
            bwMV [1][by][bx] = dbMV[1][by][bx] = all_bmv[i][j][dref][0][1];
          }
          else
          {
            dref = 0;
            fwMV [0][by][bx] = dfMV[0][by][bx] = all_mv [i][j][dref][0][0];
            fwMV [1][by][bx] = dfMV[1][by][bx] = all_mv [i][j][dref][0][1];
            bwMV [0][by][bx] = dbMV[0][by][bx] = all_bmv[i][j][0][0][0];
            bwMV [1][by][bx] = dbMV[1][by][bx] = all_bmv[i][j][0][0][1];
          }
        }
      }
    }
  }
}

/*
*************************************************************************
* Function:R-D Cost for a macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: RDCost_for_macroblocks (double lambda, int_32_t mode, double*  min_rdcost)
{
  int_32_t         i, j;
  int_32_t         rate=0, distortion=0;
  double      rdcost;
  Macroblock  *currMB       = &img->mb_data[img->current_mb_nr];
  int_32_t         bframe   = (img->type==B_IMG);
  int_32_t         tmp_cc;
  int_32_t         use_of_cc     =  (img->type!=INTRA_IMG);
  int_32_t         cc_rate, dummy;
  byte        **imgY_orig   = imgY_org;
  byte        ***imgUV_orig = imgUV_org;
  int_32_t         pix_y         = img->pix_y;
  int_32_t         pix_c_y       = img->pix_c_y;
  int_32_t         pix_x         = img->pix_x;
  int_32_t         pix_c_x       = img->pix_c_x;
  int_32_t         rate_tmp,rate_total=0, distortion_blk, distortion_total = 0;
  int_32_t         block8x8,block_x,block_y;
  int_32_t         lum_bits;
  int_32_t         cp_table[2][6] = {{1,1,1,1,1,1},{0,0,1,1,0,1}};
  //=====  SET REFERENCE FRAMES AND BLOCK MODES
  SetModesAndRefframeForBlocks (mode);
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  if (mode==I4MB)
  {
    currMB->cbp = Mode_Decision_for_AVSIntraMacroblock (lambda, &dummy);
  }
  else
  {
    LumaResidualCoding ();
  }

  //Rate control
  if (input->RCEnable == 1)
  {
	  for(j=0; j<16; j++)
		  memcpy(pred[j], img->mpr[j], 16*sizeof(int_16_t));
  }

  dummy = 0;
  ChromaResidualCoding (&dummy);

  //=====   GET DISTORTION
  // LUMA

  for(block8x8=0; block8x8<4;block8x8++)
  {
    block_y = (block8x8/2)*8;
    block_x = (block8x8%2)*8;

    distortion_blk =0;

    for (j=0; j<8; j++)
    {
      for (i=0; i<8; i++)
      {
        distortion_blk += img->quad [imgY_orig[block_y+j+pix_y][block_x+i+pix_x] - imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]];
      }
    }
    block_y = (block8x8/2)*4;
    block_x = (block8x8%2)*4;

    // CHROMA
    for (j=0; j<4; j++)
    {
      for (i=0; i<4; i++)
      {
        distortion_blk += img->quad [imgUV_orig[0][block_y+j+pix_c_y][block_x+i+pix_c_x] - imgUV[0][img->pix_c_y+block_y+j][img->pix_c_x+block_x+i]];
        distortion_blk += img->quad [imgUV_orig[1][block_y+j+pix_c_y][block_x+i+pix_c_x] - imgUV[1][img->pix_c_y+block_y+j][img->pix_c_x+block_x+i]];
      }
    }
    distortion_total += distortion_blk;
  }

  //=====   S T O R E   C O D I N G   S T A T E   =====
  store_coding_state (cs_cm);
  //=====   GET RATE
  //----- macroblock header -----
  if (use_of_cc)
  {
    if (currMB->mb_type!=0 || (bframe && currMB->cbp!=0))
    {
      // cod counter and macroblock mode are written ==> do not consider code counter
      tmp_cc = img->cod_counter;
      rate_total   = writeMBHeader (1);
      ue_linfo (tmp_cc, dummy, &cc_rate, &dummy);
      rate_total  -= cc_rate;
      img->cod_counter = tmp_cc;
    }
    else
    {
      // cod counter is just increased  ==> get additional rate
      ue_linfo (img->cod_counter+1, dummy, &rate_total,    &dummy);
      ue_linfo (img->cod_counter,   dummy, &cc_rate, &dummy);
      rate_total -= cc_rate;
    }
  }
  else
  {
    rate_total = writeMBHeader (1);
  }

  if (mode)
  {
    //----- motion information -----
    storeMotionInfo (0);
    writeReferenceIndex(&rate_total);
    writeMVD(&rate_total);
    writeCBPandDqp(&rate_total);
  }

  if (mode || (bframe && (currMB->cbp!=0)))
  {
    for(block8x8=0; block8x8<6;block8x8++)
      if(cp_table[0][block8x8])
      {
        rate_tmp = writeBlockCoeff (block8x8);

        if(block8x8 < 2)
          rate_total +=rate_tmp;
        else if(block8x8< 4)
          rate_total +=rate_tmp;
        if(block8x8 == 4)
          rate_total +=rate_tmp;
        if(block8x8 == 5)
          rate_total +=rate_tmp;

        if(block8x8 < 4)
          lum_bits = rate_total;
      }

  }

  //=====   R E S T O R E   C O D I N G   S T A T E   =====
  reset_coding_state (cs_cm);
  rdcost = (double)distortion_total + lambda * (double)rate_total;

  if (rdcost >= *min_rdcost)
  {
    return 0;
  }
  //=====   U P D A T E   M I N I M U M   C O S T   =====
  *min_rdcost = rdcost;
  return 1;
}

/*
*************************************************************************
* Function:Store macroblock parameters
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: store_macroblock_parameters (int_32_t mode)
{
  int_32_t  i, j, k,l;
  int_16_t ****i4p, ***i3p;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int_32_t        bframe   = (img->type==B_IMG);
  int_32_t        **frefar = ((img->type==B_IMG)? fw_refFrArr : refFrArr);
  int_32_t        **brefar = bw_refFrArr;
  int_32_t        block_x  = img->block_x;
  int_32_t        pix_x    = img->pix_x;
  int_32_t        pix_c_x  = img->pix_c_x;
  int_32_t        block_y  = img->block_y;
  int_32_t        pix_y    = img->pix_y;
  int_32_t        pix_c_y  = img->pix_c_y;
  int_32_t        temp_chroma;

  //--- store best mode ---
  best_mode = mode;
  best_c_imode = currMB->c_ipred_mode;

  best_weight_flag = img->mbweightflag ;

  for (i=0; i<4; i++)
  {
    b8mode[i]   = currMB->b8mode[i];
    b8pdir[i]   = currMB->b8pdir[i];
  }

  //--- reconstructed blocks ----
  for (j=0; j<16; j++)
  {
    for (i=0; i<16; i++)
    {
      rec_mbY[j][i] = imgY[img->pix_y+j][img->pix_x+i];
    }
  }
  for (j=0; j<8; j++)
  {
    for (i=0; i<8; i++)
    {
      rec_mbU[j][i] = imgUV[0][img->pix_c_y+j][img->pix_c_x+i];
      rec_mbV[j][i] = imgUV[1][img->pix_c_y+j][img->pix_c_x+i];
    }
  }
  if(best_mode == I4MB)
  {
    for(i=0; i<4; i++)
    {
      best_intra_pred_modes_tmp[i]=currMB->intra_pred_modes[i];
    }
    for(j=0;j<2;j++)
    {
      for(i=0;i<2;i++)
      {
        best_ipredmode_tmp[i][j]=img->ipredmode[ 1 + img->block8_x + i ][ 1 + img->block8_y + j ];
      }
    }

    for(j=0;j<16;j++)
    {
      for(i=0;i<16;i++)
      {
        best_mpr_tmp[j][i]=img->mpr[j][i];
      }
    }
  }
  //--- coeff, cbp, kac ---
  if (mode || bframe)
  {
    i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
    i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
    cbp     = currMB->cbp;
    cbp_blk = currMB->cbp_blk;
    for(j=0;j<4;j++)
    {
      for(k=0;k<2;k++)
      {
        for(l=0;l<17;l++)
        {
          temp_chroma = chromacofAC4x4[0][j][k][l];
          chromacofAC4x4[0][j][k][l] = img->chromacofAC[0][j][k][l];
          img->chromacofAC[0][j][k][l] = temp_chroma;
          temp_chroma = chromacofAC4x4[1][j][k][l];
          chromacofAC4x4[1][j][k][l] = img->chromacofAC[1][j][k][l];
          img->chromacofAC[1][j][k][l] = temp_chroma;
        }
      }
    }
  }
  else
  {
    cbp = cbp_blk = 0;
  }

  for(j=0;j<2;j++)
  {
    for (i=0; i<2; i++)
    {
      frefframe[j][i] = frefar[(block_y>>1)+j  ][(block_x>>1) +i];
      if (bframe)
      {
        brefframe[j][i] = brefar[(block_y>>1) +j ][(block_x>>1)+i];
      }
    }
  }
}

/*
*************************************************************************
* Function:Set stored macroblock parameters
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: set_stored_macroblock_parameters ()
{
  int_32_t  i, j, k;
  int_16_t ****i4p, ***i3p,l;
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  int_32_t     mode     = best_mode;
  int_32_t     bframe   = (img->type==B_IMG);
  int_32_t     **frefar = ((img->type==B_IMG)? fw_refFrArr : refFrArr);
  int_32_t     **brefar = bw_refFrArr;
  int_32_t     **ipredmodes = img->ipredmode;
  int_32_t     block_x = img->block_x;
  int_32_t     block_y = img->block_y;
  int_32_t     cp_table[2][6] = {{1,1,1,1,1,1},{0,0,1,1,0,1}};
  int_32_t     temp_chroma;
  //===== reconstruction values =====
  for (j=0; j<16; j++)
  {
    for (i=0; i<16; i++)
    {
      imgY[img->pix_y+j][img->pix_x+i] = rec_mbY[j][i];
    }
  }

  for (j=0; j<8; j++)
  {
    for (i=0; i<8; i++)
    {
      imgUV[0][img->pix_c_y+j][img->pix_c_x+i] = rec_mbU[j][i];
      imgUV[1][img->pix_c_y+j][img->pix_c_x+i] = rec_mbV[j][i];
    }
  }

  //===== coefficients and cbp =====
  i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
  i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
  currMB->cbp      = cbp;
  currMB->cbp_blk = cbp_blk;
  //==== macroblock type ====
  currMB->mb_type = mode;

  img->mbweightflag = best_weight_flag ;

  for(j=0;j<4;j++)
  {
    for(k=0;k<2;k++)
    {
      for(l=0;l<17;l++)
      {
        temp_chroma = chromacofAC4x4[0][j][k][l];
        chromacofAC4x4[0][j][k][l] = img->chromacofAC[0][j][k][l];
        img->chromacofAC[0][j][k][l] = temp_chroma;

        temp_chroma = chromacofAC4x4[1][j][k][l];
        chromacofAC4x4[1][j][k][l] = img->chromacofAC[1][j][k][l];
        img->chromacofAC[1][j][k][l] = temp_chroma;
      }
    }
  }
  for (i=0; i<4; i++)
  {
    currMB->b8mode[i]   = b8mode[i];
    currMB->b8pdir[i]   = b8pdir[i];
  }
  //==== reference frames =====
  for (j=0; j<2; j++)
  {
    for (i=0; i<2; i++)
    {
      frefar[(block_y>>1)+j][(block_x>>1)+i] = frefframe[j][i];
    }
  }

  if (bframe)
  {
    for (j=0; j<2; j++)
    {
      for (i=0; i<2; i++)
      {
        brefar[(block_y>>1)+j][(block_x>>1)+i] = brefframe[j][i];
      }
    }
  }

  if(currMB->mb_type == I4MB)
  {
    for(i=0; i<4; i++)
    {
      currMB->intra_pred_modes[i]=best_intra_pred_modes_tmp[i];
    }

    for(j=0;j<2;j++)
    {
      for(i=0;i<2;i++)
      {
        img->ipredmode[ 1 + (img->mb_x<<1) + i ][ 1 + (img->mb_y<<1) + j ]=best_ipredmode_tmp[i][j];
      }
    }

    for(j=0;j<16;j++)
    {
      for(i=0;i<16;i++)
      {
        img->mpr[j][i]=best_mpr_tmp[j][i];
      }
    }
  }

  //==== intra prediction modes ====
  currMB->c_ipred_mode = best_c_imode;

  if (mode==P8x8)
  {
    for (j=img->block8_y+1; j<img->block8_y+3; j++)
    {
      for ( i=img->block8_x+1; i<img->block8_x+3; i++)
      {
        k = 2*(j-img->block8_y-1)+i-img->block8_x-1;
        ipredmodes[i][j] = b8_ipredmode[k];
        currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
      }
    }
  }
  else if (mode!=I4MB)
  {
    for (j=img->block8_y+1; j<img->block8_y+3; j++)
    {
      for ( i=img->block8_x+1; i<img->block8_x+3; i++)
      {
        k = 2*(j-img->block8_y-1)+i-img->block8_x-1;
        ipredmodes    [i][j] = -1;
        currMB->intra_pred_modes[k] = DC_PRED;
      }
    }
  }
  //==== motion vectors =====
  SetMotionVectorsMB (currMB, bframe);
  storeMotionInfo(0);
}

/*
*************************************************************************
* Function:Set reference frames and motion vectors
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::SetRefAndMotionVectors (int_32_t block, int_32_t mode, int_32_t ref, int_32_t bw_ref,int_32_t pdir)
{
  int_32_t i, j;
  int_32_t     bframe       = (img->type==B_IMG);
  int_32_t     pmode        = (mode==1||mode==2||mode==3?mode:4);
  int_32_t     j0           = (block/2);
  int_32_t     i0           = (block%2);
  int_32_t     j1           = j0 + (input->blc_size[pmode][1]>>3);
  int_32_t     i1           = i0 + (input->blc_size[pmode][0]>>3);
  int_32_t**   frefArr      = (bframe ? fw_refFrArr : refFrArr);
  int_32_t**   brefArr      = bw_refFrArr;
  int_32_t***  bmvArr       = tmp_bwMV;
  int_32_t***** all_bmv_arr = img->all_bmv;
  int_32_t   block_x        = img->block_x>>1;
  int_32_t   block_y        = img->block_y>>1;

  int_32_t***  fmvArr  = (bframe ? tmp_fwMV    : tmp_mv);
  int_32_t***** all_mv_arr = (pdir==2 && mode != 0)?img->all_omv:img->all_mv; //mode != 0 added by xyji 7.11

  if ((pdir==0 || pdir==2) && (mode!=IBLOCK && mode!=0))
  {
    for (j=j0; j<j1; j++)
    {
      for (i=i0; i<i1; i++)
      {
        fmvArr[0][block_y+j][block_x+i+4] = all_mv_arr[i][j][ref][mode][0];
        fmvArr[1][block_y+j][block_x+i+4] = all_mv_arr[i][j][ref][mode][1];
        frefArr  [block_y+j][block_x+i  ] = ref;
      }
    }
  }
  else
  {
    if(!mode&&bframe)
    {
      for (j=j0; j<j0+1; j++)
      {
        for (i=i0; i<i0+1; i++)
        {
          if(img->picture_structure)
            ref = 0;
          if(ref==-1)
          {
            fmvArr[0][block_y+j][block_x+i+4] = 0;
            fmvArr[1][block_y+j][block_x+i+4] = 0;
            frefArr  [block_y+j][block_x+i  ] = -1;
          }
          else
          {
            fmvArr[0][block_y+j][block_x+i+4] = all_mv_arr[i][j][ref][mode][0];
            fmvArr[1][block_y+j][block_x+i+4] = all_mv_arr[i][j][ref][mode][1];
            if(img->picture_structure)
              ref = 0;
            frefArr  [block_y+j][block_x+i  ] = ref;
          }
        }
      }
    }
    else
      for (j=j0; j<j0+1; j++)
        for (i=i0; i<i0+1; i++)
        {
          fmvArr[0][block_y+j][block_x+i+4] = 0;
          fmvArr[1][block_y+j][block_x+i+4] = 0;
          frefArr  [block_y+j][block_x+i  ] = -1;
        }
  }

  if ((pdir==1 || pdir==2) && (mode!=IBLOCK && mode!=0))
  {
    for (j=j0; j<j0+1; j++)
      for (i=i0; i<i0+1; i++)
      {
        if(pdir == 2)
        {
          {
            int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,refframe,delta_PB;
            refframe = ref;
            delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
            delta_P = (delta_P + 512)%512;
            if(img->picture_structure)
              TRp = (refframe+1)*delta_P;
            else
            {
              TRp = delta_P;
            }
            delta_PB = 2*(picture_distance - img->imgtr_last_P_frm);
            TRp = (TRp + 512)%512;
            delta_PB = (delta_PB + 512)%512;
            if(!img->picture_structure)
            {
              if(img->current_mb_nr_fld < img->total_number_mb) //top field
                DistanceIndexFw =  refframe == 0 ? delta_PB-1:delta_PB;
              else
                DistanceIndexFw =  refframe == 0 ? delta_PB:delta_PB+1;
            }
            else
              DistanceIndexFw = delta_PB;
            DistanceIndexBw    = (TRp - DistanceIndexFw+512)%512;
            bmvArr[0][block_y+j][block_x+i+4] = - ((all_mv_arr[i][j][ref][mode][0]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
            bmvArr[1][block_y+j][block_x+i+4] = - ((all_mv_arr[i][j][ref][mode][1]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
          }
        }
        else
        {
          bmvArr[0][block_y+j][block_x+i+4] = all_bmv_arr[i][j][bw_ref][mode][0];
          bmvArr[1][block_y+j][block_x+i+4] = all_bmv_arr[i][j][bw_ref][mode][1];
        }
        brefArr  [block_y+j][block_x+i  ] = bw_ref;/*lgp*13*/
      }
  }
  else if (bframe)
  {
    if(!mode)
    {
      for (j=j0; j<j0+1; j++)
        for (i=i0; i<i0+1; i++)
        {
          {
            if(!img->picture_structure)
            {
              if(refFrArr[block_y+j][block_x+i] >= 0)
              {
                bmvArr[0][block_y+j][block_x+i+4] = all_bmv_arr[i][j][bw_ref][mode][0];
                bmvArr[1][block_y+j][block_x+i+4] = all_bmv_arr[i][j][bw_ref][mode][1];
              }
              else//intra prediction
              {
                bmvArr[0][block_y+j][block_x+i+4] = all_bmv_arr[i][j][1][mode][0];
                bmvArr[1][block_y+j][block_x+i+4] = all_bmv_arr[i][j][1][mode][1];
              }

              if(img->current_mb_nr_fld < img->total_number_mb)
                brefArr  [block_y+j][block_x+i  ] = 1;
              else
                brefArr  [block_y+j][block_x+i  ] = 0;
              if(refFrArr[block_y+j][block_x+i] < 0)
                brefArr  [block_y+j][block_x+i  ] = 1;
            }
            else
            {
              bmvArr[0][block_y+j][block_x+i+4] = all_bmv_arr[i][j][0][mode][0];
              bmvArr[1][block_y+j][block_x+i+4] = all_bmv_arr[i][j][0][mode][1];
              ref = refFrArr[block_y+j][block_x+i];
              if(img->picture_structure)
                brefArr  [block_y+j][block_x+i  ] = 0;
              else
                brefArr  [block_y+j][block_x+i  ] = 0;
            }
          }
        }
    }
    else
    {
      for (j=j0; j<j0+1; j++)
      {
        for (i=i0; i<i0+1; i++)
        {
          bmvArr[0][block_y+j][block_x+i+4] = 0;
          bmvArr[1][block_y+j][block_x+i+4] = 0;
          brefArr  [block_y+j][block_x+i  ] = -1;
        }
      }
    }
  }
}
void c_avs_enc::encode_one_intra_macroblock_rdo()
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  double      qp, lambda_mode, min_rdcost;
  int_32_t         i, j;
  int_32_t         max_mcost=(1<<30);
  int_32_t         mb_available_up, mb_available_left, mb_available_up_left;
  //===== SET LAGRANGE PARAMETERS =====
  qp = (double)img->qp - SHIFT_QP;
  if (input->successive_Bframe>0)
    lambda_mode   = 0.68 * pow (2, qp/4.0) * (img->type==B_IMG? max(2.00,min(4.00,(qp / 8.0))):1.0);
  else
    lambda_mode   = 0.85 * pow (2, qp/4.0) * (img->type==B_IMG? 4.0:1.0);
  min_rdcost = 1e30;
  img->NoResidueDirect = 0;
  // pre compute all new chroma intra prediction modes
  IntraChromaPrediction8x8(&mb_available_up, &mb_available_left, &mb_available_up_left);
  for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=PLANE_8; currMB->c_ipred_mode++)
    {
    // bypass if c_ipred_mode is not allowed
    if ((currMB->c_ipred_mode==VERT_PRED_8 && !mb_available_up) ||
      (currMB->c_ipred_mode==HOR_PRED_8 && !mb_available_left) ||
      (currMB->c_ipred_mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
      continue;
    //===== GET BEST MACROBLOCK MODE =====
    SetModesAndRefframeForBlocks (I4MB);
    if (RDCost_for_macroblocks (lambda_mode, I4MB, &min_rdcost))
      {
      if (input->RCEnable == 1)
      {
        //Rate control
        for (j=0; j<16; j++)
          {
          for(i=0; i<16; i++)
            {
            diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
            }
          }
      }
      store_macroblock_parameters (I4MB);
      }
    }
  set_stored_macroblock_parameters ();
}

void c_avs_enc::encode_one_intra_macroblock_not_rdo()
  {
  Macroblock* currMB      = &img->mb_data[img->current_mb_nr];
  double      lambda_mode;
  int_32_t         cost, dummy;
  int_32_t         max_mcost=(1<<30);
  lambda_mode = QP2QUANT[max(0,img->qp-SHIFT_QP)];
  // reset chroma intra predictor to default
  currMB->c_ipred_mode = DC_PRED_8;
  //cost = (1<<20);            //xzhao 20081106
  cost=0;                      //xzhao 20081106
  /*currMB->cbp = */
  Mode_Decision_for_AVSIntraMacroblock_not_rdo (lambda_mode, &cost);
  //===== set parameters for chosen mode =====
  SetModesAndRefframeForBlocks (I4MB);
  // pecompute all chroma intra prediction modes
  IntraChromaPrediction8x8(NULL, NULL, NULL);
  dummy = 0;
  img->NoResidueDirect = 0;
  ChromaResidualCoding (&dummy);
  if (img->current_mb_nr==0)
    intras=0;
  }
void c_avs_enc::encode_one_inter_macroblock_rdo()
{
  TLS static const int_32_t  mb_mode_table[7]  = {0, 1, 2, 3, P8x8, I16MB, I4MB}; // DO NOT CHANGE ORDER !!!
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  double   qp, lambda_mode, lambda_motion, min_rdcost, rdcost = 0, max_rdcost=1e30;
  int_32_t valid[MAXMODE];
  int_32_t block, index, mode, ref, i, j;
  int_32_t lambda_motion_factor;
  int_32_t fw_mcost, mcost, max_mcost=(1<<30),min_cost,cost8x8,cost;
  int_32_t curr_cbp_blk, cnt_nonz = 0, best_cnt_nonz = 0, best_fw_ref = 0, best_bw_ref = 0, best_pdir = 0;
  int_32_t write_ref   = input->no_multpred>1;
  int_32_t max_ref     = img->nb_references;
  int_32_t **ipredmodes = img->ipredmode;
  int_32_t block_x   = img->block_x;
  int_32_t block_y   = img->block_y;
  int_32_t ***tmpmvs = tmp_mv;
  int_32_t *****allmvs = img->all_mv;
  int_32_t **refar     = refFrArr;
  int_32_t bid_best_fw_ref = 0, bid_best_bw_ref = 0;
  int_32_t mb_available_up, mb_available_left, mb_available_up_left;
  int_32_t adjust_ref = 0;
  max_ref = img->nb_references;
  if (img->number % input->intra_period == 1)//for close GOP
  {
    max_ref = 1;
  }

#ifdef FastME
  decide_intrabk_SAD();
#endif

  //===== SET VALID MODES =====
  valid[I4MB]   = 1;
  valid[I16MB]  = 0;
  valid[0]      = 1;
  valid[1]      = (input->InterSearch16x16);
  valid[2]      = (input->InterSearch16x8);
  valid[3]      = (input->InterSearch8x16);
  valid[4]      = (input->InterSearch8x8);
  valid[5]      = 0;
  valid[6]      = 0;
  valid[7]      = 0;
  valid[P8x8]   = valid[4];

  //===== SET LAGRANGE PARAMETERS =====

  qp = (double)img->qp - SHIFT_QP;

  if (input->successive_Bframe>0)
    lambda_mode   = 0.68 * pow (2, qp/4.0);
  else
    lambda_mode   = 0.85 * pow (2, qp/4.0);

  lambda_motion = sqrt(lambda_mode);

  lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);
  // reset chroma intra predictor to default
  currMB->c_ipred_mode = DC_PRED_8;

  //===== set direct motion vectors =====
#ifdef FastME

  //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
  cost8x8 = 1<<20;
  min_cost=cost8x8;
  best_mode=P8x8;

  mode = INTER16x16;
  if (valid[mode])
  {
    cost = 0;
    block = 0;
    PartitionMotionSearch (mode, block, lambda_motion);
    //--- get cost and reference frame for forward prediction ---
    fw_mcost=max_mcost;
    for (ref=0; ref<max_ref-adjust_ref; ref++)
    {
      mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
      mcost += motion_cost[mode][ref+1][block];
      if (mcost < fw_mcost)
      {
        fw_mcost    = mcost;
        best_fw_ref = ref;
      }
    }
    best_bw_ref = 0;
    best_pdir  = 0;
    cost += fw_mcost;
    //----- set reference frame and direction parameters -----
    best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = best_fw_ref;
    best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
    best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;
    best8x8symref [1][0][0] = best8x8symref [1][1][0] = best8x8symref [1][2][0] = best8x8symref [1][3][0] = bid_best_fw_ref;
    best8x8symref [1][0][1] = best8x8symref [1][1][1] = best8x8symref [1][2][1] = best8x8symref [1][3][1] = bid_best_bw_ref;
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER16x8;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);
      //--- get cost and reference frame for forward prediction ---
      fw_mcost=max_mcost;
      for (ref=0; ref<max_ref; ref++)
      {
        mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
        mcost += motion_cost[mode][ref+1][block];
        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }
      best_bw_ref = 0;
      best_pdir  = 0;
      cost      += fw_mcost;

      //----- set reference frame and direction parameters -----
      best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
      best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
      best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
      best8x8symref   [2][2*block][0] = best8x8symref   [2][2*block+1][0] = bid_best_fw_ref;
      best8x8symref   [2][2*block][1] = best8x8symref   [2][2*block+1][1] = bid_best_bw_ref;

      if (block==0)
        SetRefAndMotionVectors (block, mode, best_fw_ref, best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER8x16;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);
      //--- get cost and reference frame for forward prediction ---
      fw_mcost=max_mcost;
      for (ref=0; ref<max_ref-adjust_ref; ref++)
      {
        mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
        mcost += motion_cost[mode][ref+1][block];
        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }
      best_bw_ref = 0;
      best_pdir  = 0;
      cost      += fw_mcost;

      //----- set reference frame and direction parameters -----
      best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
      best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
      best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
      best8x8symref  [3][block][0] = best8x8symref   [3][block+2][0] = bid_best_fw_ref;
      best8x8symref  [3][block][1] = best8x8symref   [3][block+2][1] = bid_best_bw_ref;

      if (block==0)
        SetRefAndMotionVectors (block, mode, best_fw_ref, best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

#endif

  if (valid[P8x8])
  {
    cbp8x8 = cnt_nonz_8x8 = 0;
    mode = 4;
    best_pdir = 0;
    //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
    for (block=0; block<4; block++)
    {
      //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
      curr_cbp_blk = 0;
      //--- motion estimation for all reference frames ---
      PartitionMotionSearch (mode, block, lambda_motion);
      //--- get cost and reference frame for forward prediction ---
      fw_mcost=max_mcost;
      for (ref=0; ref<max_ref; ref++)
      {
        mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
        mcost += motion_cost[mode][ref+1][block];
        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }

      //--- set variables if best mode has changed ---
        best8x8mode           [block] = mode;
        best8x8pdir     [P8x8][block] = best_pdir;
        best8x8ref      [P8x8][block] = best_fw_ref;
        //--- store number of nonzero coefficients ---
        best_cnt_nonz  = cnt_nonz;  

      //----- set cbp and count of nonzero coefficients ---
      if (best_cnt_nonz)
      {
        cbp8x8        |= (1<<block);
        cnt_nonz_8x8  += best_cnt_nonz;
      }
      if (block<3)
      {
        //===== set motion vectors and reference frames (prediction) =====
        SetRefAndMotionVectors (block, mode, best8x8ref[P8x8][block], best8x8bwref[P8x8][block],best8x8pdir[P8x8][block]);
       } // if (block<3)

    }
  }
#ifndef FastME
  //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====

  mode = INTER16x16; //INTER16x16
  if (valid[mode])
  {
    block = 0;
    PartitionMotionSearch (mode, block, lambda_motion);
    fw_mcost=max_mcost;
    for (ref=0; ref<max_ref; ref++)
    {
      mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
      mcost += motion_cost[mode][ref+1][block];

      if (mcost < fw_mcost)
      {
        fw_mcost    = mcost;
        best_fw_ref = ref;
      }
    }
    best_pdir  = 0;

    best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = best_fw_ref;
    best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
    best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;

  }//INTER16x16

  mode = INTER16x8; //INTER16x8
  if (valid[mode])
  {
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);
      fw_mcost=max_mcost;
      for (ref=0; ref<max_ref; ref++)
      {
        mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
        mcost += motion_cost[mode][ref+1][block];

        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }
      best_pdir  = 0;

      best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
      best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
      best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;

      if (block==0)
        SetRefAndMotionVectors(block, mode, best_fw_ref, best_bw_ref, best_pdir);

    }
  }//INTER16x8

  mode = INTER8x16; //INTER8x16
  if (valid[mode])
  {
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);
      fw_mcost=max_mcost;
      for (ref=0; ref<max_ref; ref++)
      {
        mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0);
        mcost += motion_cost[mode][ref+1][block];

        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }

      best_pdir  = 0;

      best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
      best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
      best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;

      if (block==0)
        SetRefAndMotionVectors (block, mode, best_fw_ref, best_bw_ref, best_pdir);

    }
  }//INTER8x16

#endif
  // Find a motion vector for the Skip mode
  FindSkipModeMotionVector ();
  min_rdcost = max_rdcost;
  //===== GET BEST MACROBLOCK MODE =====
  img->NoResidueDirect = 0;
  if (valid[I4MB])
  {
    IntraChromaPrediction8x8(&mb_available_up, &mb_available_left, &mb_available_up_left);
  }
  for (index=0; index<7; index++)
  {
    mode = mb_mode_table[index];
    //--- for INTER16x16 check all prediction directions ---
    if (valid[mode])
    {
      // bypass if c_ipred_mode not used

      if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
      {
        if (input->RCEnable == 1)
        {
          //Rate control
          if(mode == P8x8)
          {
            for (j=0; j<16; j++)
            {
              for(i=0; i<16; i++)
              {
                diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - mpr8x8[j][i];
              }
            }
          }
          else
          {
            for (j=0; j<16; j++)
            {
              for(i=0; i<16; i++)
              {
                diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
              }
            }
          }
        }
        store_macroblock_parameters (mode);
      }
    }
  }
  if(input->RCEnable == 1)
  {
    img->MADofMB[img->current_mb_nr] = calc_MAD();
    if(input->basicunit<img->total_number_mb)
    {
      img->TotalMADBasicUnit +=img->MADofMB[img->current_mb_nr];
      if ((cbp!=0 || best_mode==I16MB))
        currMB->prev_cbp = 1;
      else
      {
        img->qp -= currMB->delta_qp;
        currMB->delta_qp = 0;
        currMB->qp = img->qp;
        currMB->prev_cbp = 0;
      }
    }
  }
  set_stored_macroblock_parameters ();

  if (img->current_mb_nr==0)
    intras=0;
  if (IS_INTRA(currMB))
    intras++;
  //for intra block consideration
#ifdef FastME
  skip_intrabk_SAD(best_mode, max_ref);
#endif
}

void c_avs_enc::encode_one_inter_macroblock_not_rdo()
{
	TLS static const int_32_t  mb_mode_table[7]  = {0, 1, 2, 3, P8x8, I16MB, I4MB}; // DO NOT CHANGE ORDER !!!
	const int_32_t  num_blocks_in_mb[7] = {0,1,2,2,4,0,0};

	int_32_t it,blk0,blk1,block,loop;  /* iterator */
	int_32_t mode,dummy,lambda_motion_factor;
	int_32_t ref,max_ref,best_fw_ref,best_bw_ref = 0,best_pdir = 0;
	int_32_t min_cost = (1 << 30), temp_cost, fw_cost,mcost;
	double lambda_mode, lambda_motion;
	Macroblock * currMB = &img->mb_data[img->current_mb_nr];

	int_32_t valid_mode[MAXMODE] = {0};
	valid_mode[I4MB] = 1;
	valid_mode[I16MB] = 0;
	valid_mode[0] = 1; //skip mode ?
	valid_mode[1] = input->InterSearch16x16;
	valid_mode[2] = input->InterSearch16x8;
	valid_mode[3] = input->InterSearch8x16;
	valid_mode[4] = input->InterSearch8x8;
	valid_mode[P8x8] = valid_mode[4];
	for (loop=0; loop<5; loop++)
	{
		only_motion_cost[loop][0] = 99999;
		only_motion_cost[loop][1] = 99999;
		only_motion_cost[loop][2] = 99999;
		only_motion_cost[loop][3] = 99999;
	}
#ifdef FastME
	decide_intrabk_SAD();
#endif
	lambda_mode = lambda_motion = QP2QUANT[max(0,img->qp-SHIFT_QP)];
	lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);
	max_ref = img->nb_references;
	if (img->number % input->intra_period == 1)//for close GOP
	{
		max_ref = 1;
	}
	currMB->c_ipred_mode = DC_PRED_8;
	if (img->current_mb_nr == 166)
	{
		loop = 1;
	}
	for (loop = 0; loop < 7; loop++)
	{
		if (valid_mode[mode = mb_mode_table[loop]])
		{
			temp_cost = 0;
			switch (mode)
			{
			case 0:      /* skip mode */
				SetMotionVectorPredictor(img->mv[0][0][0][1],refFrArr,tmp_mv,0,0,0,16,16,0);
				FindSkipModeMotionVector();
				temp_cost = Get_Skip_CostMB(imgY_org,0,img->pix_x,img->pix_y,1,img->mv[0][0][0][1][0],img->mv[0][0][0][1][1],&img->all_mv[0][0][0][0][0],&img->all_mv[0][0][0][0][1],1,1,99999,0);
				break;

			case I4MB:   /* intra mode */
				currMB->c_ipred_mode = DC_PRED_8;
				Mode_Decision_for_AVSIntraMacroblock_not_rdo(lambda_mode,&temp_cost);
				break;
			
			default:    /* inter16*16/16*8/8*16/8*8 mode */

				for (block=0; block<num_blocks_in_mb[loop]; block++)
				{
					/* partition motion search for all refrence frames */
					PartitionMotionSearch((mode<4?mode:4), block,lambda_motion);
					fw_cost = 1<<30;
					for (ref=0; ref<max_ref; ref++)
					{
						mcost = (int_32_t) (2*lambda_motion * min(ref,1));
						mcost += motion_cost[(mode<4?mode:4)][ref+1][block];
						if (mcost < fw_cost)
						{
							fw_cost = mcost;
							best_fw_ref = ref;
						}
					}
					temp_cost += fw_cost;
					blk0 = (mode==1) ? 0 :((mode==2) ? (block<<1) : block);
					blk1 = (mode==1) ? 4 :((mode==2) ? (blk0+2)   : ((mode==3) ?(blk0+3) : blk0+1));
					for (it=blk0; it<blk1; it+=((mode==3) ? 2 : 1))
					{
						best8x8pdir[mode][it] = best_pdir;
						best8x8ref[mode][it] = best_fw_ref, best8x8bwref[mode][it] = best_bw_ref;
					}
					if (mode == P8x8)
					{
						best8x8mode[block] = 4;
						temp_cost += (REF_COST (lambda_motion_factor, B8Mode2Value ((mode<4?mode:4), best_pdir)) - 1);
					}
					if (block < num_blocks_in_mb[loop]-1) /* store the motion info */
					{
						SetRefAndMotionVectors(block,(mode<4?mode:4),best_fw_ref,best_bw_ref,best_pdir);
					}
				}
				break;
			}
			if (temp_cost < min_cost)
			{
				min_cost  = temp_cost;
				best_mode = mode;
			}
#ifdef _FAST_MODE_DECISION_
			if ((mode == SKIP && min_cost < _SKIP_MODE_COST_THRESHOLD_) ||
			    (mode == 2 && abs(only_motion_cost[INTER16x16][0]-(only_motion_cost[INTER16x8][0]+only_motion_cost[INTER16x8][1])) < _COST_THRESHOLD_) ||
			    (mode == 3 && abs((only_motion_cost[INTER8x16][0]+only_motion_cost[INTER8x16][1])-(only_motion_cost[INTER16x8][0]+only_motion_cost[INTER16x8][1])) < _COST_THRESHOLD_))
			{
				break;
			}
#endif
		}
	}

	SetModesAndRefframeForBlocks(best_mode);
	storeMotionInfo(0);
    if (best_mode == I4MB)
	{
		IntraChromaPrediction8x8(NULL,NULL,NULL);
		dummy = 0;
		img->NoResidueDirect = 0;
		ChromaResidualCoding(&dummy);
	}
	else
	{
		img->reconflag = 1;
		LumaResidualCoding();
		dummy = 0;
		img->NoResidueDirect = 0;
		ChromaResidualCoding(&dummy);
	}
	SetMotionVectorsMB(currMB,0);
	if (img->current_mb_nr == 0) intras = 0;
	if (currMB->mb_type == I4MB)
		intras++;
#ifdef FastME
	skip_intrabk_SAD(best_mode, max_ref);
#endif
}


void c_avs_enc::encode_one_b_frame_macroblock_rdo_fast()
  {
  TLS static const int_32_t  b8_mode_table[2]  = {0, 4};         // DO NOT CHANGE ORDER !!!
  TLS static const int_32_t  mb_mode_table[7]  = {0, 1, 2, 3, P8x8}; // DO NOT CHANGE ORDER !!!

  Macroblock* currMB      = &img->mb_data[img->current_mb_nr];
  double      qp, lambda_mode, lambda_motion, min_rdcost, rdcost = 0, max_rdcost=1e30;
  int_32_t         valid[MAXMODE];
  int_32_t         block, index, mode, i0, j0, ref, i, j, k, ctr16x16;
  int_32_t         lambda_motion_factor;
  int_32_t         fw_mcost, bw_mcost, bid_mcost, mcost, max_mcost=(1<<30);
  int_32_t         curr_cbp_blk, cnt_nonz = 0, best_cnt_nonz = 0, best_fw_ref = 0, best_bw_ref = 0, best_pdir = 0;
  int_32_t         cost, min_cost, cost8x8, cost_direct=0, have_direct=0;
  int_32_t         write_ref   = input->no_multpred>1;
  int_32_t         max_ref     = img->nb_references;
  int_32_t         adjust_ref;
  int_32_t         **ipredmodes = img->ipredmode;
  int_32_t         block_x   = img->block_x;
  int_32_t         block_y   = img->block_y;
  int_32_t         ***tmpmvs = tmp_mv;
  int_32_t         *****allmvs = img->all_mv;
  int_32_t         **refar     = refFrArr;
  int_32_t         bid_best_fw_ref = 0, bid_best_bw_ref = 0;
  int_32_t         min_cost_8x8;

  max_ref = min (img->nb_references, img->buf_cycle);
  adjust_ref = 0;
  if(!img->picture_structure)
    {
    if (max_ref > 2)
      max_ref = 2;
    if(img->old_type != img->type)
      max_ref = 1;
    }
  else
    {
    if (max_ref > 1)
      max_ref = 1;
    }
  //===== SET VALID MODES =====
  valid[I4MB]   = 0;
  valid[I16MB]  = 0;
  valid[0]      = 1;
  valid[1]      = (input->InterSearch16x16);
  valid[2]      = (input->InterSearch16x8);
  valid[3]      = (input->InterSearch8x16);
  valid[4]      = (input->InterSearch8x8);
  valid[5]      = 0;
  valid[6]      = 0;
  valid[7]      = 0;
  valid[P8x8]   = valid[4];

  //===== SET LAGRANGE PARAMETERS =====

  qp = (double)img->qp - SHIFT_QP;

  if (input->successive_Bframe>0)
    lambda_mode   = 0.68 * pow (2, qp/4.0) * (max(2.00,min(4.00,(qp / 8.0))));
  else
    lambda_mode   = 0.85 * pow (2, qp/4.0) * 4.0;

  lambda_motion = sqrt(lambda_mode);

  lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);
  // reset chroma intra predictor to default
  currMB->c_ipred_mode = DC_PRED_8;

  //===== set direct motion vectors =====
  Get_IP_direct();
  if (valid[P8x8])
    {
    cost8x8 = 0;
    min_cost = (1<<20);
    //===== store coding state of macroblock =====
    store_coding_state (cs_mb);
    //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
    cbp8x8 = cbp_blk8x8 = cnt_nonz_8x8 = 0;
    for (block=0; block<4; block++)
      {
      //--- set coordinates ---
      j0 = ((block/2)<<3);    //j1 = (j0>>2);
      i0 = ((block%2)<<3);    //i1 = (i0>>2);
      //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
      min_cost_8x8=(1<<20);
      min_rdcost=1e30;

      //check DIRECT mode rdcost
      mode = 0;
      if (valid[mode])
        {
        //for B frame
        best_pdir = 2;
        curr_cbp_blk = 0;
        //--- store coding state before coding with current mode ---
        store_coding_state (cs_cm);

        //--- get and check rate-distortion cost ---
        rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,block, mode, best_pdir, bid_best_fw_ref,bid_best_bw_ref);

        //--- set variables if best mode has changed ---
        if (rdcost < min_rdcost)
          {
          min_cost_8x8                  = cost;
          min_rdcost                    = rdcost;
          best8x8mode           [block] = mode;
          best8x8pdir     [P8x8][block] = best_pdir;
          best8x8ref      [P8x8][block] = best_fw_ref;
          best8x8bwref    [P8x8][block] = best_bw_ref;
          best8x8symref[P8x8][block][0] = bid_best_fw_ref;
          best8x8symref[P8x8][block][1] = bid_best_bw_ref;
          //--- store number of nonzero coefficients ---
          best_cnt_nonz  = cnt_nonz;

          //--- store block cbp ---
          cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
          cbp_blk8x8    |= curr_cbp_blk;

          //--- store coefficients ---
          for (k=0; k< 4; k++)
            for (j=0; j< 2; j++)
              for (i=0; i<65; i++)
                cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

          //--- store reconstruction and prediction ---
          for (j=0; j<8; j++)
            for (i=i0; i<i0+8; i++)
              {
              rec_mbY8x8[j0+j][i] = imgY[img->pix_y+j0+j][img->pix_x+i];
              mpr8x8    [j+j0][i] = (byte)img->mpr[j+j0][i];
              }

            //--- store coding state ---
            store_coding_state (cs_b8);
          } // if (rdcost <= min_rdcost)
        //--- re-set coding state as it was before coding with current mode was performed ---
        reset_coding_state (cs_cm);
        }

      mode = 4;
      if (valid[mode])
        {
        curr_cbp_blk = 0;
        //--- motion estimation for all reference frames ---
        PartitionMotionSearch (mode, block, lambda_motion);
        //--- get cost and reference frame for forward prediction ---
        fw_mcost=max_mcost;
        for (ref=0; ref<max_ref-adjust_ref; ref++)
          {
          mcost = motion_cost[mode][img->picture_structure ? ref+1 : ref+2][block];
          if (mcost < fw_mcost)
            {
            fw_mcost    = mcost;
            best_fw_ref = ref;
            }
          }
        PartitionMotionSearch_bid (mode, block, lambda_motion);
        best_bw_ref = 0;
        bw_mcost   = motion_cost[mode][0][block];
        for (bw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
          {
          mcost = motion_cost[mode][ref][block];
          if (mcost < bw_mcost)
            {
            bw_mcost    = mcost;
            best_bw_ref = ref;
            }
          }
        for (bid_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
          {
          mcost = motion_cost_bid[mode][img->picture_structure ? ref+1 : ref+2][block];
          if (mcost < bid_mcost)
            {
            bid_mcost    = mcost;
            bid_best_fw_ref = ref;
            bid_best_bw_ref = ref;
            }
          }
        //--- get prediction direction ----
        if(fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
          {
          best_pdir = 0;
          cost = fw_mcost;
          best_bw_ref = 0;
          }
        else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
          {
          best_pdir = 1;
          cost = bw_mcost;
          best_fw_ref = 0;
          }
        else
          {
          best_pdir = 2;
          cost = bid_mcost;
          }
        //--- store coding state before coding with current mode ---
        store_coding_state (cs_cm);

        //--- get and check rate-distortion cost ---
        rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,block, mode, best_pdir,
          best_pdir==2?bid_best_fw_ref:best_fw_ref,best_pdir==2?bid_best_bw_ref:best_bw_ref);

        //--- set variables if best mode has changed ---
        if (rdcost < min_rdcost)
          {
          min_cost_8x8                  = cost;
          min_rdcost                    = rdcost;
          best8x8mode           [block] = mode;
          best8x8pdir     [P8x8][block] = best_pdir;
          best8x8ref      [P8x8][block] = best_fw_ref;
          best8x8bwref    [P8x8][block] = best_bw_ref;
          best8x8symref[P8x8][block][0] = bid_best_fw_ref;
          best8x8symref[P8x8][block][1] = bid_best_bw_ref;
          //--- store number of nonzero coefficients ---
          best_cnt_nonz  = cnt_nonz;

          //--- store block cbp ---
          cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
          cbp_blk8x8    |= curr_cbp_blk;

          //--- store coefficients ---
          for (k=0; k< 4; k++)
            for (j=0; j< 2; j++)
              for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

          //--- store reconstruction and prediction ---
          for (j=0; j<8; j++)
            for (i=i0; i<i0+8; i++)
              {
              rec_mbY8x8[j0+j][i] = imgY[img->pix_y+j0+j][img->pix_x+i];
              mpr8x8    [j+j0][i] = (byte)img->mpr[j+j0][i];
              }

            //--- store coding state ---
            store_coding_state (cs_b8);
          } // if (rdcost <= min_rdcost)
        //--- re-set coding state as it was before coding with current mode was performed ---
        reset_coding_state (cs_cm);
        } // if (valid[mode=b8_mode_table[index]])*/

      cost8x8 += min_cost_8x8;

      //----- set cbp and count of nonzero coefficients ---
      if (best_cnt_nonz)
        {
        cbp8x8        |= (1<<block);
        cnt_nonz_8x8  += best_cnt_nonz;
        }

      mode=best8x8mode[block];

      if (block<3)
        {
        //===== set motion vectors and reference frames (prediction) =====
        SetRefAndMotionVectors (block, mode,best8x8pdir[P8x8][block]==2?best8x8symref[P8x8][block][0]:
          best8x8ref[P8x8][block],best8x8pdir[P8x8][block]==2?best8x8symref[P8x8][block][1]:
          best8x8bwref[P8x8][block],best8x8pdir[P8x8][block]);
      //===== re-set reconstructed block =====
      j0   = 8*(block/2);
      i0   = 8*(block%2);
      for (j=0; j<8; j++)
        {
        for (i=i0; i<i0+8; i++)
          {
          imgY[img->pix_y+j0+j][img->pix_x+i] = rec_mbY8x8[j0+j][i];
          }
        }
        } // if (block<3)

      //===== set the coding state after current block =====
      reset_coding_state (cs_b8);
      } // for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)

    //--- re-set coding state (as it was before 8x8 block coding) ---
    reset_coding_state (cs_mb);
    //Rate control
    for (j=0; j<16; j++)
      {
      for(i=0; i<16; i++)
        {
        diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mpr[j][i];
        }
      }
    if(cost8x8<min_cost)
      {
      best_mode = P8x8;
      min_cost = best_mode;
      }
    }
  else // if (valid[P8x8])
    {
    cost8x8 = (1<<20);
    }
  //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
  min_cost=cost8x8;
  best_mode=P8x8;

  mode = INTER16x16;
  if (valid[mode])
    {
    best_pdir = pMbInfo[img->current_mb_nr].pdir;
    if (best_pdir != 2)
    {
      PartitionMotionSearch (mode, block, lambda_motion);
    }
    else
    {
      PartitionMotionSearch_bid (mode, block, lambda_motion);
    }
    //----- set reference frame and direction parameters -----
    best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = 0;
    best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
    best8x8bwref[1][0] = best8x8bwref[1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = 0;
    }

  mode = INTER16x8;
  if (valid[mode])
    {
    cost = 0;
    for (block=0; block<2; block++)
      {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, 0) : 0);

      mcost += motion_cost[mode][1][block];
      fw_mcost    = mcost;
      best_fw_ref = 0;

      //--- get cost for bidirectional prediction ---
      PartitionMotionSearch_bid (mode, block, lambda_motion);

      best_bw_ref = 0;
      bw_mcost   = motion_cost[mode][0][block];
      bid_mcost = motion_cost_bid[mode][best_fw_ref+1][block];
      //--- get prediction direction ----
      if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
        best_pdir = 0;
        best_bw_ref = 0;
        cost += fw_mcost;
        }
      else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
        best_pdir = 1;
        cost += bw_mcost;
        best_fw_ref = 0;
        }
      else
        {
        best_pdir = 2;
        cost += bid_mcost;
        }
      //----- set reference frame and direction parameters -----
      best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
      best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
      best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors(block, mode, best_pdir==2? bid_best_fw_ref: best_fw_ref,
        best_pdir==2? bid_best_bw_ref: best_bw_ref, best_pdir);
      }
    if (cost < min_cost)
      {
      best_mode = mode;
      min_cost  = cost;
      }
    }

  mode = INTER8x16;
  if (valid[mode])
    {
    cost = 0;
    for (block=0; block<2; block++)
      {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, 0) : 0);

      mcost += motion_cost[mode][1][block];
      fw_mcost    = mcost;
      best_fw_ref = 0;

      //--- get cost for bidirectional prediction ---
      PartitionMotionSearch_bid (mode, block, lambda_motion);

      best_bw_ref = 0;
      bw_mcost   = motion_cost[mode][0][block];
      bid_mcost = motion_cost_bid[mode][best_fw_ref+1][block];
      //--- get prediction direction ----
      if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
        best_pdir = 0;
        best_bw_ref = 0;
        cost += fw_mcost;
        }
      else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
        best_pdir = 1;
        cost += bw_mcost;
        best_fw_ref = 0;
        }
      else
        {
        best_pdir = 2;
        cost += bid_mcost;
        }
      //----- set reference frame and direction parameters -----
      best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
      best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
      best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors(block, mode, best_pdir==2? bid_best_fw_ref: best_fw_ref,
        best_pdir==2? bid_best_bw_ref: best_bw_ref, best_pdir);
      }
    if (cost < min_cost)
      {
      best_mode = mode;
      min_cost  = cost;
      }
    }


  min_rdcost = max_rdcost;

  //===== GET BEST MACROBLOCK MODE =====
  for (ctr16x16=0, index=0; index<7; index++)
    {
    mode = mb_mode_table[index];
    img->NoResidueDirect = 0;
    if (valid[mode])
      {
      // bypass if c_ipred_mode not used
      SetModesAndRefframeForBlocks (mode);

      if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
        {
        //Rate control
        if(mode == P8x8)
          {
          for (j=0; j<16; j++)
            {
            for(i=0; i<16; i++)
              {
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - mpr8x8[j][i];
              }
            }
          }
        else
          {
          for (j=0; j<16; j++)
            {
            for(i=0; i<16; i++)
              {
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
              }
            }
          }
        store_macroblock_parameters (mode);
        }

      }

    if (mode == 0 && currMB->cbp && (currMB->cbp&15) != 15)
      {
      img->NoResidueDirect = 1;
      if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
        {
        //Rate control
        for (j=0; j<16; j++)
          for(i=0; i<16; i++)
            diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
        store_macroblock_parameters (mode);
        }
      }
    }

  // Rate control
  set_stored_macroblock_parameters ();

  if (img->current_mb_nr==0)
    intras=0;
  //for intra block consideration
#ifdef FastME
  skip_intrabk_SAD(best_mode, max_ref-adjust_ref);
#endif
  }

void c_avs_enc::encode_one_b_frame_macroblock_rdo()
{
  TLS static const int_32_t  b8_mode_table[2]  = {0, 4};         // DO NOT CHANGE ORDER !!!
  TLS static const int_32_t  mb_mode_table[7]  = {0, 1, 2, 3, P8x8}; // DO NOT CHANGE ORDER !!!

  Macroblock* currMB      = &img->mb_data[img->current_mb_nr];
  double   qp, lambda_mode, lambda_motion, min_rdcost, rdcost = 0, max_rdcost=1e30;
  int_32_t valid[MAXMODE];
  int_32_t block, index, mode, i0, j0, i, j, k, ctr16x16;
  int_32_t lambda_motion_factor;
  int_32_t fw_mcost, bw_mcost, bid_mcost, mcost, max_mcost=(1<<30);
  int_32_t curr_cbp_blk, cnt_nonz = 0, best_cnt_nonz = 0, best_fw_ref = 0, best_bw_ref = 0, best_pdir = 0;
  int_32_t cost, min_cost, cost8x8, cost_direct=0, have_direct=0;
  int_32_t write_ref   = input->no_multpred>1;
  int_32_t max_ref     = img->nb_references;
  int_32_t **ipredmodes = img->ipredmode;
  int_32_t block_x   = img->block_x;
  int_32_t block_y   = img->block_y;
  int_32_t ***tmpmvs = tmp_mv;
  int_32_t *****allmvs = img->all_mv;
  int_32_t **refar     = refFrArr;
  int_32_t bid_best_fw_ref = 0, bid_best_bw_ref = 0;
  int_32_t adjust_ref;
  int_32_t mb_available_up, mb_available_left, mb_available_up_left;
  max_ref = min (img->nb_references, img->buf_cycle);
  if(!img->picture_structure)
  {
	  adjust_ref = 0;
    if (max_ref > 2)
      max_ref = 2;
    if(img->old_type != img->type)
      max_ref = 1;
  }
  else
  {
	  adjust_ref = 0;
    if (max_ref > 1)
      max_ref = 1;
  }
#ifdef FastME
  decide_intrabk_SAD();
#endif
  //===== SET VALID MODES =====
  valid[I4MB]   = 1;
  valid[I16MB]  = 0;
  valid[0]      = 1;
  valid[1]      = (input->InterSearch16x16);
  valid[2]      = (input->InterSearch16x8);
  valid[3]      = (input->InterSearch8x16);
  valid[4]      = (input->InterSearch8x8);
  valid[5]      = 0;
  valid[6]      = 0;
  valid[7]      = 0;
  valid[P8x8]   = valid[4];

  //===== SET LAGRANGE PARAMETERS =====

  qp = (double)img->qp - SHIFT_QP;

  lambda_mode   = 0.68 * pow (2, qp/4.0) * (max(2.00,min(4.00,(qp / 8.0))));
  lambda_motion = sqrt(lambda_mode);

  lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);
  // reset chroma intra predictor to default
  currMB->c_ipred_mode = DC_PRED_8;
  memset(best8x8symref, 0, MAXMODE*4*2*sizeof(int_32_t));
  //===== set direct motion vectors =====
  Get_IP_direct();

#ifdef FastME
  int ref;
  cost8x8 = 1<<20;
  //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
  min_cost=cost8x8;
  best_mode=P8x8;

  mode = INTER16x16;
  if (valid[mode])
  {
    cost=0;
    block=0;
    PartitionMotionSearch (mode, block, lambda_motion);
    //--- get cost and reference frame for forward prediction ---
    fw_mcost=max_mcost;
    for (ref=0; ref<max_ref-adjust_ref; ref++)
    {
      mcost = motion_cost[mode][img->picture_structure ? ref+1 : ref+2][block];
      if (mcost < fw_mcost)
      {
        fw_mcost    = mcost;
        best_fw_ref = ref;
      }
    }
    best_bw_ref = 0;

    if (img->number % input->intra_period == 0)  //for close GOP
    {
      best_pdir = 0;
      cost += fw_mcost;
    }
    else
    {
      //--- get cost for bidirectional prediction ---
      PartitionMotionSearch_bid (mode, block, lambda_motion);
      best_bw_ref = 0;
      bw_mcost   = motion_cost[mode][0][block];
      for (bw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
      {
        mcost = motion_cost[mode][ref][block];
        if (mcost < bw_mcost)
        {
          bw_mcost    = mcost;
          best_bw_ref = ref;
        }
      }

      for (bid_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
      {
        mcost = motion_cost_bid[mode][img->picture_structure ? ref+1 : ref+2][block];
        if (mcost < bid_mcost)
        {
          bid_mcost    = mcost;
          bid_best_fw_ref = ref;
          bid_best_bw_ref = ref;
        }
      }
      //--- get prediction direction ----
      if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
      {
        best_pdir = 0;
        best_bw_ref = 0;
        cost += fw_mcost;
      }
      else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
      {
        best_pdir = 1;
        cost += bw_mcost;
        best_fw_ref = 0;
      }
      else
      {
        best_pdir = 2;
        cost += bid_mcost;
      }
    }
    //----- set reference frame and direction parameters -----
    best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = best_fw_ref;
    best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
    best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;
    best8x8symref [1][0][0] = best8x8symref [1][1][0] = best8x8symref [1][2][0] = best8x8symref [1][3][0] = bid_best_fw_ref;
    best8x8symref [1][0][1] = best8x8symref [1][1][1] = best8x8symref [1][2][1] = best8x8symref [1][3][1] = bid_best_bw_ref;

    //--- set reference frames and motion vectors ---
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER16x8;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      for (fw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
      {
        mcost = motion_cost[mode][img->picture_structure ? ref+1 : ref+2][block];
        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }
      best_bw_ref = 0;

      if (img->number % input->intra_period == 0)  //for close GOP
      {
        best_pdir = 0;
        cost += fw_mcost;
      }
      else
      {
        //--- get cost for bidirectional prediction ---
        PartitionMotionSearch_bid (mode, block, lambda_motion);
        best_bw_ref = 0;
        bw_mcost   = motion_cost[mode][0][block];
        for (bw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
        {
          mcost = motion_cost[mode][ref][block];
          if (mcost < bw_mcost)
          {
            bw_mcost    = mcost;
            best_bw_ref = ref;
          }
        }

        for (bid_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
        {
          mcost = motion_cost_bid[mode][img->picture_structure ? ref+1 : ref+2][block];
          if (mcost < bid_mcost)
          {
            bid_mcost    = mcost;
            bid_best_fw_ref = ref;
            bid_best_bw_ref = ref;
          }
        }
        //--- get prediction direction ----
        if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
          best_pdir = 0;
          best_bw_ref = 0;
          cost += fw_mcost;
        }
        else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
          best_pdir = 1;
          cost += bw_mcost;
          best_fw_ref = 0;
        }
        else
        {
          best_pdir = 2;
          cost += bid_mcost;
        }
      }
      //----- set reference frame and direction parameters -----
      best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
      best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
      best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
      best8x8symref   [2][2*block][0] = best8x8symref   [2][2*block+1][0] = bid_best_fw_ref;
      best8x8symref   [2][2*block][1] = best8x8symref   [2][2*block+1][1] = bid_best_bw_ref;

      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors (block, mode, best_pdir==2?bid_best_fw_ref:best_fw_ref,
        best_pdir==2?bid_best_bw_ref:best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER8x16;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      for (fw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
      {
        mcost = motion_cost[mode][img->picture_structure ? ref+1 : ref+2][block];
        if (mcost < fw_mcost)
        {
          fw_mcost    = mcost;
          best_fw_ref = ref;
        }
      }
      best_bw_ref = 0;
      if (img->number % input->intra_period == 0)  //for close GOP
      {
        best_pdir = 0;
        cost += fw_mcost;
      }
      else
      {
        //--- get cost for bidirectional prediction ---
        PartitionMotionSearch_bid (mode, block, lambda_motion);
        best_bw_ref = 0;
        bw_mcost   = motion_cost[mode][0][block];
        for (bw_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
        {
          mcost = motion_cost[mode][ref][block];
          if (mcost < bw_mcost)
          {
            bw_mcost    = mcost;
            best_bw_ref = ref;
          }
        }

        for (bid_mcost=max_mcost, ref=0; ref<max_ref-adjust_ref; ref++)
        {
          mcost = motion_cost_bid[mode][img->picture_structure ? ref+1 : ref+2][block];
          if (mcost < bid_mcost)
          {
            bid_mcost    = mcost;
            bid_best_fw_ref = ref;
            bid_best_bw_ref = ref;
          }
        }
        //--- get prediction direction ----
        if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
          best_pdir = 0;
          best_bw_ref = 0;
          cost += fw_mcost;
        }
        else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
          best_pdir = 1;
          cost += bw_mcost;
          best_fw_ref = 0;
        }
        else
        {
          best_pdir = 2;
          cost += bid_mcost;
        }
      }
      //----- set reference frame and direction parameters -----
      best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
      best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
      best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
      best8x8symref  [3][block][0] = best8x8symref   [3][block+2][0] = bid_best_fw_ref;
      best8x8symref  [3][block][1] = best8x8symref   [3][block+2][1] = bid_best_bw_ref;

      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors (block, mode, best_pdir==2?bid_best_fw_ref:best_fw_ref,
        best_pdir==2?bid_best_bw_ref:best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }
#endif

  if (valid[P8x8])
  {
    cost8x8 = 0;
    min_cost = (1<<20);
    //===== store coding state of macroblock =====
    store_coding_state (cs_mb);
    //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
    cbp8x8 = cbp_blk8x8 = cnt_nonz_8x8 = 0;
    for (block=0; block<4; block++)
    {
      //--- set coordinates ---
      j0 = ((block/2)<<3); i0 = ((block%2)<<3);
      //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
      min_rdcost=1e30;
      //check DIRECT mode rdcost
      mode = 0;
      if (valid[mode])
      {
        //for B frame
        best_pdir = 2;
        curr_cbp_blk = 0;
        //--- store coding state before coding with current mode ---
        store_coding_state (cs_cm);
        //--- get and check rate-distortion cost ---
        rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,block, mode, best_pdir, bid_best_fw_ref,bid_best_bw_ref);
        //--- set variables if best mode has changed ---
        if (rdcost < min_rdcost)
        {
          min_rdcost                    = rdcost;
          best8x8mode           [block] = mode;
          best8x8pdir     [P8x8][block] = best_pdir;
          best8x8ref      [P8x8][block] = best_fw_ref;
          best8x8bwref    [P8x8][block] = best_bw_ref;
          best8x8symref[P8x8][block][0] = bid_best_fw_ref;
          best8x8symref[P8x8][block][1] = bid_best_bw_ref;
          //--- store number of nonzero coefficients ---
          best_cnt_nonz  = cnt_nonz;

          //--- store block cbp ---
          cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
          cbp_blk8x8    |= curr_cbp_blk;

          //--- store coefficients ---
          for (k=0; k< 4; k++)
            for (j=0; j< 2; j++)
              for (i=0; i<65; i++)
                cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

          //--- store reconstruction and prediction ---
          for (j=j0; j<j0+8; j++)
          {
            memcpy(&mpr8x8[j][i0], &img->mpr[j][i0], 8*sizeof(int_16_t));
            memcpy(&rec_mbY8x8[j][i0], &imgY[img->pix_y+j][img->pix_x+i0], 8*sizeof(byte));
          }
          //--- store coding state ---
          store_coding_state (cs_b8);
        } // if (rdcost <= min_rdcost)
        //--- re-set coding state as it was before coding with current mode was performed ---
        reset_coding_state (cs_cm);
      }

      mode = 4;
      if (valid[mode])
      {
        curr_cbp_blk = 0;
        //forward
        //--- motion estimation for all reference frames ---
        PartitionMotionSearch (mode, block, lambda_motion);
        //--- get cost and reference frame for forward prediction ---
        fw_mcost = motion_cost[mode][1][block];
        bw_mcost = motion_cost[mode][0][block];
        best_fw_ref = 0;
        //backward
        best_bw_ref = 0;
        if (img->number % input->intra_period == 0)  //for close GOP
        {
          best_pdir = 0;
          cost = fw_mcost;
        }
        else
        {
          PartitionMotionSearch_bid (mode, block, lambda_motion);
          bid_mcost = motion_cost_bid[mode][1][block];
          bid_best_fw_ref = bid_best_bw_ref = 0;
          //--- get prediction direction ----
          if(fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
          {
            best_pdir = 0;
            cost = fw_mcost;
            best_bw_ref = 0;
          }
          else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
          {
            best_pdir = 1;
            cost = bw_mcost;
            best_fw_ref = 0;
          }
          else
          {
            best_pdir = 2;
            cost = bid_mcost;
          }
        }
        //--- store coding state before coding with current mode ---
        store_coding_state (cs_cm);

        //--- get and check rate-distortion cost ---
        rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,block, mode, best_pdir, best_pdir==2?bid_best_fw_ref:best_fw_ref,best_pdir==2?bid_best_bw_ref:best_bw_ref);

        //--- set variables if best mode has changed ---
        if (rdcost < min_rdcost)
        {
          min_rdcost                    = rdcost;
          best8x8mode           [block] = mode;
          best8x8pdir     [P8x8][block] = best_pdir;
          best8x8ref      [P8x8][block] = best_fw_ref;
          best8x8bwref    [P8x8][block] = best_bw_ref;
          best8x8symref[P8x8][block][0] = bid_best_fw_ref;
          best8x8symref[P8x8][block][1] = bid_best_bw_ref;
          //--- store number of nonzero coefficients ---
          best_cnt_nonz  = cnt_nonz;

          //--- store block cbp ---
          cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
          cbp_blk8x8    |= curr_cbp_blk;

          //--- store coefficients ---
          for (k=0; k< 4; k++)
            for (j=0; j< 2; j++)
              for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

          //--- store reconstruction and prediction ---
          for (j=j0; j<j0+8; j++)
          {
            memcpy(&mpr8x8[j][i0], &img->mpr[j][i0], 8*sizeof(int_16_t));
            memcpy(&rec_mbY8x8[j][i0], &imgY[img->pix_y+j][img->pix_x+i0], 8*sizeof(byte));
          }
          //--- store coding state ---
          store_coding_state (cs_b8);
        } // if (rdcost <= min_rdcost)
        //--- re-set coding state as it was before coding with current mode was performed ---
        reset_coding_state (cs_cm);
      } // if (valid[mode=b8_mode_table[index]])*/

      cost8x8 += (int_32_t)min_rdcost;

      //----- set cbp and count of nonzero coefficients ---
      if (best_cnt_nonz)
      {
        cbp8x8        |= (1<<block);
        cnt_nonz_8x8  += best_cnt_nonz;
      }

      mode=best8x8mode[block];

      if (block<3)
      {
        //===== set motion vectors and reference frames (prediction) =====
        SetRefAndMotionVectors (block, mode,best8x8pdir[P8x8][block]==2?best8x8symref[P8x8][block][0]:best8x8ref[P8x8][block],best8x8pdir[P8x8][block]==2?best8x8symref[P8x8][block][1]:best8x8bwref[P8x8][block],best8x8pdir[P8x8][block]);
        //===== re-set reconstructed block =====
        j0   = 8*(block/2);
        i0   = 8*(block%2);
        for (j=j0; j<j0+8; j++)
        {
          memcpy(&imgY[img->pix_y+j][img->pix_x+i0], &rec_mbY8x8[j][i0], 8*sizeof(byte));
        }
      } // if (block<3)

      //===== set the coding state after current block =====
      reset_coding_state (cs_b8);
    } // for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)

    //--- re-set coding state (as it was before 8x8 block coding) ---
    reset_coding_state (cs_mb);
    //Rate control
    if (input->RCEnable == 1)
    {
      for (j=0; j<16; j++)
      {
        for(i=0; i<16; i++)
        {
          diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mpr[j][i];
        }
      }
    }
    if(cost8x8<min_cost)
    {
      best_mode = P8x8;
      min_cost = cost8x8;
    }
  }
  else // if (valid[P8x8])
  {
    cost8x8 = (1<<20);
  }
#ifndef FastME
  //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
  min_cost=cost8x8;
  best_mode=P8x8;

  mode = INTER16x16;
  if (valid[mode])
  {
    cost = 0;
    block = 0;
    PartitionMotionSearch (mode, block, lambda_motion);
    //--- get cost and reference frame for forward prediction ---
    mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, 0) : 0);
    mcost += motion_cost[mode][1][block];
    fw_mcost    = mcost;
    best_fw_ref = 0;
    if (img->number % input->intra_period == 0)  //for close GOP
    {
      best_pdir = 0;
      cost += fw_mcost;
      best_bw_ref = 0;
    }
    else
    {
      //--- get cost for bidirectional prediction ---
      PartitionMotionSearch_bid (mode, block, lambda_motion);

      best_bw_ref = 0;
      bw_mcost   = motion_cost[mode][0][block];
      bid_mcost = motion_cost_bid[mode][best_fw_ref+1][block];
      //--- get prediction direction ----
      if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
      {
        best_pdir = 0;
        best_bw_ref = 0;
        cost += fw_mcost;
      }
      else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
      {
        best_pdir = 1;
        cost += bw_mcost;
        best_fw_ref = 0;
      }
      else
      {
        best_pdir = 2;
        cost += bid_mcost;
      }
    }
    //----- set reference frame and direction parameters -----
    best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = best_fw_ref;
    best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
    best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;

    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER16x8;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      fw_mcost = motion_cost[mode][1][block];
      best_fw_ref = 0;
      if (img->number % input->intra_period == 0)  //for close GOP
      {
        best_pdir = 0;
        cost += fw_mcost;
        best_bw_ref = 0;
      }
      else
      {
        //--- get cost for bidirectional prediction ---
        PartitionMotionSearch_bid (mode, block, lambda_motion);
        best_bw_ref = 0;
        bw_mcost  = motion_cost[mode][0][block];
        bid_mcost = motion_cost_bid[mode][best_fw_ref+1][block];
        //--- get prediction direction ----
        if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
          best_pdir = 0;
          best_bw_ref = 0;
          cost += fw_mcost;
        }
        else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
          best_pdir = 1;
          cost += bw_mcost;
          best_fw_ref = 0;
        }
        else
        {
          best_pdir = 2;
          cost += bid_mcost;
        }
      }
      //----- set reference frame and direction parameters -----
      best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
      best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
      best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors(block, mode, best_pdir==2? bid_best_fw_ref: best_fw_ref,
        best_pdir==2? bid_best_bw_ref: best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }

  mode = INTER8x16;
  if (valid[mode])
  {
    cost = 0;
    for (block=0; block<2; block++)
    {
      PartitionMotionSearch (mode, block, lambda_motion);

      //--- get cost and reference frame for forward prediction ---
      mcost  = (write_ref ? REF_COST_FWD (lambda_motion_factor, 0) : 0);

      mcost += motion_cost[mode][1][block];
      fw_mcost    = mcost;
      best_fw_ref = 0;
      if (img->number % input->intra_period == 0)  //for close GOP
      {
        best_pdir = 0;
        cost += fw_mcost;
        best_bw_ref = 0;
      }
      else
      {
        //--- get cost for bidirectional prediction ---
        PartitionMotionSearch_bid (mode, block, lambda_motion);

        best_bw_ref = 0;
        bw_mcost   = motion_cost[mode][0][block];
        bid_mcost = motion_cost_bid[mode][best_fw_ref+1][block];
        //--- get prediction direction ----
        if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
        {
          best_pdir = 0;
          best_bw_ref = 0;
          cost += fw_mcost;
        }
        else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
        {
          best_pdir = 1;
          cost += bw_mcost;
          best_fw_ref = 0;
        }
        else
        {
          best_pdir = 2;
          cost += bid_mcost;
        }
      }
      //----- set reference frame and direction parameters -----
      best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
      best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
      best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
      //--- set reference frames and motion vectors ---
      if (block==0)
        SetRefAndMotionVectors(block, mode, best_pdir==2? bid_best_fw_ref: best_fw_ref,
        best_pdir==2? bid_best_bw_ref: best_bw_ref, best_pdir);
    }
    if (cost < min_cost)
    {
      best_mode = mode;
      min_cost  = cost;
    }
  }
#endif

  min_rdcost = max_rdcost;

  //===== GET BEST MACROBLOCK MODE =====
  if (valid[I4MB])
  {
    IntraChromaPrediction8x8(&mb_available_up, &mb_available_left, &mb_available_up_left);
  }
  for (ctr16x16=0, index=0; index<7; index++)
  {
    mode = mb_mode_table[index];
    //--- for INTER16x16 check all prediction directions ---
    if (mode==1)
    {
      best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = ctr16x16;
      if (img->number % input->intra_period != 0)  //for close GOP
      {
        if (ctr16x16 < 2) index--;
        ctr16x16++;
      }
    }
    img->NoResidueDirect = 0;
    if (valid[mode])
    {
      // bypass if c_ipred_mode not used
      SetModesAndRefframeForBlocks (mode);

      if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
      {
        if (input->RCEnable == 1)
        {
          if(mode == P8x8)
          {
            for (j=0; j<16; j++)
            {
              for(i=0; i<16; i++)
              {
                diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - mpr8x8[j][i];
              }
            }
          }
          else
          {
            for (j=0; j<16; j++)
            {
              for(i=0; i<16; i++)
              {
                diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
              }
            }
          }
        }
        store_macroblock_parameters (mode);
      }
    }

    if (mode == 0 && currMB->cbp && (currMB->cbp&15) != 15)
    {
      img->NoResidueDirect = 1;
      if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
      {
        if (input->RCEnable == 1)
        {
          //Rate control
          for (j=0; j<16; j++)
            for(i=0; i<16; i++)
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i] - pred[j][i];
        }
        store_macroblock_parameters (mode);
      }
    }
  }
  set_stored_macroblock_parameters ();

  if (img->current_mb_nr==0)
    intras=0;
  //for intra block consideration
#ifdef FastME
  skip_intrabk_SAD(best_mode, max_ref-adjust_ref);
#endif
}

void c_avs_enc::encode_one_b_frame_macroblock_not_rdo()
{
	TLS static const int_32_t mb_mode_table[7] = {0,1,2,3,P8x8,I16MB,I4MB};
	TLS static const int_32_t num_blocks_in_mb[7] = {0,1,2,2,4,0,0};

	int_32_t it,blk0,blk1,block,loop;
	int_32_t mode,sub_mode,best_mode,dummy;
	int_32_t min_cost,temp_cost,fw_cost,bw_cost,bid_cost,temp_cost_8x8[2];
	int_32_t max_ref, adjust_ref,best_fw_ref = 0,best_bw_ref = 0, best_bid_fw_ref = 0, best_bid_bw_ref = 0, best_pdir;
	int_32_t lambda_motion_factor;
	double lambda_motion,lambda_mode;
	Macroblock * currMB = &img->mb_data[img->current_mb_nr];
	bool valid[MAXMODE] ={0};
	
	valid[I4MB] = 1;
	valid[I16MB] = 0;
	valid[0] = 1;
	valid[1] = input->InterSearch16x16;
	valid[2] = input->InterSearch16x8;
	valid[3] = input->InterSearch8x16;
	valid[4] = input->InterSearch8x8;
	valid[P8x8] = valid[4];
#ifdef FastME
	decide_intrabk_SAD();
#endif	
	max_ref = min (img->nb_references, img->buf_cycle);
	if(!img->picture_structure)
	{
		adjust_ref = 0;
		if (max_ref > 2)
			max_ref = 2;
		if(img->old_type != img->type)
			max_ref = 1;
	}
	else
	{
		adjust_ref = 0;
		if (max_ref > 1)
			max_ref = 1;
	}
	lambda_mode = lambda_motion = QP2QUANT[max(0,img->qp-SHIFT_QP)];
	lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);

	currMB->c_ipred_mode = DC_PRED_8;

	for (it=0; it<5; it++)
	{
		only_motion_cost[it][0] = 99999;
		only_motion_cost[it][1] = 99999;
		only_motion_cost[it][2] = 99999;
		only_motion_cost[it][3] = 99999;
	}

	min_cost = 1 << 30;

	for (loop=0; loop<7; loop++)
	{
		if (valid[mode = mb_mode_table[loop]])
		{
			temp_cost = 0;
			switch (mode)
			{
			case DIRECT:
				Get_IP_direct();
				temp_cost = Get_Direct_CostMB(lambda_mode);
				temp_cost -= (int_32_t)floor(16*lambda_motion+0.4999);
				break;

			case I4MB:
				Mode_Decision_for_AVSIntraMacroblock_not_rdo(lambda_mode,&temp_cost);
				break;

			default:
				for (block=0; block<num_blocks_in_mb[loop]; block++)
				{
					sub_mode = mode < 4 ? mode : 4;
					PartitionMotionSearch(sub_mode,block,lambda_motion);
					PartitionMotionSearch_bid(sub_mode,block,lambda_motion);
					fw_cost = motion_cost[sub_mode][1][block];
					bw_cost = motion_cost[sub_mode][0][block];
					bid_cost = motion_cost_bid[sub_mode][1][block];
					best_pdir = ((fw_cost<=bw_cost && fw_cost<=bid_cost) ? 0 :  ((bw_cost<=fw_cost && bw_cost<=bid_cost) ? 1 : 2));
					blk0 = (mode==1) ? 0 :((mode==2) ? (block<<1) : block);
					blk1 = (mode==1) ? 4 :((mode==2) ? (blk0+2)   : ((mode==3) ?(blk0+3) : blk0+1));
					for (it=blk0; it<blk1; it+=((mode==3) ? 2 : 1))
					{
						best8x8pdir[mode][it] = best_pdir;
						best8x8ref[mode][it] = best_fw_ref, best8x8bwref[mode][it] = best_bw_ref;
						best8x8symref[mode][it][0] = best_bid_fw_ref, best8x8symref[mode][it][1] = best_bid_bw_ref;
					}
					if (mode != P8x8)
					{
						temp_cost += ((best_pdir == 0) ? fw_cost : ((best_pdir==1) ? bw_cost : bid_cost));
					}
					else
					{
						temp_cost_8x8[0] = 1 << 30;
						if (valid[DIRECT])
						{
							temp_cost_8x8[0] = Get_Direct_Cost8x8(block,lambda_mode);
							temp_cost_8x8[0] += (REF_COST (lambda_motion_factor, B8Mode2Value (0, 2)) - 1);
						}
						temp_cost_8x8[1] = ((best_pdir == 0) ? fw_cost : ((best_pdir==1) ? bw_cost : bid_cost));
						temp_cost_8x8[1] += (REF_COST (lambda_motion_factor, B8Mode2Value (4, best_pdir)) - 1);
						if (temp_cost_8x8[0] < temp_cost_8x8[1])
						{
							temp_cost += temp_cost_8x8[0];
							best8x8mode[block] = DIRECT;
							best8x8pdir[mode][block] = 2;
						}
						else
						{
							temp_cost += temp_cost_8x8[1];
							best8x8mode[block] = 4;
							best8x8pdir[mode][block] = best_pdir;
						}
						sub_mode  = best8x8mode[block];
						best_pdir = best8x8pdir[mode][block];
					}
					if (block < num_blocks_in_mb[loop]-1)
					{
						SetRefAndMotionVectors(block,mode==P8x8?sub_mode:mode,best_pdir==2 ? best_bid_fw_ref : best_fw_ref,
							best_pdir==2? best_bid_bw_ref : best_bw_ref, best_pdir);
					}
				}
				break;
			}
			if (temp_cost < min_cost)
			{
				min_cost = temp_cost;
				best_mode = mode;
			}
#ifdef _FAST_MODE_DECISION_
#define _DIRECT_THRESHOLD_COST_ 200
			if ((mode == DIRECT && min_cost < _DIRECT_THRESHOLD_COST_) ||
				(mode == 2 && abs(only_motion_cost[INTER16x16][0]-(only_motion_cost[INTER16x8][0]+only_motion_cost[INTER16x8][1])) < _COST_THRESHOLD_) ||
				(mode == 3 && abs((only_motion_cost[INTER8x16][0]+only_motion_cost[INTER8x16][1])-(only_motion_cost[INTER16x8][0]+only_motion_cost[INTER16x8][1])) < _COST_THRESHOLD_))
			{
				break;
			}
#endif
		}
	}
	SetModesAndRefframeForBlocks(best_mode);
	storeMotionInfo(0);
	if (best_mode == I4MB)
	{
		IntraChromaPrediction8x8(NULL,NULL,NULL);
		dummy = 0;
		img->NoResidueDirect = 0;
		ChromaResidualCoding(&dummy);
	}
	else
	{
		img->reconflag = 1;
		LumaResidualCoding();
		dummy = 0;
		img->NoResidueDirect = 0;
		ChromaResidualCoding(&dummy);
	}
	SetMotionVectorsMB(currMB,1);
	if (img->current_mb_nr == 0)
		intras = 0;
	if (best_mode == I4MB)
		intras++;
#ifdef FastME
	skip_intrabk_SAD(best_mode, max_ref-adjust_ref);
#endif  
}

int_32_t  c_avs_enc::Mode_Decision_for_AVSIntraMacroblock_not_rdo (double lambda,  int_32_t* total_cost)
{
  int_32_t  cbp=0, b8;
  int_32_t  min_cost;

  int_32_t ipmode,best_ipmode,cost,loc_cbp,loc_cbp_blk;
  int_32_t nonzero;
  int_32_t block_x;
  int_32_t block_y;
  int_32_t pic_pix_x;
  int_32_t pic_pix_y;
  int_32_t pic_block_x;
  int_32_t pic_block_y;
  __declspec(align(16)) int_16_t tmp_block_88[8][8];
  int_32_t upMode;
  int_32_t leftMode;
  int_32_t mostProbableMode;
  Macroblock *currMB;
  int_32_t MBRowSize = img->width / MB_BLOCK_SIZE;
  int_32_t loc_cost = (int_32_t)floor(6.0 * lambda + 0.4999);//,best_AVS_cost;
  byte  *d1, *d2, *d3, *d4;
  int_16_t *d5, *d6;

  for (min_cost=0, b8=0; b8<4; b8++)
  {
    currMB=img->mb_data+img->current_mb_nr;

    block_x    =8*(b8&1);
    block_y    =8*(b8>>1);
    pic_pix_x  =img->pix_x+block_x;
    pic_pix_y  =img->pix_y+block_y;
    pic_block_x=pic_pix_x>>3;
    pic_block_y=pic_pix_y>>3;
    //min_rdcost =1e30;
    min_cost = (1<<20);

    upMode           = img->ipredmode[pic_block_x+1][pic_block_y  ];
    leftMode         = img->ipredmode[pic_block_x  ][pic_block_y+1];
    mostProbableMode = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

    //===== INTRA PREDICTION FOR 4x4 BLOCK =====
    intrapred_luma_AVS(pic_pix_x,pic_pix_y);
    cost = 0;
    //===== LOOP OVER ALL INTRA PREDICTION MODES =====
    for (ipmode=0;ipmode<NO_INTRA_PMODE;ipmode++)
    {
      if(img->available_intra_mode[ipmode] == 1)
      {
        d1 = &imgY_org[pic_pix_y][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+1][pic_pix_x];
        d3 = &img->mprr[ipmode][0][0];
        d4 = &img->mprr[ipmode][1][0];
        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm0, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm0, xmm7

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm6, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm6, xmm7

          psadbw      xmm0, xmm6
        }

        d1 = &imgY_org[pic_pix_y+2][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+3][pic_pix_x];
        d3 = &img->mprr[ipmode][2][0];
        d4 = &img->mprr[ipmode][3][0];
        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm1, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm1, xmm7

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm6, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm6, xmm7

          psadbw      xmm1, xmm6
        }

        d1 = &imgY_org[pic_pix_y+4][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+5][pic_pix_x];
        d3 = &img->mprr[ipmode][4][0];
        d4 = &img->mprr[ipmode][5][0];
        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm2, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm2, xmm7

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm6, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm6, xmm7

          psadbw      xmm2, xmm6
        }

        d1 = &imgY_org[pic_pix_y+6][pic_pix_x];
        d2 = &imgY_org[pic_pix_y+7][pic_pix_x];
        d3 = &img->mprr[ipmode][6][0];
        d4 = &img->mprr[ipmode][7][0];
        __asm
        {
          mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
          movdqu    xmm3, xmmword ptr [esi]
          mov      esi,  dword ptr [d2]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm3, xmm7

          mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
          movdqu    xmm6, xmmword ptr [esi]
          mov      esi,  dword ptr [d4]
          movdqu    xmm7, xmmword ptr [esi]
          punpcklqdq  xmm6, xmm7

          psadbw      xmm3, xmm6
          paddw       xmm0, xmm1
          paddw       xmm0, xmm2
          paddw       xmm0, xmm3
          movdqa    xmm1, xmm0
          punpckhqdq  xmm0, xmm0
          paddw       xmm0, xmm1
          pextrw      eax,  xmm0, 0
          mov      cost, eax
        }
        cost += (ipmode == mostProbableMode) ? 0 : (int_32_t)floor(3 * lambda );

        if (cost < min_cost)
        {
          best_ipmode = ipmode;
          min_cost    = cost;
        }
      }
    }
    //===== set intra mode prediction =====
    img->ipredmode[pic_block_x+1][pic_block_y+1] = best_ipmode;
//#define _OUTPUT_TRACE_
#ifdef _OUTPUT_TRACE_
        FILE *pf_trace = NULL;
        int i, j;
        if (frame_no ==2 && img->current_mb_nr == 25 && b8==0)
          {
          pf_trace = fopen("enc_trace.txt", "a");
      fprintf(pf_trace, "best_ipmode: %5d\n", best_ipmode);
          for (j=0; j<8; j++)
            {
            for (i=0; i<8; i++)
              {
              //fprintf(pf_trace, "%5d ", curr_blk[j][i]);
          fprintf(pf_trace, "%5d ", img->mpr[block_y +j][block_x +i]);
              }
            fprintf(pf_trace, "\n");
            }
          fprintf(pf_trace, "\n");
          fclose(pf_trace);
          //exit(0);
          }
#endif
    currMB->intra_pred_modes[b8] = mostProbableMode == best_ipmode ? -1 : best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1;
    d1 = &imgY_org[pic_pix_y][pic_pix_x];
    d2 = &imgY_org[pic_pix_y+1][pic_pix_x];
    d3 = &img->mprr[best_ipmode][0][0];
    d4 = &img->mprr[best_ipmode][1][0];
    d5 = &img->mpr[block_y][block_x];
    d6 = &img->mpr[block_y+1][block_x];
    /*if (img->current_mb_nr == 26 && frame_no == 0 && b8==1)
    {
      int x, y;
      for (y=0; y<8; y++)
      {
        for (x=0; x<8; x++)
        {
          printf("%4d ", img->mprr[best_ipmode][y][x]);
        }
        printf("\n");
      }
      printf("\n");
      exit(0);
    }  */
    __asm
    {
      mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
      movdqu    xmm0, xmmword ptr [esi]
      mov      esi,  dword ptr [d2]
      movdqu    xmm1, xmmword ptr [esi]

      mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
      movdqu    xmm2, xmmword ptr [esi]
      mov      esi,  dword ptr [d4]
      movdqu    xmm3, xmmword ptr [esi]

      pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
      psubw       xmm0, xmm2
      psubw       xmm1, xmm3
      movdqa      tmp_block_88, xmm0
      movdqa      tmp_block_88+16, xmm1
      mov      esi,  dword ptr [d5]
      movdqa      xmmword ptr [esi],   xmm2
      mov      esi,  dword ptr [d6]
      movdqa      xmmword ptr [esi],   xmm3
    }

    d1 = &imgY_org[pic_pix_y+2][pic_pix_x];
    d2 = &imgY_org[pic_pix_y+3][pic_pix_x];
    d3 = &img->mprr[best_ipmode][2][0];
    d4 = &img->mprr[best_ipmode][3][0];
    d5 = &img->mpr[block_y+2][block_x];
    d6 = &img->mpr[block_y+3][block_x];
    __asm
    {
      mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
      movdqu    xmm0, xmmword ptr [esi]
      mov      esi,  dword ptr [d2]
      movdqu    xmm1, xmmword ptr [esi]

      mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
      movdqu    xmm2, xmmword ptr [esi]
      mov      esi,  dword ptr [d4]
      movdqu    xmm3, xmmword ptr [esi]

      pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
      psubw       xmm0, xmm2
      psubw       xmm1, xmm3
      movdqa      tmp_block_88+32, xmm0
      movdqa      tmp_block_88+48, xmm1
      mov      esi,  dword ptr [d5]
      movdqa      xmmword ptr [esi],   xmm2
      mov      esi,  dword ptr [d6]
      movdqa      xmmword ptr [esi],   xmm3
    }

    d1 = &imgY_org[pic_pix_y+4][pic_pix_x];
    d2 = &imgY_org[pic_pix_y+5][pic_pix_x];
    d3 = &img->mprr[best_ipmode][4][0];
    d4 = &img->mprr[best_ipmode][5][0];
    d5 = &img->mpr[block_y+4][block_x];
    d6 = &img->mpr[block_y+5][block_x];
    __asm
    {
      mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
      movdqu    xmm0, xmmword ptr [esi]
      mov      esi,  dword ptr [d2]
      movdqu    xmm1, xmmword ptr [esi]

      mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
      movdqu    xmm2, xmmword ptr [esi]
      mov      esi,  dword ptr [d4]
      movdqu    xmm3, xmmword ptr [esi]

      pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
      psubw       xmm0, xmm2
      psubw       xmm1, xmm3
      movdqa      tmp_block_88+64, xmm0
      movdqa      tmp_block_88+80, xmm1
      mov      esi,  dword ptr [d5]
      movdqa      xmmword ptr [esi],   xmm2
      mov      esi,  dword ptr [d6]
      movdqa      xmmword ptr [esi],   xmm3
    }

    d1 = &imgY_org[pic_pix_y+6][pic_pix_x];
    d2 = &imgY_org[pic_pix_y+7][pic_pix_x];
    d3 = &img->mprr[best_ipmode][6][0];
    d4 = &img->mprr[best_ipmode][7][0];
    d5 = &img->mpr[block_y+6][block_x];
    d6 = &img->mpr[block_y+7][block_x];
    __asm
    {
      mov      esi,  dword ptr [d1]  //read in imgY_org[pic_pix_y][pic_pix_x]
      movdqu    xmm0, xmmword ptr [esi]
      mov      esi,  dword ptr [d2]
      movdqu    xmm1, xmmword ptr [esi]

      mov      esi,  dword ptr [d3]  //read in img->mprr[ipmode][][]
      movdqu    xmm2, xmmword ptr [esi]
      mov      esi,  dword ptr [d4]
      movdqu    xmm3, xmmword ptr [esi]

      pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
      psubw       xmm0, xmm2
      psubw       xmm1, xmm3
      movdqa      tmp_block_88+96, xmm0
      movdqa      tmp_block_88+112, xmm1
      mov      esi,  dword ptr [d5]
      movdqa      xmmword ptr [esi],   xmm2
      mov      esi,  dword ptr [d6]
      movdqa      xmmword ptr [esi],   xmm3
    }
    loc_cbp = loc_cbp_blk = 0;
    avs_dct_sse(tmp_block_88);
    // xzhao {
    //scanquant_B8   (img->qp,4,b8,tmp_block_88,0,&loc_cbp,&loc_cbp_blk); // '|4' indicate intra for quantization
    scanquant_B8_recon   (img->qp,4,b8,tmp_block_88,0,&loc_cbp,&loc_cbp_blk); // '|4' indicate intra for quantization
    // xzhao }
    //scanquant returns a SCR value!.....
    nonzero=(loc_cbp!=0);
    currMB->cbp|=loc_cbp;
    currMB->cbp_blk|=loc_cbp_blk;
    // xzhao { 2007.7.24
    min_cost += loc_cost;
    // xzhao }

    if (nonzero)
    {
      cbp |= (1<<b8);
    }
    *total_cost += min_cost;
  }
  return cbp;
}
