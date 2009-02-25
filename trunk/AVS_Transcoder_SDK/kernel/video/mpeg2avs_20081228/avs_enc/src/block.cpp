/*$T block.cpp GC 1.140 10/28/07 20:29:46 */


/*$6
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>


#include "global.h"
#include "const_data.h"
#include "sse_header.h"

TLS byte clip0d[16];
TLS byte clip255d[16];

/*
 =======================================================================================================================
 =======================================================================================================================
 */
__inline void c_avs_enc::avs_const_initialize_block()
{
  int i=0;
  for(i=0;i<16;i=i+2)
    {
    clip0d[i]=clip255d[i]=0;
    clip0d[i+1]=128;
    clip255d[i+1]=127;
    }
  clip0q = _mm_loadu_si128((const __m128i *)clip0d);
  clip255q = _mm_loadu_si128((const __m128i *)clip255d);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
__inline __m128i c_avs_enc::avs_clip_0_255_w(__m128i xmm0)
{
  xmm0 = _mm_adds_epi16(xmm0, clip255q);
  xmm0 = _mm_subs_epi16(xmm0, clip255q);
  xmm0 = _mm_adds_epi16(xmm0, clip0q);
  xmm0 = _mm_subs_epi16(xmm0, clip0q);
  return xmm0;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_16_t c_avs_enc::scanquant_B8
(
  int_32_t  qp,
  int_32_t  mode,
  int_32_t  block8x8,
  int_16_t  curr_blk[B8_SIZE][B8_SIZE],
  int_32_t  scrFlag,
  int_32_t  *cbp,
  int_32_t  *cbp_blk
)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_16_t  coeff_cost = 0;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  if(img->reconflag)
  {
    coeff_cost = scanquant_B8_recon(qp, mode, block8x8, curr_blk, scrFlag, cbp, cbp_blk);
  }
  else
  {
    coeff_cost = scanquant_B8_cost(qp, mode, block8x8, curr_blk, scrFlag, cbp, cbp_blk);
  }

  return coeff_cost;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_16_t c_avs_enc::scanquant_B8_recon
(
  int_32_t  qp,
  int_32_t  mode,
  int_32_t  block8x8,
  int_16_t  curr_blk[B8_SIZE][B8_SIZE],
  int_32_t  scrFlag,
  int_32_t  *cbp,
  int_32_t  *cbp_blk
)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_16_t  run;
  int_16_t  xx, yy;
  int_16_t  icoef, ipos;
  int_16_t  b8_y = (block8x8 / 2) << 3;
  int_16_t  b8_x = (block8x8 % 2) << 3;
  int_16_t  coeff_cost = 0;
  int_16_t  *ACLevel;
  int_16_t  *ACRun;
  int_16_t  curr_val;
  /* xzhao { 2007.10.09 */
  __m64    tmp;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*
   * xzhao } ;
   * if(img->current_mb_nr==263 && frame_no==5) ;
   * { ;
   * printf("block8x8: %d\n",block8x8);
   * * for(jj=0;
   * jj<8;
   * jj++) ;
   * { ;
   * for(ii=0;
   * ii<8;
   * ii++) ;
   * { ;
   * printf("%3d ",curr_blk[jj][ii]);
   * * } ;
   * printf("\n");
   * * } ;
   * } ;
   * Quantization ;
   * quant_B8(qp, mode, curr_blk);
   */
  avs_quant_sse(qp, mode, curr_blk);
  avs_const_initialize_block();
  /* General block information */
  ACLevel = img->cofAC[block8x8][0][0];
  ACRun = img->cofAC[block8x8][0][1];

  for(xx = 0; xx < 64; xx = xx + 8)
  {
    *(__m128i *) (ACRun + xx) = _mm_setzero_si128();
    *(__m128i *) (ACLevel + xx) = _mm_setzero_si128();
  }

  ACRun[64] = ACLevel[64] = 0;
  run = -1;
  ipos = 0;

  for(icoef = 0; icoef < 64; icoef++)
  {
    run++;
    tmp = *(__m64 *) (AVS_SCAN[img->picture_structure][icoef]);
    xx = _mm_extract_pi16(tmp, 0);
    yy = _mm_extract_pi16(tmp, 2);

    curr_val = curr_blk[yy][xx];

    if(curr_val != 0)
    {
      ACLevel[ipos] = curr_val;
      ACRun[ipos] = run;

      if(scrFlag && absm(ACLevel[ipos]) == 1)
        coeff_cost += (scrFlag == 1) ? 1 : AVS_COEFF_COST[run];
      else
        coeff_cost = MAX_VALUE1; /* block has to be saved */

      run = -1;
      ipos++;
    }
  }
  {
    avs_dequant_sse(qp, curr_blk);
  }

  if(ipos > 0)  /* there are coefficients */
    if(block8x8 <= 3)
      (*cbp) |= (1 << block8x8);
    else
      (*cbp) |= (1 << (block8x8));
  {
    avs_idct_sse(curr_blk);
    if(block8x8 <= 3)
    {
      /*~~~~~~~~~~~~~~~*/
      __m128i xmm0, xmm1;
      __m128i curr_val_w;
      __m64  curr_val_b;
      /*~~~~~~~~~~~~~~~*/
        {
//#define _OUTPUT_TRACE_
#ifdef _OUTPUT_TRACE_
        FILE *pf_trace = NULL;
        int i, j;
        if (frame_no ==2 && img->current_mb_nr == 25 && block8x8==0)
          {
          pf_trace = fopen("enc_trace.txt", "a");
          for (j=0; j<8; j++)
            {
            for (i=0; i<8; i++)
              {
              //fprintf(pf_trace, "%5d ", curr_blk[j][i]);
          fprintf(pf_trace, "%5d ", img->mpr[b8_y +j][b8_x +i]);
              }
            fprintf(pf_trace, "\n");
            }
          fprintf(pf_trace, "\n");
          fclose(pf_trace);
          //exit(0);
          }
#endif
        }
      for(yy = 0; yy < 8; yy++)
      {
        xmm0 = _mm_loadu_si128((const __m128i *) (img->mpr[b8_y + yy] + b8_x));
        xmm1 = _mm_loadu_si128((const __m128i *) curr_blk[yy]);
        curr_val_w = _mm_add_epi16(xmm0, xmm1);
        curr_val_w = avs_clip_0_255_w(curr_val_w);
        _mm_storeu_si128((__m128i *) img->m7[yy], curr_val_w);
        curr_val_w = _mm_packus_epi16(curr_val_w, curr_val_w);

        /*
         * mm_storeu_si128(test_blk[yy], curr_val_w);
         */
        curr_val_b = _mm_movepi64_pi64(curr_val_w);
        *(__m64 *) (imgY[img->pix_y + b8_y + yy] + img->pix_x + b8_x) = curr_val_b;

        /*
         * (__m64 *)(rst_blk[yy])=curr_val_b;
         */
      }
//#define _OUTPUT_TRACE_
#ifdef _OUTPUT_TRACE_
      FILE *pf_trace = NULL;
      int i, j;
      if (frame_no ==2 && img->current_mb_nr == 25 && block8x8==0)
        {
        pf_trace = fopen("enc_trace.txt", "a");
        for (j=0; j<8; j++)
          {
          for (i=0; i<8; i++)
            {
            fprintf(pf_trace, "%5d ", imgY[img->pix_y + b8_y + j][img->pix_x+b8_x+i]);
            }
          fprintf(pf_trace, "\n");
          }
        fprintf(pf_trace, "\n");
        fclose(pf_trace);
        exit(0);
        }
#endif
      _mm_empty();
    }
    else
    {
      /*~~~~~~~~~~~~~~~*/
      __m128i xmm0, xmm1;
      __m128i curr_val_w;
      __m64  curr_val_b;
      /*~~~~~~~~~~~~~~~*/

      for(yy = 0; yy < 8; yy++)
      {
        xmm0 = _mm_loadu_si128((const __m128i *) img->mpr[yy]);
        xmm1 = _mm_loadu_si128((const __m128i *) curr_blk[yy]);
        curr_val_w = _mm_add_epi16(xmm0, xmm1);
        curr_val_w = avs_clip_0_255_w(curr_val_w);
        _mm_storeu_si128((__m128i *) img->m7[yy], curr_val_w);
        curr_val_w = _mm_packus_epi16(curr_val_w, curr_val_w);

        /*
         * mm_storeu_si128(test_blk[yy], curr_val_w);
         */
        curr_val_b = _mm_movepi64_pi64(curr_val_w);
        *(__m64 *) (imgUV[block8x8 - 4][img->pix_c_y + yy] + img->pix_c_x) = curr_val_b;

        /*
         * (__m64 *)(rst_blk[yy])=curr_val_b;
         */
      }

      _mm_empty();
    }
  }

  return coeff_cost;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_16_t c_avs_enc::scanquant_B8_cost
(
  int_32_t  qp,
  int_32_t  mode,
  int_32_t  block8x8,
  int_16_t  curr_blk[B8_SIZE][B8_SIZE],
  int_32_t  scrFlag,
  int_32_t  *cbp,
  int_32_t  *cbp_blk
)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_16_t  run;
  int_16_t  xx, yy;
  int_16_t  icoef, ipos;
  int_16_t  b8_y = (block8x8 / 2) << 3;
  int_16_t  b8_x = (block8x8 % 2) << 3;
  int_16_t  coeff_cost = 0;
  int_16_t  *ACLevel;
  int_16_t  *ACRun;
  int_16_t  curr_val;
  /* xzhao { 2007.10.09 */
  __m64    tmp;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*
   * xzhao } ;
   * Quantization ;
   * quant_B8(qp, mode, curr_blk);
   */
  avs_quant_sse(qp, mode, curr_blk);

  /* General block information */
  ACLevel = img->cofAC[block8x8][0][0];
  ACRun = img->cofAC[block8x8][0][1];

  for(xx = 0; xx < 64; xx = xx + 8)
  {
    *(__m128i *) (ACRun + xx) = _mm_setzero_si128();
    *(__m128i *) (ACLevel + xx) = _mm_setzero_si128();
  }

  ACRun[64] = ACLevel[64] = 0;
  run = -1;
  ipos = 0;

  for(icoef = 0; icoef < 64; icoef++)
  {
    run++;
    tmp = *(__m64 *) (AVS_SCAN[img->picture_structure][icoef]);
    xx = _mm_extract_pi16(tmp, 0);
    yy = _mm_extract_pi16(tmp, 2);

    curr_val = curr_blk[yy][xx];

    if(curr_val != 0)
    {
      ACLevel[ipos] = curr_val;
      ACRun[ipos] = run;

      if(scrFlag && absm(ACLevel[ipos]) == 1)
        coeff_cost += (scrFlag == 1) ? 1 : AVS_COEFF_COST[run];
      else
        coeff_cost = MAX_VALUE1; /* block has to be saved */

      run = -1;
      ipos++;
    }
  }

  _mm_empty();
  if(ipos > 0)  /* there are coefficients */
    if(block8x8 <= 3)
      (*cbp) |= (1 << block8x8);
    else
      (*cbp) |= (1 << (block8x8));

  return coeff_cost;
}

/*
 =======================================================================================================================
        Function: Calculate SAD or SATD for a prediction error block of size iSizeX x iSizeY. Input: Output: Return:
        Attention:
 =======================================================================================================================
 */
int_32_t c_avs_enc::find_sad_8x8
(
  int_32_t  iMode,
  int_32_t  iSizeX,
  int_32_t  iSizeY,
  int_32_t  iOffX,
  int_32_t  iOffY,
  int_32_t  m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE]
)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j;
  int_32_t  bmode;
  int_32_t  ishift = 0;
  int_32_t  sad = 0;
  /*~~~~~~~~~~~~~~~~~~~~~~~*/

  assert((iSizeX + iOffX) <= MB_BLOCK_SIZE);
  assert((iSizeY + iOffY) <= MB_BLOCK_SIZE);

  /* m7[y,j,line][x,i,pixel] */
  switch(iMode)
  {
  case 0:      /* SAD */
    for(j = iOffY; j < iOffY + iSizeY; j++)
      for(i = iOffX; i < iOffX + iSizeX; i++) sad += absm(m7[j][i]);
    break;

  case 1:      /* SATD */
    bmode = iSizeX + iSizeY;
    if(bmode < 24)  /* 8x8 */
    {
      sad = sad_hadamard(iSizeY, iSizeX, iOffY, iOffX, m7);  /* Attention: sad_hadamard() is X/Y
                     * flipped */
      ishift = 2;
      sad = (sad + (1 << (ishift - 1))) >> ishift;
    }
    else      /* 8x16-16x16 */
    {
      switch(bmode)
      {
      case 24:  /* 16x8 8x16 */
        sad = sad_hadamard(8, 8, iOffY, iOffX, m7);
        sad += sad_hadamard
        (
          8,
          8,
          iOffY + ((iSizeY == 16) ? 8 : 0),
          iOffX + ((iSizeX == 16) ? 8 : 0),
          m7
        );
        ishift = 2;
        break;

      case 32:  /* 16x16 */
        sad = sad_hadamard(8, 8, 0, 0, m7);
        sad += sad_hadamard(8, 8, 8, 0, m7);
        sad += sad_hadamard(8, 8, 0, 8, m7);
        sad += sad_hadamard(8, 8, 8, 8, m7);
        ishift = 2;
        break;

      default:
        assert(0 == 1);
      }

      sad = (sad + (1 << (ishift - 1))) >> ishift;
    }
    break;

  default:
    assert(0 == 1);    /* more switches may be added here later */
  }

  return sad;
}

/*
 =======================================================================================================================
        Function: calculates the SAD of the Hadamard transformed block of size iSizeX*iSizeY. Block may have an offset
        of (iOffX,iOffY). If offset!=0 then iSizeX/Y has to be <=8. Input: Output: Return: Attention:
 =======================================================================================================================
 */
int_32_t c_avs_enc::sad_hadamard
(
  int_32_t  iSizeX,
  int_32_t  iSizeY,
  int_32_t  iOffX,
  int_32_t  iOffY,
  int_32_t  m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE]
)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, ii;
  int_32_t  m1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int_32_t  m2[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int_32_t  m3[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  int_32_t  sad = 0;
  int_32_t  iy[MB_BLOCK_SIZE] =
  {
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY,
    iOffY
  };
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /* in this routine, cols are j,y and rows are i,x */
  assert(((iOffX == 0) || (iSizeX <= 8)) && ((iOffY == 0) || (iSizeY <= 8)));

  for(j = 1; j < iSizeY; j++) iy[j] += j;

  /* vertical transform */
  if(iSizeY == 4)
  {
    for(i = 0; i < iSizeX; i++)
    {
      ii = i + iOffX;

      m1[i][0] = m7[ii][iy[0]] + m7[ii][iy[3]];
      m1[i][1] = m7[ii][iy[1]] + m7[ii][iy[2]];
      m1[i][2] = m7[ii][iy[1]] - m7[ii][iy[2]];
      m1[i][3] = m7[ii][iy[0]] - m7[ii][iy[3]];

      m3[i][0] = m1[i][0] + m1[i][1];
      m3[i][1] = m1[i][0] - m1[i][1];
      m3[i][2] = m1[i][2] + m1[i][3];
      m3[i][3] = m1[i][3] - m1[i][2];
    }
  }
  else
  {
    for(i = 0; i < iSizeX; i++)
    {
      ii = i + iOffX;

      m1[i][0] = m7[ii][iy[0]] + m7[ii][iy[4]];
      m1[i][1] = m7[ii][iy[1]] + m7[ii][iy[5]];
      m1[i][2] = m7[ii][iy[2]] + m7[ii][iy[6]];
      m1[i][3] = m7[ii][iy[3]] + m7[ii][iy[7]];
      m1[i][4] = m7[ii][iy[0]] - m7[ii][iy[4]];
      m1[i][5] = m7[ii][iy[1]] - m7[ii][iy[5]];
      m1[i][6] = m7[ii][iy[2]] - m7[ii][iy[6]];
      m1[i][7] = m7[ii][iy[3]] - m7[ii][iy[7]];

      m2[i][0] = m1[i][0] + m1[i][2];
      m2[i][1] = m1[i][1] + m1[i][3];
      m2[i][2] = m1[i][0] - m1[i][2];
      m2[i][3] = m1[i][1] - m1[i][3];
      m2[i][4] = m1[i][4] + m1[i][6];
      m2[i][5] = m1[i][5] + m1[i][7];
      m2[i][6] = m1[i][4] - m1[i][6];
      m2[i][7] = m1[i][5] - m1[i][7];

      m3[i][0] = m2[i][0] + m2[i][1];
      m3[i][1] = m2[i][0] - m2[i][1];
      m3[i][2] = m2[i][2] + m2[i][3];
      m3[i][3] = m2[i][2] - m2[i][3];
      m3[i][4] = m2[i][4] + m2[i][5];
      m3[i][5] = m2[i][4] - m2[i][5];
      m3[i][6] = m2[i][6] + m2[i][7];
      m3[i][7] = m2[i][6] - m2[i][7];
    }
  }

  /* horizontal transform */
  if(iSizeX == 4)
  {
    for(j = 0; j < iSizeY; j++)
    {
      m1[0][j] = m3[0][j] + m3[3][j];
      m1[1][j] = m3[1][j] + m3[2][j];
      m1[2][j] = m3[1][j] - m3[2][j];
      m1[3][j] = m3[0][j] - m3[3][j];

      m2[0][j] = m1[0][j] + m1[1][j];
      m2[1][j] = m1[0][j] - m1[1][j];
      m2[2][j] = m1[2][j] + m1[3][j];
      m2[3][j] = m1[3][j] - m1[2][j];

      for(i = 0; i < iSizeX; i++) sad += absm(m2[i][j]);
    }
  }
  else
  {
    for(j = 0; j < iSizeY; j++)
    {
      m2[0][j] = m3[0][j] + m3[4][j];
      m2[1][j] = m3[1][j] + m3[5][j];
      m2[2][j] = m3[2][j] + m3[6][j];
      m2[3][j] = m3[3][j] + m3[7][j];
      m2[4][j] = m3[0][j] - m3[4][j];
      m2[5][j] = m3[1][j] - m3[5][j];
      m2[6][j] = m3[2][j] - m3[6][j];
      m2[7][j] = m3[3][j] - m3[7][j];

      m1[0][j] = m2[0][j] + m2[2][j];
      m1[1][j] = m2[1][j] + m2[3][j];
      m1[2][j] = m2[0][j] - m2[2][j];
      m1[3][j] = m2[1][j] - m2[3][j];
      m1[4][j] = m2[4][j] + m2[6][j];
      m1[5][j] = m2[5][j] + m2[7][j];
      m1[6][j] = m2[4][j] - m2[6][j];
      m1[7][j] = m2[5][j] - m2[7][j];

      m2[0][j] = m1[0][j] + m1[1][j];
      m2[1][j] = m1[0][j] - m1[1][j];
      m2[2][j] = m1[2][j] + m1[3][j];
      m2[3][j] = m1[2][j] - m1[3][j];
      m2[4][j] = m1[4][j] + m1[5][j];
      m2[5][j] = m1[4][j] - m1[5][j];
      m2[6][j] = m1[6][j] + m1[7][j];
      m2[7][j] = m1[6][j] - m1[7][j];

      for(i = 0; i < iSizeX; i++) sad += absm(m2[i][j]);
    }
  }

  return(sad);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_32_t c_avs_enc::writeLumaCoeffAVS_B8(int_32_t b8, int_32_t intra)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  no_bits = 0;
  int_32_t  mb_nr = img->current_mb_nr;
  Macroblock  *currMB = &img->mb_data[mb_nr];
  const int_32_t  cbp = currMB->cbp;
  int_32_t  *bitCount = currMB->bitcounter;
  SyntaxElement  *currSE = &img->MB_SyntaxElements[currMB->currSEnr];
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  const char(*AVS_2DVLC_table_intra)[26][27];
  const char(*AVS_2DVLC_table_inter)[26][27];

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  inumblk;  /* number of blocks per CBP */
  int_32_t  inumcoeff;  /* number of coeffs per block */
  int_32_t  icoef;    /* current coefficient */
  int_32_t  ipos;    /* current position in cof_AVS */
  int_32_t  run, level, abs_level;
  int_16_t  *ACLevel;
  int_16_t  *ACRun;
  int_32_t  symbol2D;
  int_32_t  escape_level_diff;
  int_32_t  tablenum;
  int_32_t  incVlc_intra[7] = { 0, 1, 2, 4, 7, 10, 3000 };
  int_32_t  incVlc_inter[7] = { 0, 1, 2, 3, 6, 9, 3000 };
  int_32_t  blk_table[4] = { 0, -1, -1, -1 };
  int_32_t  blkmode2ctx[4] = { LUMA_8x8, LUMA_8x4, LUMA_4x8, LUMA_4x4 };
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  inumblk = 1;
  inumcoeff = 65;      /* all positions + EOB */

  AVS_2DVLC_table_intra = AVS_2DVLC_INTRA;
  AVS_2DVLC_table_inter = AVS_2DVLC_INTER;
  if(cbp & (1 << b8))
  {
    /* code all symbols */
    level = 1;    /* get inside loop */
    ipos = 0;
    ACLevel = img->cofAC[b8][0][0];
    ACRun = img->cofAC[b8][0][1];

    for(icoef = 0; icoef < inumcoeff; icoef++)    /* count coeffs */
      if(!ACLevel[icoef]) break;
    tablenum = 0;
    icoef--;
    if(intra)
    {
      for(icoef; icoef >= 0; icoef--)
      {
        level = ACLevel[icoef];
        abs_level = abs(level);
        run = ACRun[icoef];
        symbol2D = CODE2D_ESCAPE_SYMBOL;  /* symbol for out-of-table */
        if(level > -27 && level < 27 && run < 26)
        {
          if(tablenum == 0)
            symbol2D = AVS_2DVLC_table_intra[tablenum][run][abs_level - 1];
          else
            symbol2D = AVS_2DVLC_table_intra[tablenum][run][abs_level];
          if(symbol2D >= 0 && level < 0) symbol2D++;
          if(symbol2D < 0)
            symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
        }
        else
        {
          symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
        }

        currSE->type = SE_LUM_AC_INTER;
        currSE->value1 = symbol2D;
        currSE->value2 = 0;
        currSE->golomb_grad = VLC_Golomb_Order[0][tablenum][0];
        currSE->golomb_maxlevels = VLC_Golomb_Order[0][tablenum][1];

        writeSyntaxElement_GOLOMB(currSE, currBitStream);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        no_bits += currSE->len;
        currSE++;
        currMB->currSEnr++;

        if(symbol2D >= CODE2D_ESCAPE_SYMBOL)
        {
          currSE->type = SE_LUM_AC_INTER;
          currSE->golomb_grad = 1;
          currSE->golomb_maxlevels = 11;
          escape_level_diff = abs_level - ((run > MaxRun[0][tablenum]) ? 1 : RefAbsLevel[tablenum][run]);
          currSE->value1 = escape_level_diff;

          writeSyntaxElement_GOLOMB(currSE, currBitStream);
          bitCount[BITS_COEFF_Y_MB] += currSE->len;
          no_bits += currSE->len;
          currSE++;
          currMB->currSEnr++;
        }

        if(abs_level > incVlc_intra[tablenum])
        {
          if(abs_level <= 2)
            tablenum = abs_level;
          else if(abs_level <= 4)
            tablenum = 3;
          else if(abs_level <= 7)
            tablenum = 4;
          else if(abs_level <= 10)
            tablenum = 5;
          else
            tablenum = 6;
        }
      }

      /* coding EOB */
      abs_level = level = 0;
      run = 0;
      symbol2D = AVS_2DVLC_table_intra[tablenum][run][abs_level];
      currSE->type = SE_LUM_AC_INTER;
      currSE->value1 = symbol2D;
      currSE->value2 = 0;
      currSE->golomb_grad = VLC_Golomb_Order[0][tablenum][0];
      currSE->golomb_maxlevels = VLC_Golomb_Order[0][tablenum][1];

      writeSyntaxElement_GOLOMB(currSE, currBitStream);
      bitCount[BITS_COEFF_Y_MB] += currSE->len;
      no_bits += currSE->len;
      currSE++;  /* proceed to next SE */
      currMB->currSEnr++;
    }
    else      /* !intra */
    {
      for(icoef; icoef >= 0; icoef--)
      {
        level = ACLevel[icoef];
        abs_level = abs(level);
        run = ACRun[icoef];
        symbol2D = CODE2D_ESCAPE_SYMBOL;  /* symbol for out-of-table */
        if(level > -27 && level < 27 && run < 26)
        {
          if(tablenum == 0)
            symbol2D = AVS_2DVLC_table_inter[tablenum][run][abs_level - 1];
          else
            symbol2D = AVS_2DVLC_table_inter[tablenum][run][abs_level];
          if(symbol2D >= 0 && level < 0) symbol2D++;
          if(symbol2D < 0)
            symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
        }
        else
        {
          symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
        }

        /*
         * currSE->type = SE_LUM_AC_INTER;
         */
        currSE->value1 = symbol2D;
        currSE->value2 = 0;

        currSE->golomb_grad = VLC_Golomb_Order[1][tablenum][0];
        currSE->golomb_maxlevels = VLC_Golomb_Order[1][tablenum][1];

        writeSyntaxElement_GOLOMB(currSE, currBitStream);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        no_bits += currSE->len;
        currSE++;
        currMB->currSEnr++;

        if(symbol2D >= CODE2D_ESCAPE_SYMBOL)
        {
          /*
           * currSE->type = SE_LUM_AC_INTER;
           */
          currSE->golomb_grad = 0;
          currSE->golomb_maxlevels = 11;
          escape_level_diff = abs_level - ((run > MaxRun[1][tablenum]) ? 1 : RefAbsLevel[tablenum + 7][run]);
          currSE->value1 = escape_level_diff;

          writeSyntaxElement_GOLOMB(currSE, currBitStream);
          bitCount[BITS_COEFF_Y_MB] += currSE->len;
          no_bits += currSE->len;
          currSE++;      /* proceed to next SE */
          currMB->currSEnr++;

          /* end */
        }

        if(abs_level > incVlc_inter[tablenum])
        {
          if(abs_level <= 3)
            tablenum = abs_level;
          else if(abs_level <= 6)
            tablenum = 4;
          else if(abs_level <= 9)
            tablenum = 5;
          else
            tablenum = 6;
        }
      }

      abs_level = level = 0;
      run = 0;
      symbol2D = AVS_2DVLC_table_inter[tablenum][run][abs_level];

      /*
       * currSE->type = SE_LUM_AC_INTER;
       */
      currSE->value1 = symbol2D;
      currSE->value2 = 0;
      currSE->golomb_grad = VLC_Golomb_Order[1][tablenum][0];
      currSE->golomb_maxlevels = VLC_Golomb_Order[1][tablenum][1];

      writeSyntaxElement_GOLOMB(currSE, currBitStream);
      bitCount[BITS_COEFF_Y_MB] += currSE->len;
      no_bits += currSE->len;
      currSE++;  /* proceed to next SE */
      currMB->currSEnr++;
    }      /* inter */
  }  /* if ( cbp & (1<<b8) ) */

  return no_bits;
}

/*
 =======================================================================================================================
        Function:Write chroma coefficients of one 8x8 block Input: Output: Return: Attention:
 =======================================================================================================================
 */
int_32_t c_avs_enc::writeChromaCoeffAVS_B8(int_32_t b8)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  no_bits = 0;
  int_32_t  mb_nr = img->current_mb_nr;
  Macroblock  *currMB = &img->mb_data[mb_nr];
  const int_32_t  cbp = currMB->cbp;
  int_32_t  *bitCount = currMB->bitcounter;
  SyntaxElement  *currSE = &img->MB_SyntaxElements[currMB->currSEnr];
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  const char(*AVS_2DVLC_table_chroma)[26][27];

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t    inumblk;  /* number of blocks per CBP */
  int_32_t    inumcoeff;  /* number of coeffs per block */
  int_32_t    icoef;    /* current coefficient */
  int_32_t    ipos;    /* current position in cof_AVS */
  int_32_t    run, level, abs_level;
  int_16_t    *ACLevel;
  int_16_t    *ACRun;
  int_32_t    tablenum;
  int_32_t    symbol2D;
  int_32_t    escape_level_diff;
  TLS static const int_32_t  incVlc_chroma[5] = { 0, 1, 2, 4, 3000 };
  TLS static const int_32_t  blk_table[4] = { 0, -1, -1, -1 };
  TLS static const int_32_t  blkmode2ctx[4] = { LUMA_8x8, LUMA_8x4, LUMA_4x8, LUMA_4x4 };
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  inumblk = 1;
  inumcoeff = 65;    /* all positions + EOB */

  AVS_2DVLC_table_chroma = AVS_2DVLC_CHROMA;
  if(cbp & (1 << b8))
  {
    /* code all symbols */
    level = 1;  /* get inside loop */
    ipos = 0;
    ACLevel = img->cofAC[b8][0][0];
    ACRun = img->cofAC[b8][0][1];

    for(icoef = 0; icoef < inumcoeff; icoef++)  /* count coeffs */
    {
      if(!ACLevel[icoef]) break;
    }

    tablenum = 0;
    icoef--;
    for(icoef; icoef >= 0; icoef--)
    {
      level = ACLevel[icoef];
      abs_level = abs(level);
      run = ACRun[icoef];
      symbol2D = CODE2D_ESCAPE_SYMBOL;  /* symbol for out-of-table */
      if(level > -27 && level < 27 && run < 26)
      {
        if(tablenum == 0)
          symbol2D = AVS_2DVLC_table_chroma[tablenum][run][abs_level - 1];
        else
          symbol2D = AVS_2DVLC_table_chroma[tablenum][run][abs_level];
        if(symbol2D >= 0 && level < 0) symbol2D++;
        if(symbol2D < 0)
          symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
      }
      else
      {
        symbol2D = (CODE2D_ESCAPE_SYMBOL + (run << 1) + ((level > 0) ? 1 : 0));
      }

      currSE->type = SE_LUM_AC_INTER;
      currSE->value1 = symbol2D;
      currSE->value2 = 0;
      currSE->golomb_grad = VLC_Golomb_Order[2][tablenum][0];
      currSE->golomb_maxlevels = VLC_Golomb_Order[2][tablenum][1];

      writeSyntaxElement_GOLOMB(currSE, currBitStream);
      bitCount[BITS_COEFF_UV_MB] += currSE->len;
      no_bits += currSE->len;
      currSE++;
      currMB->currSEnr++;
      if(symbol2D >= CODE2D_ESCAPE_SYMBOL)
      {
        currSE->type = SE_LUM_AC_INTER;
        currSE->golomb_grad = 0;
        currSE->golomb_maxlevels = 11;
        escape_level_diff = abs_level - ((run > MaxRun[2][tablenum]) ? 1 : RefAbsLevel[tablenum + 14][run]);
        currSE->value1 = escape_level_diff;

        writeSyntaxElement_GOLOMB(currSE, currBitStream);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        no_bits += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }

      if(abs_level > incVlc_chroma[tablenum])
      {
        if(abs_level <= 2)
          tablenum = abs_level;
        else if(abs_level <= 4)
          tablenum = 3;
        else
          tablenum = 4;
      }
    }

    /* write eob */
    abs_level = level = 0;
    run = 0;
    symbol2D = AVS_2DVLC_table_chroma[tablenum][run][abs_level];
    currSE->type = SE_LUM_AC_INTER;
    currSE->value1 = symbol2D;
    currSE->value2 = 0;
    currSE->golomb_grad = VLC_Golomb_Order[2][tablenum][0];
    currSE->golomb_maxlevels = VLC_Golomb_Order[2][tablenum][1];

    writeSyntaxElement_GOLOMB(currSE, currBitStream);
    bitCount[BITS_COEFF_UV_MB] += currSE->len;
    no_bits += currSE->len;
    currSE++;  /* proceed to next SE */
    currMB->currSEnr++;
  }      /* if ( cbp & (1<<b8) ) */

  return no_bits;
}

/*
 =======================================================================================================================
        Function: Make Intra prediction for all 9 modes for 8*8, img_x and img_y are pixels offsets in the picture.
        Input: Output: Return: Attention:
 =======================================================================================================================
 */
//xzhao 20081106
void c_avs_enc::intrapred_luma_AVS(int_32_t img_x, int_32_t img_y)
{
  unsigned char edgepixels[40];
  int x,y,last_pix,new_pix;
  int b8_x,b8_y;
  int i;
  int block_available_up,block_available_up_right;
  int block_available_left,block_available_left_down;
  int bs_x=8;
  int bs_y=8;
  int MBRowSize = img->width / MB_BLOCK_SIZE;/*lgp*/
  int off_up=1,incr_y=1,off_y=0;/*lgp*/
  const int mb_nr    = img->current_mb_nr;    /*oliver*/
  Macroblock *currMB = &img->mb_data[mb_nr];  /*oliver*/
  int mb_left_available = (img_x >= MB_BLOCK_SIZE)?currMB->slice_nr == img->mb_data[mb_nr-1].slice_nr:0;  /*oliver*/
  int mb_up_available = (img_y >= MB_BLOCK_SIZE)?currMB->slice_nr == img->mb_data[mb_nr-MBRowSize].slice_nr:0;  /*oliver*/
  int mb_up_right_available;
  int mb_left_down_available;


  if((img->mb_y==0)||(img->mb_x==img->width/MB_BLOCK_SIZE-1))
  mb_up_right_available =1;
  else if((img_y-img->pix_y)>0)
  mb_up_right_available =(img_x-img->pix_x)>0? (currMB->slice_nr == img->mb_data[mb_nr+1].slice_nr):1;  /*oliver*/
   else
  mb_up_right_available =((img_x-img->pix_x)>0? (currMB->slice_nr == img->mb_data[mb_nr-MBRowSize+1].slice_nr):(currMB->slice_nr== img->mb_data[mb_nr-MBRowSize].slice_nr));  /*oliver*/


  if((img->mb_x==0)||(img->mb_y==img->height/MB_BLOCK_SIZE-1))
  mb_left_down_available = 1;
  else if(img_x-img->pix_x>0)
  mb_left_down_available =(img_y-img->pix_y)>0? (currMB->slice_nr == img->mb_data[mb_nr+MBRowSize].slice_nr):1;  /*oliver*/
  else
  mb_left_down_available =((img_y-img->pix_y)>0? (currMB->slice_nr == img->mb_data[mb_nr+MBRowSize-1].slice_nr):(currMB->slice_nr == img->mb_data[mb_nr-1].slice_nr));  /*oliver*/
  b8_x=img_x>>3;
  b8_y=img_y>>3;

  //check block up

  block_available_up=( b8_y-1>=0 && (mb_up_available||(img_y/8)%2));

  //check block up right
  block_available_up_right=( b8_x+1<(img->width>>3) && b8_y-1>=0 && mb_up_right_available);

  //check block left
  block_available_left=( b8_x - 1 >=0 && (mb_left_available || (img_x/8)%2));

  //check block left down
  block_available_left_down=( b8_x - 1>=0 && b8_y + 1 < (img->height>>3) && mb_left_down_available);
  if( ((img_x/8)%2)  && ((img_y/8)%2))
    block_available_up_right=0;

  if( ((img_x/8)%2)  || ((img_y/8)%2))
    block_available_left_down=0;

  //get prediciton pixels
  if(block_available_up)
  {
    for(x=0;x<bs_x;x++)
      EP[x+1]=imgY[img_y-/*1*/off_up/*lgp*/][img_x+x];

    if(block_available_up_right)
    {
      for(x=0;x<bs_x;x++)
        EP[1+x+bs_x]=imgY[img_y-/*1*/off_up/*lgp*/][img_x+bs_x+x];
      for(;x<bs_y;x++)
        EP[1+x+bs_x]=imgY[img_y-/*1*/off_up/*lgp*/][img_x+bs_x+bs_x-1];
    }
    else
    {
      for(x=0;x<bs_y;x++)
        EP[1+x+bs_x]=EP[bs_x];
    }

    for(;x<bs_y+2;x++)
      EP[1+x+bs_x]=EP[bs_x+x];

    EP[0]=imgY[img_y-/*1*/off_up/*lgp*/][img_x];
  }
  if(block_available_left)
  {
    for(y=0;y<bs_y;y++)
      EP[-1-y]=imgY[img_y+/*y*/incr_y*y-off_y][img_x-1];

    if(block_available_left_down)
    {
      for(y=0;y<bs_y;y++)
        EP[-1-y-bs_y]=imgY[img_y+bs_y+y][img_x-1];
      for(;y<bs_x;y++)
        EP[-1-y-bs_y]=imgY[img_y+bs_y+bs_y-1][img_x-1];
    }
    else
    {
      for(y=0;y<bs_x;y++)
        EP[-1-y-bs_y]=EP[-bs_y];
    }

    for(;y<bs_x+2;y++)
      EP[-1-y-bs_y]=EP[-y-bs_y];

    EP[0]=imgY[img_y-off_y/*lgp*/][img_x-1];
  }
  if(block_available_up&&block_available_left)
    EP[0]=imgY[img_y-/*1*/off_y-off_up/*lgp*/][img_x-1];

  //lowpass (Those emlements that are not needed will not disturb)
  last_pix=EP[-(bs_x+bs_y)];
  for(i=-(bs_x+bs_y);i<=(bs_x+bs_y);i++)
  {
    new_pix=( last_pix + (EP[i]<<1) + EP[i+1] + 2 )>>2;
    last_pix=EP[i];
    EP[i]=(unsigned char)new_pix;
  }
  //xzhao 20081106
  for(i=0;i<9;i++)
    img->available_intra_mode[i] = -1;       //初始化为-1，为了在Mode_Decision_for_AVS_IntraBlocks中判断当前的mode是否可用

  // 0 DC
  if(!block_available_up && !block_available_left)
  {
    img->available_intra_mode[DC_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DC_PRED][y][x]=128UL;
  }

  if(block_available_up && !block_available_left)
  {
    img->available_intra_mode[DC_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DC_PRED][y][x]=EP[1+x];
  }

  if(!block_available_up && block_available_left)
  {
    img->available_intra_mode[DC_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DC_PRED][y][x]=EP[-1-y];
  }


  if(block_available_up && block_available_left)
  {
    img->available_intra_mode[DC_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DC_PRED][y][x]=(EP[1+x]+EP[-1-y])>>1;
  }

  // 1 vertical
  if(block_available_up)
  {
    img->available_intra_mode[VERT_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[VERT_PRED][y][x]=imgY[img_y-1][img_x+x];
  }

  // 2 horizontal
  if(block_available_left)
  {
    img->available_intra_mode[HOR_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[HOR_PRED][y][x]=imgY[img_y+y][img_x-1];
  }

  // 3 down-right
  if(block_available_left&&block_available_up)
  {
    img->available_intra_mode[DOWN_RIGHT_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DOWN_RIGHT_PRED][y][x]=EP[x-y];
  }

  // 4 up-right bidirectional
  if(block_available_left&&block_available_up)
  {
    img->available_intra_mode[DOWN_LEFT_PRED] = 1;
    for(y=0UL;y<bs_y;y++)
      for(x=0UL;x<bs_x;x++)
        img->mprr[DOWN_LEFT_PRED][y][x]=(EP[2+x+y]+EP[-2-(x+y)])>>1;
  }


}

//void c_avs_enc::intrapred_luma_AVS(int_32_t img_x, int_32_t img_y)
//{
//  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//  byte        edgepixels[40];
//  int_32_t  x, y, last_pix, new_pix;
//  int_32_t  b8_x, b8_y;
//  int_32_t  i;
//  int_32_t  block_available_up, block_available_up_right;
//  int_32_t  block_available_left, block_available_left_down;
//  int_32_t  bs_x = 8;
//  int_32_t  bs_y = 8;
//  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//
//  b8_x = img_x >> 3;
//  b8_y = img_y >> 3;
//
//  /* check block up */
//  block_available_up = (b8_y - 1 >= 0 && img->ipredmode[1 + b8_x][b8_y] >= 0);
//
//  /* check block up right */
//  block_available_up_right =
//    (
//      b8_x +
//      1 < (img->width >> 3)
//    &&  b8_y - 1 >= 0
//    &&  img->ipredmode[b8_x + 2][b8_y] >= 0
//    );
//
//  /* check block left */
//  block_available_left = (b8_x - 1 >= 0 && img->ipredmode[b8_x][1 + b8_y] >= 0);
//
//  /* check block left down */
//  block_available_left_down =
//    (
//      b8_x -
//      1 >= 0
//    &&  b8_y + 1 < (img->height >> 3)
//    &&  img->ipredmode[b8_x][b8_y + 2] >= 0
//    );
//
//  /*
//   * xzhao { 2007.7.24 ;
//   * if(((img_x/8)%2) && ((img_y/8)%2))
//   */
//  if((img_x & 8) && (img_y & 8)) block_available_up_right = 0;
//
//  /* if( ((img_x/8)%2) || ((img_y/8)%2)) */
//  if((img_x & 8) || (img_y & 8)) block_available_left_down = 0;
//
//  /*
//   * xzhao } ;
//   * get prediction pixels
//   */
//  if(block_available_up)
//  {
//    /* xzhao { 2007.7.29 */
//    EP[1] = imgY[img_y - 1][img_x];
//    EP[2] = imgY[img_y - 1][img_x + 1];
//    EP[3] = imgY[img_y - 1][img_x + 2];
//    EP[4] = imgY[img_y - 1][img_x + 3];
//    EP[5] = imgY[img_y - 1][img_x + 4];
//    EP[6] = imgY[img_y - 1][img_x + 5];
//    EP[7] = imgY[img_y - 1][img_x + 6];
//    EP[8] = imgY[img_y - 1][img_x + 7];
//
//    if(block_available_up_right)
//    {
//      /* xzhao { 2007.7.29 */
//      EP[9] = imgY[img_y - 1][img_x + 8];
//      EP[10] = imgY[img_y - 1][img_x + 9];
//      EP[11] = imgY[img_y - 1][img_x + 10];
//      EP[12] = imgY[img_y - 1][img_x + 11];
//      EP[13] = imgY[img_y - 1][img_x + 12];
//      EP[14] = imgY[img_y - 1][img_x + 13];
//      EP[15] = imgY[img_y - 1][img_x + 14];
//      EP[16] = imgY[img_y - 1][img_x + 15];
//    }
//    else
//    {
//      /* xzhao { 2007.7.29 */
//      EP[9] = EP[8];
//      EP[10] = EP[8];
//      EP[11] = EP[8];
//      EP[12] = EP[8];
//      EP[13] = EP[8];
//      EP[14] = EP[8];
//      EP[15] = EP[8];
//      EP[16] = EP[8];
//    }
//
//    EP[17] = EP[16];
//    EP[18] = EP[17];
//
//    EP[0] = imgY[img_y - 1][img_x];
//  }
//
//  if(block_available_left)
//  {
//    /* xzhao { 2007.7.29 */
//    EP[-1] = imgY[img_y][img_x - 1];
//    EP[-2] = imgY[img_y + 1][img_x - 1];
//    EP[-3] = imgY[img_y + 2][img_x - 1];
//    EP[-4] = imgY[img_y + 3][img_x - 1];
//    EP[-5] = imgY[img_y + 4][img_x - 1];
//    EP[-6] = imgY[img_y + 5][img_x - 1];
//    EP[-7] = imgY[img_y + 6][img_x - 1];
//    EP[-8] = imgY[img_y + 7][img_x - 1];
//
//    if(block_available_left_down)
//    {
//      /* xzhao { 2007.7.29 */
//      EP[-9] = imgY[img_y + 8][img_x - 1];
//      EP[-10] = imgY[img_y + 9][img_x - 1];
//      EP[-11] = imgY[img_y + 10][img_x - 1];
//      EP[-12] = imgY[img_y + 11][img_x - 1];
//      EP[-13] = imgY[img_y + 12][img_x - 1];
//      EP[-14] = imgY[img_y + 13][img_x - 1];
//      EP[-15] = imgY[img_y + 14][img_x - 1];
//      EP[-16] = imgY[img_y + 15][img_x - 1];
//    }
//    else
//    {
//      EP[-9] = EP[-8];
//      EP[-10] = EP[-8];
//      EP[-11] = EP[-8];
//      EP[-12] = EP[-8];
//      EP[-13] = EP[-8];
//      EP[-14] = EP[-8];
//      EP[-15] = EP[-8];
//      EP[-16] = EP[-8];
//    }
//
//    EP[-17] = EP[-16];
//    EP[-18] = EP[-17];
//
//    EP[0] = imgY[img_y][img_x - 1];
//  }
//
//  if(block_available_up && block_available_left) EP[0] = imgY[img_y - 1][img_x - 1];
//
//  /*
//   * xzhao { 2007.7.29 ;
//   * low pass (Those elements that are not needed will not disturb)
//   */
//  last_pix = EP[-(bs_x + bs_y)];
//
//  for(i = -(bs_x + bs_y); i <= (bs_x + bs_y); i++)
//  {
//    new_pix = (last_pix + (EP[i] << 1) + EP[i + 1] + 2) >> 2;
//    last_pix = EP[i];
//    EP[i] = (unsigned char) new_pix;
//  }
//
//  for(i = 0; i < 5; i++)
//    img->available_intra_mode[i] = -1;       //初始化为-1，为了在Mode_Decision_for_AVS_IntraBlocks中判断当前的mode是否可用
//
//  /* 2 DC */
//  if(!block_available_up && !block_available_left)
//  {
//    img->available_intra_mode[DC_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[DC_PRED][y][x] = 128UL;
//  }
//
//  if(block_available_up && !block_available_left)
//  {
//    img->available_intra_mode[DC_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[DC_PRED][y][x] = EP[1 + x];
//  }
//
//  if(!block_available_up && block_available_left)
//  {
//    img->available_intra_mode[DC_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[DC_PRED][y][x] = EP[-1 - y];
//  }
//
//  if(block_available_up && block_available_left)
//  {
//    img->available_intra_mode[DC_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[DC_PRED][y][x] = (EP[1 + x] + EP[-1 - y]) >> 1;
//  }
//
//  /* 0 vertical */
//  if(block_available_up)
//  {
//    img->available_intra_mode[VERT_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[VERT_PRED][y][x] = imgY[img_y - 1][img_x + x];
//  }
//
//  /* 1 horizontal */
//  if(block_available_left)
//  {
//    img->available_intra_mode[HOR_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[HOR_PRED][y][x] = imgY[img_y + y][img_x - 1];
//  }
//
//  /* 4 down-right */
//  if(block_available_left && block_available_up)
//  {
//    img->available_intra_mode[DOWN_RIGHT_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++) img->mprr[DOWN_RIGHT_PRED][y][x] = EP[x - y];
//  }
//
//  /* 3 up-right bidirectional */
//  if(block_available_left && block_available_up)
//  {
//    img->available_intra_mode[DOWN_LEFT_PRED] = 1;
//    for(y = 0UL; y < bs_y; y++)
//      for(x = 0UL; x < bs_x; x++)
//        img->mprr[DOWN_LEFT_PRED][y][x] = (EP[2 + x + y] + EP[-2 - (x + y)]) >> 1;
//  }
//  /* xzhao } */
//}