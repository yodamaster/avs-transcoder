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
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "global.h"
#include "const_data.h"
#include "sse_header.h"
/*
*************************************************************************
* Function:Update the coordinates for the next macroblock to be processed
* Input:mb: MB address in scan order
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::set_MB_parameters (int_32_t mb)
{
  img->current_mb_nr = mb;
  img->mb_x = mb % img->img_width_in_mb;
  img->mb_y = mb / img->img_width_in_mb;

  // Define vertical positions
  img->block8_y= img->mb_y<<1;  //img->mb_y * BLOCK_SIZE/2;
  img->block_y = img->mb_y<<2;  //img->mb_y * BLOCK_SIZE;      // vertical luma block position
  img->pix_y   = img->mb_y<<4;  //img->mb_y * MB_BLOCK_SIZE;   // vertical luma macroblock position
  img->pix_c_y = img->mb_y<<3;  //img->mb_y * MB_BLOCK_SIZE/2; // vertical chroma macroblock position

  // Define horizontal positions
  img->block8_x  = img->mb_x<<1;  //img->mb_x * BLOCK_SIZE/2;
  img->block_x   = img->mb_x<<2;  //img->mb_x * BLOCK_SIZE;      // luma block
  img->pix_x     = img->mb_x<<4;  //img->mb_x * MB_BLOCK_SIZE;   // luma pixel
  img->block_c_x = img->mb_x<<1;  //img->mb_x * BLOCK_SIZE/2;    // chroma block
  img->pix_c_x   = img->mb_x<<3;  //img->mb_x * MB_BLOCK_SIZE/2; // chroma pixel
}

int_32_t clip1a(int_32_t a)
{
  return ((a)>255?255:((a)<0?0:(a)));
}

/*
*************************************************************************
* Function:Update the coordinates and statistics parameter for the
next macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::proceed2nextMacroblock()
{
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int_32_t*        bitCount = currMB->bitcounter;
#if TRACE
  int_32_t i;
  if (p_trace)
  {
    fprintf(p_trace, "\n*********** Pic: %i (I/P) MB: %i Slice: %i **********\n\n", frame_no, img->current_mb_nr, img->current_slice_nr);
    // Write out the tracestring for each symbol
    for (i=0; i<currMB->currSEnr; i++)
    {
      trace2out(&(img->MB_SyntaxElements[i]));
    }
  }
#endif

  // Update the statistics
  stat->bit_use_mb_type[img->type]      += bitCount[BITS_MB_MODE];
  stat->bit_use_coeffY[img->type]       += bitCount[BITS_COEFF_Y_MB] ;
  stat->tmp_bit_use_cbp[img->type]      += bitCount[BITS_CBP_MB];
  stat->bit_use_coeffC [img->type]      += bitCount[BITS_COEFF_UV_MB];
  stat->bit_use_delta_quant[img->type]  += bitCount[BITS_DELTA_QUANT_MB];
  switch(img->type)
  {
  case INTRA_IMG:
    ++stat->mode_use_intra[currMB->mb_type];
    break;
  case INTER_IMG:
    ++stat->mode_use_inter[0][currMB->mb_type];
    stat->bit_use_mode_inter[0][currMB->mb_type]+= bitCount[BITS_INTER_MB];
    ++stat->quant0;
    stat->quant1 += img->qp;      // to find average quant for inter frames
    break;
  case B_IMG:
    stat->bit_use_mode_inter[1][currMB->mb_type]+= bitCount[BITS_INTER_MB];
    ++stat->mode_use_inter[1][currMB->mb_type];
    break;
  default:
    printf("unknow img type:%d\n", frame_no);
    break;
  }
  //set the intra prediction mode to -1 if the mb_type is inter
  if (currMB->mb_type != I4MB)
  {
    //img->ipredmode[img->block_x+1][img->block_y+1] is currMB
    img->ipredmode[img->block8_x+1][img->block8_y+1] = -1;
    img->ipredmode[img->block8_x+2][img->block8_y+1] = -1;
    img->ipredmode[img->block8_x+1][img->block8_y+2] = -1;
    img->ipredmode[img->block8_x+2][img->block8_y+2] = -1;
  }
}

void c_avs_enc::start_macroblock()
{
  int_32_t j,k;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  //Rate control
  int_32_t predict_error,dq;
  int_32_t DELTA_QP,DELTA_QP2;
  int_32_t QP,QP2;
  const int_32_t intra_dc_pred_modes[4]  = {DC_PRED, DC_PRED, DC_PRED, DC_PRED};
  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  currMB->slice_nr = img->current_slice_nr;

  // Rate control
  if(input->RCEnable == 1)
  {
    /*frame layer rate control*/
    if(input->basicunit==img->total_number_mb)
    {
      currMB->delta_qp = 0;
      currMB->qp       = img->qp;
    }
    else /*basic unit layer rate control*/
    {
      /*each I or B frame has only one QP*/
      if((img->type==INTRA_IMG)||(img->type==B_IMG))
      {
        currMB->delta_qp = 0;
        currMB->qp       = img->qp;
      }
      else if(img->type==INTER_IMG)
      {
        if (img->current_mb_nr == 0) //first macroblock
        {
          // Initialize delta qp change from last macroblock. Feature may be used for future rate control
          currMB->delta_qp = 0;
          currMB->qp       = img->qp;
          DELTA_QP = DELTA_QP2 = currMB->delta_qp;
          QP = QP2 = currMB->qp;
        }
        else
        {
          if (img->mb_data[img->current_mb_nr-1].prev_cbp == 1)
          {
            currMB->delta_qp = 0;
            currMB->qp       = img->qp;
          }
          else
          {
            currMB->qp = img->mb_data[img->current_mb_nr-1].prev_qp;
            currMB->delta_qp = currMB->qp - img->mb_data[img->current_mb_nr-1].qp;
            img->qp = currMB->qp;
          }
          DELTA_QP = DELTA_QP2 = currMB->delta_qp;
          QP = QP2 = currMB->qp;
        }

        /*compute the quantization parameter for each basic unit of P frame*/
        if((img->NumberofCodedMacroBlocks>0) && (img->NumberofCodedMacroBlocks%img->BasicUnit==0))
        {
          /*frame coding*/
          if(input->InterlaceCodingOption==0)
          {
            updateRCModel();
            img->BasicUnitQP = updateQuantizationParameter(img->TopFieldFlag);
          }
          /*adaptive field/frame coding*/
          else if((input->InterlaceCodingOption==2) && (img->IFLAG==0))
          {
            updateRCModel();
            img->BasicUnitQP = updateQuantizationParameter(img->TopFieldFlag);
          }
          /*field coding*/
          else if((input->InterlaceCodingOption==1) && (img->IFLAG==0))
          {
            updateRCModel();
            img->BasicUnitQP = updateQuantizationParameter(img->TopFieldFlag);
          }
        }

        if(img->current_mb_nr==0)
          img->BasicUnitQP=img->qp;

        currMB->predict_qp=img->BasicUnitQP;

        if(currMB->predict_qp>currMB->qp+25)
          currMB->predict_qp=currMB->qp+25;
        else if(currMB->predict_qp<currMB->qp-26)
          currMB->predict_qp=currMB->qp-26;

        currMB->prev_qp = currMB->predict_qp;

        dq = currMB->delta_qp + currMB->predict_qp-currMB->qp;
        if(dq < -26)
        {
          dq = -26;
          predict_error = dq-currMB->delta_qp;
          img->qp = img->qp+predict_error;
          currMB->delta_qp = -26;
        }
        else if(dq > 25)
        {
          dq = 25;
          predict_error = dq - currMB->delta_qp;
          img->qp = img->qp + predict_error;
          currMB->delta_qp = 25;
        }
        else
        {
          currMB->delta_qp = dq;
          predict_error=currMB->predict_qp-currMB->qp;
          img->qp = currMB->predict_qp;
        }
        currMB->qp =  img->qp;
        currMB->predict_error=predict_error;
      }
    }
  }
  else
  {
    currMB->delta_qp = 0;
    currMB->qp       = img->qp;       // needed in loop filter (even if constant QP is used)
  }
  // If MB is next to a slice boundary, mark neighboring blocks unavailable for prediction
  CheckAvailabilityOfNeighbors();

  currMB->mb_type   = 0;
  memset(currMB->mvd, 0, 16*sizeof(int_32_t));
  if(img->type != B_IMG)
  {
    for (k=0; k < 2; k++)
      for (j=0; j < BLOCK_MULTIPLE; j++)
        memset(&tmp_mv[k][img->block_y+j][img->block_x+4], 0, BLOCK_MULTIPLE*sizeof(int_32_t));
  }
  else
  {
    memcpy(currMB->intra_pred_modes, intra_dc_pred_modes, sizeof(intra_dc_pred_modes));
  }
  // Reset syntax element entries in MB struct
  currMB->mb_type   = 0;
  currMB->cbp_blk   = 0;
  currMB->cbp       = 0;
  currMB->mb_field  = 0;
  currMB->cbp_bits  = 0;
  currMB->c_ipred_mode = DC_PRED_8;
  currMB->lf_disable         = input->loop_filter_disable;
  currMB->lf_alpha_c0_offset = input->alpha_c_offset;
  currMB->lf_beta_offset     = input->beta_offset;
  // Initialize counter for MB symbols
  currMB->currSEnr=0;
  // Initialize bitcounters for this macroblock
  if(img->current_mb_nr == 0) // No slice header to account for
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }
  else if (currMB->slice_nr == img->mb_data[img->current_mb_nr-1].slice_nr)
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }

  currMB->bitcounter[BITS_MB_MODE] = 0;
  currMB->bitcounter[BITS_COEFF_Y_MB] = 0;
  currMB->bitcounter[BITS_INTER_MB] = 0;
  currMB->bitcounter[BITS_CBP_MB] = 0;
  currMB->bitcounter[BITS_DELTA_QUANT_MB] = 0;
  currMB->bitcounter[BITS_COEFF_UV_MB] = 0;
  currMB->bitcounter[BITS_MVD]         = 0;
}

void c_avs_enc::terminate_macroblock(myboolean *end_of_picture)
{
  TLS static int_32_t skip = myfalse;
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t rlc_bits=0;
  int_32_t mb_width = img->width/16;
  int_32_t slice_mb = input->slice_row_nr*mb_width;

  img->coded_mb_nr++;
  if(input->slice_row_nr && (img->coded_mb_nr != img->total_number_mb))
  {
    if(img->coded_mb_nr%slice_mb == 0 )
    {
      img->mb_data[img->current_mb_nr+1].slice_nr = img->mb_data[img->current_mb_nr].slice_nr+1;
      img->mb_no_currSliceLastMB =  min(img->mb_no_currSliceLastMB + slice_mb, img->total_number_mb) ;
    }
    else
    {
      img->mb_data[img->current_mb_nr+1].slice_nr = img->mb_data[img->current_mb_nr].slice_nr;
    }
  }

  if (img->coded_mb_nr == img->total_number_mb) // maximum number of MBs reached
  {
    *end_of_picture = mytrue;
    img->coded_mb_nr= 0;
  }

  if(*end_of_picture == mytrue && img->cod_counter)
  {
    currSE->value1 = img->cod_counter;
    currSE->mapping = &c_avs_enc::ue_linfo;
    currSE->type = SE_MBTYPE;
    writeSyntaxElement_UVLC(currSE, currBitStream);
    rlc_bits=currSE->len;
    currMB->bitcounter[BITS_MB_MODE] += rlc_bits;
    img->cod_counter = 0;
  }
}


void c_avs_enc::CheckAvailabilityOfNeighbors()
{
  int_32_t i,j;
  int_32_t mb_nr    = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int_32_t remove_prediction;

  img->block_x /= 2;
  img->block_y /= 2;
  // mark all neighbors as unavailable
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      img->mb_data[mb_nr].mb_available[i][j] = NULL;
    }
  }
  img->mb_data[mb_nr].mb_available[1][1] = currMB; // current MB

  // Check MB to the left
  if(img->pix_x >= MB_BLOCK_SIZE)
  {
    remove_prediction = (currMB->slice_nr != img->mb_data[mb_nr-1].slice_nr);
    // upper blocks
    if (remove_prediction)
    {
      img->ipredmode[img->block_x][img->block_y+1] = -1;
      img->ipredmode[img->block_x][img->block_y+2] = -1;
      // lower blocks 下面这些可否去掉?因为是8x8的块
      img->ipredmode[img->block_x][img->block_y+3] = -1;
      img->ipredmode[img->block_x][img->block_y+4] = -1;
    }
    else
    {
      currMB->mb_available[1][0] = &(img->mb_data[mb_nr-1]);
    }
  }

  // Check MB above
  if(img->pix_y >= MB_BLOCK_SIZE)
  {
    remove_prediction = (currMB->slice_nr != img->mb_data[mb_nr-img->img_width_in_mb].slice_nr);
    // left blocks
    if (remove_prediction)
    {
      img->ipredmode[img->block_x+1][img->block_y] = -1;
      img->ipredmode[img->block_x+2][img->block_y] = -1;
    }
    else
    {
      currMB->mb_available[0][1]=&(img->mb_data[mb_nr-img->img_width_in_mb]);
    }
  }

  // Check MB left above
  if(img->pix_x >= MB_BLOCK_SIZE && img->pix_y >= MB_BLOCK_SIZE )
  {
    remove_prediction = (currMB->slice_nr != img->mb_data[mb_nr-img->img_width_in_mb-1].slice_nr);
    if (remove_prediction)
    {
      img->ipredmode[img->block_x][img->block_y] = -1;
    }
    else
    {
      currMB->mb_available[0][0] = &(img->mb_data[mb_nr-img->img_width_in_mb-1]);
    }
  }

  // Check MB right above
  if(img->pix_y >= MB_BLOCK_SIZE && img->pix_x < (img->width-MB_BLOCK_SIZE ))
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-img->img_width_in_mb+1].slice_nr)
    {
      currMB->mb_available[0][2]=&(img->mb_data[mb_nr-img->img_width_in_mb+1]);
    }
  }
  img->block_x *= 2;
  img->block_y *= 2;
}

/*
*************************************************************************
* Function:Predict one component of a 4x4 Luma block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
void c_avs_enc::OneComponentLumaPrediction4x4 (int_32_t*   mpred,      //  --> array of prediction values (row by row)
                                               int_32_t    pic_pix_x,  // <--  absolute horizontal coordinate of 4x4 block
                                               int_32_t    pic_pix_y,  // <--  absolute vertical   coordinate of 4x4 block
                                               int_32_t*   mv,         // <--  motion vector
                                               int_32_t    ref)        // <--  reference frame (0.. / -1:backward)
{
  int_32_t incr=0;
  pel_t** ref_pic;

  int_32_t     pix_add = 4;
  int_32_t     j0      = (pic_pix_y << 2) + mv[1], j1=j0+pix_add, j2=j1+pix_add, j3=j2+pix_add;
  int_32_t     i0      = (pic_pix_x << 2) + mv[0], i1=i0+pix_add, i2=i1+pix_add, i3=i2+pix_add;
  int_32_t     i,j;

  pel_t (c_avs_enc::*get_pel) (pel_t**, int_32_t, int_32_t) = &c_avs_enc::UMVPelY_14;
  incr = 1;
  if(!img->picture_structure) //field coding
  {
    incr =2;
  }

  ref_pic  = img->type==B_IMG? (byte**)mref [ref+incr] : (byte**)mref [ref];

  for(j=0;j<4;j++)
  {
    for(i=0;i<4;i++)
    {
      *mpred++ = (this->*get_pel) (ref_pic, j0+pix_add*j, i0+pix_add*i);
    }
  }
}

/*
*************************************************************************
* Function:Predict one 4x4 Luma block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void
c_avs_enc::LumaPrediction4x4 (int_32_t  block_x,    // <--  relative horizontal block coordinate of 4x4 block
                              int_32_t  block_y,    // <--  relative vertical   block coordinate of 4x4 block
                              int_32_t  fw_mode,    // <--  forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                              int_32_t  bw_mode,    // <--  backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                              int_32_t  fw_ref,     // <--  reference frame for forward prediction (-1: Intra4x4 pred. with fw_mode)
                              int_32_t  bw_ref  )
{
  TLS static int_32_t fw_pred[16];
  TLS static int_32_t bw_pred[16];
  int_32_t  i, j;
  int_32_t  block_x4  = block_x+4;
  int_32_t  block_y4  = block_y+4;
  int_32_t  pic_pix_x = img->pix_x + block_x;
  int_32_t  pic_pix_y = img->pix_y + block_y;
  int_32_t  bx        = block_x >> 3;
  int_32_t  by        = block_y >> 3;
  int_32_t* fpred     = fw_pred;
  int_32_t* bpred     = bw_pred;
  int_32_t  direct    = (fw_mode == 0 && bw_mode == 0 && (img->type==B_IMG));
  int_32_t  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type!=B_IMG));
  int_32_t scale_frame;
  int_32_t  *****fmv_array = (fw_mode && bw_mode)?img->all_omv:img->all_mv;    // For MB level frame/field coding
  int_32_t  *****bmv_array = img->all_bmv;                                     // For MB level frame/field coding
  int_32_t fw_lum_scale , fw_lum_shift ;
  int_32_t bw_lum_scale , bw_lum_shift ;
  int_32_t bw_ref_num ;
  int_32_t fw_ref_num ;

  int_32_t structshift ;

  if(img->picture_structure)
  {
    structshift = 1 ;
  }
  else
  {
    structshift = 0 ;
  }

  if (direct)
  {
    fw_ref= 0;
    bw_ref= 0;
    if(!img->picture_structure) //field coding
    {
      scale_frame = refFrArr[img->block8_y + by][img->block8_x+bx];
      bw_ref = 0;
      if (img->current_mb_nr_fld < img->total_number_mb) //top field
      {
        fw_ref =  scale_frame >= 1 ? 1 : 0;
        bw_ref = scale_frame >= 0 ?1 :0;
      }
      else
      {
        fw_ref =(scale_frame == 1 || scale_frame < 0)? 0 : 1;
        bw_ref = 0;
      }
      if(scale_frame < 0)
        bw_ref = 1;
    }
  }

  direct_mode = direct;

  if (fw_mode || (direct && (fw_ref !=-1) ) || skipped)
  {
    if(!img->picture_structure)  // !! field
      OneComponentLumaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, fmv_array [bx][by][fw_ref][fw_mode], fw_ref);
    else  // !! frame
      OneComponentLumaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, fmv_array [bx][by][fw_ref][fw_mode], (direct ?0:fw_ref));
  }
  if (bw_mode || (direct && (bw_ref !=-1)) ||(direct && !img->picture_structure && (bw_ref ==-1)))
  {
    {
      int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,refframe,delta_PB;
      int_32_t mv[2];
      refframe = fw_ref;
      delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
      delta_P = (delta_P + 512) % 512;
      if(img->picture_structure)
        TRp = (refframe+1)*delta_P;  //the latest backward reference
      else
      {
        TRp = delta_P;//refframe == 0 ? delta_P-1 : delta_P+1;
      }
      delta_PB = 2*(img->tr - img->imgtr_last_P_frm);
      delta_PB = (delta_PB + 512)%512;
      TRp = (TRp+512) % 512;
      if(!img->picture_structure)
      {
        if(img->current_mb_nr_fld < img->total_number_mb) //top field
          DistanceIndexFw =  refframe == 0 ? delta_PB-1:delta_PB;
        else
          DistanceIndexFw =  refframe == 0 ? delta_PB:delta_PB+1;
      }
      else
        DistanceIndexFw = delta_PB;
      DistanceIndexBw    = TRp - DistanceIndexFw;
      mv[0] = - ((fmv_array[bx][by][fw_ref][fw_mode][0]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
      mv[1] = - ((fmv_array[bx][by][fw_ref][fw_mode][1]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
      if(fw_mode && bw_mode)
        OneComponentLumaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, mv,(!img->picture_structure) ? (-2+bw_ref) : (-1-bw_ref));
      else
        OneComponentLumaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, bmv_array[bx][by][bw_ref][bw_mode],(!img->picture_structure) ? (-2+bw_ref) : (-1-bw_ref));
    }
  }

  if (direct || (fw_mode && bw_mode))
  {
    bw_ref_num = bw_ref ;
    fw_ref_num = fw_ref + structshift;
    fw_lum_scale = img->lum_scale[fw_ref_num];
    fw_lum_shift = img->lum_shift[fw_ref_num];
    bw_lum_scale = img->lum_scale[bw_ref_num];
    bw_lum_shift = img->lum_shift[bw_ref_num];

    for(j=block_y; j<block_y4; j++)
    {
      for (i=block_x; i<block_x4; i++)
      {
        img->mpr[j][i] = (*fpred + *bpred + 1) / 2;
        if(img->LumVarFlag == 1)
          img->mpr_weight[j][i] = (((((*fpred)*fw_lum_scale + 16)>>5) + fw_lum_shift) + ((((*bpred)*bw_lum_scale + 16)>>5) + bw_lum_shift) + 1) / 2;
        fpred++;
        bpred++;
      }
    }
  }
  else if (fw_mode || skipped)
  {
    if(img->type==B_IMG)
    {
      fw_lum_scale = img->lum_scale[fw_ref + structshift];
      fw_lum_shift = img->lum_shift[fw_ref + structshift];
    }
    else
    {
      fw_lum_scale = img->lum_scale[fw_ref];
      fw_lum_shift = img->lum_shift[fw_ref];
    }

    for(j=block_y; j<block_y4; j++)
    {
      for (i=block_x; i<block_x4; i++)
      {
        img->mpr[j][i] = *fpred;
        if((!skipped) && (img->LumVarFlag == 1))
        {
          img->mpr_weight[j][i] = ((((*fpred)*fw_lum_scale + 16)>>5) + fw_lum_shift);
        }
        fpred++;
      }
    }
  }
  else
  {
    bw_ref_num = bw_ref ;
    bw_lum_scale = img->lum_scale[bw_ref_num];
    bw_lum_shift = img->lum_shift[bw_ref_num];
    for(j=block_y; j<block_y4; j++)
    {
      for (i=block_x; i<block_x4; i++)
      {
        img->mpr[j][i] = *bpred;
        if(img->LumVarFlag == 1)
        {
          img->mpr_weight[j][i] = ((((*bpred)*bw_lum_scale + 16)>>5) + bw_lum_shift);
        }
        bpred++;
      }
    }
  }
}


/*
*************************************************************************
* Function:Residual Coding of an 8x8 Luma block (not for intra)
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t                                       //  ==> coefficient cost
c_avs_enc::LumaResidualCoding8x8 (int_32_t  *cbp,         //  --> cbp (updated according to processed 8x8 luminance block)
								  int_32_t  *cbp_blk,     //  --> block cbp (updated according to processed 8x8 luminance block)
								  int_32_t  block8x8,     // <--  block number of 8x8 block
								  int_32_t  fw_mode,      // <--  forward  prediction mode (1-7, 0=DIRECT)
								  int_32_t  bw_mode,      // <--  backward prediction mode (1-7, 0=DIRECT)
								  int_32_t  fw_refframe,  // <--  reference frame for forward prediction
								  int_32_t  bw_refframe   // <--  reference frame for backward prediction
								  )
{
	
	int_32_t coeff_cost = 0;
	int_32_t cbp_mask   = 1 << block8x8;
	int_32_t scrFlag = 0;                // 0=noSCR, 1=strongSCR, 2=jmSCR
	int_32_t mv[2],fw_ref,it;
	int_32_t h_fw,v_fw,h_bw,v_bw,h,v,blkh = block8x8 % 2, blkv = block8x8 / 2;
	int_32_t pix_h = img->pix_x, pix_v = img->pix_y,mb_h = blkh  << 3, mb_v = blkv << 3;
	int_32_t width4  = ((img->width + (IMG_PAD_SIZE << 1) - 1) << 2) - 32;
	int_32_t height4 = ((img->height+ (IMG_PAD_SIZE << 1) - 1) << 2) - 32;
	int_32_t skipped = (fw_mode==0 && bw_mode==0 && img->type != B_IMG);
    int_32_t direct = (fw_mode==0 && bw_mode==0 && img->type == B_IMG);
	int_32_t *****fw_mv = ((fw_mode && bw_mode) ? img->all_omv : img->all_mv);
	int_32_t *****bw_mv = ((fw_mode && bw_mode) ? img->all_bw_omv : img->all_bmv);
	byte     **fw_ref_pic, ** bw_ref_pic, **imgY_original = imgY_org;
	int_32_t dir    = ((direct || (fw_mode && bw_mode)) ? 2 : ((fw_mode || skipped) ? 0 : 1));
	byte     *d1,*d2,*d3,*d4,*d5,*d6;
	int_16_t *d7,*d8;
	__declspec(align(16)) int_16_t curr_blk[B8_SIZE][B8_SIZE],*pos;
	__m128i xmm0, zero128={0};
	scrFlag = (img->type == B_IMG ? 1 : scrFlag);
	fw_ref = (img->type == B_IMG ? 1 : fw_refframe);
	if (direct)
	{
		fw_refframe = bw_refframe = 0;
	}

	switch (dir)
	{
	case 0:  /* forward */
		mv[0] = fw_mv[blkh][blkv][fw_refframe][fw_mode][0];
		mv[1] = fw_mv[blkh][blkv][fw_refframe][fw_mode][1];
		h_fw = ((pix_h + mb_h + IMG_PAD_SIZE) << 2) + mv[0];
		v_fw = ((pix_v + mb_v + IMG_PAD_SIZE) << 2) + mv[1];
		if (h_fw<0)
			h_fw &= 3;
		else if (h_fw>width4)
			h_fw = width4 + (h_fw & 3);
		if (v_fw<0)
			v_fw &= 3;
		else if (v_fw>height4)
			v_fw = height4 + (v_fw & 3);
		h = h_fw & 3, v = v_fw & 3;
		h_fw >>= 2, v_fw >>= 2;

		fw_ref_pic = mref[fw_ref][v][h];
		for (pos = (int_16_t *)curr_blk , it = 0; it<8; it+=2)
		{
			d1 = &imgY_original[pix_v + mb_v + it][pix_h + mb_h];
			d2 = &imgY_original[pix_v + mb_v + it + 1][pix_h + mb_h];
			d3 = &fw_ref_pic[v_fw + it][h_fw];
			d4 = &fw_ref_pic[v_fw + it + 1][h_fw];
			d7 = &img->mpr[mb_v+it][mb_h];
			d8 = &img->mpr[mb_v+it+1][mb_h];
			_asm
			{
				mov    esi, dword ptr [d1]
				movdqu xmm0, xmmword ptr[esi]
				mov    esi, dword ptr[d2]
				movdqu xmm1, xmmword ptr[esi]
				mov    esi, dword ptr[d3]
				movdqu xmm2, xmmword ptr[esi]
				mov    esi, dword ptr[d4]
				movdqu xmm3, xmmword ptr[esi]

				pxor xmm7, xmm7
				punpcklbw xmm0, xmm7
				punpcklbw xmm1, xmm7
				punpcklbw xmm2, xmm7
				punpcklbw xmm3, xmm7
				
				psubw xmm0, xmm2
				psubw xmm1, xmm3
				
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm0
				add pos, 16
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm1
                add pos, 16
				mov esi, dword ptr [d7]
				movdqa xmmword ptr[esi], xmm2
				mov esi, dword ptr [d8]
				movdqa xmmword ptr[esi], xmm3
			}
		}
		break;

	case 1:  /* backward */
		mv[0] = bw_mv[blkh][blkv][bw_refframe][bw_mode][0];
		mv[1] = bw_mv[blkh][blkv][bw_refframe][bw_mode][1];
		h_bw = ((pix_h + mb_h + IMG_PAD_SIZE) << 2) + mv[0];
		v_bw = ((pix_v + mb_v + IMG_PAD_SIZE) << 2) + mv[1];
		if (h_bw<0)
			h_bw &= 3;
		else if (h_bw>width4)
			h_bw = width4 + (h_bw & 3);
		if (v_bw<0)
			v_bw &= 3;
		else if (v_bw>height4)
			v_bw = height4 + (v_bw & 3);
		h = h_bw & 3, v = v_bw & 3;
		h_bw >>= 2, v_bw >>= 2;
		bw_ref_pic = mref[bw_refframe][v][h];
		for (pos = (int_16_t *)curr_blk, it = 0; it<8; it+=2)
		{
			d1 = &imgY_original[pix_v + mb_v + it][pix_h + mb_h];
			d2 = &imgY_original[pix_v + mb_v + it + 1][pix_h + mb_h];
			d3 = &bw_ref_pic[v_bw + it][h_bw];
			d4 = &bw_ref_pic[v_bw + it + 1][h_bw];
			d7 = &img->mpr[mb_v + it][mb_h];
			d8 = &img->mpr[mb_v + it + 1][mb_h];
			_asm
			{
				mov esi, dword ptr [d1]
				movdqu xmm0, xmmword ptr [esi]
				mov esi, dword ptr [d2]
				movdqu xmm1, xmmword ptr [esi]
				mov esi, dword ptr [d3]
				movdqu xmm2, xmmword ptr [esi]
				mov esi, dword ptr [d4]
				movdqu xmm3, xmmword ptr [esi]

				pxor xmm7, xmm7
				punpcklbw xmm0, xmm7
				punpcklbw xmm1, xmm7
				punpcklbw xmm2, xmm7
				punpcklbw xmm3, xmm7
				
				psubw xmm0, xmm2
				psubw xmm1, xmm3
				
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm0
				add pos, 16
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm1
				add pos, 16
				mov esi, dword ptr [d7]
				movdqa xmmword ptr [esi], xmm2
				mov esi, dword ptr [d8]
				movdqa xmmword ptr [esi], xmm3
 			}
		}
		break;

	case 2:  /* bi-direction */
        /* forward  */
		mv[0] = fw_mv[blkh][blkv][0][fw_mode][0];
		mv[1] = fw_mv[blkh][blkv][0][fw_mode][1];
		h_fw = ((pix_h + mb_h + IMG_PAD_SIZE) << 2) + mv[0];
		v_fw = ((pix_v + mb_v + IMG_PAD_SIZE) << 2) + mv[1];
		if (h_fw<0) 
			h_fw &= 3;
		else if (h_fw>width4) 
			h_fw = width4 + (h_fw & 3);
		if (v_fw<0) 
			v_fw &= 3;
		else if (v_fw>height4) 
			v_fw = height4 + (v_fw & 3);
		h = h_fw & 3, v = v_fw & 3;
		h_fw >>= 2, v_fw >>= 2;
		fw_ref_pic = mref[1][v][h];
		/* backward  */
		mv[0] = bw_mv[blkh][blkv][0][bw_mode][0];
		mv[1] = bw_mv[blkh][blkv][0][bw_mode][1];
		h_bw = ((pix_h + mb_h + IMG_PAD_SIZE) << 2) + mv[0];
		v_bw = ((pix_v + mb_v + IMG_PAD_SIZE) << 2) + mv[1];
		if (h_bw<0) 
			h_bw &= 3;
		else if (h_bw>width4) 
			h_bw = width4 + (h_bw & 3);
		if (v_bw<0) 
			v_bw &= 3;
		else if (v_bw>height4) 
			v_bw = height4 + (v_bw & 3);
		h = h_bw & 3, v = v_bw & 3;
		h_bw >>= 2, v_bw >>= 2;
		bw_ref_pic = mref[0][v][h];

		for (pos = (int_16_t*)curr_blk, it = 0; it<8; it+=2)
		{
			d1 = &imgY_original[pix_v + mb_v + it][pix_h + mb_h];
			d2 = &imgY_original[pix_v + mb_v + it + 1][pix_h + mb_h];
			d3 = &fw_ref_pic[v_fw + it][h_fw];
			d4 = &fw_ref_pic[v_fw + it + 1][h_fw];
			d5 = &bw_ref_pic[v_bw + it][h_bw];
			d6 = &bw_ref_pic[v_bw + it + 1][h_bw];
			d7 = &img->mpr[mb_v + it][mb_h];
			d8 = &img->mpr[mb_v + it + 1][mb_h];

			_asm
			{
				mov esi, dword ptr [d1]
				movdqu xmm0, xmmword ptr [esi]
				mov esi, dword ptr [d2]
				movdqu xmm1, xmmword ptr [esi]
				mov esi, dword ptr [d3]
				movdqu xmm2, xmmword ptr [esi]
				mov esi, dword ptr [d4]
				movdqu xmm3, xmmword ptr [esi]
				mov esi, dword ptr [d5]
				movdqu xmm4, xmmword ptr [esi]
				mov esi, dword ptr [d6]
				movdqu xmm5, xmmword ptr [esi]

				pxor xmm7, xmm7
				punpcklbw xmm0, xmm7
				punpcklbw xmm1, xmm7
				punpcklbw xmm2, xmm7
				punpcklbw xmm3, xmm7
				punpcklbw xmm4, xmm7
				punpcklbw xmm5, xmm7

				pavgw xmm2, xmm4
				pavgw xmm3, xmm5
				psubw xmm0, xmm2
				psubw xmm1, xmm3
				
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm0
				add pos, 16
				mov esi, dword ptr [pos]
				movdqa xmmword ptr [esi], xmm1
				add pos, 16

				mov esi, dword ptr [d7]
				movdqa xmmword ptr [esi], xmm2
				mov esi, dword ptr [d8]
				movdqa xmmword ptr [esi], xmm3
			}
		}
		break;
	default:
		break;
	}
	
	if (!skipped)
	{
		avs_dct_sse(curr_blk);
		coeff_cost = scanquant_B8   (img->qp, 0, block8x8, curr_blk, scrFlag, cbp, cbp_blk);
	}
  /*
  The purpose of the action below is to prevent that single or 'expensive' coefficients are coded.
  With 4x4 transform there is larger chance that a single coefficient in a 8x8 or 16x16 block may be nonzero.
  A single small (level=1) coefficient in a 8x8 block will cost: 3 or more bits for the coefficient,
  4 bits for EOBs for the 4x4 blocks,possibly also more bits for CBP.  Hence the total 'cost' of that single
  coefficient will typically be 10-12 bits which in a RD consideration is too much to justify the distortion improvement.
  The action below is to watch such 'single' coefficients and set the reconstructed block equal to the prediction according
  to a given criterium.  The action is taken only for inter luma blocks.

  Notice that this is a pure encoder issue and hence does not have any implication on the standard.
  coeff_cost is a parameter set in dct_luma() and accumulated for each 8x8 block.  If level=1 for a coefficient,
  coeff_cost is increased by a number depending on RUN for that coefficient.The numbers are (see also dct_luma()): 3,2,2,1,1,1,0,0,...
  when RUN equals 0,1,2,3,4,5,6, etc.
  If level >1 coeff_cost is increased by 9 (or any number above 3). The threshold is set to 3. This means for example:
  1: If there is one coefficient with (RUN,level)=(0,1) in a 8x8 block this coefficient is discarded.
  2: If there are two coefficients with (RUN,level)=(1,1) and (4,1) the coefficients are also discarded
  sum_cnt_nonz is the accumulation of coeff_cost over a whole macro block.  If sum_cnt_nonz is 5 or less for the whole MB,
  all nonzero coefficients are discarded for the MB and the reconstructed block is set equal to the prediction.
  */

  if (!skipped && coeff_cost <= _LUMA_COEFF_COST_)
  {
    coeff_cost  = 0;
    (*cbp)     &=  (63 - cbp_mask);
    (*cbp_blk) &= ~(51 << (4*block8x8-2*(block8x8%2)));
    for (v=mb_v; v<mb_v+8; v++)
    {
      //imgY[img->pix_y+j][img->pix_x+i] = img->mpr[j][i];
      h=mb_h;
      xmm0 = _mm_load_si128((__m128i*)(img->mpr[v]+h));
      xmm0 = _mm_packus_epi16(xmm0,zero128);
      _mm_storel_epi64((__m128i *)(imgY[img->pix_y+v]+img->pix_x+h),xmm0);
    }
  }
  return coeff_cost;	
}

void
c_avs_enc::LumaPrediction (int_32_t  *cbp,         //  --> cbp (updated according to processed 8x8 luminance block)
                           int_32_t  *cbp_blk,     //  --> block cbp (updated according to processed 8x8 luminance block)
                           int_32_t  block8x8,     // <--  block number of 8x8 block
                           int_32_t  fw_mode,      // <--  forward  prediction mode (1-7, 0=DIRECT)
                           int_32_t  bw_mode,      // <--  backward prediction mode (1-7, 0=DIRECT)
                           int_32_t  fw_refframe,  // <--  reference frame for forward prediction
                           int_32_t  bw_refframe   // <--  reference frame for backward prediction
                           )
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int_32_t    block_y, block_x, pic_pix_y, pic_pix_x, cbp_blk_mask;
  int_32_t    coeff_cost = 0;
  int_32_t    mb_y       = (block8x8 / 2) << 3;
  int_32_t    mb_x       = (block8x8 % 2) << 3;
  int_32_t    cbp_mask   = 1 << block8x8;
  int_32_t    bxx, byy;                   // indexing curr_blk
  int_32_t    scrFlag = 0;                // 0=noSCR, 1=strongSCR, 2=jmSCR
  byte** imgY_original = imgY_org;
  int_32_t  pix_x    = img->pix_x;
  int_32_t  pix_y    = img->pix_y;
  int_32_t  skipped  = (fw_mode == 0 && bw_mode == 0 && (img->type!=B_IMG));

  if (img->type==B_IMG)
    scrFlag = 1;
  //===== loop over 4x4 blocks =====
  for (byy=0, block_y=mb_y; block_y<mb_y+8; byy+=4, block_y+=4)
  {
    pic_pix_y = pix_y + mb_y;
    for (bxx=0, block_x=mb_x; block_x<mb_x+8; bxx+=4, block_x+=4)
    {
      pic_pix_x = pix_x + block_x;
      cbp_blk_mask = (block_x>>2) + block_y;
      //===== prediction of 4x4 block =====
      LumaPrediction4x4 (block_x, block_y, fw_mode, bw_mode, fw_refframe, bw_refframe);
    }
  }
  return ;
}

/*
*************************************************************************
* Function:Set mode parameters and reference frames for an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::SetModesAndRefframe (int_32_t b8, int_32_t* fw_mode, int_32_t* bw_mode, int_32_t* fw_ref, int_32_t* bw_ref)
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int_32_t         j      = (b8/2);
  int_32_t         i      = (b8%2);
  int_32_t**     frefarr = refFrArr;       // For MB level field/frame coding
  int_32_t**     fw_refarr = fw_refFrArr;  // For MB level field/frame coding
  int_32_t**     bw_refarr = bw_refFrArr;  // For MB level field/frame coding
  int_32_t       block_x = img->block_x;
  int_32_t       block_y = img->block_y;   // For MB level field/frame coding

  *fw_mode = *bw_mode = *fw_ref = *bw_ref = -1;

  if (img->type != B_IMG)
  {
    *fw_ref = frefarr[(block_y>>1)+j][(block_x>>1)+i];
    *bw_ref = 0;
    *bw_mode  = 0;
    *fw_mode  = currMB->b8mode[b8];
  }
  else
  {
    if (currMB->b8pdir[b8]==-1)
    {
      *fw_ref   = -1;
      *bw_ref   = -1;
      *fw_mode  =  0;
      *bw_mode  =  0;
    }
    else if (currMB->b8pdir[b8]==0)
    {
      *fw_ref   = fw_refarr[(block_y>>1)+j][(block_x>>1)+i];
      *bw_ref   = 0;
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = 0;
    }
    else if (currMB->b8pdir[b8]==1)
    {
      *fw_ref   = 0;
      *bw_ref   = bw_refarr[(block_y>>1)+j][(block_x>>1)+i];
      *fw_mode  = 0;
      *bw_mode  = currMB->b8mode[b8];
    }
    else
    {
      *fw_ref   = fw_refarr[(block_y>>1)+j][(block_x>>1)+i];
      *bw_ref   = bw_refarr[(block_y>>1)+j][(block_x>>1)+i];
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = currMB->b8mode[b8];

      if (currMB->b8mode[b8]==0) // direct
      {
        if (img->type==B_IMG)
        {
          *fw_ref = 0;// max(0,frefarr[(block_y>>1)+j][(block_x>>1)+i]);
          *bw_ref = 0;
        }
        else
        {
          *fw_ref = max(0,frefarr[(block_y>>1)+j][(block_x>>1)+i]);
          *bw_ref = 0;
        }
      }
    }
  }
}

/*
*************************************************************************
* Function:Residual Coding of a Luma macroblock (not for intra)
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::LumaResidualCoding ()
{
  int_32_t j,block8x8;
  int_32_t fw_mode, bw_mode, refframe;
  int_32_t sum_cnt_nonz;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int_32_t skipped;
  byte** imgY_original = imgY_org;
  int_32_t bw_ref;

  __m128i xmm0,zero128={0};

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;
  sum_cnt_nonz    = 0 ;

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;
  sum_cnt_nonz    = 0 ;

  for (block8x8=0; block8x8<4; block8x8++)
  {
    SetModesAndRefframe (block8x8, &fw_mode, &bw_mode, &refframe, &bw_ref);
    skipped       = (fw_mode == 0 && bw_mode == 0 && (img->type!=B_IMG));
    sum_cnt_nonz += LumaResidualCoding8x8 (&(currMB->cbp), &(currMB->cbp_blk), block8x8,fw_mode, bw_mode, refframe, bw_ref);
  }

  if (sum_cnt_nonz <= 5)
  {
    currMB->cbp     &= 0xfffff0 ;
    currMB->cbp_blk &= 0xff0000 ;
    for (j=0; j < MB_BLOCK_SIZE; j++)
    {
      //imgY[img->pix_y+j][img->pix_x+i]=img->mpr[j][i];
      xmm0 = _mm_load_si128((__m128i*)(img->mpr[j]+0));
      xmm0 = _mm_packus_epi16(xmm0,zero128);
      _mm_storel_epi64((__m128i *)(imgY[img->pix_y+j]+img->pix_x+0),xmm0);

      xmm0 = _mm_load_si128((__m128i*)(img->mpr[j]+8));
      xmm0 = _mm_packus_epi16(xmm0,zero128);
      _mm_storel_epi64((__m128i *)(imgY[img->pix_y+j]+img->pix_x+8),xmm0);
    }
  }
}

/*
*************************************************************************
* Function: Predict one component of a chroma 4x4 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::OneComponentChromaPrediction4x4 (int_16_t* mpred, int_32_t pix_c_x, int_32_t pix_c_y, int_32_t***** mv, int_32_t ref, int_32_t      blocktype,  int_32_t  uv, int_32_t directforword)
{
  int_32_t     i, j, ii, jj, if0, if1, jf0, jf1;
  int_32_t     incr;
  int_32_t*    mvb;
  int_32_t     refframe  = (ref<0 ?      0 :    ref);
  pel_t** refimage;
  int_32_t     je        = pix_c_y + 4;
  int_32_t     ie        = pix_c_x + 4;
  int_32_t     img_pic_c_x = img->pix_c_x;
  int_32_t     img_pic_c_y = img->pix_c_y;
  // int_32_t     scale   = 1;


  __m64 cmp0={0};
  __declspec(align(16))  int_16_t cph[4]={img->width_cr -1,img->height_cr -1,img->width_cr -1,img->height_cr -1};
  __declspec(align(16))  int_16_t cpa[4];



  incr = 1;
  if(img->type==B_IMG && !img->picture_structure) // field coding
  {
    refframe = ref<0 ? ref+2 : ref;
    incr = 2;
  }
  ref = (img->type==B_IMG) ? ref+incr : ref;

  refimage  = mcef [ref][uv];
  for (j=pix_c_y; j<je; j++)
    for (i=pix_c_x; i<ie; i++)
      //for (i=pix_c_x; i<ie; i++)
    {
      mvb  = mv [(i-img_pic_c_x)>>2][(j-img_pic_c_y)>>2][refframe][blocktype];
      ii   = (i<<3) + mvb[0];
      jj   = (j<<3) + mvb[1];
      cpa [0]=ii>>3;
      cpa [1]=jj>>3;
      cpa [2]=(ii+7)>>3;
      cpa [3]=(jj+7)>>3;
      _asm
      {
        lea         eax, cph
          movq        mm0,mmword ptr [eax]
        lea         eax, cpa
          movq        mm1,mmword ptr [eax]
        pminsw      mm1,mm0
          lea         ecx, cmp0
          movq        mm0,mmword ptr [ecx]
        pmaxsw      mm1,mm0
          movq        mmword ptr [eax],mm1
      }
      _mm_empty();
      //ii0  = max (0, min (img->width_cr -1, ii>>3     ));
      //jj0  = max (0, min (img->height_cr-1, jj>>3     ));    // For MB level field/frame -- scale chroma height by 2
      //ii1  = max (0, min (img->width_cr -1, (ii+7)>>3));
      //jj1  = max (0, min (img->height_cr-1, (jj+7)>>3));
      if1  = (ii&7);  if0 = 8-if1;
      jf1  = (jj&7);  jf0 = 8-jf1;
      *mpred++ = (if0 * jf0 * refimage[cpa [1]][cpa [0]] + if1 * jf0 * refimage[cpa [1]][cpa [2]] +if0 * jf1 * refimage[cpa [3]][cpa [0]] +if1 * jf1 * refimage[cpa [3]][cpa [2]] + 32) / 64;
    }

}


/*
*************************************************************************
* Function:Predict one component of a chroma 4x4 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: OneComponentChromaPrediction4x4_dir (int_16_t*     mpred,      //  --> array to store prediction values
                                                      int_32_t      pix_c_x,    // <--  horizontal pixel coordinate of 4x4 block
                                                      int_32_t      pix_c_y,    // <--  vertical   pixel coordinate of 4x4 block
                                                      int_32_t***** mv,         // <--  motion vector array
                                                      int_32_t      ref,        // <--  reference frame parameter (0.../ -1: backward)
                                                      int_32_t      blocktype,  // <--  block type
                                                      int_32_t      uv,
                                                      int_32_t    refframe
                                                      )
{
  int_32_t     i, j, ii, jj, if0, if1, jf0, jf1;
  int_32_t     incr;
  int_32_t*    mvb;
  pel_t** refimage;
  int_32_t     je        = pix_c_y + 4;
  int_32_t     ie        = pix_c_x + 4;
  int_32_t     f1        = 8 , f2=f1-1, f3=f1*f1, f4=f3>>1;
  int_32_t     s1        = 3;
  int_32_t     img_pic_c_y = img->pix_c_y;
  __m64 cmp0={0};
  __declspec(align(16))  int_16_t cph[4]={img->width_cr -1,img->height_cr -1,img->width_cr -1,img->height_cr -1};
  __declspec(align(16))  int_16_t cpa[4];
  //int_32_t   scale   = 1;
  // int_32_t field_mode;

  incr = 1;
  if(img->type==B_IMG && !img->picture_structure)
    incr = 2;

  ref = (img->type==B_IMG) ? ref+incr : ref;

  // field_mode = (!img->picture_structure);
  refimage  = mcef [ref][uv];
  for (j=pix_c_y; j<je; j++)
  {
    for (i=pix_c_x; i<ie; i++)
    {
      mvb  = mv [(i-img->pix_c_x)>>2][(j-img_pic_c_y)>>2][refframe][blocktype];
      {
        int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,delta_PB;
        delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
        delta_P = (delta_P + 512) % 512;
        if(img->picture_structure)
          TRp = (refframe+1)*delta_P;  //the lates backward reference
        else
        {
          TRp = delta_P;//refframe == 0 ? delta_P-1 : delta_P+1;
        }
        delta_PB = 2*(img->tr - img->imgtr_last_P_frm);
        delta_PB = (delta_PB + 512)%512;
        TRp = (TRp+512) % 512;
        if(!img->picture_structure)
        {
          if(img->current_mb_nr_fld < img->total_number_mb) //top field
            DistanceIndexFw =  refframe == 0 ? delta_PB-1:delta_PB;
          else
            DistanceIndexFw =  refframe == 0 ? delta_PB:delta_PB+1;
        }
        else
          DistanceIndexFw = delta_PB;
        //DistanceIndexBw    = TRp - DistanceIndexFw;
        DistanceIndexBw    = (TRp - DistanceIndexFw+512)%512;
        //ii   = (i<<3) - ((mvb[0]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
        //jj   = (j<<3) - ((mvb[1]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
        ii   = (i<<3) - ((mvb[0]*DistanceIndexBw*(512/DistanceIndexFw)+256)>>9);
        jj   = (j<<3) - ((mvb[1]*DistanceIndexBw*(512/DistanceIndexFw)+256)>>9);

      }

      cpa [0]=ii>>3;
      cpa [1]=jj>>3;
      cpa [2]=(ii+7)>>3;
      cpa [3]=(jj+7)>>3;
      _asm
      {
        lea         eax, cph
          movq        mm0,mmword ptr [eax]
        lea         eax, cpa
          movq        mm1,mmword ptr [eax]
        pminsw      mm1,mm0
          lea         ecx, cmp0
          movq        mm0,mmword ptr [ecx]
        pmaxsw      mm1,mm0
          movq        mmword ptr [eax],mm1
      }
      _mm_empty();

      if1  = (ii&7);  if0 = 8-if1;
      jf1  = (jj&7);  jf0 = 8-jf1;

      *mpred++ = (if0 * jf0 * refimage[cpa [1]][cpa [0]] + if1 * jf0 * refimage[cpa [1]][cpa [2]] +if0 * jf1 * refimage[cpa [3]][cpa [0]] +if1 * jf1 * refimage[cpa [3]][cpa [2]] + 32) / 64;
    }
  }
}


/*
*************************************************************************
* Function:Predict an intra chroma 4x4 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: IntraChromaPrediction4x4 (int_32_t  uv,       // <-- colour component
                                           int_32_t  block_x,  // <-- relative horizontal block coordinate of 4x4 block
                                           int_32_t  block_y)  // <-- relative vertical   block coordinate of 4x4 block
{
  int_32_t mode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  int_32_t i, j;
  //===== prediction =====
  for (j=block_y; j<block_y+4; j++)
  {
    for (i=block_x; i<block_x+4; i++)
    {
      img->mpr[j][i] = img->mprr_c[uv][mode][j][i];
    }
  }
}

/*
*************************************************************************
* Function:Predict one chroma 4x4 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/


void c_avs_enc::ChromaPrediction4x4 (int_32_t  uv,           // <-- colour component
                                     int_32_t  block_x,      // <-- relative horizontal block coordinate of 4x4 block
                                     int_32_t  block_y,      // <-- relative vertical   block coordinate of 4x4 block
                                     int_32_t  fw_mode,      // <-- forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                                     int_32_t  bw_mode,      // <-- backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                                     int_32_t  fw_ref_frame, // <-- reference frame for forward prediction (if (<0) -> intra prediction)
                                     int_32_t  bw_ref_frame) // <-- reference frame for backward prediction
{
  __declspec(align(16)) int_16_t fw_pred[16];
  __declspec(align(16)) int_16_t bw_pred[16];


  int_32_t  block_x4  = block_x+4;
  int_32_t  block_y4  = block_y+4;
  int_32_t  pic_pix_x = img->pix_c_x + block_x;
  int_32_t  pic_pix_y = img->pix_c_y + block_y;
  int_16_t   *fpred     = fw_pred;
  int_16_t   *bpred     = bw_pred;
  int_32_t  by = block_y >>2;
  int_32_t  bx = block_x >>2;

  int_32_t  direct    = (fw_mode == 0 && bw_mode == 0 && (img->type==B_IMG));

  int_32_t  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type!=B_IMG));
  int_32_t  *****fmv_array = (fw_mode && bw_mode)?img->all_omv:img->all_mv;    // For MB level frame/field coding

  int_32_t***** bmv_array = img->all_bmv;
  int_32_t fw_ref_idx, bw_ref_idx;
  int_32_t scale_frame;

  int_32_t directforward = (img->type==B_IMG && fw_mode==0);
  int_16_t *d1, *d2, *d3, *d4;
  //  byte add1b[16]={1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0};
  //  __m128i *m1, *m2, *m3;
  //  __m128i add1=_mm_loadu_si128(add1b);

  //===== INTRA PREDICTION =====
  if (fw_ref_frame < 0)
  {
    IntraChromaPrediction4x4 (uv, block_x, block_y);
    return;
  }

  if (direct)
  {
    fw_ref_frame= 0;
    bw_ref_frame= 0;
    if(!img->picture_structure)
    {
      scale_frame = refFrArr[img->block8_y + by][img->block8_x+bx];
      bw_ref_frame = 0;
      if (img->current_mb_nr_fld < img->total_number_mb)
      {
        fw_ref_frame =  scale_frame >= 1 ? 1 : 0;
        bw_ref_frame = scale_frame >= 0 ? 1 : 0;
      }
      else
      {
        fw_ref_frame = (scale_frame == 1 || scale_frame < 0)? 0 : 1;
        bw_ref_frame = 0;
      }

      if(scale_frame < 0)
        bw_ref_frame = 1;
    }
  }

  fw_ref_idx = fw_ref_frame;
  bw_ref_idx = bw_ref_frame;

  //===== INTER PREDICTION =====
  if (fw_mode || (direct && (fw_ref_frame!=-1)) || skipped)
  {
    if(!img->picture_structure)
      OneComponentChromaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, fmv_array , fw_ref_frame, fw_mode, uv, (directforward ?1 :0));
    else
      OneComponentChromaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, fmv_array , (directforward ?0 :fw_ref_frame), fw_mode, uv, (directforward ?1 :0));
  }
  if (bw_mode || (direct && (bw_ref_frame!=-1)) ||(direct && !img->picture_structure && (bw_ref_frame ==-1)))
  {
    if(fw_mode && bw_mode)
      OneComponentChromaPrediction4x4_dir (bw_pred, pic_pix_x, pic_pix_y, fmv_array, (!img->picture_structure) ? -2+bw_ref_frame : -1, bw_mode, uv,fw_ref_frame);
    else
      OneComponentChromaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, bmv_array, (!img->picture_structure) ? -2+bw_ref_frame : -1, bw_mode, uv,0);
  }

  d1 = &img->mpr[block_y][block_x];
  d2 = &img->mpr[block_y+1][block_x];
  d3 = &img->mpr[block_y+2][block_x];
  d4 = &img->mpr[block_y+3][block_x];

  if (direct || (fw_mode && bw_mode))
  {
    /*  bw_ref_num = bw_ref_frame + 1 ;
    fw_ref_num = fw_ref_frame + 1 ;
    fw_scale = img->chroma_scale[fw_ref_num];
    fw_shift = img->chroma_shift[fw_ref_num];
    bw_scale = img->chroma_scale[bw_ref_num];
    bw_shift = img->chroma_shift[bw_ref_num];
    for (j=block_y; j<block_y4; j++)
    {
    for (i=block_x; i<block_x4; i++)
    {
    img->mpr[j][i] = (*fpred + *bpred + 1) / 2;
    img->mpr_weight[j][i] = (((((*fpred)*fw_scale)>>5) + fw_shift) + ((((*bpred)*bw_scale)>>5) + bw_shift) + 1) / 2;
    fpred++ ;
    bpred++ ;
    }
    }
    */
    /*    m1 = (__m128i*) &img->mpr[block_x][block_y];
    m2 = (__m128i*) fpred;
    m3 = (__m128i*) bpred;
    *m1=_mm_add_epi32(*m2,*m3);
    *m1=_mm_add_epi32(*m1,add1);
    *m1=_mm_srli_epi32 (*m1, 1);

    m1 = (__m128i*) &img->mpr[block_x+1][block_y];
    m2 ++;
    *m1=_mm_add_epi32(*m2,*m3);
    *m1=_mm_add_epi32(*m1,add1);
    *m1=_mm_srli_epi32 (*m1, 1);
    m1 = (__m128i*) &img->mpr[block_x+2][block_y];
    m2 ++;
    *m1=_mm_add_epi32(*m2,*m3);
    *m1=_mm_add_epi32(*m1,add1);
    *m1=_mm_srli_epi32 (*m1, 1);
    m1 = (__m128i*) &img->mpr[block_x+3][block_y];
    m2 ++;
    *m1=_mm_add_epi32(*m2,*m3);
    *m1=_mm_add_epi32(*m1,add1);
    *m1=_mm_srli_epi32 (*m1, 1);
    _mm_empty();
    */
    __asm
    {
      movq    mm0, fw_pred
        movq    mm1, fw_pred+8
        movq    mm2, fw_pred+16
        movq    mm3, fw_pred+24
        movq    mm4, bw_pred
        movq    mm5, bw_pred+8
        movq    mm6, bw_pred+16
        movq    mm7, bw_pred+24
        pavgw       mm0, mm4
        pavgw       mm1, mm5
        pavgw       mm2, mm6
        pavgw       mm3, mm7
        mov      esi,  dword ptr [d1]
      movq    xmmword ptr [esi],  mm0
        mov      esi,  dword ptr [d2]
      movq        xmmword ptr [esi],  mm1
        mov      esi,  dword ptr [d3]
      movq    xmmword ptr [esi],  mm2
        mov      esi,  dword ptr [d4]
      movq        xmmword ptr [esi],  mm3
        emms
    }

  }
  else if (fw_mode || skipped)
  {
    /*  fw_scale = img->chroma_scale[fw_ref_frame];
    fw_shift = img->chroma_shift[fw_ref_frame];
    for (j=block_y; j<block_y4; j++)
    {
    for (i=block_x; i<block_x4; i++)
    {
    img->mpr[j][i] = *fpred;
    img->mpr_weight[j][i] = ((((*fpred)*fw_scale)>>5) + fw_shift);
    fpred++;
    }
    }
    */
    /*    // m1 = (__m128i*) &img->mpr[block_x][block_y];
    *(__m128i*) &img->mpr[block_x][block_y] =*(__m128i*) fpred;
    //*m1=*(__m128i*) fpred;
    fpred+=4;
    //m1 = (__m128i*) &img->mpr[block_x+1][block_y];
    //*m1=*(__m128i*) fpred;
    *(__m128i*) &img->mpr[block_x+1][block_y] =*(__m128i*) fpred;
    fpred+=4;
    //m1 = (__m128i*) &img->mpr[block_x+2][block_y];
    //*m1=*(__m128i*) fpred;
    *(__m128i*) &img->mpr[block_x+2][block_y] =*(__m128i*) fpred;
    fpred+=4;
    //m1 = (__m128i*) &img->mpr[block_x+3][block_y];
    //*m1=*(__m128i*) fpred;
    *(__m128i*) &img->mpr[block_x+3][block_y] =*(__m128i*) fpred;
    _mm_empty();
    */
    __asm
    {
      movq    mm0, fw_pred
        movq    mm1, fw_pred+8
        movq    mm2, fw_pred+16
        movq    mm3, fw_pred+24

        mov      esi,  dword ptr [d1]
      movq    xmmword ptr [esi],  mm0
        mov      esi,  dword ptr [d2]
      movq        xmmword ptr [esi],  mm1
        mov      esi,  dword ptr [d3]
      movq    xmmword ptr [esi],  mm2
        mov      esi,  dword ptr [d4]
      movq        xmmword ptr [esi],  mm3
        emms
    }
    //    b=1;
  }
  else
  {
    /*  bw_ref_num = bw_ref_frame + 1;
    bw_scale = img->chroma_scale[bw_ref_num];
    bw_shift = img->chroma_shift[bw_ref_num];
    for (j=block_y; j<block_y4; j++)
    {
    for (i=block_x; i<block_x4; i++)
    {
    img->mpr[j][i] = *bpred;
    img->mpr_weight[j][i] = ((((*bpred)*bw_scale)>>5) + bw_shift);
    bpred++;
    }
    }
    */

    /*    m1 = (__m128i*) &img->mpr[block_x][block_y];
    m2 = (__m128i*) bpred;
    *m1 =*m2;
    m1 = (__m128i*) &img->mpr[block_x+1][block_y];
    m2 ++;
    *m1 =*m2;
    m1 = (__m128i*) &img->mpr[block_x+2][block_y];
    m2 ++;
    *m1 =*m2;
    m1 = (__m128i*) &img->mpr[block_x+3][block_y];
    m2 ++;
    *m1 =*m2;
    _mm_empty();
    */
    __asm
    {

      movq    mm0, bw_pred
        movq    mm1, bw_pred+8
        movq    mm2, bw_pred+16
        movq    mm3, bw_pred+24

        mov      esi,  dword ptr [d1]
      movq    xmmword ptr [esi],  mm0
        mov      esi,  dword ptr [d2]
      movq        xmmword ptr [esi],  mm1
        mov      esi,  dword ptr [d3]
      movq    xmmword ptr [esi],  mm2
        mov      esi,  dword ptr [d4]
      movq        xmmword ptr [esi],  mm3
        emms
    }
  }


}

/*
*************************************************************************
* Function:Chroma residual coding for an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: ChromaResidualCoding (int_32_t* cr_cbp)
{
  __declspec(align(16)) int_16_t tmp_block_88[8][8];
  int_32_t   uv, block8, block_y, block_x, j, i;
  int_32_t   fw_mode, bw_mode, refframe;
  int_32_t   skipped = (img->mb_data[img->current_mb_nr].mb_type == 0 && img->type == INTER_IMG);
  int_32_t   bw_ref;
  int_32_t   tmp_cbp_blk;
  int_32_t   cbpc=0;
  int_32_t   mode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  *cr_cbp=0;
  for (uv=0; uv<2; uv++)
  {
    //===== prediction of chrominance blocks ===d==
    block8 = 0;
    for (block_y=0; block_y<8; block_y+=4)
    {
      for (block_x=0; block_x<8; block_x+=4)
      {
        SetModesAndRefframe (block8, &fw_mode, &bw_mode, &refframe, &bw_ref);
        //fw_mode = 11; bw_mode = 0; refframe = -1; bw_ref = 0;
        ChromaPrediction4x4 (uv, block_x, block_y, fw_mode, bw_mode, refframe, bw_ref);
        if (skipped||img->NoResidueDirect)
        {
          for (j=0; j<4; j++)
          {
            for (i=block_x; i<block_x+4; i++)
            {
              imgUV[uv][img->pix_c_y+block_y+j][img->pix_c_x+i] = (byte)img->mpr[block_y+j][i];
            }
          }

        }
        else
        {
          for (j=0; j<4; j++)
          {
            for (i=block_x; i<block_x+4; i++)
            {
              if (refframe<0)
                img->m7[block_y+j][i] = imgUV_org[uv][img->pix_c_y+block_y+j][img->pix_c_x+i] - img->mprr_c[uv][mode][block_y+j][i];
              else
                img->m7[block_y+j][i] = imgUV_org[uv][img->pix_c_y+block_y+j][img->pix_c_x+i] - img->mpr[block_y+j][i];
            }
          }
        }
        block8++;
      }
    }
    //===== DCT, Quantization, inverse Quantization, IDCT, and Reconstruction =====

    if (!skipped && !img->NoResidueDirect)
    {
      for (j=0;j<8; j++)
      {
        for (i=0;i<8; i++)
        {
          tmp_block_88[j][i] = img->m7[j][i];
        }
      }
      avs_dct_sse(tmp_block_88);
      scanquant_B8_recon(QP_SCALE_CR[img->qp], 4, 4+uv, tmp_block_88, 0, cr_cbp, &tmp_cbp_blk);
    }
  }
  //===== update currMB->cbp =====
  img->mb_data[img->current_mb_nr].cbp += (*cr_cbp);
}

/*
*************************************************************************
* Function:Predict an intra chroma 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: IntraChromaPrediction8x8 (int_32_t *mb_up, int_32_t *mb_left, int_32_t*mb_up_left)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  unsigned char edgepixu[40];
#define EPU (edgepixu+20)
  unsigned char edgepixv[40];
#define EPV (edgepixv+20)
  int_32_t last_pix,new_pix;
  int_32_t bs_x=8;
  int_32_t bs_y=8;
  int_32_t x, y;
  int_32_t i, j, k;
  pel_t** image;
  int_32_t     block_x, block_y, b4;
  int_32_t     img_cx            = img->pix_c_x;
  int_32_t     img_cy            = img->pix_c_y;
  int_32_t     img_cx_1          = img->pix_c_x-1;
  int_32_t     img_cx_4          = img->pix_c_x+4;
  int_32_t     img_cy_1          = img->pix_c_y-1;
  int_32_t     img_cy_4          = img->pix_c_y+4;
  int_32_t     b8_x              = img->pix_c_x>>2;
  int_32_t     b8_y              = img->pix_c_y>>2;
  int_32_t     mb_nr             = img->current_mb_nr;
  int_32_t     mb_width          = img->width>>4;
  byte    stuffUV;
  //int_32_t     mb_available_up   = (img_cy/BLOCK_SIZE == 0 || (img_cy/BLOCK_SIZE >0 && img->ipredmode[1+b8_x][1+b8_y-1]<0)) ? 0 : 1;
  //int_32_t     mb_available_left = (img_cx/BLOCK_SIZE == 0 || (img_cx/BLOCK_SIZE >0 && img->ipredmode[1+b8_x - 1][1+b8_y]<0)) ? 0 : 1;
  //int_32_t     mb_available_up_left = (img_cx/BLOCK_SIZE == 0 || img_cy/BLOCK_SIZE == 0 || (img_cy/BLOCK_SIZE >0 && img->ipredmode[1+b8_x][1+b8_y-1]<0) ||
  //  (img_cx/BLOCK_SIZE >0 && img->ipredmode[1+b8_x - 1][1+b8_y]<0)) ? 0 : 1;
  int  mb_available_up_right=((img_cy==0)||(b8_x>=(img->width_cr/BLOCK_SIZE-2))) ? 0 : (img->mb_data[img->current_mb_nr].slice_nr == img->mb_data[img->current_mb_nr-mb_width+1].slice_nr);
  int mb_available_left_down=((img_cx==0)||(b8_y>=(img->height_cr/BLOCK_SIZE-2))) ? 0 : (img->mb_data[img->current_mb_nr].slice_nr == img->mb_data[img->current_mb_nr+mb_width-1].slice_nr);
  int mb_available_up   = (img_cy == 0) ? 0 : (img->mb_data[img->current_mb_nr].slice_nr == img->mb_data[img->current_mb_nr-mb_width].slice_nr);
  int mb_available_left = (img_cx == 0) ? 0 : (img->mb_data[img->current_mb_nr].slice_nr == img->mb_data[img->current_mb_nr-1].slice_nr);
  int mb_available_up_left = (img_cx/BLOCK_SIZE == 0 || img_cy/BLOCK_SIZE == 0 ) ? 0 : (img->mb_data[img->current_mb_nr].slice_nr == img->mb_data[img->current_mb_nr-mb_width-1].slice_nr);

  int_16_t     ih,iv;
  int_16_t     ib,ic,iaa;
  int_16_t     uv;
  int_16_t     hline[8], vline[8];

  int_32_t     mode;
  int_32_t     best_mode = DC_PRED_8;         //just an initilaization here, should always be overwritten
  int_32_t     cost;
  int_32_t     min_cost;
  /*__declspec(align(16))*/ int_16_t   diff[16];
  int_32_t pred_plane1;

  if (mb_up)
    *mb_up = mb_available_up;
  if (mb_left)
    *mb_left = mb_available_left;
  if( mb_up_left )
    *mb_up_left = mb_available_up_left;
  // compute all chroma intra prediction modes for both U and V
  uv=0;
  if(mb_available_up)
  {
    stuffUV=imgUV[0][img_cy-1][img_cx+bs_x-1];
    for(x=0;x<bs_x;x++)
    {
      EPU[x+1]=imgUV[0][img_cy-1][img_cx+x];
      /*EPU[1+x+bs_x]=stuffUV;*/
    }
    if(mb_available_up_right)
    {
      for(x=0;x<bs_y;x++)
        EPU[1+x+bs_x]=imgUV[uv][img_cy-1][img_cx+bs_x+x];
    }
    else
    {
      for(x=0;x<bs_y;x++)
        EPU[1+x+bs_x]=EPU[bs_x];  //bs_x=8; EPU[9~16]=r[8]
    }

    EPU[0]=imgUV[uv][img_cy-1][img_cx];
  }
  if(mb_available_left)
  {
    stuffUV=imgUV[0][img_cy+bs_y-1][img_cx-1];
    for(y=0;y<bs_y;y++)
    {
      EPU[-1-y]=imgUV[0][img_cy+y][img_cx-1];
      EPU[-1-y-bs_y]=stuffUV;
    }
    EPU[0]=imgUV[uv][img_cy][img_cx-1];
  }
  if(mb_available_up&&mb_available_left)
    EPU[0]=imgUV[uv][img_cy-1][img_cx-1];

  //lowpass (Those emlements that are not needed will not disturb)
  last_pix=EPU[-8];
  for(i=-8;i<=8;i++)
  {
    new_pix=( last_pix + (EPU[i]<<1) + EPU[i+1] + 2 )>>2;
    last_pix=EPU[i];
    EPU[i]=(unsigned char)new_pix;
  }

  uv=1;
  if(mb_available_up)
  {
    stuffUV=imgUV[1][img_cy-1][img_cx+bs_x-1];
    for(x=0;x<bs_x;x++)
    {
      EPV[x+1]=imgUV[1][img_cy-1][img_cx+x];
    }
    if(mb_available_up_right){
      for(x=0;x<bs_y;x++)
        EPV[1+x+bs_x]=imgUV[uv][img_cy-1][img_cx+bs_x+x];
    }
    else{
      for(x=0;x<bs_y;x++)
        EPV[1+x+bs_x]=EPV[bs_x];  //bs_x=8; EPV[9~16]=r[8]
    }

    EPV[0]=imgUV[uv][img_cy-1][img_cx];
  }
  if(mb_available_left)
  {
    stuffUV=imgUV[1][img_cy+bs_y-1][img_cx-1];
    for(y=0;y<bs_y;y++)
    {
      EPV[-1-y]=imgUV[uv][img_cy+y][img_cx-1];
      EPV[-1-y-bs_y]=stuffUV;
    }
    EPV[0]=imgUV[uv][img_cy][img_cx-1];
  }
  if(mb_available_up&&mb_available_left)
    EPV[0]=imgUV[uv][img_cy-1][img_cx-1];

  //lowpass (Those emlements that are not needed will not disturb)
  last_pix=EPV[-8];
  for(i=-8;i<=8;i++)
  {
    new_pix=( last_pix + (EPV[i]<<1) + EPV[i+1] + 2 )>>2;
    last_pix=EPV[i];
    EPV[i]=(unsigned char)new_pix;
  }

  // compute all chroma intra prediction modes for both U and V
  if(!mb_available_up && !mb_available_left)
  {
    for (uv=0; uv<2; uv++)
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
          img->mprr_c[uv][DC_PRED_8][j][i] = 128;

  }


  if(mb_available_up && !mb_available_left)
  {
    for (j=0; j<8; j++)
    {
      for (i=0; i<8; i++)
      {
        img->mprr_c[0][DC_PRED_8][j][i] = EPU[1+i];
        img->mprr_c[1][DC_PRED_8][j][i] = EPV[1+i];
      }
    }
  }
  if(!mb_available_up && mb_available_left)
  {
    for (j=0; j<8; j++)
    {
      for (i=0; i<8; i++)
      {
        img->mprr_c[0][DC_PRED_8][j][i] = EPU[-1-j];
        img->mprr_c[1][DC_PRED_8][j][i] = EPV[-1-j];
      }
    }
  }
  if(mb_available_up && mb_available_left)
  {
    for (j=0; j<8; j++)
    {
      for (i=0; i<8; i++)
      {
        img->mprr_c[0][DC_PRED_8][j][i] = (EPU[1+i]+EPU[-1-j])>>1;
        img->mprr_c[1][DC_PRED_8][j][i] = (EPV[1+i]+EPV[-1-j])>>1;
      }
    }
  }
  // vertical prediction
  for (uv=0; uv<2; uv++)
  {
    image = imgUV[uv];
    // vertical prediction
    if (mb_available_up)
    {
      for (i=0; i<8; i++)
      {
        hline[i] = image[img_cy_1][img_cx+i];
        for (j=0; j<8; j++)
          img->mprr_c[uv][VERT_PRED_8][j][i] = hline[i];
      }
    }
    // horizontal prediction
    if (mb_available_left)
    {
      for (j=0; j<8; j++)
      {
        vline[j] = image[img_cy+j][img_cx_1];
        for (i=0; i<8; i++)
          img->mprr_c[uv][HOR_PRED_8][j][i] = vline[j];
      }
    }
    if (mb_available_up_left)
    {
      ih = (hline[7] - image[img_cy_1][img_cx_1])<<2;
      iv = (vline[7] - image[img_cy_1][img_cx_1])<<2;
      for (i=1;i<4;i++)
      {
        ih += i*(hline[3+i] - hline[3-i]);
        iv += i*(vline[3+i] - vline[3-i]);
      }
      ib=(((ih+1)<<4)+ih)>>5;
      ic=(((iv+1)<<4)+iv)>>5;

      //iaa=16*(hline[7]+vline[7]);
      iaa=(hline[7]+vline[7])<<4;

      for (j=0; j<8; j++)
      {
        pred_plane1=(j-3)*ic;
        for (i=0; i<8; i++)
        {
          img->mprr_c[uv][PLANE_8][j][i]=max(0,min(255,(iaa+(i-3)*ib +pred_plane1 + 16)>>5));// store plane prediction
        }
      }
    }
  }

  if (!input->rdopt) // the rd-opt part does not work correctly (see encode_one_macroblock)
  {                       // since ipredmodes could be overwritten => encoder-decoder-mismatches
    // pick lowest cost prediction mode
    min_cost = 1<<20;
    for (mode=DC_PRED_8; mode<=PLANE_8; mode++)
    {
      if ((mode==VERT_PRED_8 && !mb_available_up) ||
        (mode==HOR_PRED_8 && !mb_available_left) ||
        (mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
        continue;

      cost = 0;
      for (uv=0; uv<2; uv++)
      {
        image = imgUV_org[uv];
        for (b4=0,block_y=0; block_y<8; block_y+=4)
        {
          for (block_x=0; block_x<8; block_x+=4,b4++)
          {
            for (k=0,j=block_y; j<block_y+4; j++)
            {
              for (i=block_x; i<block_x+4; i++,k++)
              {
                diff[k] = image[img_cy+j][img_cx+i] - img->mprr_c[uv][mode][j][i];
              }
            }
            cost += SATD(diff, input->hadamard);
          }
        }
      }

      if (cost < min_cost)
      {
        best_mode = mode;
        min_cost = cost;
      }
    }
    currMB->c_ipred_mode = best_mode;
  }
}


int_32_t c_avs_enc::SubMBType2Value (Macroblock* currMB,int_32_t layer)
{
  int_32_t mbtype;

  if (img->type!=B_IMG)
  {
    if (currMB->mb_type==4)
      currMB->mb_type = P8x8;

    if (currMB->mb_type==I4MB)
      return (img->type==INTRA_IMG ? 0 : 5);

    else if (currMB->mb_type==P8x8)
    {           //delete by xfwang 2004.7.29
      //if (ZeroRef (currMB))
      //  return 5;
      //else
      return 4;
    }
    else
      return currMB->mb_type;
  }
  else
  {
    if(currMB->mb_type==4)
      currMB->mb_type = P8x8;

    mbtype = currMB->mb_type;

    if      (mbtype==0)
      return 0;
    else if (mbtype==I4MB)
      return 5;
    else if (mbtype==P8x8)
      return 4;
    else if (mbtype==2)
      return 1 + currMB->b8pdir[2*layer];
    else return 0;

  }
}


/*
*************************************************************************
* Function:Converts macroblock type to coding value
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: MBType2Value (Macroblock* currMB)
{
  TLS static const int_32_t dir1offset[3]    =  { 1,  2, 3};
  TLS static const int_32_t dir2offset[3][3] =
  {
    { 0,  4,  8},   // 1. block forward
    { 6,  2, 10},   // 1. block backward
    {12, 14, 16}
  };  // 1. block bi-directional
  int_32_t mbtype, pdir0, pdir1;

  if (img->type != B_IMG)
  {
    if (currMB->mb_type==4)
      currMB->mb_type = P8x8;
    if (currMB->mb_type==I4MB)
    {
      return (img->type==INTRA_IMG ? 0 : 5)+NCBP[currMB->cbp][0];
    }
    else if (currMB->mb_type==P8x8)
    {
      return 4;
    }
    else
    {
      return currMB->mb_type;
    }
  }
  else
  {
    if(currMB->mb_type==4)
      currMB->mb_type = P8x8;

    mbtype = currMB->mb_type;
    pdir0  = currMB->b8pdir[0];
    pdir1  = currMB->b8pdir[3];

    if (mbtype==0)
      return 0;
    else if (mbtype==I4MB)
      return 23+NCBP[currMB->cbp][0];// qhg;
    else if (mbtype==P8x8)
      return 22;
    else if (mbtype==1)
      return dir1offset[pdir0];
    else if (mbtype==2)
      return 4 + dir2offset[pdir0][pdir1];
    else
      return 5 + dir2offset[pdir0][pdir1];
  }
}

/*
*************************************************************************
* Function:Writes intra prediction modes for an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeIntra4x4Modes(int_32_t only_this_block)
{
  int_32_t i;
  int_32_t block8x8;
  int_32_t rate;
  int_32_t ipred_array[16],cont_array[16],ipred_number;
  Macroblock    *currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t           *bitCount   = currMB->bitcounter;

  ipred_number=0;

  for(block8x8=0;block8x8<4;block8x8++)
  {
    if( currMB->b8mode[block8x8]==IBLOCK && (only_this_block<0||only_this_block==block8x8) )
    {
      ipred_array[ipred_number]=currMB->intra_pred_modes[block8x8];
      cont_array[ipred_number]=block8x8;
      ipred_number++;
    }
  }

  rate=0;

  for(i=0;i<ipred_number;i++)
  {
    currMB->IntraChromaPredModeFlag = 1;
    currSE->context = cont_array[i];
    currSE->value1  = ipred_array[i];

#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d %d",currSE->value1,currSE->context);
#endif

    /*--- set symbol type and function pointers ---*/
    currSE->type = SE_INTRAPREDMODE;

    /*--- encode and update rate ---*/
    writeSyntaxElement_Intra4x4PredictionMode(currSE, currBitStream);
    bitCount[BITS_COEFF_Y_MB]+=currSE->len;
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }

  return rate;
}

/*
*************************************************************************
* Function:Converts 8x8 block tyoe to coding value
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: B8Mode2Value (int_32_t b8mode, int_32_t b8pdir)
{
  TLS static const int_32_t b8start[8] = {0,0,0,0, 1, 4, 5, 10};
  TLS static const int_32_t b8inc  [8] = {0,0,0,0, 1, 2, 2, 1};

  if (img->type!=B_IMG)
  {
    return (b8mode-4);
  }
  else
  {
    return b8start[b8mode] + b8inc[b8mode] * b8pdir;
  }

}

/*
*************************************************************************
* Function: Codes macroblock header
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeMBHeader (int_32_t rdopt)  // GB CHROMA !!!!!!!!
{
  int_32_t             i;
  int_32_t             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];
  SyntaxElement *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*            bitCount  = currMB->bitcounter;
  int_32_t             no_bits   = 0;
  int_32_t             mb_x      = img->mb_x;
  int_32_t             mb_y      = img->mb_y;
  int_32_t             skip      = currMB->mb_type ? 0:((img->type==B_IMG) ? !currMB->cbp:1);
  //the mb header include skip_mb_number, mb_type
  currMB->IntraChromaPredModeFlag = IS_INTRA(currMB);
  currMB->mb_field = img->field_mode;
  switch(img->type)
  {
  case INTRA_IMG:
    break;
  case INTER_IMG:
    if (input->skip_mode_flag)
    {
      if (currMB->mb_type != 0)
      {
        //===== Run Length Coding: Non-Skipped macorblock =====
        //写入跳过宏块
        currSE->value1  = img->cod_counter;
        currSE->mapping = &c_avs_enc::ue_linfo;
        currSE->type    = SE_MBTYPE;
        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        // Reset cod counter
        img->cod_counter = 0;
        // Put out mb mode
        currSE->value1 = MBType2Value (currMB);
        currSE->value1--;
        currSE->mapping = &c_avs_enc::ue_linfo;
        currSE->type    = SE_MBTYPE;
        writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
        if (img->type==B_IMG)
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else
          snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
      else
      {
        //Run Length Coding: Skipped macroblock
        img->cod_counter++;
        if (img->current_mb_nr == img->mb_no_currSliceLastMB )
        {
          // Put out run
          currSE->value1  = img->cod_counter;
          currSE->mapping = &c_avs_enc::ue_linfo;
          currSE->type    = SE_MBTYPE;
          writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
          // Reset cod counter
          img->cod_counter = 0;
        }
      }
    }
    else
    {
      // Put out mb mode
      currSE->value1  = MBType2Value (currMB);
      currSE->value1--;
      if (currMB->mb_type == 0)
        currSE->value1 = 0;
      else
        currSE->value1++;
      currSE->mapping = &c_avs_enc::ue_linfo;
      currSE->type    = SE_MBTYPE;
      writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
      if (img->type==B_IMG)
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
      else
        snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    break;
  case B_IMG:
    if (input->skip_mode_flag)
    {
      if(currMB->mb_type != 0 || currMB->cbp != 0)
      {
        //===== Run Length Coding: Non-Skipped macorblock =====
        //写入跳过宏块
        currSE->value1  = img->cod_counter;
        currSE->mapping = &c_avs_enc::ue_linfo;
        currSE->type    = SE_MBTYPE;
        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
        // Reset cod counter
        img->cod_counter = 0;
        // Put out mb mode
        currSE->value1  = MBType2Value (currMB);
        currSE->mapping = &c_avs_enc::ue_linfo;
        currSE->type    = SE_MBTYPE;
        writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
        if (img->type==B_IMG)
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else
          snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
      else
      {
        //Run Length Coding: Skipped macroblock
        img->cod_counter++;
        if (img->current_mb_nr == img->mb_no_currSliceLastMB )
        {
          // Put out run
          currSE->value1  = img->cod_counter;
          currSE->mapping = &c_avs_enc::ue_linfo;
          currSE->type    = SE_MBTYPE;
          writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
          // Reset cod counter
          img->cod_counter = 0;
        }
      }
    }
    else
    {
      // Put out mb mode
      if (currMB->mb_type == 0 && currMB->cbp == 0)
      {
        currSE->value1 = 0;
      }
      else
      {
        currSE->value1  = MBType2Value (currMB);
        currSE->value1++;
      }
      currSE->mapping = &c_avs_enc::ue_linfo;
      currSE->type    = SE_MBTYPE;
      writeSyntaxElement_UVLC( currSE, currBitStream);
#if TRACE
      if (img->type==B_IMG)
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
      else
        snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    if (IS_P8x8 (currMB))
    {
      for (i=0; i<4; i++)
      {
        //mb_part_type is fix length coding(fix length equal 2)
        currSE->value1  =  currSE->bitpattern = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
        currSE->type    = SE_MBTYPE;
        currSE->len     = 2;
        writeUVLC2buffer(currSE, currBitStream);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "8x8 mode/pdir(%2d) = %3d/%d", i,currMB->b8mode[i],currMB->b8pdir[i]);
#endif
        bitCount[BITS_MB_MODE]+= currSE->len;
        no_bits               += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
    break;
  }
  if(!currMB->IntraChromaPredModeFlag && !rdopt) //GB CHROMA !!!!!
    currMB->c_ipred_mode = DC_PRED_8; //setting c_ipred_mode to default is not the right place here
  if(IS_INTRA(currMB))
  {
    for (i=0; i<4; i++)
    {
      currSE->context = i;
      currSE->value1  = currMB->intra_pred_modes[i];
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d %d",currSE->value1,currSE->context);
#endif
      /*--- set symbol type and function pointers ---*/
      currSE->type = SE_INTRAPREDMODE;
      /*--- encode and update rate ---*/
      writeSyntaxElement_Intra4x4PredictionMode(currSE, currBitStream);
      bitCount[BITS_MB_MODE]+=currSE->len;
      no_bits += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    currSE->mapping = &c_avs_enc::ue_linfo;
    currSE->value1 = currMB->c_ipred_mode;
    currSE->type = SE_INTRAPREDMODE;
    writeSyntaxElement_UVLC(currSE, currBitStream);
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Chroma intra pred mode %d", currMB->c_ipred_mode);
#endif
    currSE++;
    currMB->currSEnr++;
  }
  return no_bits;
}

/*
*************************************************************************
* Function:Write chroma intra prediction mode.
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeChromaIntraPredMode()
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*            bitCount  = currMB->bitcounter;
  int_32_t             rate      = 0;

  //===== BITS FOR CHROMA INTRA PREDICTION MODES
  currSE->mapping = &c_avs_enc::ue_linfo;

  currSE->value1 = currMB->c_ipred_mode;
  currSE->type = SE_INTRAPREDMODE;
  writeSyntaxElement_UVLC(currSE, currBitStream);
  bitCount[BITS_COEFF_UV_MB] += currSE->len;
  rate                    += currSE->len;

#if TRACE
  snprintf(currSE->tracestring, TRACESTRING_SIZE, "Chroma intra pred mode");
#endif

  currSE++;
  currMB->currSEnr++;

  return rate;
}
/*
*************************************************************************
* Function:Writes motion vectors of an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc:: writeMotionVector8x8 (int_32_t i0, int_32_t  j0, int_32_t  i1, int_32_t  j1, int_32_t  refframe, int_32_t  dmv_flag, int_32_t  fwd_flag, int_32_t  mv_mode)
{
  int_32_t            i, j, k, l, m;
  int_32_t            curr_mvd;
  int_32_t            bwflag      = ((refframe<0 || (!fwd_flag))?1:0);
  int_32_t            rate        = 0;
  int_32_t            step_h      = input->blc_size[mv_mode][0] >> 3;
  int_32_t            step_v      = input->blc_size[mv_mode][1] >> 3;
  Macroblock*    currMB      = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE      = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*           bitCount    = currMB->bitcounter;
  int_32_t            refindex    = (refframe<0 ? 0 : refframe);
  int_32_t*****       all_mv      = (fwd_flag ? img->all_mv : img->all_bmv);
  int_32_t*****       pred_mv     = ((img->type!=B_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));

  for (j=j0; j<j1; j+=step_v)
  {
    for (i=i0; i<i1; i+=step_h)
    {
      for (k=0; k<2; k++)
      {
        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];
        //--- store (oversampled) mvd ---
        for (l=0; l < step_v; l++)
        {
          for (m=0; m < step_h; m++)
          {
            currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;
          }
        }
        currSE->value1 = curr_mvd;
        currSE->type   =  SE_MVD;
        currSE->mapping = &c_avs_enc::se_linfo;
        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        if (fwd_flag)
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "FMVD(%d) = %3d  (org_mv %3d pred_mv %3d) %d",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k],currSE->value2);
        }
        else
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "BMVD(%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k]);
        }
#endif

        bitCount[BITS_INTER_MB] += currSE->len;
        rate                    += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
  }
  return rate;
}

/*
*************************************************************************
* Function:Writes motion vectors of an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeMotionVector8x8_bid (int_32_t  i0,int_32_t  j0, int_32_t  i1,int_32_t  j1, int_32_t  refframe, int_32_t  dmv_flag, int_32_t  fwd_flag, int_32_t  mv_mode, int_32_t  pdir)
{
  int_32_t            i, j, k, l, m;
  int_32_t            curr_mvd;
  int_32_t            bwflag     = ((refframe<0 || (!fwd_flag))?1:0);
  int_32_t            rate       = 0;
  int_32_t            step_h     = input->blc_size[mv_mode][0] >> 3;
  int_32_t            step_v     = input->blc_size[mv_mode][1] >> 3;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*           bitCount   = currMB->bitcounter;
  int_32_t            refindex   = (refframe<0 ? 0 : refframe);
  int_32_t*****       all_mv     = (fwd_flag ? img->all_mv : img->all_bmv);
  int_32_t*****       pred_mv    = ((img->type!=B_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));
  all_mv = img->all_omv;
  pred_mv = img->omv;
  for (j=j0; j<j1; j+=step_v)
    for (i=i0; i<i1; i+=step_h)
    {
      if(img->nb_references>1 && img->type!=B_IMG)
      {
        currSE->value1 = refindex;
        currSE->type   = SE_REFFRAME;//Attention

        currSE->mapping = &c_avs_enc::se_linfo;

        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "Reference Index = %d", refindex);
        }
#endif
        bitCount[BITS_INTER_MB] += currSE->len;
        rate                    += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }

      for (k=0; k<2; k++)
      {

        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];

        //--- store (oversampled) mvd ---
        for (l=0; l < step_v; l++)
          for (m=0; m < step_h; m++)    currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;

        currSE->value1 = curr_mvd;
        currSE->type   = (img->type==B_IMG ? SE_BFRAME : SE_MVD);

        currSE->mapping = &c_avs_enc::se_linfo;
        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        if (fwd_flag)
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "FMVD(%d) = %3d  (org_mv %3d pred_mv %3d) %d",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k],currSE->value2);
        }
        else
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "BMVD(%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k]);
        }
#endif
        bitCount[BITS_INTER_MB] += currSE->len;
        rate                    += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }

    return rate;
}
/*
*************************************************************************
* Function:Writes Luma Coeff of an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeLumaCoeff8x8 (int_32_t block8x8, int_32_t intra4x4mode)
{
  int_32_t  rate = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];

  rate = writeLumaCoeffAVS_B8(block8x8,intra4x4mode);

  return rate;
}

/*
*************************************************************************
* Function:Writes CBP, DQUANT, and Luma Coefficients of an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeCBPandLumaCoeff ()
{
  int_32_t             i;
  int_32_t             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t             cbp       = currMB->cbp;//((((currMB->cbp)&0x300)>>4) + ((currMB->cbp)&0xf));

  int_16_t*  DCLevel = img->cofDC[0][0];
  int_16_t*  DCRun   = img->cofDC[0][1];
  int_32_t   mb_xpos = img->mb_x;
  int_32_t   mb_ypos = img->mb_y;

  //=====   C B P   =====
  //---------------------

  //=====  L U M I N A N C E   =====
  //--------------------------------
  for (i=0; i<4; i++)
  {
    if( currMB->b8mode[i]==IBLOCK && i < 4)
    {
      currSE->context = i;
      currSE->value1  = currMB->intra_pred_modes[i];

#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d %d",currSE->value1,currSE->context);
#endif

      /*--- set symbol type and function pointers ---*/
      currSE->type = SE_INTRAPREDMODE;


      /*--- encode and update rate ---*/
      writeSyntaxElement_Intra4x4PredictionMode(currSE, currBitStream);
      bitCount[BITS_COEFF_Y_MB]+=currSE->len;
      rate += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }

    if (cbp & (1<<i))
    {
      rate += writeLumaCoeff8x8 (i, (currMB->b8mode[i]==IBLOCK));
    }
  }

  return rate;
}

int_32_t c_avs_enc:: StoreMotionVector8x8 (int_32_t i0, int_32_t j0, int_32_t i1, int_32_t j1, int_32_t refframe, int_32_t dmv_flag, int_32_t fwd_flag, int_32_t mv_mode)
{
  int_32_t            i, j, k, l, m;
  int_32_t            curr_mvd;
  int_32_t            bwflag     = ((refframe<0 || (!fwd_flag))?1:0);
  int_32_t            rate       = 0;
  int_32_t            step_h     = input->blc_size[mv_mode][0] >> 3;
  int_32_t            step_v     = input->blc_size[mv_mode][1] >> 3;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*           bitCount   = currMB->bitcounter;

  int_32_t            refindex   = (refframe<0 ? 0 : refframe);
  int_32_t*****       all_mv     = (fwd_flag ? img->all_mv : img->all_bmv);
  int_32_t*****       pred_mv    = ((img->type!=B_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));

  if (!fwd_flag) bwflag = 1;

  for (j=j0; j<j1; j+=step_v)
  {
    for (i=i0; i<i1; i+=step_h)
    {
      for (k=0; k<2; k++)
      {
        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];
        //--- store (oversampled) mvd ---
        for (l=0; l < step_v; l++)
        {
          for (m=0; m < step_h; m++)
          {
            currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;
          }
        }
      }
    }
  }
  return 0;
}

int_32_t c_avs_enc::StoreMotionVector8x8_bid (int_32_t i0, int_32_t  j0, int_32_t  i1, int_32_t  j1, int_32_t  refframe, int_32_t  dmv_flag, int_32_t  fwd_flag, int_32_t  mv_mode, int_32_t  pdir)
{
  int_32_t            i, j, k, l, m;
  int_32_t            curr_mvd;
  int_32_t            bwflag     = ((refframe<0 || (!fwd_flag))?1:0);
  int_32_t            rate       = 0;
  int_32_t            step_h     = input->blc_size[mv_mode][0] >> 3;
  int_32_t            step_v     = input->blc_size[mv_mode][1] >> 3;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*           bitCount   = currMB->bitcounter;

  int_32_t            refindex   = (refframe<0 ? 0 : refframe);
  int_32_t*****       all_mv     = (fwd_flag ? img->all_mv : img->all_bmv);
  int_32_t*****       pred_mv    = ((img->type!=B_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));
  if (!fwd_flag) bwflag = 1;

  if(pdir == 2 && mv_mode != 0)
  {
    all_mv = img->all_omv;
    pred_mv = img->omv;
  }

  for (j=j0; j<j1; j+=step_v)
  {
    for (i=i0; i<i1; i+=step_h)
    {
      for (k=0; k<2; k++)
      {
        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];

        //--- store (oversampled) mvd ---
        for (l=0; l < step_v; l++)
        {
          for (m=0; m < step_h; m++)
          {
            currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;
          }
        }
      }
    }
  }
  return 0;
}

/*
*************************************************************************
* Function:Writes motion vectors of an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc:: writeMVD8x8 (int_32_t  i0, int_32_t  j0, int_32_t  i1, int_32_t  j1, int_32_t  refframe, int_32_t  dmv_flag, int_32_t  fwd_flag, int_32_t  mv_mode)
{
  int_32_t            i, j, k;
  int_32_t            bwflag     = ((refframe<0 || (!fwd_flag))?1:0);
  int_32_t            rate       = 0;
  int_32_t            step_h     = input->blc_size[mv_mode][0] >> 3;
  int_32_t            step_v     = input->blc_size[mv_mode][1] >> 3;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*           bitCount   = currMB->bitcounter;

  int_32_t            refindex   = (refframe<0 ? 0 : refframe);
  int_32_t*****       all_mv     = (fwd_flag ? img->all_mv : img->all_bmv);
  int_32_t*****       pred_mv    = ((img->type!=B_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));

  for (j=j0; j<j1; j+=step_v)
  {
    for (i=i0; i<i1; i+=step_h)
    {
      for (k=0; k<2; k++)
      {
        currSE->value1 = currMB->mvd[bwflag][j][i][k];
        currSE->type= SE_MVD;
        currSE->mapping = &c_avs_enc::se_linfo;
        writeSyntaxElement_UVLC(currSE, currBitStream);
#if TRACE
        {
          if (fwd_flag)
          {
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "FMVD(%d) = %3d  (org_mv %3d pred_mv %3d) %d",k, currSE->value1, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k],currSE->value2);
          }
          else
          {
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "BMVD(%d) = %3d  (org_mv %3d pred_mv %3d)",k, currSE->value1, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k]);
          }
        }
#endif
        bitCount[BITS_INTER_MB] += currSE->len;
        rate                    += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
  }
  return rate;
}

void c_avs_enc:: writeCBPandDqp (int_32_t *CBPRate)
{
  int_32_t             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];

  //=====   C B P   =====
  //---------------------
  currSE->value1 = currMB->cbp; //CBP
  if (currMB->mb_type == I4MB)
  {
    currSE->mapping = &c_avs_enc::cbp_linfo_intra;
    currSE->type = SE_CBP_INTRA;
  }
  else
  {
    currSE->mapping = &c_avs_enc::cbp_linfo_inter;
    currSE->type = SE_CBP_INTER;
  }
  if(img->type==INTRA_IMG||IS_INTER(currMB)) //inter && img->type == I4MB的情况呢 可能是bug？
  {
    writeSyntaxElement_UVLC(currSE, currBitStream);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                   = currSE->len;

#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->cbp);
#endif
    currSE++;
    currMB->currSEnr++;
  }

  //=====   DQUANT   =====
  //----------------------
  if (currMB->cbp && !input->fixed_picture_qp)
  {
    currSE->value1 = currMB->delta_qp;
    currSE->mapping = &c_avs_enc::se_linfo;
    currSE->type = SE_DELTA_QUANT_INTER;
    writeSyntaxElement_UVLC(currSE, currBitStream);
    bitCount[BITS_DELTA_QUANT_MB] += currSE->len;
    rate                          += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Delta QP (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->delta_qp);
#endif
    // proceed to next SE
    currSE++;
    currMB->currSEnr++;
  }
  *CBPRate += rate;
}



/*
*************************************************************************
* Function:Writes motion info
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: storeMotionInfo (int_32_t pos)
{
  int_32_t k, j0, i0, refframe;

  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t             no_bits   = 0;

  int_32_t   bframe          = (img->type==B_IMG);
  int_32_t** refframe_array  = ((img->type==B_IMG) ? fw_refFrArr : refFrArr);
  int_32_t** bw_refframe_array  = bw_refFrArr;
  int_32_t   multframe       = (input->no_multpred>1);
  int_32_t   step_h0         = (input->blc_size[(IS_P8x8(currMB))? 4 : currMB->mb_type][0] >> 3);
  int_32_t   step_v0         = (input->blc_size[(IS_P8x8(currMB))? 4 : currMB->mb_type][1] >> 3);
  int_32_t   b8_y     = img->block8_y;
  int_32_t   b8_x     = img->block8_x;

  //===== write forward motion vectors =====
  if (IS_INTERMV (currMB))
  {
    for (j0=pos; j0<2; j0+=step_v0)
    {
      for (i0=0; i0<2; i0+=step_h0)
      {
        k=j0*2+i0;
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
        {
          refframe  = refframe_array[b8_y+j0][b8_x+i0];

          if(currMB->b8pdir[k]==2)
          {
            StoreMotionVector8x8_bid(i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 1, currMB->b8mode[k],currMB->b8pdir[k]);
          }
          else
            StoreMotionVector8x8(i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 1, currMB->b8mode[k]);

        }
      }
    }
  }

  //===== write backward motion vectors =====
  if (IS_INTERMV (currMB) && bframe)
  {
    for (j0=pos; j0<2; j0+=step_v0)
    {
      for (i0=0; i0<2; i0+=step_h0)
      {
        k=j0*2+i0;
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
        {
          refframe  = bw_refframe_array[b8_y+j0][b8_x+i0];
          StoreMotionVector8x8(i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 0, currMB->b8mode[k]);
        }
      }
    }
  }
  return no_bits;
}
int_32_t c_avs_enc:: writeFrameRef (int_32_t  mode, int_32_t  i, int_32_t  j, int_32_t  fwd_flag, int_32_t  ref)
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int_32_t*            bitCount  = currMB->bitcounter;

  int_32_t             rate      = 0;
  if (img->type != INTER_IMG && img->picture_structure)
  {
    return 0;
  }
  currSE->value1 = ref;
  currSE->bitpattern = currSE->value1;
  currSE->type=SE_REFFRAME;
  currSE->len = 1;
  writeSyntaxElement2Buf_Fixed(currSE, currBitStream);
  bitCount[BITS_INTER_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  if (fwd_flag)
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Fwd Ref frame no %d", currSE->bitpattern);
  }
  else
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Bwd Ref frame no %d", currSE->bitpattern);
  }

#endif
  currSE++;
  currMB->currSEnr++;
  return rate;
}
/*
*************************************************************************
* Function:Writes motion info
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::writeReferenceIndex (int_32_t *RefIndexRate)
{
  int_32_t k, j0, i0, refframe;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t             no_bits   = 0;
  int_32_t** refframe_array  = ((img->type==B_IMG) ? fw_refFrArr : refFrArr);
  int_32_t** bw_refframe_array  = bw_refFrArr;
  int_32_t   multframe       = (input->no_multpred>1);
  int_32_t   step_h0         ;
  int_32_t   step_v0         ;
  int_32_t   b8_y     = img->block8_y;
  int_32_t   b8_x     = img->block8_x;
  if (img->type != INTRA_IMG && img->nb_references > 1 && currMB->mb_type != I4MB)
  {
    step_h0=input->blc_size[ currMB->mb_type ][0] >> 3;
    step_v0=input->blc_size[ currMB->mb_type ][1] >> 3;
    //forward reference
    for (j0=0; j0<2; j0+=step_v0)
    {
      for (i0=0; i0<2; i0 += step_h0)
      {
        k=j0*2+i0;
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          refframe  = refframe_array[b8_y+j0][b8_x+i0];
          no_bits   = writeFrameRef (currMB->b8mode[k], i0, j0, 1, refframe);
        }
        *RefIndexRate += no_bits;
      }
    }

    //backward reference
    if(!img->picture_structure && img->type == B_IMG)
    {
      for (j0=0; j0<2; i0 += step_h0)
      {
        for (i0=0; i0<2; j0+=step_v0)
        {
          k=j0*2+i0;
          if ((currMB->b8pdir[k]==1) && currMB->b8mode[k]!=0)//has backward vector
          {
            refframe = bw_refframe_array[b8_y+j0][b8_x+i0];
            refframe = 1 - refframe;
            no_bits   = writeFrameRef (currMB->b8mode[k], i0, j0, 0, refframe);
          }
          *RefIndexRate += no_bits;
        }
      }
    }
  }
}

/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
void c_avs_enc:: writeMVD (int_32_t *MVDRate)
{
  int_32_t k, j0, i0, refframe;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t             no_bits;
  int_32_t** refframe_array  = ((img->type==B_IMG) ? fw_refFrArr : refFrArr);
  int_32_t** bw_refframe_array  = bw_refFrArr;
  int_32_t   step_h0, step_v0;
  int_32_t   b8_y     = img->block8_y;
  int_32_t   b8_x     = img->block8_x;
  if(currMB->mb_type != I4MB)
  {
    if (IS_P8x8(currMB))
    {
      step_h0 = 1;
      step_v0 = 1;
    }
    else
    {
      step_h0 = input->blc_size[ currMB->mb_type ][0] >> 3;
      step_v0 = input->blc_size[ currMB->mb_type ][1] >> 3;
    }

    for (j0=0; j0<2; j0+=step_v0)
    {
      for (i0=0; i0<2; i0+=step_h0)
      {
        k=j0*2+i0;

        no_bits = 0;
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
        {
          refframe  = refframe_array[b8_y+j0][b8_x+i0];

          if(currMB->b8pdir[k]==2 && input->InterlaceCodingOption==FRAME_CODING)//symmetry
          {
            no_bits  += writeMotionVector8x8_bid (i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 1, currMB->b8mode[k],currMB->b8pdir[k]);
          }
          else
          {
            no_bits   = writeMVD8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 1, currMB->b8mode[k]);
          }
        }
        *MVDRate += no_bits;
      }
    }

    if (img->type==B_IMG)
    {
      for (j0=0; j0<2; j0+=step_v0)
      {
        for (i0=0; i0<2; i0+=step_h0)
        {
          k=j0*2+i0;
          no_bits = 0;
          if ((currMB->b8pdir[k]==1) && currMB->b8mode[k]!=0)//has backward vector
          {
            refframe  = bw_refframe_array[b8_y+j0][b8_x+i0];
            no_bits   = writeMVD8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 0, currMB->b8mode[k]);
          }
          *MVDRate += no_bits;
        }
      }
    }
  }
}


/*
*************************************************************************
* Function:Writes CBP, DQUANT, and Luma Coefficients of an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: writeBlockCoeff (int_32_t block8x8)
{
  int_32_t             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int_32_t*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];

  //=====  L U M I N A N C E   =====
  //--------------------------------
  if (currMB->cbp & (1<<block8x8))
  {
    if(block8x8 < 4)
    {
      rate += writeLumaCoeffAVS_B8 (block8x8, (currMB->b8mode[block8x8]==IBLOCK));
    }
    else
    {
      rate +=  writeChromaCoeffAVS_B8(block8x8);
    }
  }
  return rate;
}


void c_avs_enc:: writeweightflag()
{
  int_32_t             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];
  Bitstream *bitstream      = currBitStream;

  if(!IS_INTRA(currMB))
  {
    if((img->LumVarFlag == 1)&&(img->allframeweight == 0))
    {
      u_v(1,"mb weighting flag",img->mbweightflag,bitstream);
    }
  }
  return ;
}

/*
*************************************************************************
* Function:Passes the chosen syntax elements to the NAL
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: write_one_macroblock (int_32_t eos_bit)
{
  Macroblock* currMB   = &img->mb_data[img->current_mb_nr];
  int_32_t*        bitCount = currMB->bitcounter;
  int_32_t i;
  int_32_t mb_x = img->mb_x;
  int_32_t mb_y = img->mb_y;
  int_32_t dummy;
  //--- write header ---
  writeMBHeader (0);
  //  Do nothing more if copy and inter mode
  if ((currMB->mb_type != 0) || ((img->type==B_IMG) && currMB->cbp != 0))
  {
    writeReferenceIndex(&dummy);
    writeMVD(&dummy);
    writeweightflag();
    writeCBPandDqp(&dummy);
    for (i=0; i < 6; i++)
    {
      writeBlockCoeff (i);
    }
  }

  //--- set total bit-counter ---
  bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE] + bitCount[BITS_COEFF_Y_MB] + bitCount[BITS_INTER_MB]
  + bitCount[BITS_CBP_MB] + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
  stat->bit_slice += bitCount[BITS_TOTAL_MB];

  //Rate control
  img->NumberofMBHeaderBits  = bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB] + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB];
  img->NumberofMBTextureBits = bitCount[BITS_COEFF_Y_MB]+ bitCount[BITS_COEFF_UV_MB];
  img->NumberofTextureBits  += img->NumberofMBTextureBits;
  img->NumberofHeaderBits   += img->NumberofMBHeaderBits;
  /*basic unit layer rate control*/
  if(img->BasicUnit<img->total_number_mb)
  {
    img->NumberofBasicUnitHeaderBits +=img->NumberofMBHeaderBits;
    img->NumberofBasicUnitTextureBits +=img->NumberofMBTextureBits;
  }

  //Rate control
  img->NumberofCodedMacroBlocks++;
}
