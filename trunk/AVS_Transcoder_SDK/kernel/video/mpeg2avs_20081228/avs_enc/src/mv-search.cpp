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
#include <limits.h>
#include <emmintrin.h>
#include <xmmintrin.h>

#include "global.h"
#include "transcoding_type.h"

#ifdef FastME
#include "fast_me.h"
#endif

#define MEDIAN(a,b,c)  (a + b + c - min(a, min(b, c)) - max(a, max(b, c)));

// These procedure pointers are used by motion_search() and one_eigthpel()
static pel_t (*PelY_14) (pel_t**, int_32_t, int_32_t);
static pel_t *(*PelYline_11) (pel_t *, int_32_t, int_32_t);

// Statistics, temporary
TLS int_32_t    max_mvd;
TLS int_32_t*   spiral_search_x;
TLS int_32_t*   spiral_search_y;
TLS int_32_t*   mvbits;
TLS int_32_t*   refbits;
TLS int_32_t*   byte_abs;
TLS int_32_t*** motion_cost;
TLS int_32_t  mcost_tmp, mv_x_tmp, mv_y_tmp;
TLS int_32_t*** motion_cost_bid;

/*
******************************************************************************
*  Function: calculated field or frame distance between current field(frame)
*            and the reference field(frame).
*     Input:
*    Output:
*    Return:
* Attention:
*    Author: Yulj 2004.07.14
******************************************************************************
*/
int_32_t c_avs_enc::calculate_distance(int_32_t blkref, int_32_t fw_bw )  //fw_bw>=: forward prediction.
{
  int_32_t distance=1;
  if ( img->top_bot == -1 )   // frame
  {
    if ( img->type == INTER_IMG ) // P img
    {
      if(blkref==0)
        distance = picture_distance*2 - img->imgtr_last_P_frm*2 ;
      else if(blkref==1)
        distance = picture_distance*2 - img->imgtr_last_prev_P_frm*2;
      else
      {
        assert(0); //only two reference pictures for P frame
      }
    }
    else //B_IMG
    {
      if (fw_bw >=0 ) //forward
        distance = picture_distance*2 - img->imgtr_last_P_frm*2;
      else
        distance = img->imgtr_next_P_frm*2  - img->tr*2;

    }
  }
  else  // field
  {
    if(img->type==INTER_IMG)
    {
      if(img->top_bot==0) //top field
      {
        switch ( blkref )
        {
        case 0:
          distance = picture_distance*2 - img->imgtr_last_P_frm*2 - 1 ;
          break;
        case 1:
          distance = picture_distance*2 - img->imgtr_last_P_frm*2 ;
          break;
        case 2:
          distance = picture_distance*2 - img->imgtr_last_prev_P_frm*2 - 1;
          break;
        case 3:
          distance = picture_distance*2 - img->imgtr_last_prev_P_frm*2 ;
          break;
        }
      }
      else if(img->top_bot==1) // bottom field.
      {
        switch ( blkref )
        {
        case 0:
          distance = 1 ;
          break;
        case 1:
          distance = picture_distance*2 - img->imgtr_last_P_frm*2 ;
          break;
        case 2:
          distance = picture_distance*2 - img->imgtr_last_P_frm*2 + 1;
          break;
        case 3:
          distance = picture_distance*2 - img->imgtr_last_prev_P_frm*2 ;
          break;
        }
      }
      else
      {
        printf("Error. frame picture should not run into this branch.");
        exit(-1);
      }
    }
    else if(img->type==B_IMG)
    {
      assert(blkref==0 || blkref == 1);
      if (fw_bw >= 0 ) //forward
      {
        if(img->top_bot==0) //top field
        {
          switch ( blkref )
          {
          case 0:
            distance = picture_distance*2 - img->imgtr_last_P_frm*2 - 1 ;
            break;
          case 1:
            distance = picture_distance*2 - img->imgtr_last_P_frm*2;
            break;
          }
        }
        else if(img->top_bot==1) // bottom field.
        {
          switch ( blkref )
          {
          case 0:
            distance = picture_distance*2 - img->imgtr_last_P_frm*2 ;
            break;
          case 1:
            distance = picture_distance*2 - img->imgtr_last_P_frm*2 + 1;
            break;
          }
        }
        else
        {
          printf("Error. frame picture should not run into this branch.");
          exit(-1);
        }
      }
      else // backward
      {
        if(img->top_bot==0) //top field
        {
          switch ( blkref )
          {
          case 0:
            distance = img->imgtr_next_P_frm*2 - picture_distance*2;
            break;
          case 1:
            distance = img->imgtr_next_P_frm*2 - picture_distance*2 + 1;
            break;
          }
        }
        else if(img->top_bot==1) // bottom field.
        {
          switch ( blkref )
          {
          case 0:
            distance = img->imgtr_next_P_frm*2 - picture_distance*2 -  1;
            break;
          case 1:
            distance = img->imgtr_next_P_frm*2 - picture_distance*2 ;
            break;
          }
        }
        else
        {
          printf("Error. frame picture should not run into this branch.");
          exit(-1);
        }
      }

    }
  }
  distance = (distance+512)%512;
  return distance;
}

int_32_t c_avs_enc::scale_motion_vector(int_32_t motion_vector, int_32_t currblkref, int_32_t neighbourblkref, int_32_t block_y_pos, int_32_t curr_block_y, int_32_t ref)
{
  int_32_t neighbour_coding_stage = -2;
  int_32_t current_coding_stage = -2;
  int_32_t sign = (motion_vector>0?1:-1);
  int_32_t mult_distance;
  int_32_t devide_distance;

  motion_vector = abs(motion_vector);

  if(motion_vector == 0)
    return 0;
  // The better way is to deinfe ref_frame same as index in motion search function.
  // ref_frame is different from index when it is B field back ward prediction.
  // change ref_frame to index for backwward prediction of B field.
  // ref_frame :   1  |  0  | -2  | -1
  //      index:   1  |  0  |  0  |  1
  //  direction:   f  |  f  |  b  |  b
  if (img->type == B_IMG && !img->picture_structure && ref < 0 )
  {
    currblkref = 1 - currblkref;
    neighbourblkref = 1 - neighbourblkref;
  }
  mult_distance   = calculate_distance(currblkref, ref);
  devide_distance = calculate_distance(neighbourblkref, ref);
  motion_vector = sign*((motion_vector*mult_distance*(512/devide_distance)+256)>>9);
  return motion_vector;
}

/*
*************************************************************************
* Function:setting the motion vector predictor
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::SetMotionVectorPredictor (int_32_t  pmv[2], int_32_t  **refFrArr, int_32_t  ***tmp_mv, int_32_t  ref_frame, int_32_t  mb_pix_x, int_32_t  mb_pix_y, int_32_t  blockshape_x, int_32_t  blockshape_y, int_32_t  ref)
{
  int_32_t pic_block_x          = (img->block_x>>1) + (mb_pix_x>>3);
  int_32_t pic_block_y          = (img->block_y>>1) + (mb_pix_y>>3);
  int_32_t mb_nr                = img->current_mb_nr;
  int_32_t mb_width             = img->width/16;

  int_32_t mb_available_up      = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
  int_32_t mb_available_left    = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);
  int_32_t mb_available_upleft  = (img->mb_x == 0 || img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
  int_32_t mb_available_upright = (img->mb_x >= mb_width-1 || img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);

  int_32_t block_available_up, block_available_left, block_available_upright, block_available_upleft;
  int_32_t mv_a, mv_b, mv_c, mv_d, pred_vec=0;
  int_32_t mvPredType, rFrameL, rFrameU, rFrameUR, rFrameUL;
  int_32_t hv;
  int_32_t mva[3] , mvb[3],mvc[3];
  Macroblock*     currMB = &img->mb_data[img->current_mb_nr];
  /* D B C */
  /* A X   */

  /* 1 A, B, D are set to 0 if unavailable       */
  /* 2 If C is not available it is replaced by D */

  block_available_up   = mb_available_up   || (mb_pix_y > 0);
  block_available_left = mb_available_left || (mb_pix_x > 0);

  if (mb_pix_y > 0)
  {
    if (mb_pix_x < 8)  // first column of 8x8 blocks
    {
      if (mb_pix_y==8)
      {
        if (blockshape_x == 16)
          block_available_upright = 0;
        else
          block_available_upright = 1;
      }
      else
      {
        if (mb_pix_x+blockshape_x != 8)
          block_available_upright = 1;
        else
          block_available_upright = 0;
      }
    }
    else
    {
      if (mb_pix_x+blockshape_x != 16)
        block_available_upright = 1;
      else
        block_available_upright = 0;
    }
  }
  else if (mb_pix_x+blockshape_x != MB_BLOCK_SIZE)
  {
    block_available_upright = block_available_up;
  }
  else
  {
    block_available_upright = mb_available_upright;
  }

  if (mb_pix_x > 0)
  {
    block_available_upleft = (mb_pix_y > 0 ? 1 : block_available_up);
  }
  else if (mb_pix_y > 0)
  {
    block_available_upleft = block_available_left;
  }
  else
  {
    block_available_upleft = mb_available_upleft;
  }

  mvPredType = MVPRED_MEDIAN;

  rFrameL   = block_available_left    ? refFrArr[pic_block_y]  [pic_block_x-1] : -1;
  rFrameU   = block_available_up      ? refFrArr[pic_block_y-1][pic_block_x]   : -1;
  rFrameUR  = block_available_upright ? refFrArr[pic_block_y-1][pic_block_x+blockshape_x/8] : block_available_upleft  ? refFrArr[pic_block_y-1][pic_block_x-1] : -1;
  rFrameUL  = block_available_upleft  ? refFrArr[pic_block_y-1][pic_block_x-1] : -1;

  if((rFrameL != -1)&&(rFrameU == -1)&&(rFrameUR == -1))
    mvPredType = MVPRED_L;
  else if((rFrameL == -1)&&(rFrameU != -1)&&(rFrameUR == -1))
    mvPredType = MVPRED_U;
  else if((rFrameL == -1)&&(rFrameU == -1)&&(rFrameUR != -1))
    mvPredType = MVPRED_UR;
  else if(blockshape_x == 8 && blockshape_y == 16)
  {
    if(mb_pix_x == 0)
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
    else
    {
      if(rFrameUR == ref_frame)
        mvPredType = MVPRED_UR;
    }
  }
  else if(blockshape_x == 16 && blockshape_y == 8)
  {
    if(mb_pix_y == 0)
    {
      if(rFrameU == ref_frame)
        mvPredType = MVPRED_U;
    }
    else
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
  }

  for (hv=0; hv < 2; hv++)
  {
    mva[hv] = mv_a = block_available_left    ? tmp_mv[hv][pic_block_y]  [4+pic_block_x-1]              : 0;
    mvb[hv] = mv_b = block_available_up      ? tmp_mv[hv][pic_block_y-1][4+pic_block_x]                : 0;
    mv_d    = block_available_upleft         ? tmp_mv[hv][pic_block_y-1][4+pic_block_x-1]              : 0;
    mvc[hv] = mv_c = block_available_upright ? tmp_mv[hv][pic_block_y-1][4+pic_block_x+blockshape_x/8] : mv_d;
    // mv_a, mv_b... are not scaled.
    mva[hv] = scale_motion_vector(mva[hv], ref_frame, rFrameL, pic_block_y,   pic_block_y, ref);
    mvb[hv] = scale_motion_vector(mvb[hv], ref_frame, rFrameU, pic_block_y-1, pic_block_y, ref);
    mv_d    = scale_motion_vector(mv_d,    ref_frame, rFrameUL,pic_block_y-1, pic_block_y, ref);
    mvc[hv] = block_available_upright ? scale_motion_vector(mvc[hv], ref_frame, rFrameUR, pic_block_y-1, pic_block_y, ref): mv_d;

    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(hv == 1){
        mva[2] = abs(mva[0] - mvb[0])  + abs(mva[1] - mvb[1]) ;
        mvb[2] = abs(mvb[0] - mvc[0]) + abs(mvb[1] - mvc[1]);
        mvc[2] = abs(mvc[0] - mva[0])  + abs(mvc[1] - mva[1]) ;
        pred_vec = MEDIAN(mva[2],mvb[2],mvc[2]);
        if(pred_vec == mva[2])
        {
          pmv[0] = mvc[0];
          pmv[1] = mvc[1];
        }

        else if(pred_vec == mvb[2])
        {
          pmv[0] = mva[0];
          pmv[1] = mva[1];
        }
        else
        {
          pmv[0] = mvb[0];
          pmv[1] = mvb[1];
        }
      }
      break;
    case MVPRED_L:
      pred_vec = mv_a;
      break;
    case MVPRED_U:
      pred_vec = mv_b;
      break;
    case MVPRED_UR:
      pred_vec = mv_c;
      break;
    default:
      break;
    }
    if(mvPredType != MVPRED_MEDIAN)
      pmv[hv] = pred_vec;
  }
}

/*
*************************************************************************
* Function:Initialize the motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::Init_Motion_Search_Module ()
{
  int_32_t bits, i, imin, imax, k, l;
  int_32_t search_range               = input->search_range;
  int_32_t number_of_reference_frames = img->buf_cycle;
  int_32_t max_search_points          = (2*search_range+1)*(2*search_range+1);
  int_32_t max_ref_bits               = 1 + 2 * (int_32_t)floor(log((float)(max(16,number_of_reference_frames+1))) / log((float)2) + 1e-10);
  int_32_t max_ref                    = (1<<((max_ref_bits>>1)+1))-1;
  int_32_t number_of_subpel_positions = 4 * (2*search_range+3);
  int_32_t max_mv_bits                = 3 + 2 * (int_32_t)ceil (log((float)(number_of_subpel_positions+1)) / log((float)2) + 1e-10);

  max_mvd                        = (1<<( max_mv_bits >>1)   )-1;

  //=====   CREATE ARRAYS   =====
  //-----------------------------
  if ((spiral_search_x = (int_32_t*)calloc(max_search_points, sizeof(int_32_t))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_search_x");
  if ((spiral_search_y = (int_32_t*)calloc(max_search_points, sizeof(int_32_t))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: spiral_search_y");
  if ((mvbits = (int_32_t*)calloc(2*max_mvd+1, sizeof(int_32_t))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: mvbits");
  if ((refbits = (int_32_t*)calloc(max_ref, sizeof(int_32_t))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: refbits");
  if ((byte_abs = (int_32_t*)calloc(512, sizeof(int_32_t))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: byte_abs");

  get_mem3Dint (&motion_cost, 8, 2*(img->buf_cycle+1), 4);

  get_mem3Dint (&motion_cost_bid, 8, 2*(img->buf_cycle+1), 4);

  //--- set array offsets ---
  mvbits   += max_mvd;
  byte_abs += 256;

  //=====   INIT ARRAYS   =====
  //---------------------------
  //--- init array: motion vector bits ---
  mvbits[0] = 1;

  for (bits=3; bits<=max_mv_bits; bits+=2)
  {
    imax = 1    << (bits >> 1);
    imin = imax >> 1;

    for (i = imin; i < imax; i++)
      mvbits[-i] = mvbits[i] = bits;
  }
  //--- init array: reference frame bits ---
  refbits[0] = 1;
  for (bits=3; bits<=max_ref_bits; bits+=2)
  {
    imax = (1   << ((bits >> 1) + 1)) - 1;
    imin = imax >> 1;

    for (i = imin; i < imax; i++)
      refbits[i] = bits;
  }
  //--- init array: absolute value ---
  byte_abs[0] = 0;

  for (i=1; i<256; i++)
    byte_abs[i] = byte_abs[-i] = i;
  //--- init array: search pattern ---
  spiral_search_x[0] = spiral_search_y[0] = 0;
  for (k=1, l=1; l<=max(1,search_range); l++)
  {
    for (i=-l+1; i< l; i++)
    {
      spiral_search_x[k] =  i;  spiral_search_y[k++] = -l;
      spiral_search_x[k] =  i;  spiral_search_y[k++] =  l;
    }
    for (i=-l;   i<=l; i++)
    {
      spiral_search_x[k] = -l;  spiral_search_y[k++] =  i;
      spiral_search_x[k] =  l;  spiral_search_y[k++] =  i;
    }
  }
}

/*
*************************************************************************
* Function:Free memory used by motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc:: Clear_Motion_Search_Module ()
{
  //--- correct array offset ---
  mvbits   -= max_mvd;
  byte_abs -= 256;

  //--- delete arrays ---
  free (spiral_search_x);
  free (spiral_search_y);
  free (mvbits);
  free (refbits);
  free (byte_abs);
  free_mem3Dint (motion_cost, 8);
}

/*
*************************************************************************
* Function:Full pixel block motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::FullPelBlockMotionSearch (pel_t** orig_pic, int_32_t  ref, int_32_t  pic_pix_x, int_32_t  pic_pix_y, int_32_t  blocktype, int_32_t  pred_mv_x, int_32_t  pred_mv_y, int_32_t* mv_x, int_32_t* mv_y, int_32_t  search_range, int_32_t  min_mcost, double lambda, int_32_t debug_flag)
{
  int_32_t   pos, cand_x, cand_y, y, mcost;
  pel_t *orig_line, *ref_line;
  //pel_t *(c_avs_enc::*get_ref_line)(int_32_t, pel_t*, int_32_t, int_32_t);
  pel_t*     ref_pic       = (byte*)(img->type==B_IMG? Refbuf11 [ref+(((byte***)mref==mref_fld)) +1] : Refbuf11[ref]);
  int_32_t   best_pos      = 0;                                        // position with minimum motion cost
  int_32_t   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
  int_32_t   lambda_factor = LAMBDA_FACTOR (lambda);                   // factor for determining lagragian motion cost
  int_32_t   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int_32_t   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int_32_t   blocksize_x8  = blocksize_x >> 3;  // horizontal block size in 4-pel units
  int_32_t   pred_x        = (pic_pix_x << 2) + pred_mv_x;             // predicted position x (in sub-pel units)
  int_32_t   pred_y        = (pic_pix_y << 2) + pred_mv_y;             // predicted position y (in sub-pel units)
  int_32_t   center_x      = pic_pix_x + *mv_x;                        // center position x (in pel units)
  int_32_t   center_y      = pic_pix_y + *mv_y;                        // center position y (in pel units)
  int_32_t   check_for_00  = (blocktype==1 && !input->rdopt && img->type!=B_IMG && ref==0);
  int_32_t   height        = img->height;
  int_32_t   tmp;
  int_32_t   out_flag = 0;
  //===== set function for getting reference picture lines =====
  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++)
  {
    //--- set candidate position (absolute position in pel units) ---
    cand_x = center_x + spiral_search_x[pos];
    cand_y = center_y + spiral_search_y[pos];
    //--- initialize motion cost (cost for motion vector) and check ---
    mcost = MV_COST (lambda_factor, 2, cand_x, cand_y, pred_x, pred_y);
    if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }

    if (mcost >= min_mcost)
      continue;

    tmp = cand_y * img->width;
    //--- add residual cost to motion cost ---
    for (y=0; y<blocksize_y; y++)
    {
      //ref_line  = get_ref_line (blocksize_x, ref_pic, cand_y+y, cand_x);
      if (cand_x+tmp < 0 || cand_x+tmp > img->height*img->width)
      {
        ref_line  = ref_pic + cand_x;
      }
      else
      {
        ref_line  = ref_pic + cand_x + tmp;
      }
      tmp += img->width;
      orig_line = orig_pic [y];
      /*for (x=0; x<blocksize_x8; x++)
      {
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
        mcost += byte_abs[ *orig_line++ - *ref_line++ ];
      }*/
        _asm
      {
      lea eax, orig_line
      mov         eax,dword ptr [eax]
      movq        mm0,mmword ptr [eax]
      lea eax, ref_line
      mov         eax,dword ptr [eax]
      movq        mm1,mmword ptr [eax]
      psadbw      mm1,mm0
      pextrw      eax,mm1,0
      lea ecx, mcost
      add         eax,dword ptr [ecx]
      mov         dword ptr [ecx],eax
      }

      if (blocksize_x8==2)
      {
      orig_line+=8;
      ref_line+=8;
      _asm
      {
      lea eax, orig_line
      mov         eax,dword ptr [eax]
      movq        mm0,mmword ptr [eax]
      lea eax, ref_line
      mov         eax,dword ptr [eax]
      movq        mm1,mmword ptr [eax]
      psadbw      mm1,mm0
      pextrw      eax,mm1,0
      lea ecx, mcost
      add         eax,dword ptr [ecx]
      mov         dword ptr [ecx],eax
      }
      }
      _mm_empty();
      if (mcost >= min_mcost)
      {
        break;
      }
    }

    //--- check if motion cost is less than minimum cost ---
    if (mcost < min_mcost)
    {
      best_pos  = pos;
      min_mcost = mcost;
    }
  }

  //===== set best motion vector and return minimum motion cost =====
  if (best_pos)
  {
    *mv_x += spiral_search_x[best_pos];
    *mv_y += spiral_search_y[best_pos];
  }

  return min_mcost;
}

#ifdef _THREE_STEP_MOTION_SEARCH_
/*
*************************************************************************
* Function:Full pixel block motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t                                               //  ==> minimum motion cost after search
c_avs_enc::TSSMotionSearch (pel_t**   orig_pic,     // <--  original pixel values for the AxB block
                            int_32_t       ref,          // <--  reference frame (0... or -1 (backward))
                            int_32_t       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
                            int_32_t       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
                            int_32_t       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
                            int_32_t       pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
                            int_32_t       pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
                            int_32_t*      mv_x,         // <--> in: search center (x) / out: motion vector (x) - in pel units
                            int_32_t*      mv_y,         // <--> in: search center (y) / out: motion vector (y) - in pel units
                            int_32_t       search_range, // <--  1-d search range in pel units
                            int_32_t       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
                            double    lambda,           // <--  lagrangian parameter for determining motion cost
                            int_32_t block_index
                            )
{
  int_32_t   pos, cand_x, cand_y, y, mcost;
  pel_t *orig_line, *ref_line;
  //pel_t *(c_avs_enc::*get_ref_line)(int_32_t, pel_t*, int_32_t, int_32_t);
  pel_t*  ref_pic       = (byte*)(img->type==B_IMG? Refbuf11 [ref+(((byte***)mref==mref_fld)) +1] : Refbuf11[ref]);
  int_32_t   best_pos      = 0;                                        // position with minimum motion cost
  int_32_t   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
  int_32_t   lambda_factor = LAMBDA_FACTOR (lambda);                   // factor for determining lagragian motion cost
  int_32_t   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int_32_t   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int_32_t   blocksize_x8  = blocksize_x >> 3;  // horizontal block size in 4-pel units
  int_32_t   pred_x        = (pic_pix_x << 2) + pred_mv_x;             // predicted position x (in sub-pel units)
  int_32_t   pred_y        = (pic_pix_y << 2) + pred_mv_y;             // predicted position y (in sub-pel units)
  int_32_t   center_x      = pic_pix_x + *mv_x;                        // center position x (in pel units)
  int_32_t   center_y      = pic_pix_y + *mv_y;                        // center position y (in pel units)
  int_32_t   check_for_00  = (blocktype==1 && !input->rdopt && img->type!=B_IMG && ref==0);
  int_32_t   height        = img->height;
  int_32_t   tmp;
  int_32_t   step_size;
  int_32_t   max_cand[2], min_cand[2];
  int_32_t   tmp_only_motion_cost;
  // xzhao 20080320
  __m128i    xmm0,xmm1;
  __int16    tmp_mcost;
  int_32_t   diamond_search_pattern[2][5] = {{0,1,0,-1,0},{0,0,-1,0,1}};
  int_32_t   total_pos = 0;
  int_32_t   find_flag = 0;
  //int_32_t   diamond_new_center_x, diamond_new_center_y;
  //int_32_t   diamond_new_center_x2, diamond_new_center_y2;
#define _OUTPUT_TRACE_3
#ifdef _OUTPUT_TRACE_
  FILE *pf_trace = NULL;
  if (frame_no < 14 && img->type == INTER_IMG)
  {
    pf_trace = fopen("enc_trace.txt", "a");
  }
#endif
  max_cand[0] = img->width-16/*center_x-pic_pix_x*/; //img->width  + IMG_PAD_SIZE;
  max_cand[1] = img->height-16/*center_y-pic_pix_y*/; // img->height + IMG_PAD_SIZE;
  min_cand[0] = 2;//为了提高速度，整像素mv控制着不超过图像边界
  min_cand[1] = 2;
  center_x = min(center_x, max_cand[0]);
  center_x = max(center_x, min_cand[0]);
  center_y = min(center_y, max_cand[1]);
  center_y = max(center_y, min_cand[1]);

  step_size = search_range;
  //three step motion search
  while(step_size > 1)
  {
    step_size /= 2;
    best_pos = 0;
    for (pos=0; pos<9; pos++)
    {
      tmp_only_motion_cost = 0;
      cand_x = center_x + step_size * three_step_pattern_x[pos];
      cand_y = center_y + step_size * three_step_pattern_y[pos];
      cand_x = min(cand_x, max_cand[0]);
      cand_x = max(cand_x, min_cand[0]);
      cand_y = min(cand_y, max_cand[1]);
      cand_y = max(cand_y, min_cand[1]);
      mcost = MV_COST (lambda_factor, 2, cand_x, cand_y, pred_x, pred_y);
      if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
      {
        mcost -= WEIGHTED_COST (lambda_factor, 16);
      }

      if (mcost >= min_mcost)
        continue;
      tmp = cand_y * img->width;
      //--- add residual cost to motion cost ---
      if(blocksize_x8>1)
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadu_si128((__m128i*)(orig_line));
          xmm1 = _mm_loadu_si128((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          tmp_mcost = _mm_extract_epi16(xmm0,4);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }
      else
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadl_epi64((__m128i*)(orig_line));
          xmm1 = _mm_loadl_epi64((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          if (tmp_mcost < only_motion_cost[blocktype][block_index])
          {
            only_motion_cost[blocktype][block_index] = tmp_mcost;
          }
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }

      //--- check if motion cost is less than minimum cost ---
      if (mcost < min_mcost)
      {
        best_pos  = pos;
        min_mcost = mcost;
      }
      _mm_empty();
#ifdef _OUTPUT_TRACE_1
      if (pf_trace)
      {
        fprintf(pf_trace, "pos:%4d, cand_x:%4d, cand_y:%4d, cost:%4d\n", pos, cand_x, cand_y, mcost);
      }
#endif
    }
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
    //change the search center
    if (best_pos)
    {
      center_x += step_size * three_step_pattern_x[best_pos];
      center_y += step_size * three_step_pattern_y[best_pos];
      center_x = min(center_x, max_cand[0]);
      center_x = max(center_x, min_cand[0]);
      center_y = min(center_y, max_cand[1]);
      center_y = max(center_y, min_cand[1]);
    }
    else
    {
      break;
    }
  }
#ifdef _DIAMOND_SEARCH_
  //change the search center
  center_x += step_size * three_step_pattern_x[best_pos];
  center_y += step_size * three_step_pattern_y[best_pos];
  diamond_new_center_x = pic_pix_x;
  diamond_new_center_y = pic_pix_y;
  //diamond motion search
  while (total_pos <= 20 && best_pos != 0)
  {
    best_pos = 0;
    for (pos=0; pos<5; pos++)
    {
      total_pos++;
      tmp_only_motion_cost = 0;
      cand_x = diamond_new_center_x + diamond_search_pattern[0][pos];
      cand_y = diamond_new_center_y + diamond_search_pattern[1][pos];
      cand_x = min(cand_x, max_cand[0]);
      cand_x = max(cand_x, min_cand[0]);
      cand_y = min(cand_y, max_cand[1]);
      cand_y = max(cand_y, min_cand[1]);
      mcost = MV_COST (lambda_factor, 2, cand_x, cand_y, pred_x, pred_y);
      if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
      {
        mcost -= WEIGHTED_COST (lambda_factor, 16);
      }

      if (mcost >= min_mcost)
        continue;
      tmp = cand_y * img->width;
      //--- add residual cost to motion cost ---
      if(blocksize_x8>1)
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadu_si128((__m128i*)(orig_line));
          xmm1 = _mm_loadu_si128((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          tmp_mcost = _mm_extract_epi16(xmm0,4);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }
      else
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadl_epi64((__m128i*)(orig_line));
          xmm1 = _mm_loadl_epi64((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          if (tmp_mcost < only_motion_cost[blocktype][block_index])
          {
            only_motion_cost[blocktype][block_index] = tmp_mcost;
          }
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }

      //--- check if motion cost is less than minimum cost ---
      if (mcost < min_mcost)
      {
        best_pos  = pos;
        min_mcost = mcost;
        find_flag = 1;
      }
      _mm_empty();
#ifdef _OUTPUT_TRACE_
      if (pf_trace)
      {
        fprintf(pf_trace, "pos:%4d, cand_x:%4d, cand_y:%4d, cost:%4d\n", pos, cand_x, cand_y, mcost);
      }
#endif
    }
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
    if (best_pos != 0 && find_flag == 1)
    {
      diamond_new_center_x += diamond_search_pattern[0][best_pos];
      diamond_new_center_y += diamond_search_pattern[1][best_pos];
    }
  }

  diamond_new_center_x2 = pic_pix_x + *mv_x;
  diamond_new_center_y2 = pic_pix_y + *mv_y;
  while (total_pos <= 20 && best_pos != 0)
  {
    best_pos = 0;
    for (pos=0; pos<5; pos++)
    {
      total_pos++;
      tmp_only_motion_cost = 0;
      cand_x = diamond_new_center_x + diamond_search_pattern[0][pos];
      cand_y = diamond_new_center_y + diamond_search_pattern[1][pos];
      cand_x = min(cand_x, max_cand[0]);
      cand_x = max(cand_x, min_cand[0]);
      cand_y = min(cand_y, max_cand[1]);
      cand_y = max(cand_y, min_cand[1]);
      mcost = MV_COST (lambda_factor, 2, cand_x, cand_y, pred_x, pred_y);
      if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
      {
        mcost -= WEIGHTED_COST (lambda_factor, 16);
      }

      if (mcost >= min_mcost)
        continue;
      tmp = cand_y * img->width;
      //--- add residual cost to motion cost ---
      if(blocksize_x8>1)
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadu_si128((__m128i*)(orig_line));
          xmm1 = _mm_loadu_si128((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          tmp_mcost = _mm_extract_epi16(xmm0,4);
          mcost+=tmp_mcost;
          tmp_only_motion_cost += tmp_mcost;
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }
      else
      {
        for (y=0; y<blocksize_y; y++)
        {
          ref_line  = ref_pic + cand_x + tmp;
          tmp += img->width;
          orig_line = orig_pic [y];

          xmm0 = _mm_loadl_epi64((__m128i*)(orig_line));
          xmm1 = _mm_loadl_epi64((__m128i*)(ref_line));

          xmm0 = _mm_sad_epu8(xmm0,xmm1);

          tmp_mcost = _mm_extract_epi16(xmm0,0);
          mcost+=tmp_mcost;
          if (tmp_mcost < only_motion_cost[blocktype][block_index])
          {
            only_motion_cost[blocktype][block_index] = tmp_mcost;
          }
          if (mcost >= min_mcost)
          {
            break;
          }
        }
      }

      //--- check if motion cost is less than minimum cost ---
      if (mcost < min_mcost)
      {
        best_pos  = pos;
        min_mcost = mcost;
        find_flag = 2;
      }
      _mm_empty();
#ifdef _OUTPUT_TRACE_
      if (pf_trace)
      {
        fprintf(pf_trace, "pos:%4d, cand_x:%4d, cand_y:%4d, cost:%4d\n", pos, cand_x, cand_y, mcost);
      }
#endif
    }
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
    if (best_pos != 0 && find_flag == 2)
    {
      diamond_new_center_x2 += diamond_search_pattern[0][best_pos];
      diamond_new_center_y2 += diamond_search_pattern[1][best_pos];
    }
  }

  //===== set best motion vector and return minimum motion cost =====
  if (find_flag==2)
  {
    *mv_x = diamond_new_center_x2 - pic_pix_x;
    *mv_y = diamond_new_center_y2 - pic_pix_y;
  }
  else if(find_flag == 1)
  {
    *mv_x = diamond_new_center_x - pic_pix_x;
    *mv_y = diamond_new_center_y - pic_pix_y;
  }
  else
#endif
  {
    *mv_x = center_x - pic_pix_x;
    *mv_y = center_y - pic_pix_y;
  }
#ifdef _OUTPUT_TRACE_
  if(pf_trace)
  {
    fprintf(pf_trace, "frameno:%4d, mbnr:%4d, mode:%4d, mvx:%4d, mvy:%4d, \n", frame_no, img->current_mb_nr, blocktype, *mv_x, *mv_y);
    fclose(pf_trace);
  }
#endif
  return min_mcost;
}
void c_avs_enc:: init_3_step_search()
{
  three_step_pattern_x[0] = 0;
  three_step_pattern_y[0] = 0;

  three_step_pattern_x[1] = -1;
  three_step_pattern_y[1] = 0;

  three_step_pattern_x[2] = -1;
  three_step_pattern_y[2] = -1;

  three_step_pattern_x[3] = 0;
  three_step_pattern_y[3] = -1;

  three_step_pattern_x[4] = 1;
  three_step_pattern_y[4] = -1;

  three_step_pattern_x[5] = 1;
  three_step_pattern_y[5] = 0;

  three_step_pattern_x[6] = 1;
  three_step_pattern_y[6] = 1;

  three_step_pattern_x[7] = 0;
  three_step_pattern_y[7] = 1;

  three_step_pattern_x[8] = -1;
  three_step_pattern_y[8] = -1;
}
#endif

/*
*************************************************************************
* Function:Calculate SA(T)D
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc:: SATD (int_16_t* diff, int_32_t use_hadamard)//需要修改成sse指令
{
  int_16_t k, satd = 0, m[16], dd, *d=diff;

  if (use_hadamard)
  {
    /*===== hadamard transform =====*/
    m[ 0] = d[ 0] + d[12];
    m[ 4] = d[ 4] + d[ 8];
    m[ 8] = d[ 4] - d[ 8];
    m[12] = d[ 0] - d[12];
    m[ 1] = d[ 1] + d[13];
    m[ 5] = d[ 5] + d[ 9];
    m[ 9] = d[ 5] - d[ 9];
    m[13] = d[ 1] - d[13];
    m[ 2] = d[ 2] + d[14];
    m[ 6] = d[ 6] + d[10];
    m[10] = d[ 6] - d[10];
    m[14] = d[ 2] - d[14];
    m[ 3] = d[ 3] + d[15];
    m[ 7] = d[ 7] + d[11];
    m[11] = d[ 7] - d[11];
    m[15] = d[ 3] - d[15];

    d[ 0] = m[ 0] + m[ 4];
    d[ 8] = m[ 0] - m[ 4];
    d[ 4] = m[ 8] + m[12];
    d[12] = m[12] - m[ 8];
    d[ 1] = m[ 1] + m[ 5];
    d[ 9] = m[ 1] - m[ 5];
    d[ 5] = m[ 9] + m[13];
    d[13] = m[13] - m[ 9];
    d[ 2] = m[ 2] + m[ 6];
    d[10] = m[ 2] - m[ 6];
    d[ 6] = m[10] + m[14];
    d[14] = m[14] - m[10];
    d[ 3] = m[ 3] + m[ 7];
    d[11] = m[ 3] - m[ 7];
    d[ 7] = m[11] + m[15];
    d[15] = m[15] - m[11];

    m[ 0] = d[ 0] + d[ 3];
    m[ 1] = d[ 1] + d[ 2];
    m[ 2] = d[ 1] - d[ 2];
    m[ 3] = d[ 0] - d[ 3];
    m[ 4] = d[ 4] + d[ 7];
    m[ 5] = d[ 5] + d[ 6];
    m[ 6] = d[ 5] - d[ 6];
    m[ 7] = d[ 4] - d[ 7];
    m[ 8] = d[ 8] + d[11];
    m[ 9] = d[ 9] + d[10];
    m[10] = d[ 9] - d[10];
    m[11] = d[ 8] - d[11];
    m[12] = d[12] + d[15];
    m[13] = d[13] + d[14];
    m[14] = d[13] - d[14];
    m[15] = d[12] - d[15];

    d[ 0] = m[ 0] + m[ 1];
    d[ 1] = m[ 0] - m[ 1];
    d[ 2] = m[ 2] + m[ 3];
    d[ 3] = m[ 3] - m[ 2];
    d[ 4] = m[ 4] + m[ 5];
    d[ 5] = m[ 4] - m[ 5];
    d[ 6] = m[ 6] + m[ 7];
    d[ 7] = m[ 7] - m[ 6];
    d[ 8] = m[ 8] + m[ 9];
    d[ 9] = m[ 8] - m[ 9];
    d[10] = m[10] + m[11];
    d[11] = m[11] - m[10];
    d[12] = m[12] + m[13];
    d[13] = m[12] - m[13];
    d[14] = m[14] + m[15];
    d[15] = m[15] - m[14];

    /*===== sum up =====*/
    for (dd=diff[k=0]; k<16; dd=diff[++k])
    {
      satd += (dd < 0 ? -dd : dd);
    }
    satd >>= 1;
  }
  else
  {
    /*===== sum up =====*/
    for (k = 0; k < 16; k++)
    {
      satd += byte_abs [diff [k]];
    }
  }

  return satd;
}

int_32_t                                               //  ==> minimum motion cost after search
c_avs_enc:: Get_Skip_CostMB(pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                            int_32_t       ref,           // <--  reference frame (0... or -1 (backward))
                            int_32_t       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                            int_32_t       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                            int_32_t       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                            int_32_t       pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                            int_32_t       pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                            int_32_t*      mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                            int_32_t*      mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                            int_32_t       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                            int_32_t       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
                            int_32_t       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                            double    lambda         // <--  lagrangian parameter for determining motion cost
                            )
{
  //int_32_t   diff[16], *d;
  //  char   diff[8];
  int_32_t   best_pos, mcost;
  int_32_t   y0, x0, y1, x1, ry0, rx0;
  int_32_t   cand_mv_x, cand_mv_y, mv_y0, mv_x0;
  //  pel_t *orig_line, ref_line[8];
  int   incr            = ref==-1 ? ((!img->fld_type)&&(!img->picture_structure)&&(img->type==B_IMG)) : ((byte***)mref==mref_fld)&&(img->type==B_IMG) ;
  pel_t ****ref_pic;
  byte  **ref_pic_tmp;
  int_32_t   lambda_factor   = LAMBDA_FACTOR (lambda);
  int_32_t   mv_shift        = 0;
  int_32_t   check_position0 = (blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard && !input->rdopt && img->type!=B_IMG && ref==0);
  int_32_t   blocksize_x     = input->blc_size[blocktype][0];
  int_32_t   blocksize_y     = input->blc_size[blocktype][1];
  int_32_t   pic4_pix_x      = (pic_pix_x << 2);
  int_32_t   pic4_pix_y      = (pic_pix_y << 2);
  int_32_t   max_pos_x4      = ((img->width -blocksize_x+1)<<2);
  int_32_t   max_pos_y4      = ((img->height-blocksize_y+1)<<2);
  int_32_t   min_pos2        = (input->hadamard ? 0 : 1);
  int_32_t   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  //byte  *d1, *d2, *d3, *d4, *d5, *d6, *d7, *d8,*a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8;
  int_32_t   sad = 0;
  int_32_t width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2)-32;
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2)-32;

  //__m128i   xmm_org[4],xmm_ref[4];
  __m128i   xmm_org0,xmm_org1,xmm_org2,xmm_org3;
  __m128i   xmm_ref0,xmm_ref1,xmm_ref2,xmm_ref3;

  ref_pic = img->type==B_IMG? mref [ref+incr] : mref [ref];
  //===== loop over search positions =====
  best_pos = 0;
  cand_mv_x = *mv_x;    // quarter-pel units
  cand_mv_y = *mv_y;    // quarter-pel units
  //----- set motion vector cost -----
  mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
  /*if (check_position0 && pos==0)
  {
  mcost -= WEIGHTED_COST (lambda_factor, 16);
  }*/
  for (y0=0; y0<blocksize_y ; y0+=8)
  {
    y1  = pic_pix_y + y0;
    ry0 = ((y1 + IMG_PAD_SIZE) << 2) + cand_mv_y;
    if (ry0 < 0)
      ry0 &= 3;
    else if (ry0 > height4)
      ry0 = height4 + (ry0 & 3);
    mv_y0 = ry0 % 4;
    ry0 /= 4;

    for (x0=0; x0<blocksize_x; x0+=8)
    {
      x1  = pic_pix_x + x0;
      rx0 = ((x1 + IMG_PAD_SIZE) << 2) + cand_mv_x;
      if (rx0 < 0)
        rx0 &= 3;
      else if (rx0 > width4)
        rx0 = width4 + (rx0 & 3);
      mv_x0 = rx0 % 4;
      rx0 /= 4;

      if (!input->hadamard)
      {
        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+1]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+2]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+3]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+1]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+2]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+3]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad = _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+4]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+5]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+6]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+7]+x1));

        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+4]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+5]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+6]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+7]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad += _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        //_mm_empty()
        mcost += sad;
      }
      else  //input->hadamard
      {

      }
    }
  }
  _mm_empty();
  //===== return minimum motion cost =====
  return mcost;
}
/*
*************************************************************************
* Function:Sub pixel block motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t
c_avs_enc:: SubPelBlockMotionSearch (pel_t**   orig_pic, int_32_t  ref, int_32_t  pic_pix_x, int_32_t  pic_pix_y, int_32_t  blocktype, int_32_t  pred_mv_x, int_32_t  pred_mv_y, int_32_t* mv_x, int_32_t* mv_y, int_32_t  search_pos2, int_32_t  search_pos4, int_32_t  min_mcost, double   lambda, int_32_t block_index)
{
  int_32_t   pos, best_pos, mcost, tmp_only_motion_cost;
  int_32_t   y0, x0, y1, x1, ry0, rx0;
  int_32_t   cand_mv_x, cand_mv_y, mv_y0, mv_x0;
  int_32_t   incr = ref==-1 ? ((!img->fld_type)&&(!img->picture_structure)&&(img->type==B_IMG)) : ((byte***)mref==mref_fld)&&(img->type==B_IMG) ;
  pel_t ****ref_pic;
  byte  **ref_pic_tmp;
  int_32_t   lambda_factor   = LAMBDA_FACTOR (lambda);
  int_32_t   mv_shift        = 0;
  int_32_t   check_position0 = (blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard && !input->rdopt && img->type!=B_IMG && ref==0);
  int_32_t   blocksize_x     = input->blc_size[blocktype][0];
  int_32_t   blocksize_y     = input->blc_size[blocktype][1];
  int_32_t   pic4_pix_x      = (pic_pix_x << 2);
  int_32_t   pic4_pix_y      = (pic_pix_y << 2);
  int_32_t   max_pos_x4      = ((img->width -blocksize_x+1)<<2);
  int_32_t   max_pos_y4      = ((img->height-blocksize_y+1)<<2);
  int_32_t   min_pos2        = (input->hadamard ? 0 : 1);
  int_32_t   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  int_32_t   max_mv[2], min_mv[2];
  int_32_t   sad = 0;
  int_32_t width4  = ((img->width +2*IMG_PAD_SIZE-1)<<2)-32;
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2)-32;

  __m128i   xmm_org0,xmm_org1,xmm_org2,xmm_org3;
  __m128i   xmm_ref0,xmm_ref1,xmm_ref2,xmm_ref3;
  max_mv[0] = (img->width  - pic_pix_x + 16 -1)  << 2;
  max_mv[1] = (img->height - pic_pix_y + 16 - 1) << 2;
  min_mv[0] = (-pic_pix_x-10) << 2;
  min_mv[1] = (-pic_pix_y-10) << 2;
  if (!img->picture_structure)
  {
    if (img->type==B_IMG)
    {
      incr = 2;
    }
  }
  else
  {
    if(img->type==B_IMG)
      incr = 1;
  }

  ref_pic = img->type==B_IMG? mref [ref+incr] : mref [ref];
  /*********************************
  *****                       *****
  *****  HALF-PEL REFINEMENT  *****
  *****                       *****
  *********************************/
  //===== convert search center to quarter-pel units =====
  *mv_x <<= 2;
  *mv_y <<= 2;
  //===== loop over search positions =====
  best_pos = 0;
  for (pos = min_pos2; pos < max_pos2; pos++)
  {
    tmp_only_motion_cost = 0;
    cand_mv_x = *mv_x + (spiral_search_x[pos] << 1);    // quarter-pel units
    cand_mv_y = *mv_y + (spiral_search_y[pos] << 1);    // quarter-pel units
    cand_mv_x = max(min_mv[0], cand_mv_x);
    cand_mv_x = min(max_mv[0], cand_mv_x);
    cand_mv_y = max(min_mv[1], cand_mv_y);
    cand_mv_y = min(max_mv[1], cand_mv_y);
    //----- set motion vector cost -----
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
    //----- add up SATD -----
    for (y0=0; y0<blocksize_y ; y0+=8)
    {
      y1  = pic_pix_y + y0;
      ry0 = ((y1 + IMG_PAD_SIZE) << 2) + cand_mv_y;
      if (ry0 < 0)
        ry0 &= 3;
      else if (ry0 > height4)
        ry0 = height4 + (ry0 & 3);
      mv_y0 = ry0 % 4;
      ry0 /= 4;

      for (x0=0; x0<blocksize_x; x0+=8)
      {
        x1  = pic_pix_x + x0;
        rx0 = ((x1 + IMG_PAD_SIZE) << 2) + cand_mv_x;
        if (rx0 < 0)
          rx0 &= 3;
        else if (rx0 > width4)
          rx0 = width4 + rx0 & 3;
        mv_x0 = rx0 % 4;
        rx0 /= 4;

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+1]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+2]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+3]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+1]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+2]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+3]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad = _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+4]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+5]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+6]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+7]+x1));

        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+4]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+5]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+6]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+7]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad += _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        //_mm_empty()
        mcost += sad;
        tmp_only_motion_cost += sad;
      }
    }
    _mm_empty();
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }


  if (best_pos)
  {
    *mv_x += (spiral_search_x [best_pos] << 1);
    *mv_y += (spiral_search_y [best_pos] << 1);
    cand_mv_x = max(min_mv[0], cand_mv_x);
    cand_mv_x = min(max_mv[0], cand_mv_x);
    cand_mv_y = max(min_mv[1], cand_mv_y);
    cand_mv_y = min(max_mv[1], cand_mv_y);
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
  }

  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/


  //===== loop over search positions =====
  best_pos = 0;
  for (pos = 1; pos < search_pos4; pos++)
  {
    tmp_only_motion_cost = 0;
    cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
    cand_mv_x = max(min_mv[0], cand_mv_x);
    cand_mv_x = min(max_mv[0], cand_mv_x);
    cand_mv_y = max(min_mv[1], cand_mv_y);
    cand_mv_y = min(max_mv[1], cand_mv_y);

    //----- set motion vector cost -----
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);

    //----- add up SATD -----
    for (y0=0; y0<blocksize_y ; y0+=8)
    {
      y1  = pic_pix_y + y0;
      ry0 = ((y1 + IMG_PAD_SIZE) << 2) + cand_mv_y;
      if (ry0 < 0)
        ry0 &= 3;
      else if (ry0 > height4)
        ry0 = height4 + (ry0 & 3);
      mv_y0 = ry0 % 4;
      ry0 /= 4;

      for (x0=0; x0<blocksize_x; x0+=8)
      {
        x1  = pic_pix_x + x0;
        rx0 = ((x1 + IMG_PAD_SIZE) << 2) + cand_mv_x;
        if (rx0 < 0)
          rx0 &= 3;
        else if (rx0 > width4)
          rx0 = width4 + (rx0 & 3);
        mv_x0 = rx0 % 4;
        rx0 /= 4;

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+1]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+2]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+3]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+1]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+2]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+3]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad = _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+4]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+5]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+6]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+7]+x1));

        // load ref_pic 0~3
        xmm_ref0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+4]+rx0));
        xmm_ref1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+5]+rx0));
        xmm_ref2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+6]+rx0));
        xmm_ref3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+7]+rx0));

        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_ref0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_ref1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_ref2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_ref3);

        // sum sad
        sad += _mm_extract_epi16(xmm_org0,0);
        sad += _mm_extract_epi16(xmm_org1,0);
        sad += _mm_extract_epi16(xmm_org2,0);
        sad += _mm_extract_epi16(xmm_org3,0);

        //_mm_empty();
        mcost += sad;
        tmp_only_motion_cost += sad;
      }
    }
    _mm_empty();

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }

  if (best_pos)
  {
    *mv_x += spiral_search_x [best_pos];
    *mv_y += spiral_search_y [best_pos];
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
  }
  *mv_x = max(min_mv[0], *mv_x);
  *mv_x = min(max_mv[0], *mv_x);
  *mv_y = max(min_mv[1], *mv_y);
  *mv_y = min(max_mv[1], *mv_y);

  //===== return minimum motion cost =====
  return min_mcost;
}


int_32_t c_avs_enc:: SubPelBlockMotionSearch_bid(pel_t** orig_pic, int_32_t ref, int_32_t pic_pix_x, int_32_t  pic_pix_y, int_32_t blocktype, int_32_t pred_mv_x, int_32_t pred_mv_y, int_32_t* mv_x, int_32_t* mv_y, int_32_t search_pos2, int_32_t search_pos4, int_32_t min_mcost, double lambda, int_32_t block_index)
{
  int_32_t   pos, best_pos, mcost, abort_search, tmp_only_motion_cost;
  int_32_t   y0, x0, y1, x1, ry0, rx0;
  int_32_t    ry0_bid, rx0_bid;
  int_32_t   cand_mv_x, cand_mv_y, mv_y0, mv_x0, mv_y1, mv_x1;
  int   incr = ref==-1 ? ((!img->fld_type)&&((byte***)mref==mref_fld)&&(img->type==B_IMG)) : ((byte***)mref==mref_fld)&&(img->type==B_IMG) ;
  pel_t ****ref_pic,****ref_pic_bid;
  byte  **ref_pic_tmp;
  int_32_t   lambda_factor   = LAMBDA_FACTOR (lambda);
  int_32_t   mv_shift        = 0;
  int_32_t   check_position0 = (blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard && !input->rdopt && img->type!=B_IMG && ref==0);
  int_32_t   blocksize_x     = input->blc_size[blocktype][0];
  int_32_t   blocksize_y     = input->blc_size[blocktype][1];
  int_32_t   pic4_pix_x      = (pic_pix_x << 2);
  int_32_t   pic4_pix_y      = (pic_pix_y << 2);
  int_32_t   max_pos_x4      = ((img->width -blocksize_x+1)<<2);
  int_32_t   max_pos_y4      = ((img->height-blocksize_y+1)<<2);
  int_32_t   min_pos2        = (input->hadamard ? 0 : 1);
  int_32_t   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
  int_32_t  apply_weights = 0;
  int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,refframe ,delta_PB;
  int_32_t satd;
  int_32_t tmp, cand_mv_x_tmp, cand_mv_y_tmp;
  int_32_t width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2)-32;
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2)-32;

  __m128i   xmm_org0,xmm_org1,xmm_org2,xmm_org3;
  __m128i   xmm_refa0,xmm_refa1,xmm_refa2,xmm_refa3;
  __m128i   xmm_refb0,xmm_refb1,xmm_refb2,xmm_refb3;
  int_32_t   max_mv[2], min_mv[2];
  int_32_t   max_bi_mv[2], min_bi_mv[2];
  int_32_t  target_x, target_y;
  max_mv[0] = (img->width  - pic_pix_x + 16 -1)  << 2;
  max_mv[1] = (img->height - pic_pix_y + 16 - 1) << 2;
  min_mv[0] = (-pic_pix_x-10) << 2;
  min_mv[1] = (-pic_pix_y-10) << 2;

  min_bi_mv[0] = (pic_pix_x-img->width) * 4;
  max_bi_mv[0] = pic_pix_x*4;

  min_bi_mv[1] = (pic_pix_y-img->height) * 4;
  max_bi_mv[1] = pic_pix_y*4;

  refframe = ref;
  delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
  delta_P = (delta_P + 512) % 512;
  if(img->picture_structure)
    TRp = (refframe+1)*delta_P;
  else
    TRp = delta_P;//ref == 0 ? delta_P-1 : delta_P+1;
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
  DistanceIndexBw    = (TRp - DistanceIndexFw+512)%512;
  if (!img->picture_structure)
  {
    incr = 2;
  }
  else
  {
    incr = 1;
  }

  ref_pic = mref [ref+incr];

  ref_pic_bid = mref [img->picture_structure ? 0 : ref/*2 - (ref+incr)*/];

  /*********************************
  *****                       *****
  *****  HALF-PEL REFINEMENT  *****
  *****                       *****
  *********************************/
  //===== convert search center to quarter-pel units =====
  *mv_x <<= 2;
  *mv_y <<= 2;
  //===== loop over search positions =====
  tmp = DistanceIndexBw*(512/DistanceIndexFw);
  best_pos = 0;
  for (pos = min_pos2; pos < max_pos2; pos++)
  {
    tmp_only_motion_cost = 0;
    cand_mv_x = *mv_x + (spiral_search_x[pos] << 1);    // quarter-pel units
    cand_mv_y = *mv_y + (spiral_search_y[pos] << 1);    // quarter-pel units
    cand_mv_x = max(min_mv[0], cand_mv_x);
    cand_mv_x = min(max_mv[0], cand_mv_x);
    cand_mv_y = max(min_mv[1], cand_mv_y);
    cand_mv_y = min(max_mv[1], cand_mv_y);
    target_x = (pic_pix_x*4+cand_mv_x)/4;
    target_y = (pic_pix_y*4+cand_mv_y)/4;
    cand_mv_x_tmp = cand_mv_x * tmp;
    cand_mv_x_tmp += 256;
    cand_mv_x_tmp >>= 9;
    cand_mv_y_tmp = cand_mv_y * tmp;
    cand_mv_y_tmp += 256;
    cand_mv_y_tmp >>= 9;
    target_x = (pic_pix_x*4+cand_mv_x_tmp)/4;
    target_y = (pic_pix_y*4+cand_mv_y_tmp)/4;

    //判断后向的mv，如果越界就跳过
    if (cand_mv_x_tmp > max_bi_mv[0] || cand_mv_x_tmp < max_bi_mv[0] || cand_mv_y_tmp > max_bi_mv[1] || cand_mv_y_tmp < max_bi_mv[1])
    {
      continue;
    }
    //----- set motion vector cost -----
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);
    if (check_position0 && pos==0)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }

    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=8)
    {
      y1 = pic_pix_y + y0;
      ry0 = ((y1 + IMG_PAD_SIZE) << 2) + cand_mv_y;
      ry0_bid = ((y1 + IMG_PAD_SIZE) << 2) - cand_mv_y_tmp;
      if (ry0 < 0)
        ry0 &= 3;
      else if (ry0 > height4)
        ry0 = height4 + (ry0 & 3);
      if (ry0_bid < 0)
        ry0_bid &= 3;
      else if (ry0_bid > height4)
        ry0_bid = height4 + (ry0_bid & 3);
      mv_y0 = ry0 % 4;
      mv_y1 = ry0_bid % 4;
      ry0 /= 4;
      ry0_bid /= 4;

      for (x0=0; x0<blocksize_x; x0+=8)
      {
        x1 = pic_pix_x + x0;
        rx0 = ((x1 + IMG_PAD_SIZE) << 2) + cand_mv_x;
        rx0_bid = ((x1 + IMG_PAD_SIZE) << 2) - cand_mv_x_tmp;
        if (rx0 < 0)
          rx0 &= 3;
        else if (rx0 > width4)
          rx0 = width4 + rx0 & 3;
        if (rx0_bid < 0)
          rx0_bid &= 3;
        else if (rx0_bid > width4)
          rx0_bid = width4 + rx0_bid & 3;
        mv_x0 = rx0 % 4;
        mv_x1 = rx0_bid % 4;
        rx0 /= 4;
        rx0_bid /= 4;

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+1]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+2]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+3]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pica 0~3
        xmm_refa0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0]+rx0));
        xmm_refa1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+1]+rx0));
        xmm_refa2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+2]+rx0));
        xmm_refa3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+3]+rx0));

        ref_pic_tmp = ref_pic_bid[mv_y1][mv_x1];
        // load ref_picb 0~3
        xmm_refb0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid ]+rx0_bid));
        xmm_refb1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +1]+rx0_bid));
        xmm_refb2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +2]+rx0_bid));
        xmm_refb3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +3]+rx0_bid));

        xmm_refa0 = _mm_avg_epu8(xmm_refa0,xmm_refb0);
        xmm_refa1 = _mm_avg_epu8(xmm_refa1,xmm_refb1);
        xmm_refa2 = _mm_avg_epu8(xmm_refa2,xmm_refb2);
        xmm_refa3 = _mm_avg_epu8(xmm_refa3,xmm_refb3);


        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_refa0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_refa1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_refa2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_refa3);

        // sum sad
        satd = _mm_extract_epi16(xmm_org0,0);
        satd += _mm_extract_epi16(xmm_org1,0);
        satd += _mm_extract_epi16(xmm_org2,0);
        satd += _mm_extract_epi16(xmm_org3,0);

        // load org_pic 4~7
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+4]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+5]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+6]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+7]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pica 4~7
        xmm_refa0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+4]+rx0));
        xmm_refa1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+5]+rx0));
        xmm_refa2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+6]+rx0));
        xmm_refa3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+7]+rx0));

        ref_pic_tmp = ref_pic_bid[mv_y1][mv_x1];
        // load ref_picb 4~7
        xmm_refb0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+4]+rx0_bid));
        xmm_refb1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+5]+rx0_bid));
        xmm_refb2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+6]+rx0_bid));
        xmm_refb3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+7]+rx0_bid));

        xmm_refa0 = _mm_avg_epu8(xmm_refa0,xmm_refb0);
        xmm_refa1 = _mm_avg_epu8(xmm_refa1,xmm_refb1);
        xmm_refa2 = _mm_avg_epu8(xmm_refa2,xmm_refb2);
        xmm_refa3 = _mm_avg_epu8(xmm_refa3,xmm_refb3);


        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_refa0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_refa1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_refa2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_refa3);

        // sum sad
        satd += _mm_extract_epi16(xmm_org0,0);
        satd += _mm_extract_epi16(xmm_org1,0);
        satd += _mm_extract_epi16(xmm_org2,0);
        satd += _mm_extract_epi16(xmm_org3,0);

        if (input->hadamard)
        {
          printf("\nHadmard Transform not supported yet!\n");
        }

        //if ((mcost += SATD (diff, input->hadamard)) > min_mcost)  //展开SATD
        tmp_only_motion_cost += satd;
        if ((mcost += satd) > min_mcost)
        {
          abort_search = 1;
          break;
        }
      }
    }

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    *mv_x += (spiral_search_x [best_pos] << 1);
    *mv_y += (spiral_search_y [best_pos] << 1);
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
  }


  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/
  //===== set function for getting pixel values =====
  //===== loop over search positions =====
  best_pos = 0;
  for (pos = 1; pos < search_pos4; pos++)
  {
    tmp_only_motion_cost = 0;
    cand_mv_x = *mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = *mv_y + spiral_search_y[pos];    // quarter-pel units
    cand_mv_x = max(min_mv[0], cand_mv_x);
    cand_mv_x = min(max_mv[0], cand_mv_x);
    cand_mv_y = max(min_mv[1], cand_mv_y);
    cand_mv_y = min(max_mv[1], cand_mv_y);

    cand_mv_x_tmp = cand_mv_x * tmp;
    cand_mv_x_tmp += 256;
    cand_mv_x_tmp >>= 9;
    cand_mv_y_tmp = cand_mv_y * tmp;
    cand_mv_y_tmp += 256;
    cand_mv_y_tmp >>= 9;
    if (cand_mv_x_tmp > max_bi_mv[0] || cand_mv_x_tmp < min_bi_mv[0] || cand_mv_y_tmp > max_bi_mv[1] || cand_mv_y_tmp < min_bi_mv[1])
    {
      continue;
    }
    //----- set motion vector cost -----
    mcost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);

    //----- add up SATD -----
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=8)
    {
      y1 = pic_pix_y + y0;
      ry0 = ((y1 + IMG_PAD_SIZE) << 2) + cand_mv_y;
      ry0_bid = ((y1 + IMG_PAD_SIZE) << 2) - cand_mv_y_tmp;
      if (ry0 < 0)
        ry0 &= 3;
      else if (ry0 > height4)
        ry0 = height4 + (ry0 & 3);
      if (ry0_bid < 0)
        ry0_bid &= 3;
      else if (ry0_bid > height4)
        ry0_bid = height4 + (ry0_bid & 3);
      mv_y0 = ry0 % 4;
      mv_y1 = ry0_bid % 4;
      ry0 /= 4;
      ry0_bid /= 4;

      for (x0=0; x0<blocksize_x; x0+=8)
      {
        x1 = pic_pix_x + x0;
        rx0 = ((x1 + IMG_PAD_SIZE) << 2) + cand_mv_x;
        rx0_bid = ((x1 + IMG_PAD_SIZE) << 2) - cand_mv_x_tmp;
        if (rx0 < 0)
          rx0 &= 3;
        else if (rx0 > width4)
          rx0 = width4 + rx0 & 3;
        if (rx0_bid < 0)
          rx0_bid &= 3;
        else if (rx0_bid > width4)
          rx0_bid = width4 + rx0_bid & 3;
        mv_x0 = rx0 % 4;
        mv_x1 = rx0_bid % 4;
        rx0 /= 4;
        rx0_bid /= 4;

        // load org_pic 0~3
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+1]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+2]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+3]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pica 0~3
        xmm_refa0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0]+rx0));
        xmm_refa1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+1]+rx0));
        xmm_refa2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+2]+rx0));
        xmm_refa3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+3]+rx0));

        ref_pic_tmp = ref_pic_bid[mv_y1][mv_x1];
        // load ref_picb 0~3
        xmm_refb0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid ]+rx0_bid));
        xmm_refb1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +1]+rx0_bid));
        xmm_refb2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +2]+rx0_bid));
        xmm_refb3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid +3]+rx0_bid));

        xmm_refa0 = _mm_avg_epu8(xmm_refa0,xmm_refb0);
        xmm_refa1 = _mm_avg_epu8(xmm_refa1,xmm_refb1);
        xmm_refa2 = _mm_avg_epu8(xmm_refa2,xmm_refb2);
        xmm_refa3 = _mm_avg_epu8(xmm_refa3,xmm_refb3);


        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_refa0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_refa1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_refa2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_refa3);

        // sum sad
        satd = _mm_extract_epi16(xmm_org0,0);
        satd += _mm_extract_epi16(xmm_org1,0);
        satd += _mm_extract_epi16(xmm_org2,0);
        satd += _mm_extract_epi16(xmm_org3,0);

        // load org_pic 4~7
        xmm_org0 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+4]+x1));
        xmm_org1 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+5]+x1));
        xmm_org2 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+6]+x1));
        xmm_org3 = _mm_loadl_epi64((__m128i*)(orig_pic[y1+7]+x1));

        ref_pic_tmp = ref_pic[mv_y0][mv_x0];
        // load ref_pica 4~7
        xmm_refa0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+4]+rx0));
        xmm_refa1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+5]+rx0));
        xmm_refa2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+6]+rx0));
        xmm_refa3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0+7]+rx0));

        ref_pic_tmp = ref_pic_bid[mv_y1][mv_x1];
        // load ref_picb 4~7
        xmm_refb0 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+4]+rx0_bid));
        xmm_refb1 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+5]+rx0_bid));
        xmm_refb2 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+6]+rx0_bid));
        xmm_refb3 = _mm_loadl_epi64((__m128i*)(ref_pic_tmp[ry0_bid+7]+rx0_bid));

        xmm_refa0 = _mm_avg_epu8(xmm_refa0,xmm_refb0);
        xmm_refa1 = _mm_avg_epu8(xmm_refa1,xmm_refb1);
        xmm_refa2 = _mm_avg_epu8(xmm_refa2,xmm_refb2);
        xmm_refa3 = _mm_avg_epu8(xmm_refa3,xmm_refb3);


        // sse sad
        xmm_org0 = _mm_sad_epu8(xmm_org0,xmm_refa0);
        xmm_org1 = _mm_sad_epu8(xmm_org1,xmm_refa1);
        xmm_org2 = _mm_sad_epu8(xmm_org2,xmm_refa2);
        xmm_org3 = _mm_sad_epu8(xmm_org3,xmm_refa3);

        // sum sad
        satd += _mm_extract_epi16(xmm_org0,0);
        satd += _mm_extract_epi16(xmm_org1,0);
        satd += _mm_extract_epi16(xmm_org2,0);
        satd += _mm_extract_epi16(xmm_org3,0);
        tmp_only_motion_cost += satd;
        if ((mcost += satd) > min_mcost)
        {
          abort_search = 1;
          break;
        }
      }
    }

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    *mv_x += spiral_search_x [best_pos];
    *mv_y += spiral_search_y [best_pos];
    if (tmp_only_motion_cost < only_motion_cost[blocktype][block_index])
    {
      only_motion_cost[blocktype][block_index] = tmp_only_motion_cost;
    }
  }
  *mv_x  = max(min_mv[0], *mv_x );
  *mv_x  = min(max_mv[0], *mv_x );
  *mv_y = max(min_mv[1], *mv_y);
  *mv_y = min(max_mv[1], *mv_y);

  //===== return minimum motion cost =====
  return min_mcost;
}


/*
*************************************************************************
* Function:Block motion search
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: BlockMotionSearch (int_32_t       ref,           // <--  reference frame (0... or -1 (backward))
                                        int_32_t       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                                        int_32_t       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                                        int_32_t       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                                        int_32_t       search_range,  // <--  1-d search range for integer-position search
                                        double         lambda,         // <--  Lagrangian parameter for determining motion cost
                                        int_32_t       block_index
                                        )
{
  pel_t   orig_val [256];
  pel_t  *orig_pic  [16] =
  {
    orig_val,     orig_val+ 16, orig_val+ 32, orig_val+ 48,
    orig_val+ 64, orig_val+ 80, orig_val+ 96, orig_val+112,
    orig_val+128, orig_val+144, orig_val+160, orig_val+176,
    orig_val+192, orig_val+208, orig_val+224, orig_val+240
  };

  int_32_t       pred_mv_x, pred_mv_y, mv_x, mv_y, i, j;
  int_32_t       max_value     = (1<<20);
  int_32_t       min_mcost     = max_value;
  int_32_t       refframe      = (ref==-1 ? 0 : ref);
  int_32_t*      pred_mv;
  int_32_t**     ref_array     = ((img->type!=B_IMG) ? refFrArr : ref>=0 ? fw_refFrArr : bw_refFrArr);
  int_32_t***    mv_array      = ((img->type!=B_IMG) ? tmp_mv   : ref>=0 ? tmp_fwMV    : tmp_bwMV);
  int_32_t*****  all_bmv       = img->all_bmv;
  int_32_t*****  all_mv        = (ref<0 ? img->all_bmv : img->all_mv);
  byte**    imgY_org_pic  = imgY_org;
  int_32_t       bsy           = input->blc_size[blocktype][1];
  int_32_t       bsx           = input->blc_size[blocktype][0];
  int_32_t       mb_pix_x      = pic_pix_x-img->pix_x;
  int_32_t       mb_pix_y      = pic_pix_y-img->pix_y;
  int_32_t       b8_x          = (mb_pix_x>>3);
  int_32_t       b8_y          = (mb_pix_y>>3);
  int_32_t       current_mb_nr = img->current_mb_nr;
  __m128i        xmm0;

  if (!img->picture_structure) // field coding
  {
    if (img->type==B_IMG)
    {
      refframe = ref<0 ? ref+2 : ref;
    }
  }

  pred_mv = ((img->type!=B_IMG) ? img->mv : ref>=0 ? img->p_fwMV : img->p_bwMV)[mb_pix_x>>3][mb_pix_y>>3][refframe][blocktype];

  for (j = 0; j < bsy; j++)
  {
    if(bsx>8)
    {
      xmm0 = _mm_loadu_si128((__m128i*)(imgY_org_pic[pic_pix_y+j]+pic_pix_x));
      _mm_storeu_si128((__m128i *)(orig_pic[j]),xmm0);
    }
    else
    {
      xmm0 = _mm_loadl_epi64((__m128i*)(imgY_org_pic[pic_pix_y+j]+pic_pix_x));
      _mm_storel_epi64((__m128i *)(orig_pic[j]),xmm0);
    }
  }

  //===========================================
  //=====   GET MOTION VECTOR PREDICTOR   =====
  //===========================================
  SetMotionVectorPredictor (pred_mv, ref_array, mv_array, refframe, mb_pix_x, mb_pix_y, bsx, bsy, ref);
  pred_mv_x = pred_mv[0];
  pred_mv_y = pred_mv[1];
  //==================================
  //=====   INTEGER-PEL SEARCH   =====
  //==================================
  //--- set search center ---
  mv_x = pred_mv_x>>2;
  mv_y = pred_mv_y>>2;

  if (!input->rdopt)
  {
    //--- adjust search center so that the (0,0)-vector is inside ---
    mv_x = max (-search_range, min (search_range, mv_x));
    mv_y = max (-search_range, min (search_range, mv_y));
  }

  //--- perform motion search ---
#ifdef _THREE_STEP_MOTION_SEARCH_
  min_mcost = TSSMotionSearch(orig_pic, ref, pic_pix_x, pic_pix_y, blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range, min_mcost, lambda, block_index);
#else
  min_mcost = FullPelBlockMotionSearch(orig_pic, ref, pic_pix_x, pic_pix_y, blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range, min_mcost, lambda, 0);
#endif
  if (ref == 0 && img->type == B_IMG)
    mcost_tmp = min_mcost;

  if (ref == 0)
  {
    mv_x_tmp = mv_x;
    mv_y_tmp = mv_y;
  }
  //==============================
  //=====   SUB-PEL SEARCH   =====
  //==============================
    min_mcost =  SubPelBlockMotionSearch (imgY_org_pic, ref, pic_pix_x, pic_pix_y, blocktype, pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9, min_mcost, lambda, block_index);

  if (!input->rdopt)
  {
    // Get the skip mode cost
    if (blocktype == 1 && img->type == INTER_IMG)
    {
      int_32_t cost;
      FindSkipModeMotionVector ();
      cost  = GetSkipCostMB (lambda);
      cost -= (int_32_t)floor(8*lambda+0.4999);
      if (cost < min_mcost)
      {
        min_mcost = cost;
        mv_x      = img->all_mv [0][0][0][0][0];
        mv_y      = img->all_mv [0][0][0][0][1];
      }
    }
  }

  //===============================================
  //=====   SET MV'S AND RETURN MOTION COST   =====
  //===============================================
  for (i=0; i < (bsx>>3); i++)
  {
    for (j=0; j < (bsy>>3); j++)
    {
      all_mv[b8_x+i][b8_y+j][refframe][blocktype][0] = mv_x;
      all_mv[b8_x+i][b8_y+j][refframe][blocktype][1] = mv_y;
    }
  }
  return min_mcost;
}

int_32_t c_avs_enc::GetSkipCostMB (double lambda)
{
  int_32_t block_y, block_x, pic_pix_y, pic_pix_x, x, y, x1, y1, mv[2];
  //int_16_t diff[16];
  int_32_t cost = 0, sad = 0;
  int_32_t pix_x = img->pix_x;
  int_32_t pix_y = img->pix_y;
  byte**    imgY_org_pic = imgY_org;
  byte  *d1, *d2, *d3, *d4;
  int_16_t *d5, *d6;
  int_32_t width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2)-32;
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2)-32;

  for (block_y=0; block_y<16; block_y+=8)
  {
    pic_pix_y = pix_y +block_y;
    for (block_x=0; block_x<16; block_x+=8)
    {
      pic_pix_x = pix_x + block_x;
      //根据mv计算参考帧的起始地址
      mv[0] = img->all_mv[block_x>>3][block_y>>3][0][0][0];
      mv[1] = img->all_mv[block_x>>3][block_y>>3][0][0][1];
      y  = ((pic_pix_y + IMG_PAD_SIZE) << 2) + mv[1];
      x  = ((pic_pix_x + IMG_PAD_SIZE) << 2) + mv[0];

      if (y < 0)
        y &= 3;
      else if (y > height4)
        y = height4 + (y & 3);
      if (x < 0)
        x &= 3;
      else if (x > width4)
        x = width4 + (x & 3);

      y1 = y % 4;
      x1 = x % 4;
      y /= 4;
      x /= 4;

      d1 = &imgY_org_pic[pic_pix_y][pic_pix_x];
      d2 = &imgY_org_pic[pic_pix_y + 1][pic_pix_x];
      d3 = &mref[0][y1][x1][y][x];
      d4 = &mref[0][y1][x1][y + 1][x];
      d5 = &img->mpr[block_y][block_x];
      d6 = &img->mpr[block_y + 1][block_x];

      __asm
      {
        mov      esi,  dword ptr [d1]  //read in orig
        movdqu    xmm0, xmmword ptr [esi]
        mov      esi,  dword ptr [d2]
        movdqu    xmm1, xmmword ptr [esi]

        mov      esi,  dword ptr [d3]  //read in ref_frame
        movdqu    xmm2, xmmword ptr [esi]
        mov      esi,  dword ptr [d4]
        movdqu    xmm3, xmmword ptr [esi]

        psadbw    xmm0, xmm2;    //sad
        psadbw        xmm1, xmm3;
        paddw         xmm0, xmm1;
        pextrw      eax,  xmm0, 0
          mov      sad, eax

          pxor      xmm7, xmm7        //byte -> int_16_t
          punpcklbw    xmm2, xmm7
          punpcklbw    xmm3, xmm7

          mov      esi,  dword ptr [d5]
        movdqa      xmmword ptr [esi],    xmm2
          mov      esi,  dword ptr [d6]
        movdqa      xmmword ptr [esi],    xmm3
      }
      cost += sad;

      d1 = &imgY_org_pic[pic_pix_y + 2][pic_pix_x];
      d2 = &imgY_org_pic[pic_pix_y + 3][pic_pix_x];
      d3 = &mref[0][y1][x1][y + 2][x];
      d4 = &mref[0][y1][x1][y + 3][x];
      d5 = &img->mpr[block_y + 2][block_x];
      d6 = &img->mpr[block_y + 3][block_x];

      __asm
      {
        mov      esi,  dword ptr [d1]  //read in orig
        movdqu    xmm0, xmmword ptr [esi]
        mov      esi,  dword ptr [d2]
        movdqu    xmm1, xmmword ptr [esi]

        mov      esi,  dword ptr [d3]  //read in ref_frame
        movdqu    xmm2, xmmword ptr [esi]
        mov      esi,  dword ptr [d4]
        movdqu    xmm3, xmmword ptr [esi]

        psadbw    xmm0, xmm2;    //sad
        psadbw        xmm1, xmm3;
        paddw         xmm0, xmm1;
        pextrw      eax,  xmm0, 0
          mov      sad, eax

          pxor      xmm7, xmm7        //byte -> int_16_t
          punpcklbw    xmm2, xmm7
          punpcklbw    xmm3, xmm7

          mov      esi,  dword ptr [d5]
        movdqa      xmmword ptr [esi],    xmm2
          mov      esi,  dword ptr [d6]
        movdqa      xmmword ptr [esi],    xmm3
      }
      cost += sad;

      d1 = &imgY_org_pic[pic_pix_y + 4][pic_pix_x];
      d2 = &imgY_org_pic[pic_pix_y + 5][pic_pix_x];
      d3 = &mref[0][y1][x1][y + 4][x];
      d4 = &mref[0][y1][x1][y + 5][x];
      d5 = &img->mpr[block_y + 4][block_x];
      d6 = &img->mpr[block_y + 5][block_x];

      __asm
      {
        mov      esi,  dword ptr [d1]  //read in orig
        movdqu    xmm0, xmmword ptr [esi]
        mov      esi,  dword ptr [d2]
        movdqu    xmm1, xmmword ptr [esi]

        mov      esi,  dword ptr [d3]  //read in ref_frame
        movdqu    xmm2, xmmword ptr [esi]
        mov      esi,  dword ptr [d4]
        movdqu    xmm3, xmmword ptr [esi]

        psadbw    xmm0, xmm2;    //sad
        psadbw        xmm1, xmm3;
        paddw         xmm0, xmm1;
        pextrw      eax,  xmm0, 0
          mov      sad, eax

          pxor      xmm7, xmm7        //byte -> int_16_t
          punpcklbw    xmm2, xmm7
          punpcklbw    xmm3, xmm7

          mov      esi,  dword ptr [d5]
        movdqa      xmmword ptr [esi],    xmm2
          mov      esi,  dword ptr [d6]
        movdqa      xmmword ptr [esi],    xmm3
      }
      cost += sad;

      d1 = &imgY_org_pic[pic_pix_y + 6][pic_pix_x];
      d2 = &imgY_org_pic[pic_pix_y + 7][pic_pix_x];
      d3 = &mref[0][y1][x1][y + 6][x];
      d4 = &mref[0][y1][x1][y + 7][x];
      d5 = &img->mpr[block_y + 6][block_x];
      d6 = &img->mpr[block_y + 7][block_x];
      __asm
      {
        mov      esi,  dword ptr [d1]  //read in orig
        movdqu    xmm0, xmmword ptr [esi]
        mov      esi,  dword ptr [d2]
        movdqu    xmm1, xmmword ptr [esi]

        mov      esi,  dword ptr [d3]  //read in ref_frame
        movdqu    xmm2, xmmword ptr [esi]
        mov      esi,  dword ptr [d4]
        movdqu    xmm3, xmmword ptr [esi]

        psadbw    xmm0, xmm2;    //sad
        psadbw        xmm1, xmm3;
        paddw         xmm0, xmm1;
        pextrw      eax, xmm0, 0
          mov      sad, eax

          pxor      xmm7, xmm7        //byte -> int_16_t
          punpcklbw    xmm2, xmm7
          punpcklbw    xmm3, xmm7

          mov      esi,  dword ptr [d5]
        movdqa      xmmword ptr [esi],    xmm2
          mov      esi,  dword ptr [d6]
        movdqa      xmmword ptr [esi],    xmm3
      }
      cost += sad;
    }
  }

  return cost;
}
/*
*************************************************************************
* Function:Find motion vector for the Skip mode
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::FindSkipModeMotionVector()
{
  int_32_t bx, by;
  int_32_t mb_nr = img->current_mb_nr;
  int_32_t mb_width = img->width/16;
  int_32_t mb_available_up   = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
  int_32_t mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);
  int_32_t zeroMotionAbove   = !mb_available_up  ? 1 : refFrArr[(img->block_y>>1)-1][(img->block_x>>1)] == 0 && tmp_mv[0][(img->block_y>>1)-1][4+(img->block_x>>1)]   == 0 && tmp_mv[1][(img->block_y>>1)-1][4+(img->block_x>>1)]   == 0 ? 1 : 0;
  int_32_t zeroMotionLeft    = !mb_available_left? 1 : refFrArr[(img->block_y>>1)][(img->block_x>>1)-1] == 0 && tmp_mv[0][(img->block_y>>1)  ][4+(img->block_x>>1)-1] == 0 && tmp_mv[1][(img->block_y>>1)  ][4+(img->block_x>>1)-1] == 0 ? 1 : 0;
  int_32_t mb_x = img->mb_x;
  int_32_t mb_y = img->mb_y;
  int_32_t block_x = img->block_x;
  int_32_t block_y = img->block_y;
  int_32_t **refar = refFrArr;
  int_32_t ***tmpmv = tmp_mv;
  int_32_t *****all_mv = img->all_mv;
  int_32_t *****mv  = img->mv;

  if (zeroMotionAbove || zeroMotionLeft)
  {
    for (by = 0;by < 2;by++)
    {
      for (bx = 0;bx < 2;bx++)
      {
        all_mv [bx][by][0][0][0] = 0;
        all_mv [bx][by][0][0][1] = 0;
      }
    }
  }
  else
  {
    for (by = 0;by < 2;by++)
    {
      for (bx = 0;bx < 2;bx++)
      {
        all_mv [bx][by][0][0][0] = mv[0][0][0][1][0];
        all_mv [bx][by][0][0][1] = mv[0][0][0][1][1];
      }
    }
  }
}

/*
*************************************************************************
* Function:Get cost for direct mode for an 8x8 block
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::Get_Direct_Cost8x8 (int_32_t block, double lambda)
{
  int_32_t cost  = 0, sad = 0;
  int_32_t mb_y  = (block/2)<<3;
  int_32_t mb_x  = (block%2)<<3;
  byte **imgY_original = imgY_org;
  int_32_t   y1, x1, y, x, y1_bw, x1_bw, y_bw, x_bw, mv[2];
  int_32_t   pix_x = img->pix_x;
  int_32_t   pix_y = img->pix_y;
  byte  *d1, *d2, *d3, *d4;
  int_16_t *d5, *d6;
  int_32_t width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2)-32;
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2)-32;


  mv[0] = img->all_mv[block % 2][block / 2][0][0][0];
  mv[1] = img->all_mv[block % 2][block / 2][0][0][1];
  y  = ((pix_y + mb_y + IMG_PAD_SIZE) << 2) + mv[1];
  x  = ((pix_x + mb_x + IMG_PAD_SIZE) << 2) + mv[0];

  if (y < 0)
    y &= 3;
  else if (y > height4)
    y = height4 + (y & 3);
  if (x < 0)
    x &= 3;
  else if (x > width4)
    x = width4 + (x & 3);



  y1 = y % 4;
  x1 = x % 4;
  y /= 4;
  x /= 4;

  //后向，根据mv计算参考帧的起始地址
  mv[0] = img->all_bmv[block % 2][block / 2][0][0][0];
  mv[1] = img->all_bmv[block % 2][block / 2][0][0][1];

  y_bw  = ((pix_y + mb_y + IMG_PAD_SIZE) << 2) + mv[1];
  x_bw  = ((pix_x + mb_x + IMG_PAD_SIZE) << 2) + mv[0];

  if (y_bw < 0)
    y_bw &= 3;
  else if (y_bw > height4)
    y_bw = height4 + (y_bw & 3);
  if (x_bw < 0)
    x_bw &= 3;
  else if (x_bw > width4)
    x_bw = width4 + (x_bw & 3);


  y1_bw = y_bw % 4;
  x1_bw = x_bw % 4;
  y_bw /= 4;
  x_bw /= 4;

  d3 = &mref[1][y1][x1][y][x];
  d4 = &mref[1][y1][x1][y + 1][x];

  __asm
  {
    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm2, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm3, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
  }

  d1 = &imgY_original[pix_y + mb_y][pix_x + mb_x];
  d2 = &imgY_original[pix_y + mb_y + 1][pix_x + mb_x];
  d3 = &mref[0][y1_bw][x1_bw][y_bw][x_bw];
  d4 = &mref[0][y1_bw][x1_bw][y_bw + 1][x_bw];
  d5 = &img->mpr[mb_y][mb_x];
  d6 = &img->mpr[mb_y + 1][mb_x];

  __asm
  {
    mov      esi,  dword ptr [d1]  //read in ref_frame
    movdqu    xmm0, xmmword ptr [esi]
    mov      esi,  dword ptr [d2]
    movdqu    xmm1, xmmword ptr [esi]

    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm4, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm5, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm4, xmm7
      punpcklbw  xmm5, xmm7

      pavgw       xmm2, xmm4        //(fw + bw)/2
      pavgw       xmm3, xmm5

      psubw       xmm0, xmm2        //org - pred
      psubw       xmm1, xmm3

      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm0
      pxor    xmm0, xmm7
      psrlw    xmm7, 15
      paddw       xmm0, xmm7
      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm1
      pxor    xmm1, xmm7
      psrlw    xmm7, 15
      paddw       xmm1, xmm7
      paddw    xmm0, xmm1
      pxor    xmm1, xmm1
      psadbw    xmm0, xmm1
      movdqa    xmm1, xmm0
      punpckhqdq  xmm0, xmm0
      paddw       xmm0, xmm1
      pextrw      eax,  xmm0, 0
      mov      sad,  eax

      mov      esi,  dword ptr [d5]
    movdqa      xmmword ptr [esi],    xmm2
      mov      esi,  dword ptr [d6]
    movdqa      xmmword ptr [esi],    xmm3
  }
  cost += sad;

  d3 = &mref[1][y1][x1][y + 2][x];
  d4 = &mref[1][y1][x1][y + 3][x];

  __asm
  {
    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm2, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm3, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
  }

  d1 = &imgY_original[pix_y + mb_y + 2][pix_x + mb_x];
  d2 = &imgY_original[pix_y + mb_y + 3][pix_x + mb_x];
  d3 = &mref[0][y1_bw][x1_bw][y_bw + 2][x_bw];
  d4 = &mref[0][y1_bw][x1_bw][y_bw + 3][x_bw];
  d5 = &img->mpr[mb_y + 2][mb_x];
  d6 = &img->mpr[mb_y + 3][mb_x];

  __asm
  {
    mov      esi,  dword ptr [d1]  //read in ref_frame
    movdqu    xmm0, xmmword ptr [esi]
    mov      esi,  dword ptr [d2]
    movdqu    xmm1, xmmword ptr [esi]

    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm4, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm5, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm4, xmm7
      punpcklbw  xmm5, xmm7

      pavgw       xmm2, xmm4        //(fw + bw)/2
      pavgw       xmm3, xmm5

      psubw       xmm0, xmm2        //org - pred
      psubw       xmm1, xmm3

      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm0
      pxor    xmm0, xmm7
      psrlw    xmm7, 15
      paddw       xmm0, xmm7
      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm1
      pxor    xmm1, xmm7
      psrlw    xmm7, 15
      paddw       xmm1, xmm7
      paddw    xmm0, xmm1
      pxor    xmm1, xmm1
      psadbw    xmm0, xmm1
      movdqa    xmm1, xmm0
      punpckhqdq  xmm0, xmm0
      paddw       xmm0, xmm1
      pextrw      eax,  xmm0, 0
      mov      sad,  eax

      mov      esi,  dword ptr [d5]
    movdqa      xmmword ptr [esi],    xmm2
      mov      esi,  dword ptr [d6]
    movdqa      xmmword ptr [esi],    xmm3
  }
  cost += sad;

  d3 = &mref[1][y1][x1][y + 4][x];
  d4 = &mref[1][y1][x1][y + 5][x];

  __asm
  {
    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm2, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm3, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
  }

  d1 = &imgY_original[pix_y + mb_y + 4][pix_x + mb_x];
  d2 = &imgY_original[pix_y + mb_y + 5][pix_x + mb_x];
  d3 = &mref[0][y1_bw][x1_bw][y_bw + 4][x_bw];
  d4 = &mref[0][y1_bw][x1_bw][y_bw + 5][x_bw];
  d5 = &img->mpr[mb_y + 4][mb_x];
  d6 = &img->mpr[mb_y + 5][mb_x];

  __asm
  {
    mov      esi,  dword ptr [d1]  //read in ref_frame
    movdqu    xmm0, xmmword ptr [esi]
    mov      esi,  dword ptr [d2]
    movdqu    xmm1, xmmword ptr [esi]

    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm4, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm5, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm4, xmm7
      punpcklbw  xmm5, xmm7

      pavgw       xmm2, xmm4        //(fw + bw)/2
      pavgw       xmm3, xmm5

      psubw       xmm0, xmm2        //org - pred
      psubw       xmm1, xmm3

      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm0
      pxor    xmm0, xmm7
      psrlw    xmm7, 15
      paddw       xmm0, xmm7
      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm1
      pxor    xmm1, xmm7
      psrlw    xmm7, 15
      paddw       xmm1, xmm7
      paddw    xmm0, xmm1
      pxor    xmm1, xmm1
      psadbw    xmm0, xmm1
      movdqa    xmm1, xmm0
      punpckhqdq  xmm0, xmm0
      paddw       xmm0, xmm1
      pextrw      eax,  xmm0, 0
      mov      sad,  eax

      mov      esi,  dword ptr [d5]
    movdqa      xmmword ptr [esi],    xmm2
      mov      esi,  dword ptr [d6]
    movdqa      xmmword ptr [esi],    xmm3
  }
  cost += sad;

  d3 = &mref[1][y1][x1][y + 6][x];
  d4 = &mref[1][y1][x1][y + 7][x];

  __asm
  {
    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm2, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm3, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm2, xmm7
      punpcklbw  xmm3, xmm7
  }

  d1 = &imgY_original[pix_y + mb_y + 6][pix_x + mb_x];
  d2 = &imgY_original[pix_y + mb_y + 7][pix_x + mb_x];
  d3 = &mref[0][y1_bw][x1_bw][y_bw + 6][x_bw];
  d4 = &mref[0][y1_bw][x1_bw][y_bw + 7][x_bw];
  d5 = &img->mpr[mb_y + 6][mb_x];
  d6 = &img->mpr[mb_y + 7][mb_x];

  __asm
  {
    mov      esi,  dword ptr [d1]  //read in ref_frame
    movdqu    xmm0, xmmword ptr [esi]
    mov      esi,  dword ptr [d2]
    movdqu    xmm1, xmmword ptr [esi]

    mov      esi,  dword ptr [d3]  //read in ref_frame
    movdqu    xmm4, xmmword ptr [esi]
    mov      esi,  dword ptr [d4]
    movdqu    xmm5, xmmword ptr [esi]

    pxor    xmm7, xmm7        //byte -> int_16_t
      punpcklbw  xmm0, xmm7
      punpcklbw  xmm1, xmm7
      punpcklbw  xmm4, xmm7
      punpcklbw  xmm5, xmm7

      pavgw       xmm2, xmm4        //(fw + bw)/2
      pavgw       xmm3, xmm5

      psubw       xmm0, xmm2        //org - pred
      psubw       xmm1, xmm3

      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm0
      pxor    xmm0, xmm7
      psrlw    xmm7, 15
      paddw       xmm0, xmm7
      pxor    xmm7, xmm7
      pcmpgtw    xmm7, xmm1
      pxor    xmm1, xmm7
      psrlw    xmm7, 15
      paddw       xmm1, xmm7
      paddw    xmm0, xmm1
      pxor    xmm1, xmm1
      psadbw    xmm0, xmm1
      movdqa    xmm1, xmm0
      punpckhqdq  xmm0, xmm0
      paddw       xmm0, xmm1
      pextrw      eax,  xmm0, 0
      mov      sad,  eax

      mov      esi,  dword ptr [d5]
    movdqa      xmmword ptr [esi],    xmm2
      mov      esi,  dword ptr [d6]
    movdqa      xmmword ptr [esi],    xmm3
  }
  cost += sad;


  return cost;
}

/*
*************************************************************************
* Function:Get cost for direct mode for an macroblock
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::Get_Direct_CostMB (double lambda)
{
  int_32_t i;
  int_32_t cost = 0;

  //b

  for (i=0; i<4; i++)
  {
    cost += Get_Direct_Cost8x8 (i, lambda);
  }

  return cost;
}

/*
*************************************************************************
* Function:Motion search for a partition
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::PartitionMotionSearch (int_32_t blocktype, int_32_t block8x8, double lambda)
{
  TLS static int_32_t  bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,1,0,0}, {0,1,0,1}};
  TLS static int_32_t  by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,1,0,0}, {0,0,0,0}, {0,0,1,1}};
  int_32_t   **ref_array, ***mv_array, *****all_mv;
  int_32_t   ref, refinx, refframe, v, h, mcost, search_range, i, j;
  int_32_t   pic_block_x, pic_block_y;
  int_32_t   bframe   = (img->type==B_IMG);
  int_32_t   max_ref  =  bframe ?  1 : img->nb_references;
  int_32_t   min_ref  = (bframe ? -1 : 0);

  int_32_t   step_h   = (input->blc_size[blocktype][0]>>3);
  int_32_t   step_v   = (input->blc_size[blocktype][1]>>3);
  int_32_t   block_x  = img->block_x;
  int_32_t   block_y  = img->block_y;
  static int_32_t count = 0;
  count++;
  if (max_ref > img->buf_cycle)
  {
    max_ref = img->buf_cycle;
  }
  //set the reference frame
  //now only consider the simple case: O means open-loop, C means close-loop or close-loop-me or cascade-transcoding
  //I P P P P
  //I O C O C
  //===== LOOP OVER REFERENCE FRAMES =====

  if(img->type==B_IMG)
  {
    max_ref = 1;
    //if (img->number % input->intra_period == 0)  //for close GOP
    //  min_ref = 0;
  }
  else if (img->number % input->intra_period == 1)//for close GOP
  {
    max_ref = 1;
  }

  //ref=-1的时候是后向，ref=0的时候是前向
  for (ref=min_ref; ref<max_ref; ref++)
  {
    refinx    = ref+1;
    refframe  = (ref<0 ? 0 : ref);
    //----- set search range ---
    search_range = input->search_range;
    //----- set arrays -----
    ref_array = ((img->type!=B_IMG) ? refFrArr : ref<0 ? bw_refFrArr : fw_refFrArr);
    mv_array  = ((img->type!=B_IMG) ?   tmp_mv : ref<0 ?    tmp_bwMV :    tmp_fwMV);
    all_mv    = (ref<0 ? img->all_bmv : img->all_mv);

    //----- init motion cost -----
    motion_cost[blocktype][refinx][block8x8] = 0;

    v=by0[blocktype][block8x8];
    pic_block_y = (block_y>>1) + v;

    h=bx0[blocktype][block8x8];
    pic_block_x = (block_x>>1) + h;

    //--- motion search for block ---
#ifdef FastME
    mcost = FME_BlockMotionSearch (ref, 8*pic_block_x, 8*pic_block_y, blocktype, search_range, lambda);
#else
    mcost = BlockMotionSearch (ref, 8*pic_block_x, 8*pic_block_y, blocktype, search_range, lambda, block8x8);
#endif
    motion_cost[blocktype][refinx][block8x8] += mcost;

    //--- set motion vectors and reference frame (for motion vector prediction) ---
    for (j=0; j<step_v; j++)
    {
      for (i=0; i<step_h; i++)
      {
        mv_array[0][pic_block_y+j][pic_block_x+i+4] = all_mv[h][v][refframe][blocktype][0];
        mv_array[1][pic_block_y+j][pic_block_x+i+4] = all_mv[h][v][refframe][blocktype][1];
        ref_array  [pic_block_y+j][pic_block_x+i  ] = refframe;
      }
    }
  }
}


/*
*************************************************************************
* Function:Motion search for a partition
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::PartitionMotionSearch_bid (int_32_t blocktype, int_32_t block8x8, double lambda)
{
  TLS static int_32_t  bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,1,0,0}, {0,1,0,1}};
  TLS static int_32_t  by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,1,0,0}, {0,0,0,0}, {0,0,1,1}};
  int_32_t   **ref_array, ***mv_array;
  int_32_t   ref, refinx, refframe, mcost, search_range, i, j;
  int_32_t   pic_block_x, pic_block_y;
  int_32_t   max_ref   = img->nb_references-1;
  int_32_t   parttype  = (blocktype<4?blocktype:4);
  int_32_t   step_h0   = (input->blc_size[ parttype][0]>>3);
  int_32_t   step_v0   = (input->blc_size[ parttype][1]>>3);
  int_32_t       pred_mv_x, pred_mv_y, mv_x, mv_y;
  int_32_t       mb_x;
  int_32_t       mb_y;
  int_32_t       bsx;
  int_32_t       bsy;
  int_32_t       frameref;
  int_32_t*      pred_mv;
  int_32_t       pic_pix_x;
  int_32_t       pic_pix_y;
  mb_y = by0[parttype][block8x8];
  mb_x = bx0[parttype][block8x8];
  max_ref = 1;

  if (max_ref > img->buf_cycle)
    max_ref = img->buf_cycle;
  //===== LOOP OVER REFERENCE FRAMES =====
  for (ref=0; ref<max_ref; ref++)
  {
    refinx    = ref+1;
    refframe  = (ref<0?0:ref);
    if ((byte**)mref[0] == mref_fld[0])
    {
      refinx    = ref+2;
      refframe  = (ref<0?ref+2:ref);
    }
    search_range = input->search_range;

    //----- set arrays -----
    ref_array = ref<0 ? bw_refFrArr : fw_refFrArr;
    mv_array  = ref<0 ? tmp_bwMV : tmp_fwMV;
    //----- init motion cost -----
    motion_cost_bid[blocktype][refinx][block8x8] = 0;
    //===== LOOP OVER BLOCKS =====
    pic_block_y = img->block8_y + by0[parttype][block8x8];
    pic_block_x = img->block8_x + bx0[parttype][block8x8];
    pic_pix_x = pic_block_x << 3;
    pic_pix_y = pic_block_y << 3;
    //--- motion search for block ---
#ifdef FastME
    mcost = FME_BlockMotionSearch_bid (ref, 8*pic_block_x, 8*pic_block_y, blocktype, search_range, lambda);//modify?
#else
    //---展开BlockMotionSearch_bid
    bsx       = input->blc_size[blocktype][0];
    bsy       = input->blc_size[blocktype][1];
    frameref  = (ref==-1 ? 0 : ref);
    mcost = 1<<20;
    if (!img->picture_structure) // field coding
    {
      frameref = ref<0 ? ref+2 : ref;
    }
    pred_mv = img->omv[mb_x][mb_y][frameref][blocktype];
    //===========================================
    //=====   GET MOTION VECTOR PREDICTOR   =====
    //===========================================
    SetMotionVectorPredictor (pred_mv, fw_refFrArr, tmp_fwMV, frameref, 8*mb_x, 8*mb_y, bsx, bsy, ref);
    pred_mv_x = pred_mv[0];
    pred_mv_y = pred_mv[1];
    //==================================
    //=====   INTEGER-PEL SEARCH   =====
    //==================================

    //--- set search center ---
    mv_x = pred_mv_x / 4;
    mv_y = pred_mv_y / 4;
    if (!input->rdopt)
    {
      //--- adjust search center so that the (0,0)-vector is inside ---
      if (mv_x > search_range)
      {
        mv_x = search_range;
      }
      else if (mv_x < search_range)
      {
        mv_x = -search_range;
      }
      if (mv_y > search_range)
      {
        mv_y = search_range;
      }
      else if (mv_y < -search_range)
      {
        mv_y = -search_range;
      }
    }
    //--- perform motion search ---
    mcost = mcost_tmp;
    mv_x = mv_x_tmp;
    mv_y = mv_y_tmp;
    mcost = SubPelBlockMotionSearch_bid (imgY_org, ref, pic_pix_x, pic_pix_y, blocktype,pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,mcost, lambda,block8x8);

    //===============================================
    //=====   SET MV'S AND RETURN MOTION COST   =====
    //===============================================
    for (i=0; i < step_h0; i++)
    {
      for (j=0; j < step_v0; j++)
      {
        img->all_omv[mb_x+i][mb_y+j][frameref][blocktype][0] = mv_x;
        img->all_omv[mb_x+i][mb_y+j][frameref][blocktype][1] = mv_y;
      }
    }
    //--展开BlockMotionSearch_bid结束--
    //计算后向mv
    {
      int_32_t delta_P,TRp,DistanceIndexFw,DistanceIndexBw,refframe,delta_PB;
      int_32_t mv[2];
      refframe = 0;
      delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
      delta_P = (delta_P + 512) % 512;
      TRp = (refframe+1)*delta_P;  //the latest backward reference
      TRp = (TRp+512) % 512;
      delta_PB = 2*(img->tr - img->imgtr_last_P_frm);
      delta_PB = (delta_PB + 512)%512;
      DistanceIndexFw = delta_PB;
      DistanceIndexBw = TRp - DistanceIndexFw;
      mv[0] = - ((img->all_omv[mb_x][mb_y][0][blocktype][0]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
      mv[1] = - ((img->all_omv[mb_x][mb_y][0][blocktype][1]*DistanceIndexBw*(256/DistanceIndexFw)+128)>>8);
      for (j=0; j<step_v0; j++)
      {
        for (i=0; i<step_h0; i++)
        {
          img->all_bw_omv[mb_x+i][mb_y+j][0][blocktype][0] = mv[0];
          img->all_bw_omv[mb_x+i][mb_y+j][0][blocktype][1] = mv[1];
        }
      }
    }
#endif
    motion_cost_bid[blocktype][refinx][block8x8] += mcost;
    //--- set motion vectors and reference frame (for motion vector prediction) ---
    for (j=0; j<step_v0; j++)
    {
      for (i=0; i<step_h0; i++)
      {
        mv_array[0][pic_block_y+j][pic_block_x+i+4] = img->all_omv[mb_y+j][mb_x+i][refframe][blocktype][0];
        mv_array[1][pic_block_y+j][pic_block_x+i+4] = img->all_omv[mb_y+j][mb_x+i][refframe][blocktype][1];
        ref_array  [pic_block_y+j][pic_block_x+i  ] = refframe;
      }
    }
  }
}


extern int_32_t* last_P_no;
/*********************************************
*****                                   *****
*****  Calculate Direct Motion Vectors  *****
*****                                   *****
*********************************************/

/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::Get_IP_direct()
{
  int_32_t  block_x, block_y, pic_block_x, pic_block_y;
  int_32_t  refframe, bw_ref, refP_tr, TRb, TRp, TRd, TRp1;
  int_32_t  frame_no_next_P, frame_no_B, delta_P, delta_P_scale;
  int_32_t  pix_y            = img->pix_y;
  int_32_t  **refarr         = refFrArr;
  int_32_t  ***tmpmvs        = tmp_mv;
  int_32_t  *****all_mvs     = img->all_mv;
  int_32_t  *****all_bmvs    = img->all_bmv;
  int_32_t  prev_mb_is_field = 0;
  int_32_t  mv_scale, scale_refframe;

  int_32_t  **fwrefarr         = fw_refFrArr;
  int_32_t  **bwrefarr         = bw_refFrArr;
  int_32_t  ***tmpmvfw         = tmp_fwMV;
  int_32_t  ***tmpmvbw         = tmp_bwMV;

  pic_block_x = img->block_x;
  pic_block_y = img->block_y;

  for (block_y=0; block_y<2; block_y++)
  {
    pic_block_y = (pix_y>>3) + block_y;

    for (block_x=0; block_x<2; block_x++)
    {
      pic_block_x = (img->pix_x>>3) + block_x;

      if((refframe=refarr[pic_block_y][pic_block_x]) == -1)
      {
        all_mvs [block_x][block_y][0][0][0] = 0;
        all_mvs[block_x][block_y][0][0][1]  = 0;
        all_bmvs[block_x][block_y][0][0][0] = 0;
        all_bmvs[block_x][block_y][0][0][1] = 0;
        SetMotionVectorPredictor(all_mvs [block_x][block_y][0][0],fwrefarr,tmpmvfw,0,0,0,16,16, 0);
        SetMotionVectorPredictor(all_bmvs [block_x][block_y][img->picture_structure?0:1][0],bwrefarr,tmpmvbw,img->picture_structure?0:1,0,0,16,16, -1);
      }
      else
      {
        refP_tr = nextP_tr - ((refframe+1)*img->p_interval);
        refP_tr = (refP_tr+256)%256;
        frame_no_next_P = 2*img->imgtr_next_P_frm;
        frame_no_B = 2*img->tr;
        //frame_no_B = 2*picture_distance;
        delta_P = 2*(img->imgtr_next_P_frm - img->imgtr_last_P_frm);
        delta_P = (delta_P + 512) % 512;
        delta_P_scale = 2*(img->imgtr_next_P_frm - img->imgtr_last_prev_P_frm);  // 20071009
        delta_P_scale = (delta_P_scale + 512)%512;
        if(!img->picture_structure)
        {
          if (img->current_mb_nr_fld < img->total_number_mb) //top field
            scale_refframe =   refframe == 0  ? 0 : 1;
          else
            scale_refframe =   refframe == 1  ? 1 : 2;
        }
        else
          scale_refframe = 0;

        if(!img->picture_structure)
        {
          if (img->current_mb_nr_fld < img->total_number_mb) //top field
          {
            //TRp = delta_P*(refframe/2+1)-(refframe+1)%2;  //the lates backward reference
            //TRp1 = delta_P*(scale_refframe/2+1)-(scale_refframe+1)%2;  //the lates backward reference
            TRp = refframe<2 ? delta_P-(refframe+1)%2 : delta_P_scale-(refframe+1)%2;
            TRp1 = scale_refframe<2 ? delta_P-(scale_refframe+1)%2 : delta_P_scale-(scale_refframe+1)%2;
            bw_ref = 1;
          }
          else
          {
            //TRp = 1 + delta_P*((refframe+1)/2)-refframe%2;
            //TRp1 = 1 + delta_P*((scale_refframe+1)/2)-scale_refframe%2;
            TRp  = refframe==0 ? 1 : refframe<3 ? 1 + delta_P - refframe%2 : 1 + delta_P_scale - refframe%2;
            TRp1 = scale_refframe==0 ? 1 : scale_refframe<3 ? 1 + delta_P - scale_refframe%2 : 1 + delta_P_scale - scale_refframe%2;
            bw_ref = 0;
          }
        }
        else
        {
          //TRp  = (refframe+1)*delta_P;
          //TRp1  = (scale_refframe+1)*delta_P;
          TRp  = refframe==0 ? delta_P : delta_P_scale;
          TRp1 = scale_refframe==0 ? delta_P : delta_P_scale;
        }
        TRd = frame_no_next_P - frame_no_B;
        TRb = TRp1 - TRd;

        TRp  = (TRp + 512)%512;
        TRp1 = (TRp1 + 512)%512;
        TRd  = (TRd + 512)%512;
        TRb  = (TRb + 512)%512;
        mv_scale = (TRb * 256) / TRp;       //! Note that this could be precomputed at the frame/slice level.

        refframe = 0;

        if(!img->picture_structure)
        {
          if (img->current_mb_nr_fld >= img->total_number_mb) //top field
            scale_refframe --;
          refframe = scale_refframe;
        }
        else
        {
          refframe = 0;
          bw_ref = 0;
        }

        if(tmpmvs[0][pic_block_y][pic_block_x+4] < 0)
        {
          all_mvs [block_x][block_y][refframe][0][0] = -(((16384/TRp)*(1-TRb*tmpmvs[0][pic_block_y][pic_block_x+4])-1)>>14);
          all_bmvs [block_x][block_y][bw_ref][0][0] =  ((16384/TRp)*(1-TRd*tmpmvs[0][pic_block_y][pic_block_x+4])-1)>>14;
        }
        else
        {
          all_mvs [block_x][block_y][refframe][0][0] = ((16384/TRp)*(1+TRb*tmpmvs[0][pic_block_y][pic_block_x+4])-1)>>14;
          all_bmvs [block_x][block_y][bw_ref][0][0] =  -(((16384/TRp)*(1+TRd*tmpmvs[0][pic_block_y][pic_block_x+4])-1)>>14);
        }

        if(tmpmvs[1][pic_block_y][pic_block_x+4] < 0)
        {
          all_mvs [block_x][block_y][refframe][0][1] = -(((16384/TRp)*(1-TRb*tmpmvs[1][pic_block_y][pic_block_x+4])-1)>>14);
          all_bmvs [block_x][block_y][bw_ref][0][1] =    ((16384/TRp)*(1-TRd*tmpmvs[1][pic_block_y][pic_block_x+4])-1)>>14;
        }
        else
        {
          all_mvs [block_x][block_y][refframe][0][1] = ((16384/TRp)*(1+TRb*tmpmvs[1][pic_block_y][pic_block_x+4])-1)>>14;
          all_bmvs [block_x][block_y][bw_ref][0][1] =  -(((16384/TRp)*(1+TRd*tmpmvs[1][pic_block_y][pic_block_x+4])-1)>>14);
        }
      }
    }
  }
}

/*
*************************************************************************
* Function:control the sign of a with b
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::sign(int_32_t a,int_32_t b)
{
  int_32_t x;
  x=absm(a);

  if (b >= 0)
    return x;
  else
    return -x;
}

