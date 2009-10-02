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
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <memory.h>
#include <assert.h>

#include "global.h"
#include <xmmintrin.h>
#include <emmintrin.h>

#include "windows.h"

#define TO_SAVE 4711
#define FROM_SAVE 4712
#define Clip(min,max,val) (((val)<(min))?(min):(((val)>(max))?(max):(val)))



TLS static byte clip0c[16]={0,128,0,128,0,128,0,128,0,128,0,128,0,128,0,128};
TLS static byte clip255c[16]={0,127,0,127,0,127,0,127,0,127,0,127,0,127,0,127};
TLS static byte round1c[16]={1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
TLS static byte round4c[16]={4,0,4,0,4,0,4,0,4,0,4,0,4,0,4,0};
TLS static byte round8c[16]={8,0,8,0,8,0,8,0,8,0,8,0,8,0,8,0};
TLS static byte round16c[16]={16,0,16,0,16,0,16,0,16,0,16,0,16,0,16,0};
TLS static byte round32c[16]={32,0,32,0,32,0,32,0,32,0,32,0,32,0,32,0};
TLS static byte round64c[16]={64,0,64,0,64,0,64,0,64,0,64,0,64,0,64,0};
TLS static byte round512c[16]={0,2,0,0,0,2,0,0,0,2,0,0,0,2,0,0};

#define  IClip( Min, Max, Val) (((Val)<(Min))? (Min):(((Val)>(Max))? (Max):(Val)))

/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::picture_header()
{
  int len;

  img->cod_counter = 0;

  if(img->type==INTRA_IMG)
    len = IPictureHeader(img->number);
  else
    len = PBPictureHeader();

  // Rate control
  img->NumberofHeaderBits +=len;
  if(img->BasicUnit<img->total_number_mb)
    img->NumberofBasicUnitHeaderBits +=len;

  // Update statistics
  stat->bit_slice += len;
  stat->bit_use_header[img->type] += len;

  WriteFrameFieldMBInHeader = 0;
}

int c_avs_enc::encode_one_frame ()
{
  time_t ltime1;
  time_t ltime2;
#ifdef ROI_ENABLE
  int_32_t erroCode;
#endif
#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif
  int_32_t tmp_time;
  int_32_t  bits_frm = 0, bits_fld = 0;
  float dis_frm = 0, dis_frm_y = 0, dis_frm_u = 0, dis_frm_v = 0;
  float dis_fld = 0, dis_fld_y = 0, dis_fld_u = 0, dis_fld_v = 0;
  int_32_t bits = 0;

#ifdef WIN32
  _ftime (&tstruct1);           // start time ms
#else
  ftime (&tstruct1);
#endif
  time (&ltime1);               // start time s

  init_frame();           // initial frame variables
  CalculateFrameNumber();
  ReadOneFrame ();
  CopyFrameToOldImgOrgVariables();
#ifdef ROI_ENABLE
  //detect roi
  erroCode=detect_roi(YCbCr, w, h, ROIArray);
  if (erroCode)
  {
    printf("error in detecting roi function\n");
  }
#endif
  //Rate control
  img->FieldControl=0;
  if(input->RCEnable == 1)
  {
    /*update the number of MBs in the basic unit for MB adaptive coding*/
    img->BasicUnit = input->basicunit;
    rc_init_pict(1,0,1);
    img->qp  = updateQuantizationParameter(0);
  }
  else if ( input->RCEnable == 3 )
  {
	  img->qp = (img->type == B_IMG) ? i_QPb : i_QPip;
  }
  if (input->InterlaceCodingOption == FRAME_CODING)  // !! frame coding or paff coding
  {
    put_buffer_frame ();     //initialize frame buffer
    //frame picture
    img->progressive_frame = 1;
    img->picture_structure = 1;
    img->TopFieldFlag      = 0;
    if (img->type == B_IMG)
      Bframe_ctr++;         // Bframe_ctr only used for statistics, should go to stat->

    // !!  calculate the weighting parameter
    if(img->type != INTRA_IMG && input->picture_weighting_flag == 1)
    {
      estimate_weighting_factor();
    }
    else
    {
      img->LumVarFlag = 0 ;
    }
    code_a_picture (frame_pic);
    if (img->type!=B_IMG)
    {
      Update_Picture_Buffers();
    }
    stat->bit_ctr_emulationprevention += stat->em_prev_bits_frm;
    if (img->type != B_IMG)       //all I- and P-frames
    {
#ifdef _FAST_INTERPOLATION_
      UnifiedOneForthPix_sse (imgY);
#else
      UnifiedOneForthPix_c_sse(imgY);
#endif
    }

    time (&ltime2);               // end time sec
#ifdef WIN32
    _ftime (&tstruct2);           // end time ms
#else
    ftime (&tstruct2);            // end time ms
#endif

    tmp_time = (int)((ltime2 * 1000 + tstruct2.millitm) - (ltime1 * 1000 + tstruct1.millitm));
    tot_time = tot_time + tmp_time;
  }

  writeout_picture ();
  FreeBitstream();
  find_snr ();

  GBIM_value_frm = 0;//find_GBIM(imgY);
  GBIM_value=0;
  GBIM_value += GBIM_value_frm;

#ifdef _OUTPUT_RECON_IMG_
  // Write reconstructed images
  write_reconstructed_image ();
#endif

  //Rate control
  if(input->RCEnable == 1)
  {
    bits = stat->bit_ctr-stat->bit_ctr_n;
    rc_update_pict_frame(bits);
  }
  goprate += (stat->bit_ctr-stat->bit_ctr_n);

#ifdef _DEBUG
  if (img->number == 0)
    ReportFirstframe(tmp_time);
  else
  {
    //Rate control
    if(input->RCEnable == 1)
    {
      if(input->InterlaceCodingOption == FRAME_CODING)
        bits = stat->bit_ctr-stat->bit_ctr_n;
      else
      {
        bits = stat->bit_ctr -Pprev_bits; // used for rate control update */
        Pprev_bits = stat->bit_ctr;
      }
    }

    switch (img->type)
    {
    case INTRA_IMG:
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportIntra(tmp_time);
      break;
    case B_IMG:
      stat->bit_ctr_B += stat->bit_ctr - stat->bit_ctr_n;
      ReportB(tmp_time);
      break;
    default:      // P, P_MULTPRED?
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportP(tmp_time);
    }
  }
#endif

  stat->bit_ctr_n = stat->bit_ctr;
  //Rate control
  if(input->RCEnable == 1)
  {
    rc_update_pict(bits);
    /*update the parameters of quadratic R-D model*/
    if (img->type == INTER_IMG)
    {
      if (input->InterlaceCodingOption == FRAME_CODING)
      {
        updateRCModel();
      }
      else if (img->IFLAG == 0)
      {
        updateRCModel();
      }
    }
  }
  return 1;
}



int c_avs_enc:: writeout_picture()
{
  assert (currBitStream->bits_to_go == 8);
  WriteBitstreamtoFile();
  return 0;
}


void c_avs_enc::code_a_picture (Picture *frame)
{
  stat->em_prev_bits_frm = 0;
  stat->em_prev_bits = &stat->em_prev_bits_frm;

  AllocateBitstream();
  picture_header();

  picture_data();

  frame->bits_per_picture = 8 * (currBitStream->byte_pos);
  if (input->InterlaceCodingOption != FRAME_CODING)
  {
    find_distortion ();
    frame->distortion_y = snr->snr_y;
    frame->distortion_u = snr->snr_u;
    frame->distortion_v = snr->snr_v;
  }
}


/*
*************************************************************************
* Function:Initializes the parameters for a new frame
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::init_frame ()
{
  int i, j, k;
  int prevP_no, nextP_no;
  img->top_bot = -1;
  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  img->coded_mb_nr = 0;
  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0;
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;

  refFrArr    = refFrArr_frm;
  fw_refFrArr = fw_refFrArr_frm;
  bw_refFrArr = bw_refFrArr_frm;

  if (img->type != B_IMG)
  {
    if((gframe_no%input->GopLength)>=input->GopLength - input->successive_Bframe)
    {
      img->tr = gframe_no;
    }
    else if(img->type==INTRA_IMG)
    {
      img->tr = gframe_no;
    }
    else
    {
      img->tr = gframe_no + input->successive_Bframe;
    }
    if (img->imgtr_last_P_frm < 0)
      img->imgtr_last_P_frm = 0;
    img->imgtr_last_prev_P_frm = img->imgtr_last_P_frm;
    img->imgtr_last_P_frm = img->imgtr_next_P_frm;
    img->imgtr_next_P_frm = picture_distance;

    if (img->number != 0 && input->successive_Bframe != 0)   // B pictures to encode
      nextP_tr_frm = picture_distance;

    if(!input->RCEnable)
    {
      if (img->type == INTRA_IMG)
        img->qp = input->qp0;   // set quant. parameter for I-frame
      else
      {
        img->qp = input->qpN;
      }
    }
  }
  else
  {
    img->p_interval = input->successive_Bframe + 1;
    prevP_no = (img->number - 1) * img->p_interval;
    nextP_no = (img->number) * img->p_interval;

    img->b_interval = (int) ((float) (input->successive_Bframe + 1) / (input->successive_Bframe + 1.0) + 0.49999);

    img->tr = gframe_no-1;

    if (img->tr >= nextP_no)
      img->tr = nextP_no - 1;

    if(!input->RCEnable)
      img->qp = input->qpB;

    // initialize arrays

    if(!img->picture_structure) //field coding
    {
      for (k = 0; k < 2; k++)
      {
        for (i = 0; i < img->height / BLOCK_SIZE; i++)
        {
          for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
          {
            tmp_fwMV[k][i][j] = 0;
            tmp_bwMV[k][i][j] = 0;
            dfMV[k][i][j] = 0;
            dbMV[k][i][j] = 0;
          }
        }
      }

      for (i = 0; i < img->height / B8_SIZE; i++)
      {
        for (j = 0; j < img->width / BLOCK_SIZE; j++)
        {
          fw_refFrArr[i][j] = bw_refFrArr[i][j] = -1;
        }
      }
    }
    else
    {
      for (k = 0; k < 2; k++)
      {
        for (i = 0; i < img->height / BLOCK_SIZE; i++)
        {
          for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
          {
            tmp_fwMV[k][i][j] = 0;
            tmp_bwMV[k][i][j] = 0;
            dfMV[k][i][j] = 0;
            dbMV[k][i][j] = 0;
          }
        }
      }

      for (i = 0; i < img->height / BLOCK_SIZE; i++)
      {
        for (j = 0; j < img->width / BLOCK_SIZE; j++)
        {
          fw_refFrArr[i][j] = bw_refFrArr[i][j] = -1;
        }
      }
    }
  }
  stat->bit_slice = 0;
  //for rm52j
  picture_distance = img->tr;
}


void c_avs_enc::init_field ()
{

  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  stat->bit_slice = 0;
  img->coded_mb_nr = 0;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0;
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;
  img->mb_no_currSliceLastMB = ( input->slice_row_nr != 0 )
    ? min(input->slice_row_nr * img->img_width_in_mb - 1, img->img_width_in_mb * img->img_height_in_mb - 1)
    : img->img_width_in_mb * img->img_height_in_mb - 1 ;
}

/*
*************************************************************************
* Function:Writes reconstructed image(s) to file
This can be done more elegant!
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::write_reconstructed_image ()
{
  int i, j, k;
  if (p_rec != NULL)
  {
    if (img->type != B_IMG)
    {
      if (input->successive_Bframe == 0)
      {
        for (i = 0; i < img->height; i++)
          for (j = 0; j < img->width; j++)
          {
            fputc (imgY[i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
            fputc(imgY_org[i][j], p_org_dec);
#endif
          }

          for (k = 0; k < 2; ++k)
            for (i = 0; i < img->height / 2; i++)
              for (j = 0; j < img->width / 2; j++)
              {
                fputc (imgUV[k][i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
                fputc(imgUV_org[k][i][j], p_org_dec);
#endif
              }
      }
      else if ((gframe_no%input->GopLength) == 0)
      {
        for (i = 0; i < img->height; i++)
          for (j = 0; j < img->width; j++) 
          {
            fputc (imgY[i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
            fputc(imgY_org[i][j], p_org_dec);
#endif
          }
          for (k = 0; k < 2; ++k)
            for (i = 0; i < img->height / 2; i++)
              for (j = 0; j < img->width / 2; j++)
              {
                fputc (imgUV[k][i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
                fputc(imgUV_org[k][i][j], p_org_dec);
#endif
              }
      }

      // next P picture. This is saved with recon B picture after B picture coding
      if (img->number != 0 && input->successive_Bframe != 0)
      {
        if(((gframe_no%input->GopLength)>=input->GopLength - input->successive_Bframe))
        {
          for (i = 0; i < img->height; i++)
            for (j = 0; j < img->width; j++)
            {
              fputc (imgY[i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
              fputc(imgY_org[i][j], p_org_dec);
#endif
            }

            for (k = 0; k < 2; ++k)
              for (i = 0; i < img->height / 2; i++)
                for (j = 0; j < img->width / 2; j++)
                {
                  fputc (imgUV[k][i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
                  fputc(imgUV_org[k][i][j], p_org_dec);
#endif
                }
        }
        else
        {
          for (i = 0; i < img->height; i++)
            for (j = 0; j < img->width; j++)
            {
              nextP_imgY[i][j] = imgY[i][j];
#ifdef _OUTPUT_DEC_IMG_
              org_nextP_imgY[i][j] = imgY_org[i][j];
#endif
            }

            for (k = 0; k < 2; ++k)
              for (i = 0; i < img->height / 2; i++)
                for (j = 0; j < img->width / 2; j++)
                {
                  nextP_imgUV[k][i][j] = imgUV[k][i][j];
#ifdef _OUTPUT_DEC_IMG_
                  org_nextP_imgUV[k][i][j] = imgUV_org[k][i][j];
#endif
                }
        }
      }
    }
    else
    {
      for (i = 0; i < img->height; i++)
        for (j = 0; j < img->width; j++)
        {
          fputc (imgY[i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
          fputc(imgY_org[i][j], p_org_dec);
#endif
        }

        for (k = 0; k < 2; ++k)
          for (i = 0; i < img->height / 2; i++)
            for (j = 0; j < img->width / 2; j++)
            {
              fputc (imgUV[k][i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
              fputc(imgUV_org[k][i][j], p_org_dec);
#endif
            }

            // If this is last B frame also store P frame
            if (img->b_frame_to_code == input->successive_Bframe)
            {
              // save P picture
              for (i = 0; i < img->height; i++)
                for (j = 0; j < img->width; j++)
                {
                  fputc (nextP_imgY[i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
                  fputc(org_nextP_imgY[i][j], p_org_dec);
#endif
                }

                for (k = 0; k < 2; ++k)
                  for (i = 0; i < img->height / 2; i++)
                    for (j = 0; j < img->width / 2; j++)
                    {
                      fputc (nextP_imgUV[k][i][j], p_rec);
#ifdef _OUTPUT_DEC_IMG_
                      fputc(org_nextP_imgUV[k][i][j], p_org_dec);
#endif
                    }
            }
    }
  }
  fflush(p_rec);
}

/*
*************************************************************************
* Function:Choose interpolation method depending on MV-resolution
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

__inline void c_avs_enc::avs_const_initialize()
{
  clip0    = _mm_loadu_si128((const __m128i *)clip0c);
  clip255  = _mm_loadu_si128((const __m128i *)clip255c);
  round1   = _mm_loadu_si128((const __m128i *)round1c);
  round4   = _mm_loadu_si128((const __m128i *)round4c);
  round8   = _mm_loadu_si128((const __m128i *)round8c);
  round16  = _mm_loadu_si128((const __m128i *)round16c);
  round32  = _mm_loadu_si128((const __m128i *)round32c);
  round64  = _mm_loadu_si128((const __m128i *)round64c);
  round512 = _mm_loadu_si128((const __m128i *)round512c);
}

__inline __m128i c_avs_enc::avs_combine_w2b(__m128i xmm0, __m128i xmm1)
{  
  xmm0 = _mm_packus_epi16(xmm0,xmm1);

  return xmm0;
}

__inline __m128i c_avs_enc::avs_combine_d2w(__m128i xmm0, __m128i xmm1)
{  
  xmm0 = _mm_packus_epi16(xmm0,xmm1);
  xmm1 = _mm_xor_si128(xmm1,xmm1);
  xmm0 = _mm_packus_epi16(xmm0,xmm1);
  return xmm0;
}

__inline __m128i c_avs_enc::avs_filter_halfpel_w(__m128i xmm0, __m128i xmm1, __m128i xmm2)
{
  xmm1 = _mm_add_epi16(xmm1,xmm2);
  xmm2 = _mm_slli_epi16(xmm1,2);
  xmm1 = _mm_add_epi16(xmm1,xmm2);
  
  xmm2 = _mm_sub_epi16(xmm1,xmm0);
  return xmm2;
}

__inline __m128i c_avs_enc::avs_filter_quaterpel_w(__m128i xmm0, __m128i xmm1, __m128i xmm2)
{
  xmm1 = _mm_add_epi16(xmm1,xmm2);
  xmm2 = _mm_slli_epi16(xmm1,3);
  xmm1 = _mm_sub_epi16(xmm2,xmm1);

  //xmm1 = _mm_add_epi16(xmm1,round1);
  xmm1 = _mm_srai_epi16(xmm1,1);         // bit width control

  //xmm0 = _mm_add_epi16(xmm0,round1);
  xmm0 = _mm_srai_epi16(xmm0,1);         // bit width control

  xmm2 = _mm_add_epi16(xmm1,xmm0);

  return xmm2;
}

__inline __m128i c_avs_enc::avs_filter_quaterpel_d(__m128i xmm0, __m128i xmm1, __m128i xmm2)
{
  xmm1 = _mm_add_epi32(xmm1,xmm2);
  xmm2 = _mm_slli_epi32(xmm1,3);
  xmm1 = _mm_sub_epi32(xmm2,xmm1);

  xmm2 = _mm_add_epi32(xmm1,xmm0);

  return xmm2;
}

__inline __m128i c_avs_enc::avs_clip_0_255_w(__m128i xmm0)
{
  xmm0 = _mm_adds_epi16(xmm0,clip255);
  xmm0 = _mm_subs_epi16(xmm0,clip255);
  xmm0 = _mm_adds_epi16(xmm0,clip0);
  xmm0 = _mm_subs_epi16(xmm0,clip0);
  return xmm0;
}

__inline __m128i c_avs_enc::avs_zero(__m128i xmm0)
{
  return _mm_xor_si128(xmm0,xmm0);
}


/*
*************************************************************************
* Function:Quarter Pel Intelpolation for SSE Optimization
* Output:
* Return:
* Attention:
*************************************************************************
*/
void  c_avs_enc::UnifiedOneForthPix_sse (pel_t ** imgY)
{
  int img_pad_width,img_pad_height;
  int xx,yy;
  double temp=0;

  __int16 tmpw0;
  __m128i xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9;
  __m128i qtmp0,qtmp1;

  /////////////////////////////将所有可对齐载入的数据对齐载入！！！！！！！！！
  img_pad_width  = (img->width  + (IMG_PAD_SIZE<<1));
  img_pad_height = (img->height + (IMG_PAD_SIZE<<1));

  interpolation = mref[0];

  // initialize
  avs_const_initialize();
  xmm0=_mm_setzero_si128();
  /*************************************************************
  //           Basic Quater Pel Interpolation Unit
  //
  //              ○ ==> Interger Pel Position
  //              □ ==> Half     Pel Position
  //              △ ==> Quarter  Pel Position
  //************************************************************
  //                   ○   △   □   △
  //
  //                   △   △   △   △
  //
  //                   □   △   □   △
  //
  //                   △   △   △   △
  *************************************************************/

  // Pre Stuffing for Marginal Conditions
  //        * * * *
  //        o * * o
  //        o * * o
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    // o o o o
    // x x x x
    // x x x x
    // x x x x
    tmpw0 = imgY[yy-IMG_PAD_SIZE][0];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)tmp02[yy] , xmm0);
    _mm_store_si128((__m128i*)(tmp02[yy]+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);

    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i *)interpolation[0][0][yy] , xmm0);
    _mm_store_si128((__m128i *)interpolation[0][1][yy] , xmm0);
    _mm_store_si128((__m128i *)interpolation[0][2][yy] , xmm0);
    _mm_store_si128((__m128i *)interpolation[0][3][yy] , xmm0);



    tmpw0 = imgY[yy-IMG_PAD_SIZE][img->width-1];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i *)(tmp02[yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i *)(tmp02[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i *)(interpolation[0][0][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i *)(interpolation[0][1][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i *)(interpolation[0][2][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i *)(interpolation[0][3][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
  }

  //        * o o *
  //        * * * *
  //        * * * *
  //        * o o *
  for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
  {
    // o x x x
    // o x x x
    // o x x x
    // o x x x
    xmm0 = _mm_loadu_si128((const __m128i*)(imgY[0]+xx-IMG_PAD_SIZE));
    xmm2 = _mm_unpackhi_epi8(xmm0,xmm0);
    xmm3 = _mm_unpacklo_epi8(xmm0,xmm0);
    xmm2 = _mm_srli_epi16(xmm2,8);
    xmm3 = _mm_srli_epi16(xmm3,8);
    xmm2 = _mm_slli_epi16(xmm2,3);
    xmm3 = _mm_slli_epi16(xmm3,3);

    xmm1 = _mm_loadu_si128((const __m128i*)(imgY[img->height-1]+xx-IMG_PAD_SIZE));
    xmm4 = _mm_unpackhi_epi8(xmm1,xmm1);
    xmm5 = _mm_unpacklo_epi8(xmm1,xmm1);
    xmm4 = _mm_srli_epi16(xmm4,8);
    xmm5 = _mm_srli_epi16(xmm5,8);
    xmm4 = _mm_slli_epi16(xmm4,3);
    xmm5 = _mm_slli_epi16(xmm5,3);
    for(yy=0;yy<IMG_PAD_SIZE;yy++)
    {
      _mm_store_si128((__m128i *)(interpolation[0][0][yy]+xx), xmm0);
      _mm_store_si128((__m128i *)(interpolation[1][0][yy]+xx), xmm0);
      _mm_store_si128((__m128i *)(interpolation[2][0][yy]+xx), xmm0);
      _mm_store_si128((__m128i *)(interpolation[3][0][yy]+xx), xmm0);

      _mm_store_si128((__m128i *)(tmp20[yy]+xx), xmm2);
      _mm_store_si128((__m128i *)(tmp20[yy]+xx+8), xmm3);

      _mm_store_si128((__m128i *)(interpolation[0][0][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i *)(interpolation[1][0][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i *)(interpolation[2][0][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i *)(interpolation[3][0][img_pad_height-yy-1]+xx), xmm1);

      _mm_store_si128((__m128i *)(tmp20[img_pad_height-yy-1]+xx), xmm4);
      _mm_store_si128((__m128i *)(tmp20[img_pad_height-yy-1]+xx+8), xmm5);
    }
  }

  //        o * * o
  //        * * * *
  //        * * * *
  //        o * * o
  for(yy=0; yy<IMG_PAD_SIZE; yy++)
  {
    // ALL 16 Points Need to Be Stuffed at Corner Positions
    // o o o o
    // o o o o
    // o o o o
    // o o o o
    tmpw0 = imgY[0][0];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i *)(tmp02[yy]),   xmm0);
    _mm_store_si128((__m128i *)(tmp02[yy]+8), xmm0);
    _mm_store_si128((__m128i *)(tmp20[yy]),   xmm0);
    _mm_store_si128((__m128i *)(tmp20[yy]+8), xmm0);
    xmm0  = _mm_slli_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp22[yy]),   xmm0);
    _mm_store_si128((__m128i*)(tmp22[yy]+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(interpolation[0][0][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][1][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][2][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][3][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][0][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][0][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][0][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][yy]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][yy]) , xmm0);






    tmpw0 = imgY[0][img->width-1];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(tmp02[yy]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp02[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
    _mm_store_si128((__m128i*)(tmp20[yy]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp20[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
    xmm0  = _mm_slli_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp22[yy]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp22[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(interpolation[0][0][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][1][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][2][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][3][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][0][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][0][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][0][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);






    tmpw0 = imgY[img->height-1][0];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]),   xmm0);
    _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]+8), xmm0);
    _mm_store_si128((__m128i*)(tmp20[img_pad_height-yy-1]),   xmm0);
    _mm_store_si128((__m128i*)(tmp20[img_pad_height-yy-1]+8), xmm0);
    xmm0  = _mm_slli_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]),   xmm0);
    _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(interpolation[0][0][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][1][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][2][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][3][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][0][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][0][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][0][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][img_pad_height-yy-1]) , xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][img_pad_height-yy-1]) , xmm0);





    tmpw0 = imgY[img->height-1][img->width-1];
    xmm0  = _mm_insert_epi16(xmm0,(tmpw0<<3),0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
    _mm_store_si128((__m128i*)(tmp20[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp20[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
    xmm0  = _mm_slli_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE),   xmm0);
    _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE+8), xmm0);

    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(interpolation[0][0][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][1][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][2][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[0][3][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][0][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][0][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][0][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][img_pad_height-yy-1]+img_pad_width-IMG_PAD_SIZE), xmm0);
  }

  // o x x x
  // x x x x
  // x x x x
  // x x x x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      xmm0 = _mm_loadu_si128((const __m128i*)(imgY[yy-IMG_PAD_SIZE]+xx-IMG_PAD_SIZE));
      _mm_store_si128((__m128i*)(interpolation[0][0][yy]+xx) , xmm0);

    }
  }

  // x x o x
  // x x x x
  // x x x x
  // x x x x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      // loading
      xmm0 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx-1));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);

      xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm5 = _mm_unpackhi_epi8(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi8(xmm1,xmm1);
      xmm5 = _mm_srli_epi16(xmm5,8);
      xmm1 = _mm_srli_epi16(xmm1,8);

      xmm2 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx+1));
      xmm6 = _mm_unpackhi_epi8(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi8(xmm2,xmm2);
      xmm6 = _mm_srli_epi16(xmm6,8);
      xmm2 = _mm_srli_epi16(xmm2,8);

      xmm3 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx+2));
      xmm7 = _mm_unpackhi_epi8(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi8(xmm3,xmm3);
      xmm7 = _mm_srli_epi16(xmm7,8);
      xmm3 = _mm_srli_epi16(xmm3,8);

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_halfpel_w(xmm0,xmm1,xmm2);

      _mm_store_si128((__m128i *)(tmp02[yy]+xx), xmm2);

      // rounding & cliping
      xmm2 = _mm_add_epi16(xmm2,round4);
      xmm2 = _mm_srai_epi16(xmm2,3);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_halfpel_w(xmm4,xmm5,xmm6);

      _mm_store_si128((__m128i*)(tmp02[yy]+xx+8), xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round4);
      xmm6 = _mm_srai_epi16(xmm6,3);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[0][2][yy]+xx), xmm2);
    }

  }
  // Stuffing Vertical Marginal Points for Half Pels
  //        * o o *
  //        * * * *
  //        * * * *
  //        * o o *
  for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
  {
    // x x o x
    // x x o x
    // x x o x
    // x x o x
    xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][2][IMG_PAD_SIZE]+xx));
    xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][2][img_pad_height-IMG_PAD_SIZE-1]+xx));

    xmm2 = _mm_load_si128((const __m128i*)(tmp02[IMG_PAD_SIZE]+xx));
    xmm3 = _mm_load_si128((const __m128i*)(tmp02[img_pad_height-IMG_PAD_SIZE-1]+xx));
    xmm4 = _mm_load_si128((const __m128i*)(tmp02[IMG_PAD_SIZE]+xx+8));
    xmm5 = _mm_load_si128((const __m128i*)(tmp02[img_pad_height-IMG_PAD_SIZE-1]+xx+8));
    xmm6 = _mm_slli_epi16(xmm3,3);
    xmm7 = _mm_slli_epi16(xmm5,3);
    xmm8 = _mm_slli_epi16(xmm2,3);
    xmm9 = _mm_slli_epi16(xmm4,3);
    for(yy=0; yy<IMG_PAD_SIZE; yy++)
    {
      _mm_store_si128((__m128i*)(interpolation[0][2][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[1][2][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[2][2][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[3][2][yy]+xx), xmm0);

      _mm_store_si128((__m128i*)(interpolation[0][2][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[1][2][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[2][2][img_pad_height-yy-1]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[3][2][img_pad_height-yy-1]+xx), xmm1);

      _mm_store_si128((__m128i*)(tmp02[yy]+xx), xmm2);
      _mm_store_si128((__m128i*)(tmp02[yy]+xx+8), xmm4);
      _mm_store_si128((__m128i*)(tmp22[yy]+xx), xmm8);
      _mm_store_si128((__m128i*)(tmp22[yy]+xx+8), xmm9);

      _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]+xx), xmm3);
      _mm_store_si128((__m128i*)(tmp02[img_pad_height-yy-1]+xx+8), xmm5);
      _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]+xx), xmm6);
      _mm_store_si128((__m128i*)(tmp22[img_pad_height-yy-1]+xx+8), xmm7);
    }
  }

  // x x x x
  // x x x x
  // o x o x
  // x x x x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      // x x x x
      // x x x x
      // o x x x
      // x x x x


      // loading
      xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy-1]+xx));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);

      xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm5 = _mm_unpackhi_epi8(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi8(xmm1,xmm1);
      xmm5 = _mm_srli_epi16(xmm5,8);
      xmm1 = _mm_srli_epi16(xmm1,8);

      xmm2 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy+1]+xx));
      xmm6 = _mm_unpackhi_epi8(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi8(xmm2,xmm2);
      xmm6 = _mm_srli_epi16(xmm6,8);
      xmm2 = _mm_srli_epi16(xmm2,8);

      xmm3 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy+2]+xx));
      xmm7 = _mm_unpackhi_epi8(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi8(xmm3,xmm3);
      xmm7 = _mm_srli_epi16(xmm7,8);
      xmm3 = _mm_srli_epi16(xmm3,8);

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_halfpel_w(xmm0,xmm1,xmm2);

      _mm_store_si128((__m128i*)(tmp20[yy]+xx), xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round4);
      xmm2 = _mm_srai_epi16(xmm2,3);
      xmm2 = avs_clip_0_255_w(xmm2);


      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_halfpel_w(xmm4,xmm5,xmm6);

      _mm_store_si128((__m128i*)(tmp20[yy]+xx+8), xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round4);
      xmm6 = _mm_srai_epi16(xmm6,3);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i *)(interpolation[2][0][yy]+xx), xmm2);

      // x x x x
      // x x x x
      // x x o x
      // x x x x
      // loading

      xmm0 = _mm_load_si128((const __m128i*)(tmp02[yy-1]+xx));
      xmm1 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx));
      xmm2 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx));
      xmm3 = _mm_load_si128((const __m128i*)(tmp02[yy+2]+xx));
      xmm4 = _mm_load_si128((const __m128i*)(tmp02[yy-1]+xx+8));
      xmm5 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx+8));
      xmm6 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx+8));
      xmm7 = _mm_load_si128((const __m128i*)(tmp02[yy+2]+xx+8));

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_halfpel_w(xmm0,xmm1,xmm2);

      _mm_store_si128((__m128i*)(tmp22[yy]+xx), xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round32);
      xmm2 = _mm_srai_epi16(xmm2,6);
      xmm2 = avs_clip_0_255_w(xmm2);


      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_halfpel_w(xmm4,xmm5,xmm6);

      _mm_store_si128((__m128i*)(tmp22[yy]+xx+8), xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round32);
      xmm6 = _mm_srai_epi16(xmm6,6);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[2][2][yy]+xx), xmm2);
    }
  }

  // Stuffing Horizontal Marginal Points for Half Pels
  //        * * * *
  //        o * * o
  //        o * * o
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    // x x x x
    // x x x x
    // o o o o
    // x x x x
    tmpw0 = interpolation[2][0][yy][IMG_PAD_SIZE];
    tmpw0 = tmpw0 + (tmpw0<<8);
    xmm0  = _mm_insert_epi16   (xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0, 0);
    xmm0  = _mm_shuffle_epi32  (xmm0, 0);

    _mm_store_si128((__m128i*)(interpolation[2][0][yy]), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][yy]), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][yy]), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][yy]), xmm0);

    tmpw0 = tmp20[yy][IMG_PAD_SIZE];
    xmm0  = _mm_insert_epi16   (xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32  (xmm0,0);
    _mm_store_si128((__m128i*)(tmp20[yy]),  xmm0);
    _mm_store_si128((__m128i*)(tmp20[yy]+8), xmm0);
    xmm0 = _mm_slli_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp22[yy]),  xmm0);
    _mm_store_si128((__m128i*)(tmp22[yy]+8), xmm0);

    tmpw0 = interpolation[2][0][yy][img_pad_width-IMG_PAD_SIZE-1];
    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)(interpolation[2][0][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][1][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][2][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[2][3][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);

    tmpw0 = tmp22[yy][img_pad_width-IMG_PAD_SIZE-1];
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);
    _mm_store_si128((__m128i*)(tmp22[yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(tmp22[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
    xmm0 = _mm_srai_epi16(xmm0,3);
    _mm_store_si128((__m128i*)(tmp20[yy]+img_pad_width-IMG_PAD_SIZE) , xmm0);
    _mm_store_si128((__m128i*)(tmp20[yy]+img_pad_width-IMG_PAD_SIZE+8), xmm0);
  }


  // x o x o
  // x x x x
  // x o x o
  // x x x x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      // x o x x
      // x x x x
      // x x x x
      // x x x x
      // loading

      xmm0 = _mm_loadu_si128((const __m128i*)(tmp02[yy]+xx-1));

      xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm5 = _mm_unpackhi_epi8(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi8(xmm1,xmm1);
      xmm5 = _mm_srli_epi16(xmm5,8);
      xmm1 = _mm_srli_epi16(xmm1,8);
      xmm5 = _mm_slli_epi16(xmm5,3);
      xmm1 = _mm_slli_epi16(xmm1,3);

      xmm2 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx));

      xmm3 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx+1));
      xmm7 = _mm_unpackhi_epi8(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi8(xmm3,xmm3);
      xmm7 = _mm_srli_epi16(xmm7,8);
      xmm3 = _mm_srli_epi16(xmm3,8);
      xmm7 = _mm_slli_epi16(xmm7,3);
      xmm3 = _mm_slli_epi16(xmm3,3);


      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_w(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round32);
      xmm2 = _mm_srai_epi16(xmm2,6);
      xmm2 = avs_clip_0_255_w(xmm2);

      // loading
      xmm4 = _mm_loadu_si128((const __m128i*)(tmp02[yy]+xx+8-1));
      xmm6 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx+8));

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_w(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round32);
      xmm6 = _mm_srai_epi16(xmm6,6);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[0][1][yy]+xx), xmm2);

      // x x x o
      // x x x x
      // x x x x
      // x x x x
      // loading
      xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,3);
      xmm0 = _mm_slli_epi16(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx));

      xmm2 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx+1));
      xmm6 = _mm_unpackhi_epi8(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi8(xmm2,xmm2);
      xmm6 = _mm_srli_epi16(xmm6,8);
      xmm2 = _mm_srli_epi16(xmm2,8);
      xmm6 = _mm_slli_epi16(xmm6,3);
      xmm2 = _mm_slli_epi16(xmm2,3);

      xmm3 = _mm_loadu_si128((const __m128i*)(tmp02[yy]+xx+1));



      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_w(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round32);
      xmm2 = _mm_srai_epi16(xmm2,6);
      xmm2 = avs_clip_0_255_w(xmm2);

      // loading
      xmm5 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx+8));
      xmm7 = _mm_loadu_si128((const __m128i*)(tmp02[yy]+xx+8+1));

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_w(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round32);
      xmm6 = _mm_srai_epi16(xmm6,6);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[0][3][yy]+xx), xmm2);

      // x x x x
      // x x x x
      // x o x x
      // x x x x
      //////////////////////////////////bit width in 32 bits
      //loading
      xmm0 = _mm_loadu_si128((const __m128i*)(tmp22[yy]+xx-1));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);

      xmm1 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);
      xmm5 = _mm_slli_epi32(xmm5,3);
      xmm1 = _mm_slli_epi32(xmm1,3);

      xmm2 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);

      xmm3 = _mm_loadu_si128((const __m128i*)(tmp20[yy]+xx+1));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);
      xmm7 = _mm_slli_epi32(xmm7,3);
      xmm3 = _mm_slli_epi32(xmm3,3);


      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp0 = avs_combine_d2w(xmm2,xmm6);
      // low 8 results completed!

      //loading
      xmm0 = _mm_loadu_si128((const __m128i*)(tmp22[yy]+xx+8-1));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);

      xmm1 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx+8));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);
      xmm5 = _mm_slli_epi32(xmm5,3);
      xmm1 = _mm_slli_epi32(xmm1,3);


      xmm2 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx+8));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);

      xmm3 = _mm_loadu_si128((const __m128i*)(tmp20[yy]+xx+8+1));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);
      xmm7 = _mm_slli_epi32(xmm7,3);
      xmm3 = _mm_slli_epi32(xmm3,3);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp1 = avs_combine_d2w(xmm2,xmm6);
      // high 8 results completed!

      qtmp1 = _mm_shuffle_epi32(qtmp1,78);  //[1,0,3,2]
      qtmp0 = _mm_or_si128(qtmp0,qtmp1);

      _mm_store_si128((__m128i*)(interpolation[2][1][yy]+xx), qtmp0);

      // x x x x
      // x x x x
      // x x x o
      // x x x x
      //////////////////////////////////bit width in 32 bits
      xmm0 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);
      xmm4 = _mm_slli_epi32(xmm4,3);
      xmm0 = _mm_slli_epi32(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);


      xmm2 = _mm_loadu_si128((const __m128i*)(tmp20[yy]+xx+1));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);
      xmm6 = _mm_slli_epi32(xmm6,3);
      xmm2 = _mm_slli_epi32(xmm2,3);

      xmm3 = _mm_loadu_si128((const __m128i*)(tmp22[yy]+xx+1));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp0 = avs_combine_d2w(xmm2,xmm6);
      // low 8 results completed!

      //loading
      xmm0 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx+8));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);
      xmm4 = _mm_slli_epi32(xmm4,3);
      xmm0 = _mm_slli_epi32(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx+8));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);


      xmm2 = _mm_loadu_si128((const __m128i*)(tmp20[yy]+xx+1+8));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);
      xmm6 = _mm_slli_epi32(xmm6,3);
      xmm2 = _mm_slli_epi32(xmm2,3);

      xmm3 = _mm_loadu_si128((const __m128i*)(tmp22[yy]+xx+1+8));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp1 = avs_combine_d2w(xmm2,xmm6);
      // high 8 results completed!

      qtmp1 = _mm_shuffle_epi32(qtmp1,78);  //[1,0,3,2]
      qtmp0 = _mm_or_si128(qtmp0,qtmp1);

      _mm_store_si128((__m128i*)(interpolation[2][3][yy]+xx), qtmp0);
    }

  }
  // Stuffing Vertical Marginal Points for Quater Pels
  //        * o o *
  //        * * * *
  //        * * * *
  //        * o o *
  for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
  {
    // x o x o
    // x o x o
    // x o x o
    // x o x o
    xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][1][IMG_PAD_SIZE]+xx));
    xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][3][IMG_PAD_SIZE]+xx));

    xmm2 = _mm_load_si128((const __m128i*)(interpolation[2][1][img_pad_height-IMG_PAD_SIZE-1]+xx));
    xmm3 = _mm_load_si128((const __m128i*)(interpolation[2][3][img_pad_height-IMG_PAD_SIZE-1]+xx));
    for(yy=0; yy<IMG_PAD_SIZE; yy++)
    {
      _mm_store_si128((__m128i*)(interpolation[0][1][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[1][1][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[2][1][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[3][1][yy]+xx), xmm0);
      _mm_store_si128((__m128i*)(interpolation[0][3][yy]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[1][3][yy]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[2][3][yy]+xx), xmm1);
      _mm_store_si128((__m128i*)(interpolation[3][3][yy]+xx), xmm1);

      _mm_store_si128((__m128i*)(interpolation[0][1][img_pad_height-yy-1]+xx), xmm2);
      _mm_store_si128((__m128i*)(interpolation[1][1][img_pad_height-yy-1]+xx), xmm2);
      _mm_store_si128((__m128i*)(interpolation[2][1][img_pad_height-yy-1]+xx), xmm2);
      _mm_store_si128((__m128i*)(interpolation[3][1][img_pad_height-yy-1]+xx), xmm2);
      _mm_store_si128((__m128i*)(interpolation[0][3][img_pad_height-yy-1]+xx), xmm3);
      _mm_store_si128((__m128i*)(interpolation[1][3][img_pad_height-yy-1]+xx), xmm3);
      _mm_store_si128((__m128i*)(interpolation[2][3][img_pad_height-yy-1]+xx), xmm3);
      _mm_store_si128((__m128i*)(interpolation[3][3][img_pad_height-yy-1]+xx), xmm3);
    }
  }

  // x x x x
  // o x o x
  // x x x x
  // o x o x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      // x x x x
      // o x x x
      // x x x x
      // x x x x
      xmm0 = _mm_load_si128((const __m128i*)(tmp20[yy-1]+xx));

      xmm1 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm5 = _mm_unpackhi_epi8(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi8(xmm1,xmm1);
      xmm5 = _mm_srli_epi16(xmm5,8);
      xmm1 = _mm_srli_epi16(xmm1,8);
      xmm5 = _mm_slli_epi16(xmm5,3);
      xmm1 = _mm_slli_epi16(xmm1,3);

      xmm2 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx));

      xmm3 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy+1]+xx));
      xmm7 = _mm_unpackhi_epi8(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi8(xmm3,xmm3);
      xmm7 = _mm_srli_epi16(xmm7,8);
      xmm3 = _mm_srli_epi16(xmm3,8);
      xmm7 = _mm_slli_epi16(xmm7,3);
      xmm3 = _mm_slli_epi16(xmm3,3);

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_w(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round32);
      xmm2 = _mm_srai_epi16(xmm2,6);
      xmm2 = avs_clip_0_255_w(xmm2);

      // loading
      xmm4 = _mm_load_si128((const __m128i*)(tmp20[yy-1]+xx+8));
      xmm6 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx+8));

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_w(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round32);
      xmm6 = _mm_srai_epi16(xmm6,6);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[1][0][yy]+xx), xmm2);

      // x x x x
      // x x x x
      // x x x x
      // o x x x
      xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,3);
      xmm0 = _mm_slli_epi16(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx));

      xmm2 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy+1]+xx));
      xmm6 = _mm_unpackhi_epi8(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi8(xmm2,xmm2);
      xmm6 = _mm_srli_epi16(xmm6,8);
      xmm2 = _mm_srli_epi16(xmm2,8);
      xmm6 = _mm_slli_epi16(xmm6,3);
      xmm2 = _mm_slli_epi16(xmm2,3);

      xmm3 = _mm_load_si128((const __m128i*)(tmp20[yy+1]+xx));

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_w(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi16(xmm2,round32);
      xmm2 = _mm_srai_epi16(xmm2,6);
      xmm2 = avs_clip_0_255_w(xmm2);

      // loading
      xmm5 = _mm_load_si128((const __m128i*)(tmp20[yy]+xx+8));
      xmm7 = _mm_load_si128((const __m128i*)(tmp20[yy+1]+xx+8));

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_w(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi16(xmm6,round32);
      xmm6 = _mm_srai_epi16(xmm6,6);
      xmm6 = avs_clip_0_255_w(xmm6);

      xmm2 = avs_combine_w2b(xmm2,xmm6);

      _mm_store_si128((__m128i*)(interpolation[3][0][yy]+xx), xmm2);

      // x x x x
      // x x o x
      // x x x x
      // x x x x
      xmm0 = _mm_load_si128((const __m128i*)(tmp22[yy-1]+xx));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);

      xmm1 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);
      xmm5 = _mm_slli_epi32(xmm5,3);
      xmm1 = _mm_slli_epi32(xmm1,3);

      xmm2 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);

      xmm3 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);
      xmm7 = _mm_slli_epi32(xmm7,3);
      xmm3 = _mm_slli_epi32(xmm3,3);


      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp0 = avs_combine_d2w(xmm2,xmm6);
      // low 8 results completed!

      //loading
      xmm0 = _mm_load_si128((const __m128i*)(tmp22[yy-1]+xx+8));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);

      xmm1 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx+8));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);
      xmm5 = _mm_slli_epi32(xmm5,3);
      xmm1 = _mm_slli_epi32(xmm1,3);


      xmm2 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx+8));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);

      xmm3 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx+8));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);
      xmm7 = _mm_slli_epi32(xmm7,3);
      xmm3 = _mm_slli_epi32(xmm3,3);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp1 = avs_combine_d2w(xmm2,xmm6);
      // high 8 results completed!

      qtmp1 = _mm_shuffle_epi32(qtmp1,78);  //[1,0,3,2]
      qtmp0 = _mm_or_si128(qtmp0,qtmp1);

      _mm_store_si128((__m128i*)(interpolation[1][2][yy]+xx), qtmp0);

      // x x x x
      // x x x x
      // x x x x
      // x x o x
      //////////////////////////////////bit width in 32 bits
      xmm0 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);
      xmm4 = _mm_slli_epi32(xmm4,3);
      xmm0 = _mm_slli_epi32(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);


      xmm2 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);
      xmm6 = _mm_slli_epi32(xmm6,3);
      xmm2 = _mm_slli_epi32(xmm2,3);

      xmm3 = _mm_load_si128((const __m128i*)(tmp22[yy+1]+xx));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp0 = avs_combine_d2w(xmm2,xmm6);
      // low 8 results completed!

      //loading
      xmm0 = _mm_load_si128((const __m128i*)(tmp02[yy]+xx+8));
      xmm4 = _mm_unpackhi_epi16(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi16(xmm0,xmm0);
      xmm4 = _mm_srai_epi32(xmm4,16);
      xmm0 = _mm_srai_epi32(xmm0,16);
      xmm4 = _mm_slli_epi32(xmm4,3);
      xmm0 = _mm_slli_epi32(xmm0,3);

      xmm1 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx+8));
      xmm5 = _mm_unpackhi_epi16(xmm1,xmm1);
      xmm1 = _mm_unpacklo_epi16(xmm1,xmm1);
      xmm5 = _mm_srai_epi32(xmm5,16);
      xmm1 = _mm_srai_epi32(xmm1,16);


      xmm2 = _mm_load_si128((const __m128i*)(tmp02[yy+1]+xx+8));
      xmm6 = _mm_unpackhi_epi16(xmm2,xmm2);
      xmm2 = _mm_unpacklo_epi16(xmm2,xmm2);
      xmm6 = _mm_srai_epi32(xmm6,16);
      xmm2 = _mm_srai_epi32(xmm2,16);
      xmm6 = _mm_slli_epi32(xmm6,3);
      xmm2 = _mm_slli_epi32(xmm2,3);

      xmm3 = _mm_load_si128((const __m128i*)(tmp22[yy+1]+xx+8));
      xmm7 = _mm_unpackhi_epi16(xmm3,xmm3);
      xmm3 = _mm_unpacklo_epi16(xmm3,xmm3);
      xmm7 = _mm_srai_epi32(xmm7,16);
      xmm3 = _mm_srai_epi32(xmm3,16);

      // filtering
      xmm0 = _mm_add_epi32(xmm0,xmm3);
      xmm2 = avs_filter_quaterpel_d(xmm0,xmm1,xmm2);

      // cliping
      xmm2 = _mm_add_epi32(xmm2,round512);
      xmm2 = _mm_srai_epi32(xmm2,10);
      xmm2 = avs_clip_0_255_w(xmm2);

      // filtering
      xmm4 = _mm_add_epi32(xmm4,xmm7);
      xmm6 = avs_filter_quaterpel_d(xmm4,xmm5,xmm6);

      // cliping
      xmm6 = _mm_add_epi32(xmm6,round512);
      xmm6 = _mm_srai_epi32(xmm6,10);
      xmm6 = avs_clip_0_255_w(xmm6);

      // reorder
      qtmp1 = avs_combine_d2w(xmm2,xmm6);
      // high 8 results completed!

      qtmp1 = _mm_shuffle_epi32(qtmp1,78);  //[1,0,3,2]
      qtmp0 = _mm_or_si128(qtmp0,qtmp1);

      _mm_store_si128((__m128i*)(interpolation[3][2][yy]+xx), qtmp0);
    }
  }

  // Stuffing Horizontal Marginal Points for Quater Pels
  //        * * * *
  //        o * * o
  //        o * * o
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    // x x x x
    // o o o o
    // x x x x
    // o o o o
    tmpw0 = interpolation[1][0][yy][IMG_PAD_SIZE];
    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)(interpolation[1][0][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][yy] ), xmm0);

    tmpw0 = interpolation[3][0][yy][IMG_PAD_SIZE];
    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)(interpolation[3][0][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][yy] ), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][yy] ), xmm0);

    tmpw0 = interpolation[1][2][yy][img_pad_width-IMG_PAD_SIZE-1];
    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)(interpolation[1][0][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][1][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][2][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[1][3][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);

    tmpw0 = interpolation[3][2][yy][img_pad_width-IMG_PAD_SIZE-1];
    tmpw0 = tmpw0+(tmpw0<<8);
    xmm0  = _mm_insert_epi16(xmm0,tmpw0,0);
    xmm0  = _mm_shufflelo_epi16(xmm0,0);
    xmm0  = _mm_shuffle_epi32(xmm0,0);

    _mm_store_si128((__m128i*)(interpolation[3][0][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][1][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][2][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
    _mm_store_si128((__m128i*)(interpolation[3][3][yy]+img_pad_width-IMG_PAD_SIZE), xmm0);
  }
  // x x x x
  // x o x o
  // x x x x
  // x o x o
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx=xx+16)
    {
      // loading
      xmm1 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx));
      xmm1 = _mm_srai_epi16(xmm1,1);         // bit width control

      xmm2 = _mm_load_si128((const __m128i*)(tmp22[yy]+xx+8));
      xmm2 = _mm_srai_epi16(xmm2,1);         // bit width control

      // x x x x
      // x o x x
      // x x x x
      // x x x x
      // loading
      xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy]+xx));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,5);         // shifting with bit width control
      xmm0 = _mm_slli_epi16(xmm0,5);       // shifting with bit width control

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm1);

      // cliping
      xmm0 = _mm_add_epi16(xmm0,round32);
      xmm0 = _mm_srai_epi16(xmm0,6);
      xmm0 = avs_clip_0_255_w(xmm0);

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm2);

      // cliping
      xmm4 = _mm_add_epi16(xmm4,round32);
      xmm4 = _mm_srai_epi16(xmm4,6);
      xmm4 = avs_clip_0_255_w(xmm4);

      xmm0 = avs_combine_w2b(xmm0,xmm4);

      _mm_store_si128((__m128i*)(interpolation[1][1][yy]+xx), xmm0);

      // x x x x
      // x x x o
      // x x x x
      // x x x x
      // loading

      xmm0 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy]+xx+1));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,5);         // shifting with bit width control
      xmm0 = _mm_slli_epi16(xmm0,5);       // shifting with bit width control

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm1);

      // cliping
      xmm0 = _mm_add_epi16(xmm0,round32);
      xmm0 = _mm_srai_epi16(xmm0,6);
      xmm0 = avs_clip_0_255_w(xmm0);

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm2);

      // cliping
      xmm4 = _mm_add_epi16(xmm4,round32);
      xmm4 = _mm_srai_epi16(xmm4,6);
      xmm4 = avs_clip_0_255_w(xmm4);

      xmm0 = avs_combine_w2b(xmm0,xmm4);

      _mm_store_si128((__m128i*)(interpolation[1][3][yy]+xx), xmm0);

      // x x x x
      // x x x x
      // x x x x
      // x o x x
      // loading
      xmm0 = _mm_load_si128((const __m128i*)(interpolation[0][0][yy+1]+xx));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,5);         // shifting with bit width control
      xmm0 = _mm_slli_epi16(xmm0,5);       // shifting with bit width control

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm1);

      // cliping
      xmm0 = _mm_add_epi16(xmm0,round32);
      xmm0 = _mm_srai_epi16(xmm0,6);
      xmm0 = avs_clip_0_255_w(xmm0);

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm2);

      // cliping
      xmm4 = _mm_add_epi16(xmm4,round32);
      xmm4 = _mm_srai_epi16(xmm4,6);
      xmm4 = avs_clip_0_255_w(xmm4);

      xmm0 = avs_combine_w2b(xmm0,xmm4);

      _mm_store_si128((__m128i*)(interpolation[3][1][yy]+xx), xmm0);

      // x x x x
      // x x x x
      // x x x x
      // x x x o
      // loading
      xmm0 = _mm_loadu_si128((const __m128i*)(interpolation[0][0][yy+1]+xx+1));
      xmm4 = _mm_unpackhi_epi8(xmm0,xmm0);
      xmm0 = _mm_unpacklo_epi8(xmm0,xmm0);
      xmm4 = _mm_srli_epi16(xmm4,8);
      xmm0 = _mm_srli_epi16(xmm0,8);
      xmm4 = _mm_slli_epi16(xmm4,5);         // shifting with bit width control
      xmm0 = _mm_slli_epi16(xmm0,5);       // shifting with bit width control

      // filtering
      xmm0 = _mm_add_epi16(xmm0,xmm1);

      // cliping
      xmm0 = _mm_add_epi16(xmm0,round32);
      xmm0 = _mm_srai_epi16(xmm0,6);
      xmm0 = avs_clip_0_255_w(xmm0);

      // filtering
      xmm4 = _mm_add_epi16(xmm4,xmm2);

      // cliping
      xmm4 = _mm_add_epi16(xmm4,round32);
      xmm4 = _mm_srai_epi16(xmm4,6);
      xmm4 = avs_clip_0_255_w(xmm4);

      xmm0 = avs_combine_w2b(xmm0,xmm4);

      _mm_store_si128((__m128i*)(interpolation[3][3][yy]+xx), xmm0);
    }
  }
  // test
  //posy=3;posx=3;
  /*    for(posy=0;posy<4;posy++)
  for(posx=0;posx<4;posx++)
  for(yy=0; yy<img_pad_height; yy++)
  {
  for(xx=0; xx<img_pad_width; xx++)
  if(interpolation[posy][posx][yy][xx] != mref[0][4*yy+posy][4*xx+posx])
  //if(tmp20[yy][xx]*8 != img4Y_tmp[4*yy+2][4*xx+0])
  //if(interpolation[posy][posx][yy][xx]-mref[0][4*yy+posy][4*xx+posx]>7)
  temp=temp+abs(interpolation[posy][posx][yy][xx]-mref[0][4*yy+posy][4*xx+posx]);
  }
  temp=temp/(21504*16);
  */
}
// xzhao }


/*
*************************************************************************
* Function:Upsample 4 times, store them in out4x.  Color is simply copied
* Input:srcy, srcu, srcv, out4y, out4u, out4v
* Output:
* Return:
* Attention:Side Effects_
Uses (writes) img4Y_tmp.  This should be moved to a static variable
in this module
*************************************************************************
*/
#define  IClip( Min, Max, Val) (((Val)<(Min))? (Min):(((Val)>(Max))? (Max):(Val)))

void c_avs_enc::find_snr ()
{
  int i, j;
  int diff_y, diff_u, diff_v;
  int impix;

  //  Calculate  PSNR for Y, U and V.
  //     Luma.
  impix = img->height * img->width;

  diff_y = 0;
  for (i = 0; i < img->width; ++i)
  {
    for (j = 0; j < img->height; ++j)
    {
      diff_y += img->quad[imgY_org[j][i] - imgY[j][i]];
    }
  }

  //     Chroma.
  diff_u = 0;
  diff_v = 0;

  for (i = 0; i < img->width_cr; i++)
  {
    for (j = 0; j < img->height_cr; j++)
    {
      diff_u += img->quad[imgUV_org[0][j][i] - imgUV[0][j][i]];
      diff_v += img->quad[imgUV_org[1][j][i] - imgUV[1][j][i]];
    }
  }

  //  Collecting SNR statistics
  if (diff_y != 0)
  {
    snr->snr_y = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));         // luma snr for current frame
    snr->snr_u = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));   // u chroma snr for current frame, 1/4 of luma samples
    snr->snr_v = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));   // v chroma snr for current frame, 1/4 of luma samples
  }

  if (img->number == 0)
  {
    snr->snr_y1 = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));        // keep luma snr for first frame
    snr->snr_u1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));  // keep chroma u snr for first frame
    snr->snr_v1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));  // keep chroma v snr for first frame
    snr->snr_ya = snr->snr_y1;
    snr->snr_ua = snr->snr_u1;
    snr->snr_va = snr->snr_v1;
  }
  // B pictures
  else
  {
    snr->snr_ya = (float) (snr->snr_ya * (img->number + Bframe_ctr) + snr->snr_y) / (img->number + Bframe_ctr + 1); // average snr lume for all frames inc. first
    snr->snr_ua = (float) (snr->snr_ua * (img->number + Bframe_ctr) + snr->snr_u) / (img->number + Bframe_ctr + 1); // average snr u chroma for all frames inc. first
    snr->snr_va = (float) (snr->snr_va * (img->number + Bframe_ctr) + snr->snr_v) / (img->number + Bframe_ctr + 1); // average snr v chroma for all frames inc. first
  }
}


/*
*************************************************************************
* Function:GBIM方法估计块效应
* Input: imgY
* Output:
* Return: GBIM value

*************************************************************************
*/
double c_avs_enc::find_GBIM(byte **I)
{
  int i,j,k,n;
  byte **I_trans;
  I_trans = new byte* [img->width];
  for(i=0; i<img->width; i++)
  {
    I_trans[i] = new byte[img->height];
  }

  double Weighted_Value_v=0;
  double Weighted_Value=0;
  double Weighted_Value_interpixel=0;
  double M_hor=0;
  double M_ver=0;
  double M_image=0;



  //水平块效应估计
  for (k=1;k<=img->width/8-1;k++)
  {
    for (i=0;i<img->height;i++)
    {
      double Local_mean=0;
      double Local_mean_left=0;
      double Local_mean_right=0;
      double Local_activity=0;
      double Local_activity_left=0;
      double Local_activity_right=0;
      double Weight;

      //求左右临块的亮度均值
      for(n=k*8-8;n<=k*8-1;n++)
      {
        Local_mean_left+=I[i][n]/8;
      }
      for(n=k*8;n<=k*8+7;n++)
      {
        Local_mean_right+=I[i][n]/8;
      }
      Local_mean=(Local_mean_left+Local_mean_right)/2;

      //求左右临块的亮度变化值
      for(n=k*8-8;n<=k*8-1;n++)
      {
        Local_activity_left+=(I[i][n]-Local_mean_left)*(I[i][n]-Local_mean_left)/8;
      }
      for(n=k*8;n<=k*8+7;n++)
      {
        Local_activity_right+=(I[i][n]-Local_mean_right)*(I[i][n]-Local_mean_right)/8;
      }
      Local_activity=(sqrt(Local_activity_left)+sqrt(Local_activity_right))/2;

      //求权值系数
      if (Local_mean<=81)
      {
        Weight=1.15201*log(1+(sqrt(Local_mean)/(1+Local_activity)));
      }
      else
      {
        Weight=log(1+(sqrt(255-Local_mean)/(1+Local_activity)));
      }

      Weighted_Value_v = Weight*(I[i][k*8-1]-I[i][k*8]);
      Weighted_Value += Weighted_Value_v*Weighted_Value_v;
      for(int m =1;m<=7;m++)
      {Weighted_Value_interpixel += (Weight*(I[i][k*8+m-1]-I[i][k*8+m]))*(Weight*(I[i][k*8+m-1]-I[i][k*8+m]));}

    }
  }

  M_hor=sqrt(Weighted_Value)/sqrt(Weighted_Value_interpixel)*7;
  Weighted_Value=0;
  Weighted_Value_interpixel=0;


  //图像存储矩阵转置
  for (j=0;j<img->height;j++)
  {
    for (i=0;i<img->width;i++)
    {
      I_trans[i][j]=I[j][i];
    }
  }

  //竖直块效应估计
  for (k=1;k<=img->height/8-1;k++)
  {
    for (i=0;i<img->width;i++)
    {
      double Local_mean=0;
      double Local_mean_left=0;
      double Local_mean_right=0;
      double Local_activity=0;
      double Local_activity_left=0;
      double Local_activity_right=0;
      double Weight=0;

      //求左右临块的亮度均值
      for(n=k*8-8;n<=k*8-1;n++)
      {
        Local_mean_left+=I_trans[i][n]/8;
      }
      for(n=k*8;n<=k*8+7;n++)
      {
        Local_mean_right+=I_trans[i][n]/8;
      }
      Local_mean=(Local_mean_left+Local_mean_right)/2;

      //求左右临块的亮度变化值
      for(n=k*8-8;n<=k*8-1;n++)
      {
        Local_activity_left+=(I_trans[i][n]-Local_mean_left)*(I_trans[i][n]-Local_mean_left)/8;
      }
      for(n=k*8;n<=k*8+7;n++)
      {
        Local_activity_right+=(I_trans[i][n]-Local_mean_right)*(I_trans[i][n]-Local_mean_right)/8;
      }
      Local_activity=(sqrt(Local_activity_left)+sqrt(Local_activity_right))/2;

      //求权值系数
      if (Local_mean<=81)
      {
        Weight=1.15201*log(1+(sqrt(Local_mean)/(1+Local_activity)));
      }
      else
      {
        Weight=log(1+(sqrt(255-Local_mean)/(1+Local_activity)));
      }

      Weighted_Value+=(Weight*(I_trans[i][k*8-1]-I_trans[i][k*8]))*(Weight*(I_trans[i][k*8-1]-I_trans[i][k*8]));

      for(int m=1;m<=7;m++)
      {Weighted_Value_interpixel+=(Weight*(I_trans[i][k*8+m-1]-I_trans[i][k*8+m]))*(Weight*(I_trans[i][k*8+m-1]-I_trans[i][k*8+m]));}

    }
  }
  M_ver=sqrt(Weighted_Value)/sqrt(Weighted_Value_interpixel)*7;

  for (i=0; i<img->width; i++)
  {
    delete [] I_trans[i];
  }
  delete I_trans;

  //图像块效应估计值
  M_image=(M_hor+M_ver)/2;
  return (M_image);
}

/*
*************************************************************************
* Function:Find distortion for all three components
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
void c_avs_enc::find_distortion ()
{
  int i, j;
  int diff_y, diff_u, diff_v;
  int impix;

  //  Calculate  PSNR for Y, U and V.
  //     Luma.
  impix = img->height * img->width;

  diff_y = 0;
  for (i = 0; i < img->width; ++i)
  {
    for (j = 0; j < img->height; ++j)
    {
      diff_y += img->quad[abs(imgY_org[j][i] - imgY[j][i])];
    }
  }

  //     Chroma.
  diff_u = 0;
  diff_v = 0;

  for (i = 0; i < img->width_cr; i++)
  {
    for (j = 0; j < img->height_cr; j++)
    {
      diff_u += img->quad[abs (imgUV_org[0][j][i] - imgUV[0][j][i])];
      diff_v += img->quad[abs (imgUV_org[1][j][i] - imgUV[1][j][i])];
    }
  }

  // Calculate real PSNR at find_snr_avg()
  snr->snr_y = (float) diff_y;
  snr->snr_u = (float) diff_u;
  snr->snr_v = (float) diff_v;

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

void c_avs_enc::ReportFirstframe(int tmp_time)
{
  printf ("%3d(I)  %8d %4d %7.4f %7.4f %7.4f %7.4f  %5d       %s \n",
    frame_no, stat->bit_ctr - stat->bit_ctr_n,
    img->qp, snr->snr_y, snr->snr_u, snr->snr_v, GBIM_value_frm, tmp_time, img->picture_structure ? "FRM":"FLD" );

  //Rate control
  if(input->RCEnable && input->InterlaceCodingOption != 0)
  {
    Iprev_bits = stat->bit_ctr;
  }

  stat->bitr0     = stat->bitr;
  stat->bit_ctr_0 = stat->bit_ctr;
  stat->bit_ctr   = 0;
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

void c_avs_enc::ReportIntra(int tmp_time)
{
  //FILE *file = fopen("stat.dat","at");
  //fprintf (file,"%3d(I)  %8d %4d %7.4f %7.4f %7.4f  %5d      \n",
  //  frame_no, stat->bit_ctr - stat->bit_ctr_n,
  //  img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time );
  //
  //fclose(file);

  printf ("\n%3d(I)  %8u %4d %7.4f %7.4f %7.4f %7.4f  %5d       %3s\n",
    frame_no, stat->bit_ctr - stat->bit_ctr_n,
    img->qp, snr->snr_y, snr->snr_u, snr->snr_v,GBIM_value_frm, tmp_time, img->picture_structure ? "FRM":"FLD");
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

void c_avs_enc::ReportB(int tmp_time)
{
  //FILE *file = fopen("stat.dat","at");
  //fprintf (file,"%3d(B)  %8d %4d %7.4f %7.4f %7.4f  %5d        \n",
  //  frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp,
  //  snr->snr_y, snr->snr_u, snr->snr_v, tmp_time);
  //
  //fclose(file);

  printf ("%3d(B)  %8u %4d %7.4f %7.4f %7.4f %7.4f  %5d       %3s     %3d    \n",
    frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp,
    snr->snr_y, snr->snr_u, snr->snr_v,GBIM_value_frm, tmp_time, img->picture_structure ? "FRM":"FLD",intras);

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

void c_avs_enc::ReportP(int tmp_time)
{
  //FILE *file = fopen("stat.dat","at");
  //fprintf (file,"%3d(P)  %8u %4d %7.4f %7.4f %7.4f  %5d        %3d\n",
  //  frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp, snr->snr_y,
  //  snr->snr_u, snr->snr_v, tmp_time,
  //  intras);
  //fclose(file);

  printf ("%3d(P)  %8u %4d %7.4f %7.4f %7.4f %7.4f  %5d       %3s     %3d    \n",
    frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp, snr->snr_y,
    snr->snr_u, snr->snr_v, GBIM_value_frm, tmp_time,
    img->picture_structure ? "FRM":"FLD",intras);
}

/*
*************************************************************************
* Function:Copies contents of a Source frame structure into the old-style
*    variables imgY_org_frm and imgUV_org_frm.  No other side effects
* Input:  sf the source frame the frame is to be taken from
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::CopyFrameToOldImgOrgVariables ()
{
  int x, y, xx, yy, adr;
  byte *u_buffer,*v_buffer;
  u_buffer = imgY_org_buffer + bytes_y;
  v_buffer = imgY_org_buffer + bytes_y + bytes_uv;
  for (y=0; y<img->height; y++)
  {
    for (x=0; x<img->width; x++)
    {
      imgY_org_frm [y][x] = imgY_org_buffer[y*img->width+x];
      if (y&1 && x&1)
      {
        xx=x>>1;
        yy=y>>1;
        adr=yy*img->width/2+xx;

        imgUV_org_frm[0][yy][xx] = u_buffer[adr];
        imgUV_org_frm[1][yy][xx] = v_buffer[adr];
      }
    }
  }
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    // xzhao { 2007.7.14
    for (y=0; y<img->height; y += 2)
    {
      for (x=0; x<img->width; x++)
      {
        imgY_org_top [y/2][x] = imgY_org_buffer[y*img->width    +x]; // !! Lum component for top field
        imgY_org_bot [y/2][x] = imgY_org_buffer[(y+1)*img->width+x]; // !! Lum component for bot field
        xx=x>>2;
        yy=y>>2;
        adr=yy*img->width/2+xx;
        imgUV_org_top[0][yy][xx] = u_buffer[adr];    // !! Cr and Cb component for top field
        imgUV_org_top[1][yy][xx] = v_buffer[adr];
        imgUV_org_bot[0][yy][xx] = u_buffer[adr-img->width/2]; // !! Cr and Cb component for bot field
        imgUV_org_bot[1][yy][xx] = v_buffer[adr-img->width/2];
      }
    }
  }
#ifdef ROI_ENABLE
  memcpy(YCbCr[0], imgY_org_buffer,                  bytes_y *sizeof(byte));
  memcpy(YCbCr[1], imgY_org_buffer+bytes_y,          bytes_uv*sizeof(byte));
  memcpy(YCbCr[2], imgY_org_buffer+bytes_y+bytes_uv, bytes_uv*sizeof(byte));
  w[0] = img->width;
  w[1] = w[2] = img->width_cr;
  h[0] = img->height;
  h[1] = h[2] = img->height_cr;
#endif
}

/*
*************************************************************************
* Function: Calculates the absolute frame number in the source file out
of various variables in img-> and input->
* Input:
* Output:
* Return: frame number in the file to be read
* Attention: \side effects
global variable frame_no updated -- dunno, for what this one is necessary
*************************************************************************
*/
void c_avs_enc::CalculateFrameNumber()
{
  frame_no = picture_distance;
  if (img->type == B_IMG)
  {
    // xzhao 20080329
    //frame_no = (img->number - 1) * (input->successive_Bframe + 1) + img->b_interval * img->b_frame_to_code;
    frame_no = gframe_no - 1;
  }
  else
  {
    if(img->type==INTRA_IMG)
    {
      frame_no = gframe_no;
    }
    // xzhao 20080329
    else if((gframe_no%input->GopLength)>=input->GopLength - input->successive_Bframe)
    {
      frame_no = gframe_no;
    }
    else
    {
      //frame_no = img->number * (1 + input->successive_Bframe);
      frame_no = gframe_no + input->successive_Bframe;
    }
    //frame_no = img->number * (input->successive_Bframe + 1);
  }
}


void c_avs_enc::ReadOneFrame ()
{
  int i, j;
  int stuff_height_cr = (input->img_height-input->stuff_height)/2;
  // xzhao { 2007.7.18
  if(input->img_height != input->stuff_height)
  {
    //Y
    memcpy(imgY_org_buffer, pInputImage, bytes_y);
    for(j = input->stuff_height; j<input->img_height; j++)
    {
      for(i = 0; i <input->stuff_width; i++)
      {
        imgY_org_buffer[j*input->stuff_width + i] = imgY_org_buffer[input->stuff_height*input->stuff_width - input->stuff_width +i];
      }
    }
    //U
    memcpy( imgY_org_buffer + bytes_y, pInputImage + bytes_y, bytes_uv);
    for(j = 0; j < stuff_height_cr; j++)
    {
      for(i = 0; i <input->stuff_width/2; i++)
      {
        imgY_org_buffer[bytes_y+ bytes_uv +j*input->stuff_width/2 + i] = imgY_org_buffer[bytes_y+ bytes_uv - input->stuff_width/2 +i];
      }
    }
    //V
    memcpy( imgY_org_buffer + bytes_y + bytes_uv, pInputImage + bytes_y + bytes_uv, bytes_uv);
    for(j = 0; j < stuff_height_cr; j++)
    {
      for(i = 0; i <input->stuff_width/2; i++)
      {
        imgY_org_buffer[bytes_y + bytes_uv + bytes_uv + j*input->stuff_width/2 + i] = imgY_org_buffer[bytes_y+bytes_uv+bytes_uv - input->stuff_width/2 +i];
      }
    }

  }
  else
  {
    memcpy(imgY_org_buffer, pInputImage, bytes_y+2*bytes_uv);
  } 
  // xzhao }
}


/*
*************************************************************************
* Function:point to frame coding variables
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::put_buffer_frame()
{
  int i,j;

  imgY_org  = imgY_org_frm;
  imgUV_org = imgUV_org_frm;
  tmp_mv    = tmp_mv_frm;

  //initialize ref index 1/4 pixel
  for(i=0;i<2;i++)
  {
    mref[i] = mref_frm[i];
  }

  //integer pixel for chroma
  for(i=0;i<2;i++)
  {
    for(j=0;j<2;j++)
    {
      mcef[i][j]   = ref_frm[i][j+1];
    }
  }

  //integer pixel for luma
  for(i=0;i<2;i++)
  {
    Refbuf11[i] = &ref_frm[i][0][0][0];
  }

  //current reconstructed image
  imgY  = imgY_frm  =   current_frame[0];
  imgUV = imgUV_frm =  &current_frame[1];

  refFrArr    = refFrArr_frm;
  fw_refFrArr = fw_refFrArr_frm;
  bw_refFrArr = bw_refFrArr_frm;
}


/*
*************************************************************************
* Function:update the decoder picture buffer
* Input:frame number in the bitstream and the video sequence
* Output:
* Return:
* Attention:
*************************************************************************
*/

void c_avs_enc::Update_Picture_Buffers()
{
  unsigned char ***tmp;
  byte ****tmp_y;
  int i;

  //update integer pixel reference buffer
  tmp = ref_frm[1];       //ref_frm[ref_index][yuv][height][width] ref_index = 0,1  for P frame
  ref_frm[1] = ref_frm[0];    // ref_index = 0, backward reference for B frame; 1: forward reference for B frame
  ref_frm[0] = current_frame; // current_frame: current image under reconstruction
  current_frame = tmp;

  //update luma 1/4 pixel reference buffer mref[ref_index][height][width] ref_index = 0,1 for P frame
  tmp_y = mref_frm[1];              // ref_index = 0, forward refernce for B frame ; 1: backward refernce for B frame
  mref_frm[1] = mref_frm[0];
  mref_frm[0] = tmp_y;

  //initial reference index, and for coming interpolation in mref[0]
  for(i=0;i<2;i++)
  {
    mref[i] = mref_frm[i];
  }
}



int c_avs_enc::DetectLumVar()
{
  int i , j ;
  int Histogtam_Cur[256] ;
  int Histogtam_Pre[256] ;
  int temp = 0 ;

  for( i = 0 ; i < 256 ; i++){
    Histogtam_Cur[i] = 0 ;
    Histogtam_Pre[i] = 0 ;
  }

  for(j = 0 ; j < img->height ; j++){
    for( i = 0 ; i < img->width ; i++){
      Histogtam_Cur[imgY_org[j][i]] += 1 ;
      Histogtam_Pre[Refbuf11[0][j*img->width + i]] += 1 ;
    }
  }

  for(i = 0 ; i < 256 ; i++){
    temp += abs(Histogtam_Pre[i] - Histogtam_Cur[i]);
  }

  // if(temp >= ((img->height*img->width)*2)){
  if(temp >= ((img->height*img->width)/4)){
    return 1;
  }
  else
  {
    return 0;
  }
}

void c_avs_enc::CalculateBrightnessPar(int currentblock[16][16] , int preblock[16][16] , float *c , float *d)
{
  int N = 256 ;
  int i , j ;
  int m1,m2,m3,m4,m5,m6;

  m1 = m2 = m3 = m4 = m5 = m6 = 0 ;
  for(j = 0 ; j < 16 ; j++){
    for(i = 0 ; i < 16 ; i++){
      m1 += preblock[j][i]*preblock[j][i] ;
      m2 += preblock[j][i];
      m3 += preblock[j][i];
      m4 += 1;
      m5 += preblock[j][i]*currentblock[j][i] ;
      m6 += currentblock[j][i];
    }
  }
  *c = ((float)(m4*m5 - m2*m6)) / ((float)(m1*m4 - m2*m3));
  *d = ((float)(m3*m5 - m6*m1)) / ((float)(m3*m2 - m1*m4));
  return ;
}

void c_avs_enc::CalculatePar(int refnum)
{
  int mbx , mby ;
  int currmb[16][16] ;
  int refmb[16][16] ;
  float alpha ;
  float belta ;
  int i , j ;
  int Alpha_His[256];
  int Belta_His[256];
  int max_num = 0 ;
  int max_index = -1 ;
  int belta_sum = 0 ;

  for( i = 0 ; i < 256 ; i++){
    Alpha_His[i] = 0 ;
    Belta_His[i] = 0 ;
  }

  for(mby = 0 ; mby < img->height/16 ; mby++){
    for(mbx = 0 ; mbx < img->width/16 ; mbx++){
      for( j = 0 ; j < 16 ; j++){
        for( i = 0 ; i < 16 ; i++){
          currmb[j][i] = imgY_org[mby*16+j][mbx*16+i];
          refmb [j][i] = Refbuf11[refnum][(mby*16+j)*img->width + mbx*16+i] ;
        }
      }
      CalculateBrightnessPar(currmb,refmb,&alpha,&belta);
      allalpha_lum[mby*(img->width/16)+mbx] = (int)(alpha*32);
      allbelta_lum[mby*(img->width/16)+mbx] = (int)(belta);
    }
  }

  for(i = 0 ; i < ((img->height/16)*(img->width/16)) ; i++)
  {
    if((allalpha_lum[i] < 256)&&(abs(allbelta_lum[i]) < 127))
    {
      Alpha_His[allalpha_lum[i]]++;
    }
  }

  for( i = 4 ; i < 256 ; i++) // !! 4-256 shenyanfei
  {
    if(Alpha_His[i] > max_num)
    {
      max_num = Alpha_His[i] ;
      max_index = i ;
    }
  }

  for( i = 0 ; i < ((img->height/16)*(img->width/16)) ; i++){
    if(allalpha_lum[i] == max_index){
      belta_sum += allbelta_lum[i] ;
    }
  }
  img->lum_scale[refnum] = max_index ;
  img->lum_shift[refnum] = belta_sum/max_num ;

  if(max_num > ((img->height/16)*(img->width/16) / 2))
    img->allframeweight = 1 ;
  else
    img->allframeweight = 0 ;

  img->chroma_scale[refnum] = 32 ; // !! default value
  img->chroma_shift[refnum] = 0  ; // !! default value
  return ;
}

void c_avs_enc::estimate_weighting_factor()
{
  int   bframe    = (img->type==B_IMG);
  int   max_ref   = img->nb_references;
  int   ref_num ;

  if(max_ref > img->buf_cycle)
    max_ref = img->buf_cycle;

  // !! detection luminance variation
  img->LumVarFlag = DetectLumVar();
  if(img->LumVarFlag == 1){
    for(ref_num = 0 ; ref_num < max_ref ; ref_num++){
      CalculatePar(ref_num);
    }
  }
  return;
}
void c_avs_enc::UnifiedOneForthPix_c_sse (pel_t ** imgY)
{
  int img_pad_width,img_pad_height;
  int xx,yy;
  int temp=0;
  int_16_t hh,hv;
  int qh,qv;  //qx for direction "\", qy for direction "/"




  img_pad_width  = (img->width  + (IMG_PAD_SIZE<<1));
  img_pad_height = (img->height + (IMG_PAD_SIZE<<1));

  interpolation = mref[0];

  /*************************************************************
  //           Basic Quater Pel Interpolation Unit
  //
  //              ○ ==> Interger Pel Position
  //              □ ==> Half     Pel Position
  //              △ ==> Quarter  Pel Position
  //************************************************************
  //                   ○   △   □   △
  //
  //                   △   △   △   △
  //
  //                   □   △   □   △
  //
  //                   △   △   △   △
  *************************************************************/

  /*if(frame_no==6)
  {
  printf("interpolation debug:\n");
  }*/
  // o x x x
  // x x x x
  // x x x x
  // x x x x
  //        * * * *
  //        * o o *
  //        * o o *
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
    for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx++)
      interpolation[0][0][yy][xx] = imgY[yy-IMG_PAD_SIZE][xx-IMG_PAD_SIZE];

  //        * * * *
  //        o o o o
  //        o o o o
  //        * * * *
  for(yy=IMG_PAD_SIZE; yy<img_pad_height-IMG_PAD_SIZE; yy++)
  {
    for(xx=0;xx<IMG_PAD_SIZE;xx++)
    {
      interpolation[0][0][yy][xx] = imgY[yy-IMG_PAD_SIZE][0];
      interpolation[0][0][yy][img_pad_width-xx-1] = imgY[yy-IMG_PAD_SIZE][img->width-1];
    }
  }
  //        * o o *
  //        o o o o
  //        o o o o
  //        * o o *
  for(xx=IMG_PAD_SIZE; xx<img_pad_width-IMG_PAD_SIZE; xx++)
  {
    for(yy=0;yy<IMG_PAD_SIZE;yy++)
    {
      interpolation[0][0][yy][xx] = imgY[0][xx-IMG_PAD_SIZE];
      interpolation[0][0][img_pad_height-yy-1][xx] = imgY[img->height-1][xx-IMG_PAD_SIZE];
    }
  }
  //        o o o o
  //        o o o o
  //        o o o o
  //        o o o o
  for(yy=0; yy<IMG_PAD_SIZE; yy++)
    for(xx=0; xx<IMG_PAD_SIZE; xx++)
    {
      interpolation[0][0][yy][xx]=imgY[0][0];
      interpolation[0][0][yy][img_pad_width-xx-1]=imgY[0][img->width-1];
      interpolation[0][0][img_pad_height-yy-1][xx]=imgY[img->height-1][0];
      interpolation[0][0][img_pad_height-yy-1][img_pad_width-xx-1]=imgY[img->height-1][img->width-1];
    }

    // x x o x
    // x x x x
    // x x x x
    // x x x x
    //        * o * *
    //        * o * *
    //        * o * *
    //        * o * *
    for(yy=0; yy<img_pad_height; yy++)
      for(xx=1; xx<img_pad_width-2; xx++)
      {
        hh = 5 * (interpolation[0][0][yy][xx] + interpolation[0][0][yy][xx+1])
          - (interpolation[0][0][yy][xx-1] + interpolation[0][0][yy][xx+2]);
        interpolation[0][2][yy][xx] = IClip(0, 255,((hh+4)>>3));
        tmp02[yy][xx] = hh;
      }
      //        o o o o
      //        o o o o
      //        o o o o
      //        o o o o
      for(yy=0; yy<img_pad_height; yy++)
      {
        // pos 02
        hh = 5 * (interpolation[0][0][yy][0] + interpolation[0][0][yy][1])
          - (interpolation[0][0][yy][0] + interpolation[0][0][yy][2]);
        interpolation[0][2][yy][0] = IClip(0, 255,((hh+4)>>3));
        tmp02[yy][0] = hh;

        hh = 5 * (interpolation[0][0][yy][img_pad_width-2] + interpolation[0][0][yy][img_pad_width-1])
          - (interpolation[0][0][yy][img_pad_width-3] + interpolation[0][0][yy][img_pad_width-1]);
        interpolation[0][2][yy][img_pad_width-2] = IClip(0, 255,((hh+4)>>3));
        tmp02[yy][img_pad_width-2] = hh;

        hh = 5 * (interpolation[0][0][yy][img_pad_width-1] + interpolation[0][0][yy][img_pad_width-1])
          - (interpolation[0][0][yy][img_pad_width-2] + interpolation[0][0][yy][img_pad_width-1]);
        interpolation[0][2][yy][img_pad_width-1] = IClip(0, 255,((hh+4)>>3));
        tmp02[yy][img_pad_width-1] = hh;
      }

      // x x x x
      // x x x x
      // o x o x
      // x x x x
      //        * * * *
      //        o o o o
      //        * * * *
      //        * * * *
      for(yy=1; yy<img_pad_height-2; yy++)
        for(xx=0; xx<img_pad_width; xx++)
        {

          hv = 5 * (interpolation[0][0][yy][xx] + interpolation[0][0][yy+1][xx])
            - (interpolation[0][0][yy-1][xx] + interpolation[0][0][yy+2][xx]);
          interpolation[2][0][yy][xx] = IClip(0, 255,((hv+4)>>3));
          tmp20[yy][xx] = (hv<<3);

          hv = 5 * (tmp02[yy][xx] + tmp02[yy+1][xx])
            - (tmp02[yy-1][xx] + tmp02[yy+2][xx]);
          interpolation[2][2][yy][xx] = IClip(0, 255,((hv+32)>>6));
          tmp22[yy][xx] = hv;
        }
        //        o o o o
        //        o o o o
        //        o o o o
        //        o o o o
        for(xx=0; xx<img_pad_width; xx++)
        {
          // pos 20
          hv = 5 * (interpolation[0][0][0][xx] + interpolation[0][0][1][xx])
            - (interpolation[0][0][0][xx] + interpolation[0][0][2][xx]);
          interpolation[2][0][0][xx] = IClip(0, 255,((hv+4)>>3));
          tmp20[0][xx] = (hv<<3);

          hv = 5 * (interpolation[0][0][img_pad_height-2][xx] + interpolation[0][0][img_pad_height-1][xx])
            - (interpolation[0][0][img_pad_height-3][xx] + interpolation[0][0][img_pad_height-1][xx]);
          interpolation[2][0][img_pad_height-2][xx] = IClip(0, 255,((hv+4)>>3));
          tmp20[img_pad_height-2][xx] = (hv<<3);

          hv = 5 * (interpolation[0][0][img_pad_height-1][xx] + interpolation[0][0][img_pad_height-1][xx])
            - (interpolation[0][0][img_pad_height-2][xx] + interpolation[0][0][img_pad_height-1][xx]);
          interpolation[2][0][img_pad_height-1][xx] = IClip(0, 255,((hv+4)>>3));
          tmp20[img_pad_height-1][xx] = (hv<<3);

          // pos 22
          hv = 5 * (tmp02[0][xx] + tmp02[1][xx])
            - (tmp02[0][xx] + tmp02[2][xx]);
          interpolation[2][2][0][xx] = IClip(0, 255,((hv+32)>>6));
          tmp22[0][xx] = hv;

          hv = 5 * (tmp02[img_pad_height-2][xx] + tmp02[img_pad_height-1][xx])
            - (tmp02[img_pad_height-3][xx] + tmp02[img_pad_height-1][xx]);
          interpolation[2][2][img_pad_height-2][xx] = IClip(0, 255,((hv+32)>>6));
          tmp22[img_pad_height-2][xx] = hv;

          hv = 5 * (tmp02[img_pad_height-1][xx] + tmp02[img_pad_height-1][xx])
            - (tmp02[img_pad_height-2][xx] + tmp02[img_pad_height-1][xx]);
          interpolation[2][2][img_pad_height-1][xx] = IClip(0, 255,((hv+32)>>6));
          tmp22[img_pad_height-1][xx] = hv;
        }

        // x o x o
        // x x x x
        // x o x o
        // x x x x
        for(yy=0; yy<img_pad_height; yy++)
        {
          //        o * * *
          //        o * * *
          //        o * * *
          //        o * * *
          // pos 01 & 21
          qh =   1* (interpolation[0][0][yy][0]<<3) +
            7* (interpolation[0][0][yy][0]<<3) +
            7* tmp02[yy][0] +
            1* (interpolation[0][0][yy][1]<<3);
          interpolation[0][1][yy][0] = IClip(0, 255,((qh+64)>>7));

          qh =   1* tmp20[yy][0] +
            7* tmp20[yy][0] +
            7* tmp22[yy][0] +
            1* tmp20[yy][1];
          interpolation[2][1][yy][0] = IClip(0, 255,((qh+512)>>10));

          // pos 03 & pos 23
          qh =   1* (interpolation[0][0][yy][0]<<3) +
            7* tmp02[yy][0] +
            7* (interpolation[0][0][yy][1]<<3) +
            1* tmp02[yy][1];
          interpolation[0][3][yy][0] = IClip(0, 255,((qh+64)>>7));

          qh =   1* tmp20[yy][0] +
            7* tmp22[yy][0] +
            7* tmp20[yy][1] +
            1* tmp22[yy][1];
          interpolation[2][3][yy][0] = IClip(0, 255,((qh+512)>>10));

          //        * o o *
          //        * o o *
          //        * o o *
          //        * o o *
          for(xx=1; xx<img_pad_width-1; xx++)
          {
            // pos 01
            qh =   1* tmp02[yy][xx-1] +
              7* (interpolation[0][0][yy][xx]<<3) +
              7* tmp02[yy][xx] +
              1* (interpolation[0][0][yy][xx+1]<<3);
            interpolation[0][1][yy][xx] = IClip(0, 255,((qh+64)>>7));

            // pos 21
            qh =   1* tmp22[yy][xx-1] +
              7* tmp20[yy][xx] +
              7* tmp22[yy][xx] +
              1* tmp20[yy][xx+1];
            interpolation[2][1][yy][xx] = IClip(0, 255,((qh+512)>>10));

            // pos 03
            qh =   1* (interpolation[0][0][yy][xx]<<3) +
              7* tmp02[yy][xx] +
              7* (interpolation[0][0][yy][xx+1]<<3) +
              1* tmp02[yy][xx+1];
            interpolation[0][3][yy][xx] = IClip(0, 255,((qh+64)>>7));

            // pos 23
            qh =   1* tmp20[yy][xx] +
              7* tmp22[yy][xx] +
              7* tmp20[yy][xx+1] +
              1* tmp22[yy][xx+1];
            interpolation[2][3][yy][xx] = IClip(0, 255,((qh+512)>>10));
          }

          //        * * * o
          //        * * * o
          //        * * * o
          //        * * * o
          // pos 01 & 21
          qh =   1* tmp02[yy][img_pad_width-2] +
            7* (interpolation[0][0][yy][img_pad_width-1]<<3) +
            7* tmp02[yy][img_pad_width-1] +
            1* tmp02[yy][img_pad_width-1];
          interpolation[0][1][yy][img_pad_width-1] = IClip(0, 255,((qh+64)>>7));

          qh =   1* tmp22[yy][img_pad_width-2] +
            7* tmp20[yy][img_pad_width-1] +
            7* tmp22[yy][img_pad_width-1] +
            1* tmp22[yy][img_pad_width-1];
          interpolation[2][1][yy][img_pad_width-1] = IClip(0, 255,((qh+512)>>10));

          // pos 03 & pos 23
          qh =   1* (interpolation[0][0][yy][img_pad_width-1]<<3) +
            7* tmp02[yy][img_pad_width-1] +
            7* tmp02[yy][img_pad_width-1] +
            1* tmp02[yy][img_pad_width-1];
          interpolation[0][3][yy][img_pad_width-1] = (IClip(0, 255,(qh+64)>>7));

          qh =   1* tmp20[yy][img_pad_width-1] +
            7* tmp22[yy][img_pad_width-1] +
            7* tmp22[yy][img_pad_width-1] +
            1* tmp22[yy][img_pad_width-1];
          interpolation[2][3][yy][img_pad_width-1] = IClip(0, 255,((qh+512)>>10));
        }

        // x x x x
        // o x o x
        // x x x x
        // o x o x
        //        o o o o
        //        * * * *
        //        * * * *
        //        * * * *
        for(xx=0; xx<img_pad_width; xx++)
        {
          // pos 10 & 12
          qv =   1* (interpolation[0][0][0][xx]<<6) +
            7* (interpolation[0][0][0][xx]<<6) +
            7* tmp20[0][xx] +
            1* (interpolation[0][0][1][xx]<<6);
          interpolation[1][0][0][xx] = IClip(0, 255,((qv+512)>>10));

          qv =   1* (tmp02[0][xx]<<3) +
            7* (tmp02[0][xx]<<3) +
            7* tmp22[0][xx] +
            1* (tmp02[1][xx]<<3);
          interpolation[1][2][0][xx] = IClip(0, 255,((qv+512)>>10));

          // pos 30 & pos 32
          qv =   1* (interpolation[0][0][0][xx]<<6) +
            7* tmp20[0][xx] +
            7* (interpolation[0][0][1][xx]<<6) +
            1* tmp20[1][xx];
          interpolation[3][0][0][xx] = IClip(0, 255,((qv+512)>>10));

          qv =   1* (tmp02[0][xx]<<3) +
            7* tmp22[0][xx] +
            7* (tmp02[1][xx]<<3) +
            1* tmp22[1][xx];
          interpolation[3][2][0][xx] = IClip(0, 255,((qv+512)>>10));
        }

        for(yy=1; yy<img_pad_height-1; yy++)
        {


          //        * * * *
          //        o o o o
          //        o o o o
          //        * * * *
          for(xx=0; xx<img_pad_width; xx++)
          {
            // pos 10
            qv =   1* tmp20[yy-1][xx] +
              7* (interpolation[0][0][yy][xx]<<6) +
              7* tmp20[yy][xx] +
              1* (interpolation[0][0][yy+1][xx]<<6);
            interpolation[1][0][yy][xx] = IClip(0, 255,((qv+512)>>10));

            // pos 12
            qv =   1* tmp22[yy-1][xx] +
              7* (tmp02[yy][xx]<<3) +
              7* tmp22[yy][xx] +
              1* (tmp02[yy+1][xx]<<3);
            interpolation[1][2][yy][xx] = IClip(0, 255,((qv+512)>>10));

            // pos 30
            qv =   1* (interpolation[0][0][yy][xx]<<6) +
              7* tmp20[yy][xx] +
              7* (interpolation[0][0][yy+1][xx]<<6) +
              1* tmp20[yy+1][xx];
            interpolation[3][0][yy][xx] = IClip(0, 255,((qv+512)>>10));

            // pos 32
            qv =   1* (tmp02[yy][xx]<<3) +
              7* tmp22[yy][xx] +
              7* (tmp02[yy+1][xx]<<3) +
              1* tmp22[yy+1][xx];
            interpolation[3][2][yy][xx] = IClip(0, 255,((qv+512)>>10));
          }
          //        * * * *
          //        * * * *
          //        * * * *
          //        o o o o
          for(xx=0; xx<img_pad_width; xx++)
          {
            // pos 10 & 12
            qv =   1* tmp20[img_pad_height-2][xx] +
              7* (interpolation[0][0][img_pad_height-1][xx]<<6) +
              7* tmp20[img_pad_height-1][xx] +
              1* tmp20[img_pad_height-1][xx];
            interpolation[1][0][img_pad_height-1][xx] = IClip(0, 255,((qv+512)>>10));

            qv =   1* tmp22[img_pad_height-2][xx] +
              7* (tmp02[img_pad_height-1][xx]<<3) +
              7* tmp22[img_pad_height-1][xx] +
              1* tmp22[img_pad_height-1][xx];
            interpolation[1][2][img_pad_height-1][xx] = IClip(0, 255,((qv+512)>>10));

            // pos 30 & 32
            qv =   1* (interpolation[0][0][img_pad_height-1][xx]<<6) +
              7* tmp20[img_pad_height-1][xx] +
              7* tmp20[img_pad_height-1][xx] +
              1* tmp20[img_pad_height-1][xx];
            interpolation[3][0][img_pad_height-1][xx] = IClip(0, 255,((qv+512)>>10));

            qv =   1* (tmp02[img_pad_height-1][xx]<<3) +
              7* tmp22[img_pad_height-1][xx] +
              7* tmp22[img_pad_height-1][xx] +
              1* tmp22[img_pad_height-1][xx];
            interpolation[3][2][img_pad_height-1][xx] = IClip(0, 255,((qv+512)>>10));
          }
        }
        // x x x x
        // x o x o
        // x x x x
        // x o x o
        //        o o o *
        //        o o o *
        //        o o o *
        //        * * * *
        for(yy=0; yy<img_pad_height-1; yy++)
          for(xx=0; xx<img_pad_width-1; xx++)
          {
            // "\"
            // pos 11 & 33
            interpolation[1][1][yy][xx]=IClip(0, 255,(((interpolation[0][0][yy][xx]<<6)+tmp22[yy][xx]+64)>>7));
            interpolation[3][3][yy][xx]=IClip(0, 255,(((interpolation[0][0][yy+1][xx+1]<<6)+tmp22[yy][xx]+64)>>7));

            // "/"
            // pos 13 & 31
            interpolation[1][3][yy][xx]=IClip(0, 255,(((interpolation[0][0][yy][xx+1]<<6)+tmp22[yy][xx]+64)>>7));
            interpolation[3][1][yy][xx]=IClip(0, 255,(((interpolation[0][0][yy+1][xx]<<6)+tmp22[yy][xx]+64)>>7));
          }
          //        o o o o
          //        o o o o
          //        o o o o
          //        * * * *
          for(yy=0; yy<img_pad_height-1; yy++)
          {
            // pos 11 & 31
            interpolation[1][1][yy][img_pad_width-1]=IClip(0, 255,(((interpolation[0][0][yy][img_pad_width-1]<<6)+tmp22[yy][img_pad_width-1]+64)>>7));
            interpolation[3][1][yy][img_pad_width-1]=IClip(0, 255,(((interpolation[0][0][yy+1][img_pad_width-1]<<6)+tmp22[yy][img_pad_width-1]+64)>>7));

            // pos 13 & 33
            interpolation[1][3][yy][img_pad_width-1]=interpolation[1][1][yy][img_pad_width-1];
            interpolation[3][3][yy][img_pad_width-1]=interpolation[3][1][yy][img_pad_width-1];

          }
          //        o o o o
          //        o o o o
          //        o o o o
          //        o o o *
          for(xx=0; xx<img_pad_width-1; xx++)
          {
            // pos 11 & 13
            interpolation[1][1][img_pad_height-1][xx]=IClip(0, 255,(((interpolation[0][0][img_pad_height-1][xx]<<6)+tmp22[img_pad_height-1][xx]+64)>>7));
            interpolation[1][3][img_pad_height-1][xx]=IClip(0, 255,(((interpolation[0][0][img_pad_height-1][xx+1]<<6)+tmp22[img_pad_height-1][xx]+64)>>7));

            //pos  31 & 33
            interpolation[3][1][img_pad_height-1][xx]=interpolation[1][1][img_pad_height-1][xx];
            interpolation[3][3][img_pad_height-1][xx]=interpolation[1][3][img_pad_height-1][xx];
          }
          //        o o o o
          //        o o o o
          //        o o o o
          //        o o o o
          {
            interpolation[1][1][img_pad_height-1][img_pad_width-1]=IClip(0, 255,(((interpolation[0][0][img_pad_height-1][img_pad_width-1]<<6)+tmp22[img_pad_height-1][img_pad_width-1]+64)>>7));
            interpolation[1][3][img_pad_height-1][img_pad_width-1]=interpolation[1][1][img_pad_height-1][img_pad_width-1];
            interpolation[3][1][img_pad_height-1][img_pad_width-1]=interpolation[1][1][img_pad_height-1][img_pad_width-1];
            interpolation[3][3][img_pad_height-1][img_pad_width-1]=interpolation[1][1][img_pad_height-1][img_pad_width-1];
          }
          /*  if(frame_no==6)
          {
          printf("interpolation:\n");
          printf("%3d\n",interpolation[2][1][192][367]);

          printf("%3d ",tmp22[192][img_pad_width-2]);
          printf("%3d ",tmp20[192][img_pad_width-1]);
          printf("%3d \n",tmp22[192][img_pad_width-1]);

          printf("%3d ",interpolation[2][2][192][img_pad_width-2]);
          printf("%3d ",interpolation[2][0][192][img_pad_width-1]);
          printf("%3d ",interpolation[2][2][192][img_pad_width-1]);

          for(yy=0;yy<8;yy++)
          {
          for(xx=0;xx<8;xx++)
          {
          //printf("%3d ",img->mpr[jj+mb_y][ii+mb_x]);
          //printf("%d ",interpolation[2][1][192+yy][360+xx]);

          }
          printf("\n");
          }
          }*/
          // for test
          /*for(posy=0;posy<4;posy++)
          for(posx=0;posx<4;posx++)
          for(yy=0; yy<img_pad_height; yy++)
          for(xx=0; xx<img_pad_width; xx++)
          if(interpolation[posy][posx][yy][xx] != mref[0][4*yy+posy][4*xx+posx])
          temp++;*/

}
