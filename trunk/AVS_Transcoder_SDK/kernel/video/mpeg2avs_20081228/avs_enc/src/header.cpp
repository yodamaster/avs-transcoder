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
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "defines.h"
#include "global.h"

// A little trick to avoid those horrible #if TRACE all over the source code
#if TRACE
#define SYMTRACESTRING(s) strncpy(sym->tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // do nothing
#endif

/*
*************************************************************************
* Function:
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::frametotc(int_32_t frame, int_32_t dropflag)
{
  int_32_t fps, pict, sec, minute, hour, tc;
  
  fps = (int_32_t)frame_rate;                 // modified by wangyue
  //fps = (int_32_t)(frame_rate+0.5);
  pict = frame%fps;
  frame = (frame-pict)/fps;
  sec = frame%60;
  frame = (frame-sec)/60;
  minute = frame%60;
  frame = (frame-minute)/60;
  hour = frame%24;
  tc = (dropflag<<23) | (hour<<18) | (minute<<12) | (sec<<6) | pict;
  
  return tc;
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
int_32_t c_avs_enc::IPictureHeader(int_32_t frame)
{
  Bitstream *bitstream = currBitStream;
  int_32_t len = 0;

  //rm52j
  //len += u_v(32, "video edit code", 0x01B7, bitstream);
  //len += u_v(32, "start code", 0xB7, bitstream);
  len += u_v(32, "I picture start code", 0x1B3, bitstream);

  //rate control
  // modified by wangyue
  if(input->RCEnable && img->BasicUnit==img->Frame_Total_Number_MB)
    input->fixed_picture_qp = 1;
  else
    input->fixed_picture_qp = 1;

  len += u_v(16,"bbv delay",0xffff,bitstream);

  len += u_v(1, "time_code_flag",0,bitstream);

  //rm52j marker bit
  len += u_v(1, "marker bit", 1, bitstream);
  //rm52j
  len += u_v(8,"picture distance", picture_distance, bitstream);
  //rm52c
  //len += u_v(8,"picture distance",picture_distance,bitstream);


  len+=u_v(1,"progressive frame",img->progressive_frame,bitstream);
  if(!img->progressive_frame)
  len+=u_v(1,"picture_structure",img->picture_structure,bitstream);

  len+=u_v(1,"top field first",input->top_field_first,bitstream);
  len+=u_v(1,"repeat first field",input->repeat_first_field,bitstream);
  len+=u_v(1,"fixed picture qp",input->fixed_picture_qp,bitstream);
  //rate control
  if(input->RCEnable)
    len+=u_v(6,"I picture QP",img->qp,bitstream);
  else
    {
    len+=u_v(6,"I picture QP",input->qp0,bitstream);  
    img->qp = input->qp0;
    }

  if(!img->picture_structure )
    {
    len+=u_v(1,"skip mode flag",input->skip_mode_flag,bitstream);
    }

  len+=u_v(4,"reserved bits",0,bitstream);

  len+=u_v(1,"loop filter disable",input->loop_filter_disable,bitstream);
  if (!input->loop_filter_disable)
    {
    len+=u_v(1,"loop filter parameter flag",input->loop_filter_parameter_flag,bitstream);
    if (input->loop_filter_parameter_flag)
      {
      len+=se_v("alpha offset",input->alpha_c_offset,bitstream);
      len+=se_v("beta offset",input->beta_offset,bitstream);
      }
    }

  picture_reference_flag  = 0;

  return len;
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

int_32_t c_avs_enc::PBPictureHeader()
{
  Bitstream *bitstream = currBitStream;
  int_32_t len = 0;

  if (img->type == INTER_IMG)
    picture_coding_type = 1;
  else
    picture_coding_type = 2;

  //rate control
  if(input->RCEnable&&img->BasicUnit==img->Frame_Total_Number_MB)
    input->fixed_picture_qp = 1;
  else
    input->fixed_picture_qp = 1;

  if (img->nb_references==1)
    picture_reference_flag = 1;
  else if (img->nb_references>1)
    picture_reference_flag = 0;

  len+=u_v(24,"start code prefix",1,bitstream);
  len+=u_v(8, "PB picture start code",0xB6,bitstream);
  len+=u_v(16,"bbv delay",0xffff,bitstream);
  len+=u_v(2,"picture coding type",picture_coding_type,bitstream);

  //rm52j
  len += u_v(8,"picture distance", picture_distance, bitstream);

  //rm52c
  //len+=u_v(8,"picture_distance",picture_distance,bitstream);

  len+=u_v(1,"progressive frame",img->progressive_frame,bitstream);
  if (!img->progressive_frame)
    {
    len+=u_v(1,"picture_structure",img->picture_structure,bitstream);
    if (!img->picture_structure)
      len+=u_v(1,"advanced_pred_mode_disable",img->advanced_pred_mode_disable,bitstream);
    }
  len+=u_v(1,"top field first",input->top_field_first,bitstream);
  len+=u_v(1,"repeat first field",input->repeat_first_field,bitstream);
  len+=u_v(1,"fixed qp",input->fixed_picture_qp,bitstream);
  //rate control
  if(img->type==INTER_IMG)
    {
    if(input->RCEnable)
      len+=u_v(6,"I picture QP",img->qp,bitstream);
    else
      {
      len+=u_v(6,"I picture QP",input->qpN,bitstream);
      img->qp=input->qpN;
      }
    }
  else if(img->type==B_IMG)
    {
    if(input->RCEnable)
      len+=u_v(6,"I picture QP",img->qp,bitstream);
    else
      {
      len+=u_v(6,"I picture QP",input->qpB,bitstream);
      img->qp=input->qpB;
      }
    }

  if (!(picture_coding_type == 2 && img->picture_structure==1))
    {
    len+=u_v(1,"piture reference flag",picture_reference_flag,bitstream);
    }

  //rm52j
  len += u_v(1, "no forward reference flag", 0, bitstream);     // no forward reference flag
  len += u_v(3, "reserved bits", 0, bitstream);                 // reserved bits

  //rm52c
  //len+=u_v(4,"reserved bits",0,bitstream);
  len+=u_v(1,"skip mode flag",input->skip_mode_flag, bitstream);

  len+=u_v(1,"loop filter disable",input->loop_filter_disable,bitstream);
  if (!input->loop_filter_disable)
    {
    len+=u_v(1,"loop filter parameter flag",input->loop_filter_parameter_flag,bitstream);
    if (input->loop_filter_parameter_flag)
      {
      len+=se_v("alpha offset",input->alpha_c_offset,bitstream);
      len+=se_v("beta offset",input->beta_offset,bitstream);
      }
    }
  return len;
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

int_32_t c_avs_enc::SliceHeader(int_32_t slice_nr, int_32_t slice_qp)
{
  Bitstream *bitstream = currBitStream;
  int_32_t i;
  int_32_t len = 0;
  
  len+=u_v(24,"start code prefix",1,bitstream);
  //rm52j
  /*len += u_v(8, "slice vertical position", 0, bitstream); */
  len+=u_v(8, "slice vertical position",slice_nr,bitstream);

  if(img->type != INTRA_IMG)
    {
    len += u_v(1,"picutre weighting flag",img->LumVarFlag,bitstream);
    if(img->LumVarFlag)
      {
      for(i=0;i<img->nb_references;i++)
        {
        len+=u_v(8,"luma scale",img->lum_scale[i],bitstream);
        len+=u_v(8,"luma shift",img->lum_shift[i]+127,bitstream);
        u_1 ("insert bit", 1, bitstream);
        len+=u_v(8,"chroma scale",img->chroma_scale[i],bitstream);
        len+=u_v(8,"chroma shift",img->chroma_shift[i],bitstream);
        }
      for( ; i < img->buf_cycle ; i++){
        len+=u_v(8,"luma scale",0,bitstream);
        len+=u_v(8,"luma shift",0,bitstream);//u_i
        u_1 ("insert bit", 1, bitstream);
        len+=u_v(8,"chroma scale",0,bitstream);
        len+=u_v(8,"chroma shift",0,bitstream);//u_i
        }
      len+=u_v(1,"mb weighting flag",img->allframeweight,bitstream);
      }
    }

  return len;
}