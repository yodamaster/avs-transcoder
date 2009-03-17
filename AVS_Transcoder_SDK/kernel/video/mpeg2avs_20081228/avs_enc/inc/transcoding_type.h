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

#define MAX_GOP_LENGTH 100
#ifndef _TRANSCODING_TYPE_H_
#define _TRANSCODING_TYPE_H_
#include "typedef.h"
typedef struct tag_mb_info
  {
  int_32_t  mb_type;  /* 0: skip, 1: intra, 2:intra */

  int_32_t  pdir;    /* 0: foreward, 1: backward, 2: bidirect, -1: intra */

  int_32_t  mc_type;  /* must be MC_FRAME
              * 2;
              * */
  int_32_t  mv[2][2][2];  /* mpeg2mv[mba][r][s][t] */
  } MB_INFO;
typedef struct tag_mpeg2_dec_create_t
{
  int_32_t  version;
  int_32_t  framerate;    /* [out:opt] frame rate */
  int_32_t  width;      /* [out:opt] image width */
  int_32_t  height;      /* [out:opt] image width */
  unsigned char  *pBitstream;    /* [in] n* 2048 byte buffer */
} mpeg2_dec_create_t;

typedef struct tag_mpeg2_dec_frame_t
{
  int_32_t  version;
  int_16_t    frame_type;
  void    *bitstream;    /* [in] n * 2048 bytes bitstream (read from) */
  int_32_t  length;      /* [in] bitstream length */
  unsigned char  *output;    /* [in] output image (written to) */

  MB_INFO    *MbInfo;    /* [out] mb info array */
} mpeg2_dec_frame_t;

typedef struct tag_mpeg2_dec_stats_t
{
  int_32_t  version;

  int_32_t  framerate;    /* [out:opt] frame rate */
  int_32_t  type;      /* [out] mpeg2_picture_coding_type 1:intra, 2: P, 3: B */

  int_32_t  Bitstream_Framenum;  /* [out] output frame number */  
} mpeg2_dec_stats_t;

typedef struct tag_trans_frame_t
{
  unsigned char  *pMpeg2Stream;
  int_32_t  nMpeg2Ptr;

  unsigned char  *pAvsStream;
  int_32_t  nAvsPtr;
} trans_frame_t;

typedef struct tag_avs_enc_create_t
{
  char    strConfigFile[256];

  int_32_t  frame_rate_code;  /* [in:opt] frame rate code (1: 24000/1001,2: 24,3: 25,4: 30000/1001,5:30,6: 50,7: 60000/1001,8: 60) */
  int_32_t  width;      /* [in:opt] image width */
  int_32_t  height;      /* [in:opt] image width */
} avs_enc_create_t;

/*
 -----------------------------------------------------------------------------------------------------------------------
        The structure
 -----------------------------------------------------------------------------------------------------------------------
 */
typedef struct tag_avs_enc_frame_t
{
  int_32_t  version;

  unsigned char  *input;      /* [in] input image (read from) */
  unsigned char  **inputfrm;
#ifdef  _ME_FOR_RATE_CONTROL_
  unsigned char *pre_frm;
#endif
  int_16_t  type[MAX_GOP_LENGTH];      /* [in:opt] coding type */

  int_32_t  lastFrame;    /* [in: opt] last frame */

  void    *bitstream;    /* [in:opt] bitstream ptr (written to) */
  int_32_t  length;      /* [in:opt] bitstream length (bytes) */

  int_32_t  out_flags;    /* [out] bitstream output flags */

  MB_INFO    *pEncMBInfo;    /* [in:opt] mb info for transcoding, pEncMbInfo[0.. height*width/256 */
} avs_enc_frame_t;

typedef struct tag_avs_enc_stats_t
{
  int_32_t  version;

  /* encoding parameters */
  int_32_t  type;      /* [out] coding type */

  /* bitrate */
  int_32_t  length;      /* [out] frame length */

  int_32_t  hlength;    /* [out] header length (bytes) */
  int_32_t  kblks;      /* [out] number of blocks compressed as Intra */
  int_32_t  mblks;      /* [out] number of blocks compressed as Inter */
  int_32_t  ublks;      /* [out] number of blocks marked as not_coded */
} avs_enc_stats_t;

typedef struct tag_trans_create_t
{
  int_32_t    version;
  unsigned char    *pBitstream;  /* [in] n* 2048 byte buffer */
  int_32_t    framerate;  /* [out:opt] frame rate */
  int_32_t    width;    /* [out:opt] image width */
  int_32_t    height;    /* [out:opt] image width */
  mpeg2_dec_create_t  decCreate;
  mpeg2_dec_frame_t  decFrame;
  mpeg2_dec_stats_t  decStats;
  MB_INFO      *pMbInfo;
  unsigned char    *pYUVData;

  /* encoder parameters */
  avs_enc_create_t  avsCreate;
  avs_enc_frame_t    avsFrame;
  avs_enc_stats_t    avsStats;
  trans_frame_t    *frame;
} trans_create_t;
#endif
