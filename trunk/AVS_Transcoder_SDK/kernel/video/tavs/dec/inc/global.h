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


#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <stdio.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <Windows.h>

#include "defines.h"
#include "typedef.h"
#include "transcoding_type.h"

#ifdef WIN32
#pragma warning(disable : 4996)
#endif
#ifdef ROI_ENABLE
#pragma comment(lib,"../bin/TextDet")
#endif

class  __declspec(dllexport)  c_avs_enc;

typedef int_32_t (*LocateTextRegions_t)(byte* YCbCr[3],int_32_t w[3],int_32_t h[3],byte* ROIArray);

/* ! myboolean Type */
typedef enum { myfalse, mytrue } myboolean;

typedef enum { FRAME_CODING, FIELD_CODING, PAFF_CODING } CodingType;

/* ! definition of H.26L syntax elements */
typedef enum
{
  SE_HEADER,
  SE_PTYPE,
  SE_MBTYPE,
  SE_REFFRAME,
  SE_INTRAPREDMODE,
  SE_MVD,
  SE_CBP_INTRA,
  SE_LUM_DC_INTRA,
  SE_CHR_DC_INTRA,
  SE_LUM_AC_INTRA,
  SE_CHR_AC_INTRA,
  SE_CBP_INTER,
  SE_LUM_DC_INTER,
  SE_CHR_DC_INTER,
  SE_LUM_AC_INTER,
  SE_CHR_AC_INTER,
  SE_DELTA_QUANT_INTER,
  SE_DELTA_QUANT_INTRA,
  SE_BFRAME,
  SE_EOS,
  SE_MAX_ELEMENTS /* !< number of maximum syntax elements */
} SE_type;    /* substituting the definitions in elements.h */

typedef enum
{
  BITS_HEADER,
  BITS_TOTAL_MB,
  BITS_MB_MODE,
  BITS_INTER_MB,
  BITS_MVD,  /* mvd 需要的bits */
  BITS_CBP_MB,
  BITS_COEFF_Y_MB,
  BITS_COEFF_UV_MB,
  BITS_DELTA_QUANT_MB,
  MAX_BITCOUNTER_MB
} BitCountType;

typedef enum { FRAME, TOP_FIELD, BOTTOM_FIELD } PictureType;  /* !< New enum for field processing */

/*
-----------------------------------------------------------------------------------------------------------------------
! Syntax element
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_syntaxelement
{
  int_32_t  type;      /* !< type of syntax element for data part. */
  int_32_t  value1;      /* !< numerical value of syntax element */
  int_32_t  value2;      /* !< for blocked symbols, e.g. run/level */
  int_32_t  len;      /* !< length of code */
  int_32_t  inf;      /* !< info part of UVLC code */
  uint_32_t  bitpattern;    /* !< UVLC bitpattern */
  int_32_t  context;    /* !< CABAC context */
  int_32_t  k;      /* !< CABAC context for coeff_count,uv */
  int_32_t  golomb_grad;    /* needed if type is a golomb element (AVS) */
  int_32_t  golomb_maxlevels;  /* if this is zero, do not use the golomb coding. (AVS) */

#if TRACE
#define TRACESTRING_SIZE  100    /* !< size of trace string */
  char    tracestring[TRACESTRING_SIZE];  /* !< trace string */
#endif
  /* !< for mapping of syntaxElement to UVLC */
  void_t (c_avs_enc:: *mapping) (int_32_t value1, int_32_t value2, int_32_t * len_ptr, int_32_t * info_ptr);
} SyntaxElement;

/*
-----------------------------------------------------------------------------------------------------------------------
! Macroblock
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_macroblock
{
  int_32_t    currSEnr;    /* !< number of current syntax element */
  int_32_t    slice_nr;
  int_32_t    delta_qp;
  int_32_t    qp;
  int_32_t    bitcounter[MAX_BITCOUNTER_MB];
  struct tag_macroblock  *mb_available[3][3];  /* !< pointer to neighboring MBs in a 3x3 window of current MB,
                                               * which is located at [1][1] \n NULL pointer identifies
                                               * neighboring MBs which are unavailable */

  /* some storage of macroblock syntax elements for global access */
  int_32_t    mb_type;
  int_32_t    mb_type_2;
  int_32_t    mvd[2][2][2][2];  /* !< indices correspond to [forw,backw][block_y][block_x][x,y] */
  int_32_t    intra_pred_modes[4];
  int_32_t    cbp, scbp;
  int_32_t    cbp_blk;    /* !< 1 bit set for every 4x4 block with coefs (not implemented
                          * for INTRA) */
  int_32_t    b8mode[4];
  int_32_t    b8pdir[4];
  unsigned long   cbp_bits;

  int_32_t    lf_disable;
  int_32_t    lf_alpha_c0_offset;
  int_32_t    lf_beta_offset;

  int_32_t    c_ipred_mode;    /* !< chroma intra prediction mode */
  int_32_t    IntraChromaPredModeFlag;
  int_32_t    mb_field;
  int_32_t    ****cofAC;    /* lgp*dct*modify *//* !< AC coefficients
                                                * [8x8block][4x4block][level/run][scan_pos]
                                                * */
  int_32_t    ****chromacofAC;
  int_32_t    c_ipred_mode_2;    /* !< chroma intra prediction mode */

  /* rate control */
  int_32_t    prev_cbp;
  int_32_t    prev_qp;
  int_32_t    predict_qp;
  int_32_t    predict_error;
} Macroblock;

/*
-----------------------------------------------------------------------------------------------------------------------
! Bitstream
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_bitstream
{
  int_32_t  byte_pos;    /* !< current position in bitstream;* */
  int_32_t  bits_to_go;    /* !< current bitcounter */
  byte    byte_buf;    /* !< current buffer for last written byte */
  int_32_t  stored_byte_pos;  /* !< storage for position in bitstream;* */
  int_32_t  stored_bits_to_go;  /* !< storage for bitcounter */
  byte    stored_byte_buf;  /* !< storage for buffer of last written byte */

  byte    byte_buf_skip;    /* !< current buffer for last written byte */
  int_32_t  byte_pos_skip;    /* !< storage for position in bitstream; */
  int_32_t  bits_to_go_skip;  /* !< storage for bitcounter */

  byte    *streamBuffer;    /* !< actual buffer for written bytes */
} Bitstream;
typedef struct tag_csobj
{
  Bitstream  *bitstream;
  /* syntax element number and bitcounters */
  int_32_t  currSEnr;
  int_32_t  bitcounter[MAX_BITCOUNTER_MB];

  /* elements of current macroblock */
  int_32_t  mvd[2][2][2][2];
  unsigned long  cbp_bits;
} CSobj, *CSptr;

typedef struct tag_picture
{
  int_32_t  no_slices;
  int_32_t  bits_per_picture;
  float    distortion_y;
  float    distortion_u;
  float    distortion_v;
} Picture;

typedef struct tag_copyright
{
  int_32_t  extension_id;
  int_32_t  copyright_flag;
  int_32_t  copyright_id;
  int_32_t  original_or_copy;
  int_32_t  reserved;
  int_32_t  copyright_number;
} CopyRight;

typedef struct tag_CameraParamters
{
  int_32_t  reserved;
  int_32_t  camera_id;
  int_32_t  height_of_image_device;
  int_32_t  focal_length;
  int_32_t  f_number;
  int_32_t  vertical_angle_of_view;
  int_32_t  camera_position_x;
  int_32_t  camera_position_y;
  int_32_t  camera_position_z;
  int_32_t  camera_direction_x;
  int_32_t  camera_direction_y;
  int_32_t  camera_direction_z;
  int_32_t  image_plane_vertical_x;
  int_32_t  image_plane_vertical_y;
  int_32_t  image_plane_vertical_z;
} CameraParamters;
/*
-----------------------------------------------------------------------------------------------------------------------
! SNRParameters
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_SNRParameters
{
  float  snr_y;  /* !< current Y SNR */
  float  snr_u;  /* !< current U SNR */
  float  snr_v;  /* !< current V SNR */
  float  snr_y1; /* !< SNR Y(dB) first frame */
  float  snr_u1; /* !< SNR U(dB) first frame */
  float  snr_v1; /* !< SNR V(dB) first frame */
  float  snr_ya; /* !< Average SNR Y(dB) remaining frames */
  float  snr_ua; /* !< Average SNR U(dB) remaining frames */
  float  snr_va; /* !< Average SNR V(dB) remaining frames */
} SNRParameters;

/*
-----------------------------------------------------------------------------------------------------------------------
! all input parameters
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_InputParameters
{
  int_32_t  no_frames;  /* !< number of frames to be encoded */
  int_32_t  qp0;    /* !< QP of first frame */
  int_32_t  qpN;    /* !< QP of remaining frames */
  int_32_t  thread_num;
  int_32_t  hadamard;  /* !< 0: 'normal' SAD in 1/3 pixel search. 1: use 4x4 Haphazard transform and '
                       * Sum of absolute transform difference' in 1/3 pixel search */
  int_32_t  search_range;  /* !< search range - integer pel search and 16x16 blocks. The search window is
                           * generally around the predicted vector. Max vector is 2xmcrange. For 8x8 and
                           * 4x4 block sizes the search range is 1/2 of that for 16x16 blocks. */
  int_32_t  no_multpred;  /* !< 1: prediction from the last frame only. 2: prediction from the last or
                          * second last frame etc. Maximum 5 frames */
  int_32_t  img_width;  /* !< GH: if CUSTOM image format is chosen, use this size */
  int_32_t  img_height;  /* !< GH: width and height must be a multiple of 16 pels */
  int_32_t  yuv_format;  /* !< GH: YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4,currently only 4:2:0
                         * is supported) */
  int_32_t  color_depth;  /* !< GH: YUV color depth per component in bit/pel (currently only 8 bit/pel is
                          * supported) */
  int_32_t  intra_upd;  /* !< For error robustness. 0: no special action. 1: One GOB/frame is intra
                        * coded as regular 'update'. 2: One GOB every 2 frames is intra coded etc. In
                        * connection with this intra update, restrictions is put on motion vectors to
                        * prevent errors to propagate from the past */
  int_32_t  blc_size[9][2]; /* !< array for different block sizes */
  int_32_t  infile_header;  /* !< If input file has a header set this to the length of the header */
  char    infile[100];  /* !< YUV 4:2:0 input format */
  char    outfile[100];  /* !< H.26L compressed output bitstream */
  char    ReconFile[100]; /* !< Reconstructed Pictures */
  char    TraceFile[100]; /* !< Trace Outputs */
  char   DecRecFile[100];
  int_32_t  intra_period;  /* !< Random Access period though intra */
  int_32_t  GopLength;

  /* B pictures */
  int_32_t  successive_Bframe;  /* !< number of B frames that will be used */
  int_32_t  qpB;      /* !< QP of B frames */
  int_32_t  SequenceHeaderType;

  int_32_t  InterSearch16x16;
  int_32_t  InterSearch16x8;
  int_32_t  InterSearch8x16;
  int_32_t  InterSearch8x8;
  int_32_t  rdopt;

  int_32_t  InterlaceCodingOption;

  /* AVS */
  int_32_t  aspect_ratio_information;
  int_32_t  frame_rate_code;
  float       fr;

  /*
  * int_32_t bit_rate;
  */
  int_32_t  bit_rate_lower;
  int_32_t  bit_rate_upper;
  int_32_t  picture_weighting_flag;
  int_32_t  mb_weighting_flag;

  int_32_t  bbv_buffer_size;
  int_32_t  video_format;
  int_32_t  color_description;
  int_32_t  color_primaries;
  int_32_t  transfer_characteristics;
  int_32_t  matrix_coefficients;
  int_32_t  hour;
  int_32_t  minute;
  int_32_t  second;
  int_32_t  frame_offset;
  int_32_t  profile_id;
  int_32_t  level_id;
  int_32_t  progressive_sequence;
  int_32_t  repeat_first_field;
  int_32_t  top_field_first;
  int_32_t  low_delay;
  int_32_t  chroma_format;
  int_32_t  sample_precision;
  int_32_t  video_range;
  int_32_t  stream_length_flag;
  int_32_t  picture_decoder_order_flag;
  int_32_t  frame_pred_frame_dct;
  int_32_t  progressive_frame;
  int_32_t  fixed_picture_qp;
  int_32_t  time_code_flag;
  int_32_t  display_horizontal_size;
  int_32_t  display_vertical_size;
  int_32_t  dct_adaptive_flag;
  int_32_t  slice_parameter;
  int_32_t  slice_row_nr;
  int_32_t  skip_mode_flag;
  int_32_t  loop_filter_disable;
  int_32_t  loop_filter_parameter_flag;
  int_32_t  alpha_c_offset;
  int_32_t  beta_offset;

  /* ! Rate Control on JVT standard */
  int_32_t  RCEnable;
  int_32_t  bit_rate;
  int_32_t  SeinitialQP;
  int_32_t  basicunit;
  int_32_t  channel_type;
  int_32_t  frame_rate;
  int_32_t  stuff_height;
  int_32_t  stuff_width;
} InputParameters;

/*
-----------------------------------------------------------------------------------------------------------------------
! ImageParameters
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_ImageParameters
{
  int_32_t  number;      /* !< current image number to be encoded */
  int_32_t  lindex;      /* !< next long term index to be used */
  int_32_t  max_lindex;    /* !< max long term index */
  int_32_t  nb_references;
  int_32_t  current_mb_nr;
  int_32_t  total_number_mb;
  int_32_t  current_slice_nr;
  int_32_t  type;
  int_32_t  no_multpred;    /* !< 1: prediction from the last frame only. 2: prediction from the
                            * last or second last frame etc. */
  int_32_t  qp;      /* !< quant for the current frame */
  float       framerate;

  int_32_t  width;      /* !< Number of pels */
  int_32_t  width_cr;    /* !< Number of pels chroma */
  int_32_t  height;      /* !< Number of lines */
  int_32_t  height_cr;    /* !< Number of lines chroma */
  int_32_t  mb_y;      /* !< current MB vertical */
  int_32_t  mb_x;      /* !< current MB horizontal */
  int_32_t  block_y;    /* !< current block vertical */
  int_32_t  block_x;    /* !< current block horizontal */
  int_32_t  pix_y;      /* !< current pixel vertical */
  int_32_t  pix_x;      /* !< current pixel horizontal */
  int_32_t  pix_c_y;    /* !< current pixel chroma vertical */
  int_32_t  block_c_x;    /* !< current block chroma vertical */
  int_32_t  pix_c_x;    /* !< current pixel chroma horizontal */
  int_32_t  **ipredmode;    /* !< GH
                            * ipredmode[90][74];
                            * prediction mode for inter frames */
  int_32_t  cod_counter;    /* !< Current count of number of skipped macroblocks in a row */
  int_32_t  reconflag;
  __declspec(align(16)) byte    mprr[5][16][16];
  char    available_intra_mode[5];

  int_32_t  mprr_2[5][16][16];  /* !< all 4 new intra prediction modes */
  __declspec(align(16)) int_16_t  mprr_c[2][4][8][8];  /* !< new chroma 8x8 intra prediction modes */
  int_32_t  *****mv;    /* !< motion vectors for all block types and all reference frames */
  __declspec(align(16)) int_16_t  mpr[16][16];    /* !< current best prediction mode */
  __declspec(align(16)) int_16_t  m7[16][16];    /* !< the diff pixel values between orginal image and prediction */

  int_32_t  ****chromacofAC;  /* !< AC coefficients [uv][4x4block][level/run][scan_pos] */
  int_16_t  ****cofAC;    /* !< AC coefficients [8x8block][4x4block][level/run][scan_pos] */
  int_16_t  ***cofDC;    /* !< DC coefficients [yuv][level/run][scan_pos] */

  Macroblock  *mb_data;    /* !< array containing all MBs of a whole frame */
  SyntaxElement  MB_SyntaxElements[MAX_SYMBOLS_PER_MB];  /* !< temporal storage for all chosen syntax elements
                                                         * of one MB */

  int_32_t  *quad;      /* !< Array containing square values,used for snr computation */
  int_32_t  **intra_block;

  int_32_t  tr;
  int_32_t  fld_type;    /* !< top or bottom field */
  uint_32_t  fld_flag;
  int_32_t  direct_intraP_ref[4][4];
  int_32_t  imgtr_next_P_frm;
  int_32_t  imgtr_last_P_frm;
  int_32_t  imgtr_next_P_fld;
  int_32_t  imgtr_last_P_fld;
  int_32_t  imgtr_last_prev_P_frm;  /* Lou 1016 */

  /* B pictures */
  int_32_t  b_interval;
  int_32_t  p_interval;
  int_32_t  b_frame_to_code;
  int_32_t  fw_mb_mode;
  int_32_t  bw_mb_mode;
  int_32_t  *****p_fwMV;    /* !< for MVDFW */
  int_32_t  *****p_bwMV;    /* !< for MVDBW */

  int_32_t  *****all_mv;    /* !< replaces local all_mv */
  int_32_t  *****all_bmv;    /* !< replaces local all_mv */
  int_32_t  *****all_bw_omv;  /* !< replaces local all_mv */

  int_32_t  num_ref_pic_active_fwd_minus1;
  int_32_t  num_ref_pic_active_bwd_minus1;

  int_32_t  *****mv_fld;
  int_32_t  *****p_fwMV_fld;
  int_32_t  *****p_bwMV_fld;
  int_32_t  *****all_mv_fld;
  int_32_t  *****all_bmv_fld;

  int_32_t  field_mb_y;        /* Macroblock number of a field MB */
  int_32_t  field_block_y;    /* Vertical block number for the first block of a field MB */
  int_32_t  field_pix_y;    /* Co-ordinates of current macroblock in terms of field pixels (luma) */
  int_32_t  field_pix_c_y;    /* Co-ordinates of current macroblock in terms of field pixels (chroma) */
  int_32_t  *****mv_top;    /* !< For MB level field/frame coding tools */
  int_32_t  *****mv_bot;    /* !< For MB level field/frame coding tools */
  int_32_t  *****p_fwMV_top;  /* !< For MB level field/frame coding tools */
  int_32_t  *****p_fwMV_bot;  /* !< For MB level field/frame coding tools */
  int_32_t  *****p_bwMV_top;  /* !< For MB level field/frame coding tools */
  int_32_t  *****p_bwMV_bot;  /* !< For MB level field/frame coding tools */
  int_32_t  *****all_mv_top;  /* !< For MB level field/frame coding tools */
  int_32_t  *****all_mv_bot;  /* !< For MB level field/frame coding tools */
  int_32_t  *****all_bmv_top;  /* !< For MB level field/frame coding tools */
  int_32_t  *****all_bmv_bot;  /* !< For MB level field/frame coding tools */
  int_32_t  **ipredmode_top;  /* !< For MB level field/frame coding tools */
  int_32_t  **ipredmode_bot;  /* !< For MB level field/frame coding tools */
  int_32_t  field_mode;    /* !< For MB level field/frame -- field mode on flag */
  int_32_t  top_field;    /* !< For MB level field/frame -- top field flag */

  int_32_t  buf_cycle;

  uint_32_t  frame_num;    /* frame_num for this frame */

  /* the following are sent in the slice header */
  int_32_t  NoResidueDirect;

  int_32_t  coding_stage;
  int_32_t  block8_x;
  int_32_t  block8_y;
  int_32_t  coded_mb_nr;

  int_32_t  *****omv;
  int_32_t  *****all_omv;    /* !< replaces local all_mv */
  int_32_t  *****omv_fld;
  int_32_t  *****all_omv_fld;  /* !< replaces local all_mv */

  int_32_t  current_slice_start_mb;
  int_32_t  current_slice_qp;
  int_32_t  progressive_frame;
  int_32_t  picture_structure;
  int_32_t  dropflag;
  int_32_t  advanced_pred_mode_disable;
  int_32_t  old_type;
  int_32_t  current_mb_nr_fld;

  /* !! for weighting prediction */
  int_32_t  LumVarFlag;
  int_32_t  lum_scale[4];
  int_32_t  lum_shift[4];
  int_32_t  chroma_scale[4];
  int_32_t  chroma_shift[4];
  int_32_t  allframeweight;
  int_32_t  mbweightflag;
  int_32_t  mpr_weight[16][16];
  int_32_t  top_bot;    /* -1: frame / 0: top field / 1: bottom field / Yulj 2004.07.14 */
  int_32_t  mb_no_currSliceLastMB;  /* the last MB no in current slice. Yulj 2004.07.15 */

  /* rate control */
  int_32_t  NumberofHeaderBits;
  int_32_t  NumberofTextureBits;
  int_32_t  NumberofBasicUnitHeaderBits;
  int_32_t  NumberofBasicUnitTextureBits;
  double    TotalMADBasicUnit;
  int_32_t  NumberofMBTextureBits;
  int_32_t  NumberofMBHeaderBits;
  int_32_t  NumberofCodedBFrame;
  int_32_t  NumberofCodedPFrame;
  int_32_t  NumberofGOP;
  int_32_t  TotalQpforPPicture;
  int_32_t  NumberofPPicture;
  double    MADofMB[10000];
  int_32_t  BasicUnitQP;
  int_32_t  TopFieldFlag;
  int_32_t  FieldControl;
  int_32_t  FieldFrame;
  int_32_t  IFLAG;
  int_32_t  NumberofCodedMacroBlocks;
  int_32_t  BasicUnit;
  int_32_t  bot_MB;
  int_32_t  img_width_in_mb;
  int_32_t  img_height_in_mb;
} ImageParameters;

/*
-----------------------------------------------------------------------------------------------------------------------
!< statistics
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_StatParameters
{
  int_32_t  quant0;    /* !< quant for the first frame */
  int_32_t  quant1;    /* !< average quant for the remaining frames */
  float    bitr;    /* !< bit rate for current frame, used only for output til terminal */
  float    bitr0;    /* !< stored bit rate for the first frame */
  float    bitrate;  /* !< average bit rate for the sequence except first frame */
  uint_32_t  bit_ctr;  /* !< counter for bit usage */
  uint_32_t  bit_ctr_0;  /* !< stored bit use for the first frame */
  uint_32_t  bit_ctr_n;  /* !< bit usage for the current frame */
  uint_32_t  bit_slice;  /* !< number of bits in current slice */
  uint_32_t  bit_use_mode_inter[2][MAXMODE]; /* !< statistics of bit usage */
  uint_32_t  bit_ctr_emulationprevention;  /* !< stored bits needed to prevent start code emulation */
  int_32_t  mode_use_intra[25];    /* !< Macroblock mode usage for Intra frames */
  int_32_t  mode_use_inter[2][MAXMODE];

  int_32_t  mb_use_mode[2];
  int_32_t  sequence_header;
  /* B pictures */
  int_32_t  *mode_use_Bframe;
  int_32_t  *bit_use_mode_Bframe;
  uint_32_t  bit_ctr_P;
  uint_32_t  bit_ctr_B;
  float    bitrate_P;
  float    bitrate_B;

  uint_32_t  bit_use_stuffingBits[NUM_PIC_TYPE];
  uint_32_t  bit_use_mb_type[NUM_PIC_TYPE];
  uint_32_t  bit_use_header[NUM_PIC_TYPE];
  uint_32_t  tmp_bit_use_cbp[NUM_PIC_TYPE];
  uint_32_t  bit_use_coeffY[NUM_PIC_TYPE];
  uint_32_t  bit_use_coeffC[NUM_PIC_TYPE];
  uint_32_t  bit_use_delta_quant[NUM_PIC_TYPE];

  int_32_t  em_prev_bits_frm;
  uint_32_t  em_prev_bits_fld;
  int_32_t  *em_prev_bits;
  uint_32_t  bit_ctr_parametersets;
} StatParameters;

/*
-----------------------------------------------------------------------------------------------------------------------
!< For MB level field/frame coding tools ;
!< temporary structure to store MB data for field/frame coding
-----------------------------------------------------------------------------------------------------------------------
*/
typedef struct tag_RD_DATA
{
  double    min_rdcost;
  int_32_t  tmp_mv[2][2][2];    /* to hold the motion vectors for each block */
  int_32_t  tmp_fwMV[2][2][2];    /* to hold forward motion vectors for B MB's */
  int_32_t  tmp_bwMV[2][2][2];    /* to hold backward motion vectors for B MBs */
  int_32_t  dfMV[2][2][2], dbMV[2][2][2];  /* to hold direct motion vectors for B MB's */
  int_32_t  rec_mbY[16][16];    /* hold the Y component of reconstructed MB */
  int_32_t  rec_mbU[8][8], rec_mbV[8][8];
  int_32_t  ****cofAC;
  int_32_t  ****chromacofAC;
  int_32_t  ***cofDC;
  int_32_t  mvd[2][2][2][2];
  int_32_t  mb_type;
  int_32_t  b8mode[4], b8pdir[4];
  int_32_t  frefar[2][2], brefar[2][2];
  int_32_t  **ipredmode;
  int_32_t  intra_pred_modes[4];
  int_32_t  cbp, cbp_blk;
  int_32_t  mode;
  int_32_t  *****mv, *****p_fwMV, *****p_bwMV;
  int_32_t  *****all_mv;
  int_32_t  *****all_bmv;
  int_32_t  c_ipred_mode;
} RD_DATA;
typedef struct tag_OutputStream
{
  FILE    *f;
  byte        buf[STREAM_BUF_SIZE];
  uint_32_t  uPreBytes;      /* 最近写入的3个字节，初始值是0xFFFFFFFF */
  int_32_t  iBytePosition;      /* 当前字节位置 */
  int_32_t  iBitOffset;      /* 当前位偏移，0表示最高位 */
  int_32_t  iNumOfStuffBits;    /* 已插入的填充位的个数，遇到开始码时置0 */
  int_32_t  iBitsCount;      /* 码流总位数 */
} OutputStream;

class __declspec(dllexport) c_avs_enc
{
public:
  c_avs_enc();
  ~c_avs_enc(){};

#ifdef _THREE_STEP_MOTION_SEARCH_
  int_32_t three_step_pattern_x[9], three_step_pattern_y[9];
  void_t init_3_step_search();
  int_32_t TSSMotionSearch(pel_t **orig_pic,int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t pred_mv_x,int_32_t pred_mv_y,int_32_t *mv_x,int_32_t *mv_y,int_32_t search_range,int_32_t min_mcost,double   lambda, int_32_t block_index);
#endif
  pel_t****   interpolation; // [288+16*2][352+16*2] 
  __m128i clip0;
  __m128i clip255;
  __m128i round1;
  __m128i round4;
  __m128i round8;
  __m128i round16;
  __m128i round32;
  __m128i round64;
  __m128i round512;

  Picture      *frame_pic;
  Picture      *top_pic;
  Picture      *bot_pic;

  byte      *imgY_org_buffer;    /* !< Reference luma image */
  byte      *imgY_org_buffer_restore;
  int_32_t    ***tmp_mv_frm;      /* !< motion vector buffer */
  int_32_t    **refFrArr_frm;      /* !< Array for reference frames of each block */
  byte      **imgY;      /* !< Encoded luma images */
  byte      ***imgUV;    /* !< Encoded croma images */
  byte      **imgY_org;    /* !< Reference luma image */
  byte      ***imgUV_org;    /* !< Reference croma image */
  byte      **imgY_pf;    /* !< Post filter luma image */
  byte      ***imgUV_pf;    /* !< Post filter croma image */
#ifdef _OUTPUT_DEC_IMG_
  byte **org_nextP_imgY;
  byte ***org_nextP_imgUV;
#endif
  byte      ****mref[2];
  byte      **mcef[4][2];    /* !< pix chroma */

  int_16_t    **tmp02;    /* !< for quarter pel interpolation */
  int_16_t    **tmp20;
  int_16_t    **tmp22;
  int_32_t    ***tmp_mv;    /* !< motion vector buffer */
  int_32_t    **refFrArr;    /* !< Array for reference frames of each block */

  /* B pictures */
  int_32_t    ***tmp_fwMV;
  int_32_t    ***tmp_bwMV;
  int_32_t    ***tmp_fwMV_top;  /* !< For MB level field/frame coding tools */
  int_32_t    ***tmp_fwMV_bot;  /* !< For MB level field/frame coding tools */
  int_32_t    ***tmp_bwMV_top;  /* !< For MB level field/frame coding tools */
  int_32_t    ***tmp_bwMV_bot;  /* !< For MB level field/frame coding tools */
  int_32_t    **field_mb;    /* !< For MB level field/frame coding tools */
  int_32_t    WriteFrameFieldMBInHeader;  /* ! For MB level field/frame coding tools */
  int_32_t    ***tmp_fwMV_fld;    /* !< For MB level field/frame coding tools */
  int_32_t    ***tmp_bwMV_fld;    /* !< For MB level field/frame coding tools */

  int_32_t    ***dfMV;
  int_32_t    ***dbMV;
  int_32_t    **fw_refFrArr;
  int_32_t    **bw_refFrArr;
  byte      **nextP_imgY;
  byte      ***nextP_imgUV;
  pel_t      *Refbuf11[4];      /* !< 1/1th pel (full pel) reference frame buffer */

  /*
  * global picture format dependend buffers, mem allocation in image.c (field
  * picture)
  */
  byte      **imgY_org_top;
  byte      ***imgUV_org_top;
  byte      **imgY_org_bot;
  byte      ***imgUV_org_bot;
  byte      **imgY_top;      /* !< Encoded luma images */
  byte      ***imgUV_top;      /* !< Encoded croma images */
  byte      **imgY_bot;      /* !< Encoded luma images */
  byte      ***imgUV_bot;      /* !< Encoded croma images */
  pel_t      **Refbuf11_fld;      /* !< 1/1th pel (full pel) reference frame buffer */
  int_32_t    **refFrArr_top;      /* !< Array for reference frames of each block */
  int_32_t    **refFrArr_bot;      /* !< Array for reference frames of each block */
  byte      **imgY_com;      /* !< Encoded luma images */
  byte      ***imgUV_com;      /* !< Encoded croma images */
  int_32_t    **refFrArr_fld;      /* !< Array for reference frames of each block */
  int_32_t    *parity_fld;

  /*
  * global picture format dependend buffers, mem allocation in image.c (field
  * picture)
  */
  byte      **mref_fld[4];      /* !< 1/4 pix luma */
  byte      ***mref_mbfld;      /* !< For MB level field/frame coding tools */

  /*
  * global picture format dependend buffers, mem allocation in image.c (frame
  * buffer)
  */
  byte      ****mref_frm[2];    /* !< 1/4 pix luma //[2:ref_index] */

  /*
  * B pictures ;
  * motion vector : forward, backward, direct
  */
  int_32_t    **fw_refFrArr_top;
  int_32_t    **bw_refFrArr_top;
  int_32_t    **fw_refFrArr_bot;
  int_32_t    **bw_refFrArr_bot;
  int_32_t    ***tmp_mv_top;
  int_32_t    ***tmp_mv_bot;

  int_32_t    **fwdir_refFrArr;    /* !< direct mode forward reference buffer */
  int_32_t    **bwdir_refFrArr;    /* !< direct mode backward reference buffer */

  /*
  * global picture format dependend buffers, mem allocation in image.c (frame
  * buffer)
  */
  byte      **imgY_org_frm;
  byte      ***imgUV_org_frm;
  byte      **imgY_frm;      /* !< Encoded luma images */
  byte      ***imgUV_frm;      /* !< Encoded croma images */
  int_32_t    direct_mode;

  /*
  * B pictures ;
  * motion vector : forward, backward, direct
  */
  int_32_t    **fw_refFrArr_frm;
  int_32_t    **bw_refFrArr_frm;
  int_32_t    **fw_refFrArr_fld;
  int_32_t    **bw_refFrArr_fld;
  int_32_t    ***tmp_mv_fld;      /* !< motion vector buffer */

  int_32_t    intras;      /* !< Counts the intra updates in each frame. */

  int_32_t    Bframe_ctr, frame_no, nextP_tr_frm, nextP_tr,  gframe_no; // gframe_no added by xzhao 20080324;
#ifdef _ME_FOR_RATE_CONTROL_
  int_32_t    glb_me_for_rate_control_flag;
#endif
  int_32_t        frames_GOP;//xzhao 20080331
  int_32_t    tot_time;
  int_32_t    tmp_buf_cycle;
  char      errortext[ET_SIZE];  /* !< buffer for error message for exit with error() */
  RD_DATA      *rdopt;

  int_32_t    *allalpha_lum, *allbelta_lum;

  /* files */
  FILE      *p_org_dec; //org dec img
  FILE      *p_rec;  /* recon img*/
  FILE      *p_stat;      /* !< status file for the last encoding session */
  FILE      *p_log;    /* !< SNR file */
  FILE      *p_in;    /* !< YUV */
  FILE      *p_datpart;  /* !< file to write bitlength and id of all partitions */
  FILE      *p_trace;  /* !< Trace file */
  int_32_t    *refbits;
  int_32_t    ***motion_cost;
  int_32_t    ***motion_cost_bid;
  int_32_t    ipdirect_x, ipdirect_y;

  /* reference frame buffer */
  unsigned char    **reference_frame[3][3];  /* [refnum][yuv][height][width] */
  unsigned char    **reference_field[6][3];  /* [refnum][yuv][height/2][width] */
  unsigned char    ***ref[4];      /* [refnum(4 for filed)][yuv][height(height/2)][width] */
  unsigned char    ***ref_frm[2];      /* [refnum(4 for filed)][yuv][height(height/2)][width] */
  unsigned char    ***ref_fld[6];      /* [refnum(4 for filed)][yuv][height(height/2)][width] */
  unsigned char    ***b_ref[2], ***f_ref[2];
  unsigned char    ***b_ref_frm[2], ***f_ref_frm[2];
  unsigned char    ***b_ref_fld[2], ***f_ref_fld[2];
  unsigned char    ***current_frame;    /* [yuv][height][width] */
  unsigned char    ***current_field;    /* [yuv][height/2][width] */
  MB_INFO          *pAVSMbInfo;
  unsigned char    *pInputImage;
#ifdef _ME_FOR_RATE_CONTROL_
  unsigned char    *pPreImage;
#endif

  int_32_t    bytes_y;
  int_32_t    bytes_uv;
  int_32_t    framesize_in_bytes;
  int_32_t    total_encoded_frame;
  int_32_t    current_encoded_frame;
  int_32_t    output_flag;
  InputParameters   inputs, *input;
  ImageParameters   images, *img;
  SNRParameters   snrs, *snr;
  StatParameters   stats, *stat;
  CopyRight     CopyRights, *cp;
  CameraParamters   CameraParameter, *camera;
  unsigned char    pOutBuffer[AVS_OUT_BUFFER_SIZE];
  Bitstream    *currBitStream;
  int_32_t    nOutBufPtr;
  OutputStream  ORABS, *pORABS;
  double      GBIM_value;
  double      GBIM_value_frm;
  void_t      init_img();
  void_t      report();
  void_t      LumaPrediction4x4(int_32_t, int_32_t, int_32_t, int_32_t, int_32_t, int_32_t);
  int_32_t    SATD(int_16_t *, int_32_t);

  void_t      LumaResidualCoding();
  void_t      ChromaResidualCoding(int_32_t *);
  void_t      IntraChromaPrediction8x8(int_32_t *, int_32_t *, int_32_t *);
  int_32_t    writeMBHeader(int_32_t rdopt);

  int_32_t    LumaResidualCoding8x8(int_32_t *,int_32_t *,int_32_t,int_32_t,int_32_t,int_32_t,int_32_t);
  int_32_t    writeLumaCoeff8x8(int_32_t, int_32_t);
  int_32_t    writeMotionVector8x8(int_32_t,int_32_t,int_32_t,int_32_t,int_32_t,int_32_t,int_32_t,int_32_t);
  int_32_t    writeIntra4x4Modes(int_32_t);
  int_32_t    writeChromaIntraPredMode();
  int_32_t    B8Mode2Value(int_32_t b8mode, int_32_t b8pdir);

  /* dynamic memory allocation */
  int_32_t    init_global_buffers();
  void_t      free_global_buffers();
  void_t      no_mem_exit(char *where);

  int_32_t    get_mem_mv(int_32_t ******);
  void_t      free_mem_mv(int_32_t *****);
  void_t      free_img();

  int_32_t    get_mem_ACcoeff(int_16_t *****);
  int_32_t    get_mem_DCcoeff(int_16_t ****);
  int_32_t    get_mem_ref(byte *****);
  void_t      free_mem_ACcoeff(int_16_t ****);
  void_t      free_mem_DCcoeff(int_16_t ***);
  void_t      free_mem_ref(byte ****);

  void_t      split_field_top();
  void_t      split_field_bot();

  void_t      intrapred_luma_AVS(int_32_t img_x, int_32_t img_y);

  Picture      *malloc_picture();
  void_t      free_picture(Picture *pic);
  int_32_t    encode_one_slice(Picture *pic); /* ! returns the number of MBs in the slice */

  void_t      start_macroblock();
  void_t      set_MB_parameters(int_32_t mb);

  void_t      terminate_macroblock(myboolean *end_of_picture);
  void_t      write_one_macroblock(int_32_t eos_bit);
  void_t      proceed2nextMacroblock();

  void_t      CheckAvailabilityOfNeighbors();

  void_t      free_slice_list(Picture *currPic);

#if TRACE
  void_t      trace2out(SyntaxElement *se);
#endif
  void_t      error(char *text, int_32_t code);
  int_32_t    start_sequence();
  int_32_t    terminate_sequence();
  int_32_t    writeCBPandLumaCoeff();
  int_32_t    writeChromaCoeff();

  void_t      FreeBitstream();
  void_t      AllocateBitstream();
  void_t      PatchInp();
  int_32_t    avs_enc_create();
  int_32_t    avs_enc_destroy();
  int_32_t    avs_enc_encode();
  void_t      Configure(char *av);
  int_32_t    init_global_variables();    
  int_32_t    avs_enc_frame(avs_enc_frame_t *pFrame);
  void_t      information_init();
  int_32_t    write_start_code(OutputStream *p, unsigned char code);
  void_t      CloseBitStreamFile();
  void_t      OpenBitStreamFile(char *Filename);
  int_32_t    WriteSequenceHeader();
  int_32_t    WriteSequenceDisplayExtension();
  int_32_t    WriteUserData(char *userdata);
  int_32_t    WriteSequenceEnd();
  void_t      WriteBitstreamtoFile();
  void_t      WriteSlicetoFile();
  int_32_t    WriteCopyrightExtension();
  int_32_t    WriteCameraParametersExtension();
  void_t      quant_B8(int_32_t qp, int_32_t mode, int_16_t curr_blk[8][8]);
  int_16_t    scanquant_B8(int_32_t qp,int_32_t  mode,int_32_t  b8, int_16_t  curr_blk[8][8], int_32_t  scrFlag, int_32_t  *cbp, int_32_t  *cbp_blk);
  int_16_t    scanquant_B8_recon(int_32_t  qp,int_32_t  mode,int_32_t  b8,int_16_t  curr_blk[8][8],int_32_t  scrFlag,int_32_t  *cbp,int_32_t  *cbp_blk);
  int_16_t    scanquant_B8_cost(int_32_t qp,int_32_t mode,int_32_t b8,int_16_t curr_blk[8][8],int_32_t scrFlag,int_32_t *cbp,int_32_t *cbp_bl);
  int_32_t    find_sad_8x8(int_32_t iMode, int_32_t iSizeX,int_32_t iSizeY, int_32_t iOffX,int_32_t iOffY, int_32_t m7[16][16]);
  int_32_t    sad_hadamard(int_32_t iSizeX, int_32_t iSizeY,int_32_t iOffX, int_32_t iOffY,int_32_t m7[16][16]);
  int_32_t    writeLumaCoeffAVS_B8(int_32_t b8, int_32_t intra);
  int_32_t    writeChromaCoeffAVS_B8(int_32_t b8);
  void_t      idct_transform(int_16_t *mb, int_16_t *temp);
  char      *GetConfigFileContent(char *Filename);
  void_t      ParseContent(char *buf, int_32_t bufsize);
  int_32_t    ParameterNameToMapIndex(char *s);

  /* entropy coding */
  void_t      encode_golomb_word(uint_32_t symbol, uint_32_t grad0,uint_32_t max_levels, uint_32_t *res_bits,uint_32_t *res_leu);      /* returns symbol coded. (might be cropped if max_levels is too * small) */
  void_t      encode_multilayer_golomb_word(uint_32_t  symbol, const uint_32_t *grad,const uint_32_t *max_levels, uint_32_t  *res_bits, uint_32_t  *res_len);      /* terminate using a max_levels value of 30UL. */
  uint_32_t    decode_golomb_word(const unsigned char  **buffer, uint_32_t *bitoff,uint_32_t grad0, uint_32_t max_levels);

  int_32_t    writeSyntaxElement_GOLOMB(SyntaxElement *se, Bitstream *bitstream);
  int_32_t    frametotc(int_32_t frame, int_32_t dropflag);
  int_32_t    SliceHeader(int_32_t slice_nr, int_32_t slice_qp);
  int_32_t    IPictureHeader(int_32_t frame);
  int_32_t    PBPictureHeader();
  void_t      code_a_picture(Picture *frame);
  void_t      ReadOneFrame();
  void_t      write_reconstructed_image();
  int_32_t    writeout_picture();
  int_32_t    writeout_slice();
  void_t      find_snr();
  double      find_GBIM(byte **I);
  void_t      frame_mode_buffer(int_32_t bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v);
  void_t      init_frame();
  void_t      init_field();

  void_t      top_field(Picture *pic);
  void_t      bot_field(Picture *pic);
  void_t      combine_field();

  void_t      put_buffer_frame();
  void_t      put_buffer_top();
  void_t      put_buffer_bot();
  void_t      Update_Picture_Buffers_bot_field();
  void_t      Update_Picture_Buffers_top_field();
  void_t      interpolate_frame_to_fb();

  void_t      CopyFrameToOldImgOrgVariables();

  /*
  * void_t UnifiedOneForthPix (pel_t ** imgY, pel_t ** imgU, pel_t ** imgV,pel_t ** out4Y);
  */
  void_t      UnifiedOneForthPix_sse(pel_t **imgY);
  void_t      UnifiedOneForthPix_c_sse(pel_t **imgY);
  void_t      ReportFirstframe(int_32_t tmp_time);
  void_t      ReportIntra(int_32_t tmp_time);
  void_t      ReportP(int_32_t tmp_time);
  void_t      ReportB(int_32_t tmp_time);

  void_t      CalculateFrameNumber(); /* Calculates the next frame number */

  /* !! weighting prediction */
  void_t      estimate_weighting_factor();
  void_t      find_distortion();
  int_32_t    cdecide_fld_frame (float snr_frame_Y,float snr_field_Y,int_32_t bit_field, int_32_t bit_frame,double lambda_picture);
  void_t      Update_Picture_Bufffers();
  int_32_t    DetectLumVar();
  void_t      CalculateBrightnessPar(int_32_t currentblock[16][16],int_32_t preblock[16][16],float *c,float *d);
  void_t      CalculatePar(int_32_t refnum);
  void_t      LumaPrediction(int_32_t *cbp,int_32_t *cbp_blk,int_32_t block8x8,int_32_t fw_mode,int_32_t bw_mode,int_32_t fw_refframe,int_32_t bw_refframe);
  void_t      SetModesAndRefframe(int_32_t b8,int_32_t *fw_mode,int_32_t *bw_mode,int_32_t *fw_ref,int_32_t *bw_ref);
  void_t      OneComponentChromaPrediction4x4(int_16_t *mpred,int_32_t pix_c_x,int_32_t pix_c_y,int_32_t *****mv,int_32_t ref,int_32_t blocktype,int_32_t uv,int_32_t directforword);
  void_t      OneComponentChromaPrediction4x4_dir(int_16_t *mpred,int_32_t pix_c_x,int_32_t pix_c_y,int_32_t *****mv,int_32_t ref,int_32_t blocktype,int_32_t uv,int_32_t refframe);
  void_t      IntraChromaPrediction4x4(int_32_t uv, int_32_t block_x, int_32_t block_y);
  void_t      ChromaPrediction4x4(int_32_t uv,int_32_t block_x,int_32_t block_y,int_32_t fw_mode,int_32_t bw_mode,int_32_t fw_ref_frame,int_32_t bw_ref_frame);
  int_32_t    SubMBType2Value(Macroblock *currMB, int_32_t layer);
  int_32_t    MBType2Value(Macroblock *currMB);
  int_32_t    writeMotionVector8x8_bid(int_32_t i0,int_32_t j0,int_32_t i1,int_32_t j1,int_32_t refframe,int_32_t dmv_flag,int_32_t fwd_flag,int_32_t mv_mode,int_32_t pdir);
  int_32_t    StoreMotionVector8x8(int_32_t i0,int_32_t j0,int_32_t i1,int_32_t j1,int_32_t refframe,int_32_t dmv_flag,int_32_t fwd_flag,int_32_t mv_mode);
  int_32_t    StoreMotionVector8x8_bid(int_32_t i0,int_32_t j0,int_32_t i1,int_32_t j1,int_32_t refframe,int_32_t dmv_flag,int_32_t fwd_flag,int_32_t mv_mode,int_32_t pdir);
  int_32_t    writeMVD8x8(int_32_t i0,int_32_t j0,int_32_t i1,int_32_t j1,int_32_t refframe,int_32_t dmv_flag,int_32_t fwd_flag,int_32_t mv_mode);
  void_t      writeCBPandDqp(int_32_t *CBPRate);
  int_32_t    writeBlockCoeff(int_32_t block8x8);
  int_32_t    storeMotionInfo(int_32_t pos);
  int_32_t    writeFrameRef(int_32_t mode, int_32_t i, int_32_t j, int_32_t fwd_flag, int_32_t ref);
  void_t      writeReferenceIndex(int_32_t *RefIndexRate);
  void_t      writeMVD(int_32_t *MVDRate);
  void_t      writeweightflag();
  int_32_t    get_mem2D(byte ***array2D, int_32_t rows, int_32_t columns);
  int_32_t    get_mem2Dint(int_32_t ***array2D, int_32_t rows, int_32_t columns);
  int_32_t    get_mem3D(byte ****array2D, int_32_t frames, int_32_t rows, int_32_t columns);
  int_32_t    get_mem3Dint(int_32_t ****array3D, int_32_t frames, int_32_t rows, int_32_t columns);
  int_32_t    get_mem4Dint(int_32_t *****array4D,int_32_t idx,int_32_t frames,int_32_t rows,int_32_t columns);
  int_32_t    get_mem2Dshort_int(int_16_t ***array2D, int_16_t rows, int_16_t columns);

  void_t      free_mem2D(byte **array2D);
  void_t      free_mem2Dshort_int(int_16_t **array2D);
  void_t      free_mem2Dint(int_32_t **array2D);
  void_t      free_mem3D(byte ***array2D, int_32_t frames);
  void_t      free_mem3Dint(int_32_t ***array3D, int_32_t frames);
  void_t      free_mem4Dint(int_32_t ****array4D, int_32_t idx, int_32_t frames);

  /* me */
  int_32_t        only_motion_cost[5][4];
  void_t      Init_Motion_Search_Module();
  void_t      Clear_Motion_Search_Module();
  int_32_t    FullPelBlockMotionSearch(pel_t   **orig_pic,int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t pred_mv_x,int_32_t pred_mv_y,int_32_t *mv_x,int_32_t *mv_y,int_32_t search_range,int_32_t min_mcost,double lambda, int_32_t debug_flag);
  int_32_t    SubPelBlockMotionSearch(pel_t   **orig_pic,int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t pred_mv_x,int_32_t pred_mv_y,int_32_t *mv_x,int_32_t *mv_y,int_32_t search_pos2,int_32_t search_pos4,int_32_t min_mcost,double lambda,int_32_t block_index);
  int_32_t    SubPelBlockMotionSearch_bid(pel_t  **orig_pic,int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t pred_mv_x,int_32_t pred_mv_y,int_32_t *mv_x,int_32_t *mv_y,int_32_t search_pos2,int_32_t search_pos4,int_32_t min_mcost,double lambda, int_32_t block_index);
  int_32_t    BlockMotionSearch(int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t search_range,double lambda, int_32_t block_index);
  int_32_t    GetSkipCostMB(double lambda);
  void_t      FindSkipModeMotionVector();
  void_t      SetMotionVectorPredictor(int_32_t pmv[2],int_32_t **refFrArr,int_32_t ***tmp_mv,int_32_t ref_frame,int_32_t mb_pix_x,int_32_t mb_pix_y,int_32_t blockshape_x,int_32_t blockshape_y,int_32_t ref);
  int_32_t    Get_Skip_CostMB(pel_t   **orig_pic,int_32_t ref,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t blocktype,int_32_t pred_mv_x,int_32_t pred_mv_y,int_32_t *mv_x,int_32_t *mv_y,int_32_t search_pos2,int_32_t search_pos4,int_32_t min_mcost,double lambda);
  int_32_t    Get_Direct_Cost8x8(int_32_t, double);
  int_32_t    Get_Direct_CostMB(double);
  int_32_t    scale_motion_vector(int_32_t motion_vector,int_32_t currblkref,int_32_t neighbourblkref,int_32_t block_y_pos,int_32_t curr_block_y,int_32_t ref);
  int_32_t    calculate_distance(int_32_t blkref, int_32_t fw_bw);
  void_t      PartitionMotionSearch(int_32_t, int_32_t, double);
  void_t      PartitionMotionSearch_bid(int_32_t, int_32_t, double);
  void_t      Get_IP_direct();
  int_32_t    sign(int_32_t a, int_32_t b);

  // loopfilter
  void_t DeblockFrame(ImageParameters *img, byte **imgY, byte ***imgUV);
  void_t DeblockMb(ImageParameters *img, byte **imgY, byte ***imgUV, int_32_t mb_y, int_32_t mb_x);
  void_t GetStrength(byte Strength[2],Macroblock* MbP,Macroblock* MbQ,int_32_t dir,int_32_t edge,int_32_t block_y,int_32_t block_x);
  void_t EdgeLoop(byte* SrcPtr,byte Strength[2],int_32_t QP, int_32_t dir,int_32_t width,int_32_t Chro);

  /* rate control */
  double      THETA;
  int_32_t      Switch;

  int_32_t    Iprev_bits;
  int_32_t    Pprev_bits;

  /* rate control variables */
  int_32_t    Xb;
  int_32_t    R, T_field;
  int_32_t    Np, Nb, bits_topfield, Q;
  long      T, T1;

  /* HRD consideration */
  long      UpperBound1, UpperBound2, LowerBound;
  double      InitialDelayOffset;
  double      OMEGA;

  double      Wp, Wb;
  int_32_t    TotalPFrame;
  int_32_t    DuantQp;
  int_32_t    PDuantQp;
  FILE      *BitRate;
  double      DeltaP;
  double      bit_rate;
  double      bit_rate_per_frame;
  double      GAMMAP;      /* LIZG, JVT019r1 */
  double      BETAP;      /* LIZG, JVT019r1 */
  int_32_t   goprate;    //xzhao 20081108

  int_32_t    RC_MAX_QUANT;
  int_32_t    RC_MIN_QUANT;

  double      BufferSize;    /* LIZG 25/10/2002 */
  double      GOPTargetBufferLevel;
  double      CurrentBufferFullness;  /* LIZG 25/10/2002 */
  double      TargetBufferLevel;  /* LIZG 25/10/2002 */
  double      PreviousBit_Rate;  /* LIZG 25/10/2002 */
  double      AWp;
  double      AWb;
  int_32_t    MyInitialQp;
  int_32_t    PAverageQp;

  /*
  * LIZG JVT50V2 distortion prediction model ;
  * coefficients of the prediction model
  */
  double      PreviousPictureMAD;
  double      MADPictureC1;
  double      MADPictureC2;
  double      PMADPictureC1;
  double      PMADPictureC2;

  /* LIZG JVT50V2 picture layer MAD */
  myboolean      PictureRejected[21];
  double      PPictureMAD[21];
  double      PictureMAD[21];
  double      ReferenceMAD[21];

  /* quadratic rate-distortion model */
  myboolean      m_rgRejected[21];
  double      P_frm_Qstep[21];
  double      P_frm_R[21];
  double      m_X1, m_X2; //R-Q model parameters
  int_32_t    m_Qc;
  double      m_Qstep;
  int_32_t    m_Qp;
  int_32_t    Pm_Qp;
  int_32_t    PreAveMBHeader;
  int_32_t    CurAveMBHeader;
  int_32_t    PPreHeader;
  int_32_t    PreviousQp1;
  int_32_t    PreviousQp2;

  /* basic unit layer rate control */
  int_32_t    TotalFrameQP;
  int_32_t    NumberofBasicUnit;
  int_32_t    PAveHeaderBits1;
  int_32_t    PAveHeaderBits2;
  int_32_t    PAveHeaderBits3;
  int_32_t    PAveFrameQP;
  int_32_t    TotalNumberofBasicUnit;
  int_32_t    CodedBasicUnit;
  double      MINVALUE;
  double      CurrentFrameMAD;
  double      CurrentBUMAD;
  double      TotalBUMAD;
  double      PreviousFrameMAD;
  int_32_t    m_Hp;
  int_32_t    m_windowSize;
  int_32_t    MADm_windowSize;
  int_32_t    DDquant;
  double      AverageMADPreviousFrame;
  int_32_t    TotalBasicUnitBits;
  int_32_t    QPLastPFrame;
  int_32_t    QPLastGOP;

  double      PP_frm_Qstep[20];
  double      PP_frm_R[20];
  double      Pm_X1;
  double      Pm_X2;

  /* adaptive field/frame coding */
  int_32_t    FieldQPBuffer;
  int_32_t    FrameQPBuffer;
  int_32_t    FrameAveHeaderBits;
  int_32_t    FieldAveHeaderBits;
  double      BUPFMAD[6336];    /* LIZG */
  double      BUCFMAD[6336];    /* LIZG */
  double      FCBUCFMAD[6336];
  double      FCBUPFMAD[6336];
  myboolean   GOPOverdue;

  /* compute macroblock activity for rate control */
  int_32_t    diffy[16][16];
  void_t      rc_init_seq();
  void_t      rc_init_GOP(int_32_t np, int_32_t nb);
  void_t      rc_update_pict_frame(int_32_t nbits);
  void_t      rc_init_pict(int_32_t fieldpic, int_32_t topfield, int_32_t targetcomputation);
  void_t      rc_update_pict(int_32_t nbits);
  void_t      setbitscount(int_32_t nbits);
  int_32_t    updateQuantizationParameter(int_32_t topfield); /* LIZG */
  void_t      updateRCModel();  /* LIZG */
  void_t      updateMADModel();  /* LIZG */
  myboolean      skipThisFrame();  /* LIZG */
  void_t      RCModelEstimator(int_32_t n_windowSize);  /* LIZG */
  void_t      MADModelEstimator(int_32_t n_windowSize);  /* LIZG */
  double      calc_MAD();  /* LIZG */
  double      ComputeFrameMAD();
  int_32_t    Qstep2QP(double Qstep);
  double      QP2Qstep(int_32_t QP);

  /* rdo coding state */
  int_32_t    DELTA_QP, DELTA_QP2;
  int_32_t    pred[16][16];

  /* MODULE PARAMETERS */
  int_16_t    ***cofAC4x4, ****cofAC4x4intern;
  int_32_t    best_mode;
  byte        rec_mbY[16][16], rec_mbU[8][8], rec_mbV[8][8], rec_mbY8x8[16][16];  /* reconstruction
                                                                                  * values */
  int_16_t    mpr8x8[16][16];
  int_16_t    ****cofAC, ****cofAC8x8;  /* [8x8block][4x4block][level/run][scan_pos] */
  int_16_t    ***cofDC;      /* [yuv][level/run][scan_pos] */
  int_32_t    ****chromacofAC4x4;
  int_32_t    cbp, cbp8x8, cnt_nonz_8x8;
  int_32_t    cbp_blk, cbp_blk8x8;
  int_32_t    frefframe[2][2], brefframe[2][2], b8mode[4], b8pdir[4];
  int_32_t    best8x8mode[4];      /* [block] */
  int_32_t    best8x8pdir[MAXMODE][4];  /* [mode][block] */
  int_32_t    best8x8ref[MAXMODE][4];    /* [mode][block] */
  int_32_t    b8_ipredmode[4], b8_intra_pred_modes[4];
  int_32_t    best_c_imode;
  int_32_t    best_weight_flag;    /* !! shenyanfei */

  int_32_t    best8x8bwref[MAXMODE][4];  /* [mode][block] */
  int_32_t    best8x8symref[MAXMODE][4][2];  /* [mode][block] */

  int_32_t    best_intra_pred_modes_tmp[4];
  int_32_t    best_ipredmode_tmp[2][2];
  int_32_t    best_mpr_tmp[16][16];
  int_32_t    best_dct_mode;
  CSptr      cs_mb, cs_b8, cs_cm;
  void_t      delete_coding_state(CSptr);  /* !< delete structure */
  CSptr      create_coding_state();    /* !< create structure */

  void_t      store_coding_state(CSptr);  /* !< store parameters */
  void_t      reset_coding_state(CSptr);  /* !< restore parameters */

  /* rate-distortion optimization */
  void_t      init_rdopt();
  void_t      clear_rdopt();
  double      RDCost_for_AVSIntraBlocks(int_32_t *nonzero,int_32_t b8,int_32_t ipmode,double   lambda,double   min_rdcost,int_32_t mostProbableMode);
  int_32_t    Mode_Decision_for_AVSIntraMacroblock(double lambda, int_32_t *total_cost);
  double      RDCost_for_8x8blocks(int_32_t *cnt_nonz,int_32_t *cbp_blk,double   lambda,int_32_t block,int_32_t mode,int_32_t pdir,int_32_t ref,int_32_t bwd_ref);
  void_t      SetModesAndRefframeForBlocks(int_32_t mode);
  void_t      SetCoeffAndReconstruction8x8(Macroblock *currMB);
  void_t      SetMotionVectorsMB(Macroblock *currMB, int_32_t bframe);
  int_32_t    RDCost_for_macroblocks(double lambda, int_32_t mode, double *min_rdcost);
  void_t      store_macroblock_parameters(int_32_t mode);
  void_t      set_stored_macroblock_parameters();
  void_t            (c_avs_enc:: *encode_one_macroblock) ();
  void_t      encode_one_intra_macroblock_rdo();
  void_t      encode_one_intra_macroblock_not_rdo();
  void_t      encode_one_inter_macroblock_rdo();
  void_t      encode_one_inter_macroblock_not_rdo();
#ifdef _ME_FOR_RATE_CONTROL_
  void_t      encode_one_inter_macroblock_for_rate_control();
  void_t      encode_one_b_frame_macroblock_for_rate_control();
#endif

  void_t      encode_one_b_frame_macroblock_rdo_fast();
  void_t      encode_one_b_frame_macroblock_rdo();
  void_t      encode_one_b_frame_macroblock_not_rdo();
  int_32_t    Mode_Decision_for_AVSIntraMacroblock_not_rdo(double lambda, int_32_t *total_cost);

  /* ref Buffer */
  pel_t      UMVPelY_14(pel_t **Pic, int_32_t y, int_32_t x);
  pel_t      FastPelY_14(pel_t **Pic, int_32_t y, int_32_t x);
  pel_t      *FastLineX(int_32_t dummy, pel_t *Pic, int_32_t y, int_32_t x);
  pel_t      *UMVLineX(int_32_t, pel_t *, int_32_t, int_32_t);

  /* slice */
  void_t      stuffing_byte(int_32_t n);
  int_32_t    start_slice();
  int_32_t    terminate_picture();
  void_t      picture_data();
  void_t      store_field_MV();
  /* vlc */
  int_32_t    u_1(char *tracestring, int_32_t value, Bitstream *part);
  int_32_t    se_v(char *tracestring, int_32_t value, Bitstream *part);
  int_32_t    ue_v(char *tracestring, int_32_t value, Bitstream *part);
  int_32_t    u_v(int_32_t n, char *tracestring, int_32_t value, Bitstream *part);

  int_32_t    writeSyntaxElement_UVLC(SyntaxElement *se, Bitstream *this_dataPart);
  int_32_t    writeSyntaxElement_fixed(SyntaxElement *se, Bitstream *this_dataPart);

  void_t      writeUVLC2buffer(SyntaxElement *se, Bitstream *currStream);
  int_32_t    writeSyntaxElement2Buf_Fixed(SyntaxElement *se, Bitstream *this_streamBuffer);
  int_32_t    symbol2uvlc(SyntaxElement *se);
  void_t      ue_linfo(int_32_t n, int_32_t dummy, int_32_t *len, int_32_t *info);
  void_t      se_linfo(int_32_t mvd, int_32_t dummy, int_32_t *len, int_32_t *info);
  void_t      cbp_linfo_intra(int_32_t cbp, int_32_t dummy, int_32_t *len, int_32_t *info);
  void_t      cbp_linfo_inter(int_32_t cbp, int_32_t dummy, int_32_t *len, int_32_t *info);
  int_32_t    writeSyntaxElement_Intra4x4PredictionMode(SyntaxElement *se, Bitstream *this_dataPart);
  void_t      writeVlcByteAlign(Bitstream *currStream);
  int_32_t    writeSyntaxElement_Run(SyntaxElement *se, Bitstream *bitstream);
  int_32_t    symbol2vlc(SyntaxElement *sym);
  int_32_t    encode_IP_frame(avs_enc_frame_t *pFrame);
  int_32_t    encode_B_frame(avs_enc_frame_t *pFrame);    
  void_t      OpenORABS(OutputStream *p, char *fname);
  void_t      CloseORABS(OutputStream *p);
  void_t      FlushORABS(OutputStream *p);

  /* entropy coding */
  int_32_t write_n_bit(OutputStream *p,int_32_t b,int_32_t n);
  _inline  void_t write_1_bit(OutputStream *p, int_32_t b);
  int_32_t    write_align_stuff(OutputStream *p);
  int_32_t    write_a_byte(OutputStream *p, int_32_t b);
  /* header */
  int_32_t    bbv_check_times;
  int_32_t    luma_scale;
  int_32_t    luma_shift;
  int_32_t    chroma_scale;
  int_32_t    chroma_shift;
  int_32_t    low_delay;
  int_32_t    slice_vertical_position_extension;
  int_32_t    stream_length;
  int_32_t    picture_distance;
  //rm52j
  int_32_t       pic_dist;             // picture distance relative to the I frame
  int_32_t       last_pic_dist_i;      // distance for last I picture 
  int_32_t       curr_pic_dist_i;      // distance for current I picture 

  int_32_t    picture_decoder_order;
  int_32_t    bbv_delay;
  int_32_t    picture_coding_type;
  int_32_t    picture_reference_flag;
  int_32_t    picture_reference_index;
  float           frame_rate;
  int_32_t    tc0;
  void_t      picture_header();

  /* encoding frame */
  int_32_t    encode_one_frame();

  /* quantization */
  __m128i clip0q;
  __m128i clip255q;
  __inline void_t    avs_const_initialize();
  __inline void_t       avs_const_initialize_block();
  __inline __m128i  avs_clip_0_255_w(__m128i xmm0);
  __inline __m128i  avs_combine_w2b(__m128i xmm0, __m128i xmm1);
  __inline __m128i  avs_combine_d2w(__m128i xmm0, __m128i xmm1);
  __inline __m128i  avs_filter_halfpel_w(__m128i xmm0, __m128i xmm1, __m128i xmm2);
  __inline __m128i  avs_filter_quaterpel_w(__m128i xmm0, __m128i xmm1, __m128i xmm2);
  __inline __m128i  avs_filter_quaterpel_d(__m128i xmm0, __m128i xmm1, __m128i xmm2);
  __inline __m128i  avs_zero(__m128i xmm0);
  void_t          Update_Picture_Buffers();
  void_t          OneComponentLumaPrediction4x4(int_32_t *mpred,int_32_t pic_pix_x,int_32_t pic_pix_y,int_32_t *mv,int_32_t ref);
  void_t          SetRefAndMotionVectors(int_32_t block,int_32_t mode,int_32_t ref,int_32_t bw_ref,int_32_t pdir);
  avs_enc_create_t* p_avsCreate;
  avs_enc_frame_t*  p_avs_enc_frame;
  //avs_enc_stats_t * p_stats;
  pel_t line[16];
  MB_INFO * pMbInfo;
  int_32_t dec_frm_num;
#ifdef ROI_ENABLE
  byte *ROIArray;
  byte* YCbCr[3];
  int_32_t w[3],h[3];
  LocateTextRegions_t detect_roi;
#endif
};
#endif
