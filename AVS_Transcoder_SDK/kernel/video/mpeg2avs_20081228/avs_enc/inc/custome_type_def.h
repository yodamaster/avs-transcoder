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
#ifndef _CUSTOME_TYPE_DEF_H_
#define _CUSTOME_TYPE_DEF_H_
class    c_avs_enc;
/* ! Boolean Type */
typedef enum { FALSE, TRUE } Boolean;

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
  void (c_avs_enc:: *mapping) (int_32_t value1, int_32_t value2, int_32_t * len_ptr, int_32_t * info_ptr);
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
  int_32_t  mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];
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
  /*
  * int_32_t jumpd;
  * //!< number of frames to skip in input sequence (e.g 2 takes frame 0,3,6,9...)
  */
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
  char    infile[500];  /* !< YUV 4:2:0 input format */
  char    outfile[500];  /* !< H.26L compressed output bitstream */
  char    ReconFile[500]; /* !< Reconstructed Pictures */
  char    TraceFile[500]; /* !< Trace Outputs */
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

  char    PictureTypeSequence[MAXPICTURETYPESEQUENCELEN];

  int_32_t  rdopt;

  int_32_t  InterlaceCodingOption;

  /* AVS */
  int_32_t  aspect_ratio_information;
  int_32_t  frame_rate_code;  /* xfwang 2004.7.28 */

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
  int_32_t  types;      /* !< This is for SP-Pictures, since all the syntax elements for
                * SP-Pictures are the same as P-pictures, we keep the img->type as
                * P_IMG but indicate SP-Pictures by img->types */
  int_32_t  no_multpred;    /* !< 1: prediction from the last frame only. 2: prediction from the
                  * last or second last frame etc. */
  int_32_t  qp;      /* !< quant for the current frame */
  int_32_t  framerate;

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
  int_32_t  Frame_Total_Number_MB;
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
  int_32_t  bit_ctr;  /* !< counter for bit usage */
  int_32_t  bit_ctr_0;  /* !< stored bit use for the first frame */
  int_32_t  bit_ctr_n;  /* !< bit usage for the current frame */
  int_32_t  bit_slice;  /* !< number of bits in current slice */
  int_32_t  bit_use_mode_inter[2][MAXMODE]; /* !< statistics of bit usage */
  int_32_t  bit_ctr_emulationprevention;  /* !< stored bits needed to prevent start code emulation */
  int_32_t  mode_use_intra[25];    /* !< Macroblock mode usage for Intra frames */
  int_32_t  mode_use_inter[2][MAXMODE];

  int_32_t  mb_use_mode[2];

  /* B pictures */
  int_32_t  *mode_use_Bframe;
  int_32_t  *bit_use_mode_Bframe;
  int_32_t  bit_ctr_P;
  int_32_t  bit_ctr_B;
  float    bitrate_P;
  float    bitrate_B;

  int_32_t  bit_use_stuffingBits[NUM_PIC_TYPE];
  int_32_t  bit_use_mb_type[NUM_PIC_TYPE];
  int_32_t  bit_use_header[NUM_PIC_TYPE];
  int_32_t  tmp_bit_use_cbp[NUM_PIC_TYPE];
  int_32_t  bit_use_coeffY[NUM_PIC_TYPE];
  int_32_t  bit_use_coeffC[NUM_PIC_TYPE];
  int_32_t  bit_use_delta_quant[NUM_PIC_TYPE];

  int_32_t  em_prev_bits_frm;
  int_32_t  em_prev_bits_fld;
  int_32_t  *em_prev_bits;
  int_32_t  bit_ctr_parametersets;
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
  int_32_t  mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];
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
  byte        buf[SVA_STREAM_BUF_SIZE];  /* 流缓冲区 */
  uint_32_t  uPreBytes;      /* 最近写入的3个字节，初始值是0xFFFFFFFF */
  int_32_t  iBytePosition;      /* 当前字节位置 */
  int_32_t  iBitOffset;      /* 当前位偏移，0表示最高位 */
  int_32_t  iNumOfStuffBits;    /* 已插入的填充位的个数，遇到开始码时置0 */
  int_32_t  iBitsCount;      /* 码流总位数 */
  } OutputStream;
#endif