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

#if defined WIN32
#include <conio.h>
#endif
#include <math.h>

#include "avs_enc.h"
#include "global.h"

#ifdef FastME
#include "fast_me.h"
#endif
int_32_t LocateTextRegions(byte* YCbCr[3],int_32_t w[3],int_32_t h[3],byte* ROIArray);

c_avs_enc::c_avs_enc() :
Iprev_bits(0),
Pprev_bits(0),
cofAC(NULL),
cofAC8x8(NULL),
cofDC(NULL),
cofAC4x4(NULL),
cofAC4x4intern(NULL),
chromacofAC4x4(NULL),
cs_mb(NULL),
cs_b8(NULL),
cs_cm(NULL),
THETA(1.3636),
OMEGA(0.9),
Switch(0)
{
}
/*
=======================================================================================================================
=======================================================================================================================
*/
DLL_EXPORT void avs_encoder_create(c_avs_enc *p_avs_enc)
{
  p_avs_enc->avs_enc_create();
}

DLL_EXPORT void avs_encoder_encode(c_avs_enc *p_avs_enc)
{
  p_avs_enc->avs_enc_encode();
}

/*
=======================================================================================================================
=======================================================================================================================
*/
DLL_EXPORT void avs_encoder_destropy(c_avs_enc *p_avs_enc)
{
  p_avs_enc->avs_enc_destroy();
}

/*
=======================================================================================================================
Function:Main function for encoder. Input:argc number of command line arguments argv command line arguments
Output: Return: exit code Attention:
=======================================================================================================================
*/
int_32_t c_avs_enc::avs_enc_create()
{
  p_rec = p_stat = p_log = p_datpart = p_trace = NULL;
  total_encoded_frame = -1;
  img = &images;
  input = &inputs;
  snr = &snrs;
  stat = &stats;
  cp = &CopyRights;
  camera = &CameraParameter;
  pORABS = &ORABS;
  memset(pORABS, 0, SVA_STREAM_BUF_SIZE*sizeof(byte));
  snr->snr_y  = 0;
  snr->snr_u  = 0;
  snr->snr_v  = 0;
  snr->snr_y1 = 0;
  snr->snr_u1 = 0;
  snr->snr_v1 = 0;
  snr->snr_ya = 0;
  snr->snr_ua = 0;
  snr->snr_va = 0;
  stat->quant0 = 0;
  stat->quant1 = 0;
  stat->bitr = stat->bitr0 = stat->bitrate = 0;
  stat->bit_ctr = stat->bit_ctr_0 = stat->bit_ctr_n = 0;
  stat->bit_slice = stat->bit_ctr_emulationprevention = 0;
  memset(stat->bit_use_mode_inter[0], 0, MAXMODE*sizeof(int_32_t));
  memset(stat->bit_use_mode_inter[1], 0, MAXMODE*sizeof(int_32_t));
  memset(stat->mode_use_intra, 0, 25*sizeof(int_32_t));
  memset(stat->mode_use_inter[0], 0, MAXMODE*sizeof(int_32_t));
  memset(stat->mode_use_inter[1], 0, MAXMODE*sizeof(int_32_t));
  memset(stat->mb_use_mode, 0, 2*sizeof(int_32_t));
  stat->mode_use_Bframe = stat->bit_use_mode_Bframe = stat->em_prev_bits = NULL;
  stat->bit_ctr_P = stat->bit_ctr_B = 0;
  memset(stat->bit_use_stuffingBits, 0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->bit_use_mb_type,      0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->bit_use_header,       0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->tmp_bit_use_cbp,      0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->bit_use_coeffY,       0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->bit_use_coeffC,       0, NUM_PIC_TYPE*sizeof(int_32_t));
  memset(stat->bit_use_delta_quant,  0, NUM_PIC_TYPE*sizeof(int_32_t));
  stat->em_prev_bits_frm = stat->em_prev_bits_fld = stat->bit_ctr_parametersets = 0;
  memset(snr, 0, sizeof(SNRParameters));
  //memset(cp, 0, sizeof(CopyRight));
  //memset(camera, 0, sizeof(CameraParamters));
#ifdef _THREE_STEP_MOTION_SEARCH_
  init_3_step_search();
#endif

  input->intra_period = input->GopLength/(input->successive_Bframe + 1);  //for close GOP

  PatchInp();
  init_img();

  frame_pic = malloc_picture();
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    top_pic = malloc_picture();
    bot_pic = malloc_picture();
  }

  init_rdopt();

  if(input->InterlaceCodingOption == FRAME_CODING)
    input->progressive_sequence = 1;
  else
    input->progressive_sequence = 0;

  /* allocate memory for frame buffers */
  init_global_buffers();
  Init_Motion_Search_Module();
#ifdef _DEBUG
  information_init();
#endif
  nOutBufPtr = 0;
#ifdef ROI_ENABLE
  detect_roi = LocateTextRegions;
#endif
  return 0;
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::avs_enc_destroy()
{
  if(p_rec) fclose(p_rec);
  if(p_trace) fclose(p_trace);
#ifdef _OUTPUT_DEC_IMG_
  if (p_org_dec) fclose(p_org_dec);
#endif
  Clear_Motion_Search_Module();

  /* free structure for rd-opt. mode decision */
  clear_rdopt();

  /* report everything */
  report();

  free_picture(frame_pic);

  free_global_buffers();

  /* free image mem */
  free_img();
  return 0;
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::avs_enc_encode()
{
  int_32_t i;

  init_global_variables();
  OpenBitStreamFile(input->outfile);
  p_avs_enc_frame->length = 0;
  if (img->number == 0)
  {
    total_encoded_frame = -1;
  }
  gframe_no=current_encoded_frame;
  goprate = 0;

  for (i=0; i<dec_frm_num; i++)
  {
    if (i == 0)
      p_avs_enc_frame->type[0] = AVS_TYPE_I;
    else if (inputs.successive_Bframe == 0 || i % (inputs.successive_Bframe + 1) == 1 || (i>input->GopLength-input->successive_Bframe && dec_frm_num%(inputs.successive_Bframe + 1)!=1))
      p_avs_enc_frame->type[0] = AVS_TYPE_P;
    else
      p_avs_enc_frame->type[0] = AVS_TYPE_B;

    p_avs_enc_frame->input = p_avs_enc_frame->inputfrm[i];
    avs_enc_frame(p_avs_enc_frame);
    gframe_no++;
  }
  memcpy((byte*)p_avs_enc_frame->bitstream + p_avs_enc_frame->length, pORABS->buf, pORABS->iBytePosition);
  p_avs_enc_frame->length += pORABS->iBytePosition;
 
  return 0;
}

/*
=======================================================================================================================
Function:Initializes the Image structure with appropriate parameters. Input:Input Parameters struct inp_par inp
Output:Image Parameters struct img_par *img Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::init_img()
{
  int_32_t  i;
  float    FrameRate[8] =
  {
    { 24000 / 1001 },
    { 24 },
    { 25 },
    { 30000 / 1001 },
    { 30 },
    { 50 },
    { 60000 / 1001 },
    { 60 }
  };

  img->no_multpred = input->no_multpred;
  img->buf_cycle = input->no_multpred;
  img->nb_references = 0;
  img->lindex = 0;
  img->max_lindex = 0;

  img->width  = input->img_width;
  img->height = input->img_height;
  img->width_cr = input->img_width / 2;
  img->height_cr = input->img_height / 2;
  img->reconflag = input->rdopt;

  img->framerate = input->fr;
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    img->buf_cycle *= 2;
  }

  get_mem_mv(&(img->mv));
  get_mem_mv(&(img->p_fwMV));
  get_mem_mv(&(img->p_bwMV));
  get_mem_mv(&(img->all_mv));
  get_mem_mv(&(img->all_bmv));
  get_mem_mv(&(img->all_bw_omv));
  get_mem_ACcoeff(&(img->cofAC));
  get_mem_DCcoeff(&(img->cofDC));
  get_mem_mv(&(img->omv));
  get_mem_mv(&(img->all_omv));
  get_mem_mv(&(img->omv_fld));
  get_mem_mv(&(img->all_omv_fld));
  get_mem4Dint(&(img->chromacofAC), 2, 4, 2, 17);

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    img->buf_cycle /= 2;
  }

  if((img->quad = (int_32_t *) calloc(511, sizeof(int_32_t))) == NULL)
  {
    no_mem_exit("init_img: img->quad");
  }

  img->quad += 255;

  for(i = 0; i < 256; ++i)
  {
    img->quad[i] = img->quad[-i] = i * i;
  }

  if(((img->mb_data) = (Macroblock *) calloc((img->width / MB_BLOCK_SIZE) * (img->height / MB_BLOCK_SIZE),sizeof(Macroblock))) == NULL)
  {
    no_mem_exit("init_img: img->mb_data");
  }

  for(i = 0; i < (img->width / MB_BLOCK_SIZE) * (img->height / MB_BLOCK_SIZE); i++)
  {
    get_mem4Dint(&(img->mb_data[i].cofAC), 6, 4, 2, 65);
    get_mem4Dint(&(img->mb_data[i].chromacofAC), 2, 4, 2, 17);
  }

  for(i = 0; i < (img->width / MB_BLOCK_SIZE) * (img->height / MB_BLOCK_SIZE); i++)
  {
    img->mb_data[i].slice_nr = 0;
  }

  /* allocate memory for intra pred mode buffer for each block: img->ipredmode */
  get_mem2Dint(&(img->ipredmode), img->width / B8_SIZE + 100, img->height / B8_SIZE + 100);

  /*
  * Prediction mode is set to -1 outside the frame, indicating that no prediction
  * can be made from this part
  */
  for(i = 0; i < img->width / (B8_SIZE) + 100; i++)
  {
    memset(img->ipredmode[i], -1, (img->height / (B8_SIZE) + 100)*sizeof(int_32_t));
  }

  img->img_width_in_mb  = img->width  >> 4;
  img->img_height_in_mb = img->height >> 4;
  img->mb_no_currSliceLastMB = (input->slice_row_nr != 0) ? min(input->slice_row_nr * img->img_width_in_mb - 1, img->img_width_in_mb * img->img_height_in_mb - 1) : img->img_width_in_mb * img->img_height_in_mb - 1;
  img->total_number_mb = (img->width * img->height) >> 8;

  GBIM_value = 0;
}

/*
=======================================================================================================================
Function:Free the Image structures Input:Image Parameters struct img_par *img Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_img()
{
  if(input->InterlaceCodingOption != FRAME_CODING) img->buf_cycle *= 2;

  free_mem_mv(img->mv);
  free_mem_mv(img->p_fwMV);
  free_mem_mv(img->p_bwMV);
  free_mem_mv(img->all_mv);
  free_mem_mv(img->all_bmv);
  free_mem_mv(img->omv);
  free_mem_mv(img->all_omv);
  free_mem_mv(img->all_bw_omv);
  free_mem_mv(img->omv_fld);
  free_mem_mv(img->all_omv_fld);

  if(input->InterlaceCodingOption != FRAME_CODING) img->buf_cycle /= 2;

  /* Lou Start */
  free_mem4Dint(img->chromacofAC, 2, 4);

  /* Lou End */
  free_mem_ACcoeff(img->cofAC);
  free_mem_DCcoeff(img->cofDC);
  free(img->quad - 255);
}

/*
=======================================================================================================================
Function:Allocates the picture structure along with its dependent data structures Input: Output: Return:
Pointer to a Picture Attention:
=======================================================================================================================
*/
Picture *c_avs_enc::malloc_picture()
{
  /*~~~~~~~~~*/
  Picture *pic;
  /*~~~~~~~~~*/

  if((pic = (Picture *) calloc(1, sizeof(Picture))) == NULL) no_mem_exit("malloc_picture: Picture structure");
  return pic;
}

/*
=======================================================================================================================
Function:Frees a picture Input:pic: POinter to a Picture to be freed Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_picture(Picture *pic)
{
  if(pic != NULL)
  {
    free(pic);
  }
}

/*
=======================================================================================================================
Function:Reports the gathered information to appropriate outputs Input: struct inp_par *inp, \n struct img_par
img, \n struct stat_par *stat, \n struct stat_par *stat Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::report()
{
  int_32_t  bit_use[2][2];
  int_32_t  i, j;
  int_32_t  bit_use_Bframe = 0;
  uint_32_t  total_bits;
  float    frame_rate;

  bit_use[0][0] = 1;
  bit_use[1][0] = max(1, input->no_frames - 1);

  /* Accumulate bit usage for inter and intra frames */
  bit_use[0][1] = bit_use[1][1] = 0;

  for(i = 0; i < 11; i++) bit_use[1][1] += stat->bit_use_mode_inter[0][i];

  for(j = 0; j < 2; j++)
  {
    bit_use[j][1] += stat->bit_use_header[j];
    bit_use[j][1] += stat->bit_use_mb_type[j];
    bit_use[j][1] += stat->tmp_bit_use_cbp[j];
    bit_use[j][1] += stat->bit_use_coeffY[j];
    bit_use[j][1] += stat->bit_use_coeffC[j];
    bit_use[j][1] += stat->bit_use_delta_quant[j];
    bit_use[j][1] += stat->bit_use_stuffingBits[j];
  }

  /* B pictures */
  if(Bframe_ctr != 0)
  {
    bit_use_Bframe = 0;
    for(i = 0; i < 11; i++) bit_use_Bframe += stat->bit_use_mode_inter[1][i];
    bit_use_Bframe += stat->bit_use_header[2];
    bit_use_Bframe += stat->bit_use_mb_type[2];
    bit_use_Bframe += stat->tmp_bit_use_cbp[2];
    bit_use_Bframe += stat->bit_use_coeffY[2];
    bit_use_Bframe += stat->bit_use_coeffC[2];
    bit_use_Bframe += stat->bit_use_delta_quant[2];
    bit_use_Bframe += stat->bit_use_stuffingBits[2];

    stat->bitrate_P = (stat->bit_ctr_0 + stat->bit_ctr_P) * (float) (img->framerate / (input->successive_Bframe + 1)) / input->no_frames;
    stat->bitrate_B = (stat->bit_ctr_B) * (float) (img->framerate / (input->successive_Bframe + 1)) * input->successive_Bframe / Bframe_ctr;
  }
  else
  {
    if(input->no_frames > 1)
    {
      stat->bitrate = (bit_use[0][1] + bit_use[1][1]) * (float) img->framerate / (input->no_frames * (input->successive_Bframe + 1));
    }
  }

  fprintf(stdout, "--------------------------------------------------------------------------\n");
  fprintf (stdout, " Freq. for encoded bitstream       : %1.0f\n", (float) img->framerate / (float) (input->successive_Bframe + 1));
  if(input->hadamard)
    fprintf(stdout, " Hadamard transform                : Used\n");
  else
    fprintf(stdout, " Hadamard transform                : Not used\n");

  fprintf(stdout, " Image format                      : %dx%d\n", input->img_width, input->img_height);

  if(input->intra_upd)
    fprintf(stdout, " Error robustness                  : On\n");
  else
    fprintf(stdout, " Error robustness                  : Off\n");
  fprintf(stdout, " Search range                      : %d\n", input->search_range);

  fprintf(stdout, " No of ref. frames used in P pred  : %d\n", input->no_multpred);
  if(input->successive_Bframe != 0)
    fprintf(stdout, " No of ref. frames used in B pred  : %d\n", input->no_multpred);

  fprintf(stdout, " Total encoding time for the seq.  : %.3f sec \n", tot_time * 0.001);

  /* B pictures */
  fprintf(stdout, " Sequence type                     :");

  if(input->successive_Bframe == 1)
    fprintf(stdout, " IBPBP (QP: I %d, P %d, B %d) \n", input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe == 2)
    fprintf(stdout, " IBBPBBP (QP: I %d, P %d, B %d) \n", input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe == 0)
    fprintf(stdout, " IPPP (QP: I %d, P %d) \n", input->qp0, input->qpN);

  /* report on entropy coding method */
  fprintf(stdout, " Entropy coding method             : VLC\n");

  if(input->rdopt)
    fprintf(stdout, " RD-optimized mode decision        : used\n");
  else
    fprintf(stdout, " RD-optimized mode decision        : not used\n");

  fprintf(stdout, "------------------ Average data all frames  ------------------------------\n");
  fprintf(stdout, " SNR Y(dB)                         : %5.2f\n", snr->snr_ya);
  fprintf(stdout, " SNR U(dB)                         : %5.2f\n", snr->snr_ua);
  fprintf(stdout, " SNR V(dB)                         : %5.2f\n", snr->snr_va);

  if(Bframe_ctr != 0)
  {
    total_bits = stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_B + stat->sequence_header;
    fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, B %d) \n", total_bits, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_B);
    frame_rate = (float) (img->framerate);
    stat->bitrate = ((float) total_bits * frame_rate) / ((float)(total_encoded_frame + 1));
    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate / 1000);
  }
  else
  {
    total_bits = stat->bit_ctr_P + stat->bit_ctr_0 + stat->sequence_header;
    fprintf (stdout, " Total bits                        : %u (I %5u, P %5u) \n", total_bits, stat->bit_ctr_0, stat->bit_ctr_P);
    frame_rate = (float) img->framerate / ((float) (input->successive_Bframe + 1));
    stat->bitrate = ((float) total_bits * frame_rate) / ((float) input->no_frames);
    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate / 1000);
  }

  fprintf(stdout, " Bits to avoid Startcode Emulation : %d \n", stat->bit_ctr_emulationprevention);
  fprintf(stdout, " Bits for parameter sets           : %d \n", stat->bit_ctr_parametersets);
  fprintf(stdout, " GBIM value                  : %f \n", GBIM_value / (total_encoded_frame + 1));
  fprintf(stdout, "--------------------------------------------------------------------------\n");
  fprintf(stdout, "Exit IDM_AVS_Transcoder with %s ", VERSION);
  fprintf(stdout, "\n");
}

/*
=======================================================================================================================
Function:Prints the header of the protocol. Input:struct inp_par *inp Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::information_init()
{
  printf("--------------------------------------------------------------------------\n");
  printf(" Input Video file                  : %s \n", input->infile);
  printf(" Output AVS bitstream              : %s \n", input->outfile);
  if(p_rec != NULL) printf(" Output YUV file                   : %s \n", input->ReconFile);
  printf(" Output log file                   : log.dat \n");
  printf(" Output statistics file            : stat.dat \n");
  printf("--------------------------------------------------------------------------\n");
  printf(" Frame   Bit/pic   QP   SnrY    SnrU    SnrV  GBIM    Time(ms)  FRM/FLD  IntraMBs\n");
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::init_global_buffers()
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  memory_size = 0;
  int_32_t  height_field = img->height / 2;
  int_32_t  refnum;
  int_32_t  i;
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /* 申请参考帧的内存 */
  for(i = 0; i < 2; i++)
  {
    get_mem_ref(&mref_frm[i]);
  }

  imgY_org_buffer = (byte *) _aligned_malloc(img->height * img->width * 3 / 2, 16);

  /* xzhao { 2007.7.21 */
  imgY_org_buffer_restore = imgY_org_buffer;

  /* xzhao } */
  memory_size += img->height * img->width * 3 / 2;

  /* allocate memory for reference frame buffers: imgY_org, imgUV_org */
  memory_size += get_mem2D(&imgY_org_frm, img->height, img->width);
  memory_size += get_mem3D(&imgUV_org_frm, 2, img->height_cr, img->width_cr);
#ifdef _OUTPUT_DEC_IMG_
  get_mem2D(&org_nextP_imgY, img->height, img->width);
  get_mem3D(&org_nextP_imgUV, 2, img->height_cr, img->width_cr);
#endif
  /*
  * allocate memory for temp P and B-frame motion vector buffer: tmp_mv,
  * temp_mv_block ;
  * int_32_t tmp_mv[2][72][92];
  * ([2][72][88] should be enough)
  */
  memory_size += get_mem3Dint(&tmp_mv_frm, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
#ifdef _DEBUG_B_
  memory_size += get_mem3Dint(&tmp_mv_p_frm, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
  memory_size += get_mem2Dint(&refFrArr_p_frm,  img->height / BLOCK_SIZE, img->width / BLOCK_SIZE);
#endif
  /*
  * allocate memory for reference frames of each block: refFrArr ;
  * int_32_t refFrArr[72][88];
  */
  memory_size += get_mem2Dint(&refFrArr_frm, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE);

  if(input->successive_Bframe != 0)
  {
    /*
    * allocate memory for temp B-frame motion vector buffer: fw_refFrArr, bw_refFrArr ;
    * int_32_t ...refFrArr[72][88];
    */
    memory_size += get_mem2Dint(&fw_refFrArr_frm, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE);
    memory_size += get_mem2Dint(&bw_refFrArr_frm, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE);
  }

  /*
  * allocate memory for B frame coding: nextP_imgY, nextP_imgUV ;
  * byte nextP_imgY[288][352];
  * byte nextP_imgUV[2][144][176];
  */
  memory_size += get_mem2D(&nextP_imgY, img->height, img->width);
  memory_size += get_mem3D(&nextP_imgUV, 2, img->height_cr, img->width_cr);

  if(input->successive_Bframe != 0)
  {
    /*
    * allocate memory for temp B-frame motion vector buffer: tmp_fwMV, tmp_bwMV,
    * dfMV, dbMV ;
    * int_32_t ...MV[2][72][92];
    * ([2][72][88] should be enough)
    */
    memory_size += get_mem3Dint(&tmp_fwMV, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
    memory_size += get_mem3Dint(&tmp_bwMV, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
    memory_size += get_mem3Dint(&dfMV, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
    memory_size += get_mem3Dint(&dbMV, 2, img->height / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
  }

  /*
  * allocate memory for temp quarter pel luma frame buffer: img4Y_tmp ;
  * int_32_t img4Y_tmp[576][704];
  * (previously int_32_t imgY_tmp in global.h) ;
  * memory_size += get_mem2Dint(&img4Y_tmp, img->height+2*IMG_PAD_SIZE, (img->width+2*IMG_PAD_SIZE)*4);
  * memory_size += get_mem2Dint(&img4Y_tmp, (img->height+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
  */
  memory_size += get_mem2Dshort_int
    (
    &tmp02,
    (img->height + (IMG_PAD_SIZE << 1)),
    (img->width + (IMG_PAD_SIZE << 1)) << 1
    );
  memory_size += get_mem2Dshort_int
    (
    &tmp20,
    (img->height + (IMG_PAD_SIZE << 1)),
    (img->width + (IMG_PAD_SIZE << 1)) << 1
    );
  memory_size += get_mem2Dshort_int
    (
    &tmp22,
    (img->height + (IMG_PAD_SIZE << 1)),
    (img->width + (IMG_PAD_SIZE << 1)) << 1
    );

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    /*
    * allocate memory for encoding frame buffers: imgY, imgUV ;
    * byte imgY[288][352];
    * byte imgUV[2][144][176];
    */
    memory_size += get_mem2D(&imgY_com, img->height, img->width);
    memory_size += get_mem3D(&imgUV_com, 2, img->height / 2, img->width_cr);

    memory_size += get_mem2D(&imgY_org_top, img->height / 2, img->width);
    memory_size += get_mem3D(&imgUV_org_top, 2, img->height_cr / 2, img->width_cr);

    memory_size += get_mem2D(&imgY_org_bot, img->height / 2, img->width);
    memory_size += get_mem3D(&imgUV_org_bot, 2, img->height_cr / 2, img->width_cr);

    /*
    * allocate memory for encoding frame buffers: imgY, imgUV ;
    * byte imgY[288][352];
    * byte imgUV[2][144][176];
    */
    memory_size += get_mem2D(&imgY_top, height_field, img->width);
    memory_size += get_mem3D(&imgUV_top, 2, height_field / 2, img->width_cr);
    memory_size += get_mem2D(&imgY_bot, height_field, img->width);
    memory_size += get_mem3D(&imgUV_bot, 2, height_field / 2, img->width_cr);

    if(input->successive_Bframe != 0)
    {
      /*
      * allocate memory for temp B-frame motion vector buffer: fw_refFrArr, bw_refFrArr ;
      * int_32_t ...refFrArr[72][88];
      */
      memory_size += get_mem2Dint
        (
        &fw_refFrArr_top,
        height_field / BLOCK_SIZE,
        img->width / BLOCK_SIZE
        );
      memory_size += get_mem2Dint
        (
        &bw_refFrArr_top,
        height_field / BLOCK_SIZE,
        img->width / BLOCK_SIZE
        );
      memory_size += get_mem2Dint
        (
        &fw_refFrArr_bot,
        height_field / BLOCK_SIZE,
        img->width / BLOCK_SIZE
        );
      memory_size += get_mem2Dint
        (
        &bw_refFrArr_bot,
        height_field / BLOCK_SIZE,
        img->width / BLOCK_SIZE
        );
    }

    /*
    * allocate memory for temp P and B-frame motion vector buffer: tmp_mv,
    * temp_mv_block ;
    * int_32_t tmp_mv[2][72][92];
    * ([2][72][88] should be enough)
    */
    memory_size += get_mem3Dint(&tmp_mv_top, 2, height_field / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);
    memory_size += get_mem3Dint(&tmp_mv_bot, 2, height_field / BLOCK_SIZE, img->width / BLOCK_SIZE + 4);

    /*
    * allocate memory for reference frames of each block: refFrArr ;
    * int_32_t refFrArr[72][88];
    */
    memory_size += get_mem2Dint(&refFrArr_top, height_field / BLOCK_SIZE, img->width / BLOCK_SIZE);
    memory_size += get_mem2Dint(&refFrArr_bot, height_field / BLOCK_SIZE, img->width / BLOCK_SIZE);
  }

  /* FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3 */
#ifdef FastME
  memory_size += get_mem_FME();
#endif
  for(refnum = 0; refnum < 3; refnum++)
  {
    for(i = 0; i < 3; i++)
    {
      if(i == 0)
      {
        get_mem2D(&reference_frame[refnum][i], img->height, img->width);
      }
      else
      {
        get_mem2D(&reference_frame[refnum][i], img->height_cr, img->width_cr);
      }
    }
  }

  /* forward reference frame buffer */
  ref_frm[0] = reference_frame[0];    /* reference_frame[ref_index][yuv][height][width],ref_frm[ref_index][yuv][height][width]
                                      * */
  ref_frm[1] = reference_frame[1];
  current_frame = reference_frame[2];

  /* allocate field buffer */
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    for(refnum = 0; refnum < 6; refnum++)
    {
      for(i = 0; i < 3; i++)
      {
        if(i == 0)
        {
          get_mem2D(&reference_field[refnum][i], img->height / 2, img->width);
        }
        else
        {
          get_mem2D(&reference_field[refnum][i], img->height_cr / 2, img->width_cr);
        }
      }
    }

    /* forward reference frame buffer */
    for(i = 0; i < 4; i++) ref_fld[i] = reference_field[i];
    current_field = reference_field[4];
    ref_fld[4] = reference_field[5];
  }

  /* !! */
  allalpha_lum = (int_32_t *) malloc(((img->height * img->width) / 256) * sizeof(int_32_t));
  allbelta_lum = (int_32_t *) malloc(((img->height * img->width) / 256) * sizeof(int_32_t));
#ifdef ROI_ENABLE
  ROIArray = (byte*)calloc(img->width*img->height/64, sizeof(byte));
  YCbCr[0] = (byte*)malloc(img->width*img->height*sizeof(byte));
  YCbCr[1] = (byte*)malloc(img->width_cr*img->height_cr*sizeof(byte));
  YCbCr[2] = (byte*)malloc(img->width_cr*img->height_cr*sizeof(byte));
#endif
  return(memory_size);
}

/*
=======================================================================================================================
Function:Free allocated memory of frame size related global buffers buffers are defined in global.h, allocated
memory is allocated in int_32_t get_mem4global_buffers() Input: Input Parameters struct inp_par *inp, \n Image
Parameters struct img_par *img Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_global_buffers()
{
  /*~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j;
  /*~~~~~~~~~~~~~~~~~*/

  /* xzhao { 2007.7.21 */
  _aligned_free(imgY_org_buffer_restore);

  /*
  * free(imgY_org_buffer);
  * xzhao }
  */
  free_mem2D(imgY_org_frm);
  free_mem3D(imgUV_org_frm, 2);
  free_mem3Dint(tmp_mv_frm, 2);
#ifdef _DEBUG_B_
  free_mem3Dint(tmp_mv_p_frm, 2);
  free_mem2Dint(refFrArr_p_frm);
#endif
  free_mem2Dint(refFrArr_frm);
#ifdef _OUTPUT_DEC_IMG_
  free_mem2D(org_nextP_imgY);
  free_mem3D(org_nextP_imgUV, 2);
#endif
  for(i = 0; i < 3; i++)
  {
    _aligned_free(reference_frame[0][i]);
    _aligned_free(reference_frame[1][i]);
    _aligned_free(reference_frame[2][i]);
  }

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    for(j = 0; j < 6; j++)
    {
      for(i = 0; i < 3; i++)
      {
        _aligned_free(reference_field[j][i]);
      }
    }
  }

  free_mem2D(nextP_imgY);
  free_mem3D(nextP_imgUV, 2);

  /*
  * free multiple ref frame buffers ;
  * number of reference frames increased by one for next P-frame
  */
  if(input->successive_Bframe != 0)
  {
    /* free last P-frame buffers for B-frame coding */
    free_mem3Dint(tmp_fwMV, 2);
    free_mem3Dint(tmp_bwMV, 2);
    free_mem3Dint(dfMV, 2);
    free_mem3Dint(dbMV, 2);
    free_mem2Dint(fw_refFrArr_frm);
    free_mem2Dint(bw_refFrArr_frm);
  }  /* end if B frame */

  free_mem2Dshort_int(tmp02);
  free_mem2Dshort_int(tmp20);
  free_mem2Dshort_int(tmp22);

  /*
  * free mem, allocated in init_img() ;
  * free intra pred mode buffer for blocks
  */
  free_mem2Dint(img->ipredmode);

  for(i = 0; i < 2; i++)
  {
    free_mem_ref(mref_frm[i]);
  }

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    free_mem2D(imgY_com);

    /*
    * free multiple ref frame buffers ;
    * number of reference frames increased by one for next P-frame
    */
    for(i = 0; i < 4; i++) free(mref_fld[i]);
  }

  free(img->mb_data);

  free(allbelta_lum);
  free(allalpha_lum);
#ifdef ROI_ENABLE
  for (i=0; i<3; i++)
  {
    free(YCbCr[i]);
  }
  free(ROIArray);
#endif
#ifdef FastME
  free_mem_FME();
#endif
}

/*
=======================================================================================================================
Function:Allocate memory for mv Input:Image Parameters struct img_par *img \n int_32_t mv Output: Return:
memory size in bytes Attention:
=======================================================================================================================
*/
int_32_t c_avs_enc::get_mem_mv(int_32_t ******mv)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k, l;
  /*~~~~~~~~~~~~~~~~~~~~~~~*/

  if((*mv = (int_32_t *****) calloc(2, sizeof(int_32_t ****))) == NULL) no_mem_exit("get_mem_mv: mv");
  for(i = 0; i < 2; i++)
  {
    if(((*mv)[i] = (int_32_t ****) calloc(2, sizeof(int_32_t ***))) == NULL)
      no_mem_exit("get_mem_mv: mv");
    for(j = 0; j < 2; j++)
    {
      if(((*mv)[i][j] = (int_32_t ***) calloc(img->buf_cycle, sizeof(int_32_t **))) == NULL)
        no_mem_exit("get_mem_mv: mv");
      for(k = 0; k < img->buf_cycle; k++)
      {
        if(((*mv)[i][j][k] = (int_32_t **) calloc(9, sizeof(int_32_t *))) == NULL)
          no_mem_exit("get_mem_mv: mv");
        for(l = 0; l < 9; l++)
          if(((*mv)[i][j][k][l] = (int_32_t *) calloc(2, sizeof(int_32_t))) == NULL)
            no_mem_exit("get_mem_mv: mv");
      }
    }
  }

  return 2 * 2 * img->buf_cycle * 9 * 2 * sizeof(int_32_t);
}

/*
=======================================================================================================================
Function:Free memory from mv Input:int_32_t mv Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_mem_mv(int_32_t *****mv)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k, l;
  /*~~~~~~~~~~~~~~~~~~~~~~~*/

  for(i = 0; i < 2; i++)
  {
    for(j = 0; j < 2; j++)
    {
      for(k = 0; k < img->buf_cycle; k++)
      {
        for(l = 0; l < 9; l++) free(mv[i][j][k][l]);
        free(mv[i][j][k]);
      }

      free(mv[i][j]);
    }

    free(mv[i]);
  }

  free(mv);
}

/*
=======================================================================================================================
Function:Allocate memory for AC coefficients Input: Output: Return: Attention:
=======================================================================================================================
*/
int_32_t c_avs_enc::get_mem_ACcoeff(int_16_t *****cofAC)
{
  /*~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k;
  /*~~~~~~~~~~~~~~~~~~~~*/

  /*
  * if ((*cofAC = (int_16_t )calloc (6, sizeof(int_16_t ))) == NULL) no_mem_exit
  * ("get_mem_ACcoeff: cofAC");
  */
  if((*cofAC = (int_16_t ****) _aligned_malloc(6 * sizeof(int_16_t ***), 16)) == NULL)
    no_mem_exit("get_mem_ACcoeff: cofAC");
  for(k = 0; k < 6; k++)
  {
    /*
    * if (((*cofAC)[k] = (int_16_t )calloc (4, sizeof(int_16_t**))) == NULL)
    * no_mem_exit ("get_mem_ACcoeff: cofAC");
    */
    if(((*cofAC)[k] = (int_16_t ***) _aligned_malloc(4 * sizeof(int_16_t **), 16)) == NULL)
      no_mem_exit("get_mem_ACcoeff: cofAC");
    for(j = 0; j < 4; j++)
    {
      /*
      * if (((*cofAC)[k][j] = (int_16_t**)calloc (2, sizeof(int_16_t*))) == NULL)
      * no_mem_exit ("get_mem_ACcoeff: cofAC");
      */
      if(((*cofAC)[k][j] = (int_16_t **) _aligned_malloc(2 * sizeof(int_16_t *), 16)) == NULL)
        no_mem_exit("get_mem_ACcoeff: cofAC");
      for(i = 0; i < 2; i++)
      {
        /*
        * if (((*cofAC)[k][j][i] = (int_16_t*)calloc (65, sizeof(int_16_t))) == NULL)
        * no_mem_exit ("get_mem_ACcoeff: cofAC");
        * // 18->65 for AVS
        */
        if
          (
          (
          (*cofAC)[k][j][i] = (int_16_t *) _aligned_malloc
          (
          65 * sizeof(int_16_t),
          16
          )
          ) == NULL
          ) no_mem_exit("get_mem_ACcoeff: cofAC");
      }
    }
  }

  return 6 * 4 * 2 * 65 * sizeof(int_16_t);  /* 18->65 for AVS */
}

/*
=======================================================================================================================
Function:Allocate memory for DC coefficients Input: Output: Return: Attention:
=======================================================================================================================
*/
int_32_t c_avs_enc::get_mem_DCcoeff(int_16_t ****cofDC)
{
  /*~~~~~~~~~~~~~~~~~*/
  int_32_t  j, k;
  /*~~~~~~~~~~~~~~~~~*/

  if((*cofDC = (int_16_t ***) calloc(3, sizeof(int_16_t **))) == NULL) no_mem_exit("get_mem_DCcoeff: cofDC");
  for(k = 0; k < 3; k++)
  {
    if(((*cofDC)[k] = (int_16_t **) calloc(2, sizeof(int_16_t *))) == NULL)
      no_mem_exit("get_mem_DCcoeff: cofDC");
    for(j = 0; j < 2; j++)
    {
      if(((*cofDC)[k][j] = (int_16_t *) calloc(65, sizeof(int_16_t))) == NULL)
        no_mem_exit("get_mem_DCcoeff: cofDC");  /* 18->65 for AVS */
    }
  }

  return 3 * 2 * 65 * sizeof(int_16_t);  /* 18->65 for AVS */
}

/*
=======================================================================================================================
Function:Allocate memory for reference frame Input: Output: Return: Attention:
=======================================================================================================================
*/
int_32_t c_avs_enc::get_mem_ref(byte *****ref)
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k;
  int_32_t  height = img->height + (IMG_PAD_SIZE << 1);
  int_32_t  width = img->width + (IMG_PAD_SIZE << 1);
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  if((*ref = (byte ****) _aligned_malloc(4 * sizeof(byte ***), 16)) == NULL) no_mem_exit("get_mem_ref: ref");
  for(k = 0; k < 4; k++)
  {
    if(((*ref)[k] = (byte ***) _aligned_malloc(4 * sizeof(byte **), 16)) == NULL)
      no_mem_exit("get_mem_ref: ref");
    for(j = 0; j < 4; j++)
    {
      if(((*ref)[k][j] = (byte **) _aligned_malloc(height * sizeof(byte *), 16)) == NULL)
        no_mem_exit("get_mem_ref: ref");
      for(i = 0; i < (img->height + (IMG_PAD_SIZE << 1)); i++)
      {
        if(((*ref)[k][j][i] = (byte *) _aligned_malloc(width * sizeof(byte), 16)) == NULL)
          no_mem_exit("get_mem_ref: ref");
      }
    }
  }

  return 4 * 4 * height * width * sizeof(byte);
}

/*
=======================================================================================================================
Function:Free memory of AC coefficients Input: Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_mem_ACcoeff(int_16_t ****cofAC)
{
  /*~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k;
  /*~~~~~~~~~~~~~~~~~~~~*/

  for(k = 0; k < 6; k++)
  {
    for(i = 0; i < 4; i++)
    {
      for(j = 0; j < 2; j++)
      {
        /*
        * free (cofAC[k][i][j]);
        */
        _aligned_free(cofAC[k][i][j]);
      }

      /*
      * free (cofAC[k][i]);
      */
      _aligned_free(cofAC[k][i]);
    }

    /*
    * free (cofAC[k]);
    */
    _aligned_free(cofAC[k]);
  }

  /*
  * free (cofAC);
  */
  _aligned_free(cofAC);
}

/*
=======================================================================================================================
Function:Free memory of DC coefficients Input: Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_mem_DCcoeff(int_16_t ***cofDC)
{
  /*~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j;
  /*~~~~~~~~~~~~~~~~~*/

  for(j = 0; j < 3; j++)
  {
    for(i = 0; i < 2; i++)
    {
      free(cofDC[j][i]);
    }

    free(cofDC[j]);
  }

  free(cofDC);
}

/*
=======================================================================================================================
Function:Free memory of reference frame Input: Output: Return: Attention:
=======================================================================================================================
*/
void c_avs_enc::free_mem_ref(byte ****ref)
{
  /*~~~~~~~~~~~~~~~~~~~~*/
  int_32_t  i, j, k;
  /*~~~~~~~~~~~~~~~~~~~~*/

  for(k = 0; k < 4; k++)
  {
    for(i = 0; i < 4; i++)
    {
      for(j = 0; j < (img->height + (IMG_PAD_SIZE << 1)); j++)
      {
        /*
        * free (cofAC[k][i][j]);
        */
        _aligned_free(ref[k][i][j]);
      }

      /*
      * free (cofAC[k][i]);
      */
      _aligned_free(ref[k][i]);
    }

    /*
    * free (cofAC[k]);
    */
    _aligned_free(ref[k]);
  }

  /*
  * free (cofAC);
  */
  _aligned_free(ref);
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::encode_IP_frame(avs_enc_frame_t *pFrame)
{
  TLS static int_32_t M, N, np, nb, n;
  TLS static int_32_t first_frm_flag = 1;
  int_32_t i;

  pInputImage = pFrame->input;
  pAVSMbInfo  = pFrame->pEncMBInfo;

  if(first_frm_flag == 1)
  {
    if(input->RCEnable)
    {
      rc_init_seq();
    }

#ifdef FastME
    DefineThreshold();
#endif

    /* B pictures */
    Bframe_ctr = 0;
    tot_time   = 0;  /* time for total encoding session */
    input->skip_mode_flag = 1;

    /* Write sequence header */
    if(!gframe_no)
      stat->sequence_header = start_sequence();
    first_frm_flag = 0;
    img->number = 0;
    tmp_buf_cycle = img->buf_cycle;
  }

  img->buf_cycle = tmp_buf_cycle;

  /* frame_num for this frame */
  if(pFrame->type[0] == AVS_TYPE_I)
  {
    img->type   = INTRA_IMG;
    img->number = current_encoded_frame;
    img->nb_references = 1;
  }
  else
  {
    img->type = INTER_IMG;
  }
  img->frame_num = img->number % (1 << (LOG2_MAX_FRAME_NUM_MINUS4 + 4));

  if((gframe_no%input->GopLength)>=input->GopLength - input->successive_Bframe)
  {
    picture_distance = gframe_no;
  }
  else
  {
    if(img->type==INTRA_IMG)
      picture_distance = gframe_no;
    else
      picture_distance = gframe_no+input->successive_Bframe;
  }
  if(input->RCEnable && img->type == INTRA_IMG)
  {
    M = input->successive_Bframe + 1;
    n = input->GopLength;; //Goplength和标准编码器中的IntraPeriod不同
    np = 0;
    nb = 0;

    for (i=1; i<n; i++)
    {
      if (input->successive_Bframe == 0 || i % M == 1 || (i>n-input->successive_Bframe && n%M!=1))
        np++;
      else
        nb++;
    }
    rc_init_GOP(np, nb);
  }
  encode_one_frame();  /* encode one I- or P-frame */
  if (img->number == 0)
  {
    img->nb_references = 1;
  }
  else
  {
    img->nb_references += 1;
    img->nb_references  = min(input->no_multpred, img->nb_references);
  }
  img->b_frame_to_code = 1;
  img->number++;
  memcpy((byte*)pFrame->bitstream + pFrame->length, pOutBuffer, nOutBufPtr);
  pFrame->length += nOutBufPtr;
  nOutBufPtr = 0;
  return 0;
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::encode_B_frame(avs_enc_frame_t *pFrame)
{
  pInputImage = pFrame->input;
  pAVSMbInfo  = pFrame->pEncMBInfo;

  img->number--;

  img->type = B_IMG;    /* set image type to B-frame */
  picture_coding_type = 1;

  img->frame_num++;    /* increment frame_num once for B-frames */
  img->frame_num %= (1 << (LOG2_MAX_FRAME_NUM_MINUS4 + 4));

  if(img->b_frame_to_code <= input->successive_Bframe * 2)
  {
    picture_distance = gframe_no-1;
    encode_one_frame();  /* encode one B-frame */
    img->b_frame_to_code++;
    memcpy((byte*)pFrame->bitstream + pFrame->length, pOutBuffer, nOutBufPtr);
    pFrame->length += nOutBufPtr;
    nOutBufPtr = 0;
  }

  img->number++;
  return 0;
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::avs_enc_frame(avs_enc_frame_t *pFrame)
{
  total_encoded_frame++;
  if(pFrame->type[0] == AVS_TYPE_B)
  {
    encode_B_frame(pFrame);
  }
  else
  {
    encode_IP_frame(pFrame);
  }

  return 0;
}

/*
=======================================================================================================================
=======================================================================================================================
*/
int_32_t c_avs_enc::init_global_variables()
{
  bytes_y = input->stuff_width * input->stuff_height;
  bytes_uv = bytes_y >> 2;
  framesize_in_bytes = bytes_y + (bytes_uv << 1);

  return 0;
}
