#ifndef _TRANSCODER_H_
#define _TRANSCODER_H_

//#define DLL_EXPORT __declspec(dllexport) //标志着这个这个函数将成为对外的接口

#include "transcoding_type.h"
#include "global.h"

//class c_mpeg2_dec
//  {
//  public:
//    mpeg2_dec_create_t *p_mpeg2_dec_create;
//    mpeg2_dec_frame_t  *p_mpeg2_dec_frame;
//    mpeg2_dec_stats_t  *p_mpeg2_dec_stats;
//    trans_frame_t      *p_trans_frame;    
//  };

//
//DLL_EXPORT void transcoder_create  (c_avs_enc** p_avs_enc, c_mpeg2_dec* p_mpeg2_dec, int thread_num);
//DLL_EXPORT void transcoder_decode  (c_mpeg2_dec* p_mpeg2_dec);
//DLL_EXPORT void transcoder_encode  (c_avs_enc* p_avs_enc);
//DLL_EXPORT void transcoder_destropy();
void avs_encoder_create  (c_avs_enc* p_avs_enc);
void avs_encoder_encode  (c_avs_enc* p_avs_enc);
void avs_encoder_destropy(c_avs_enc* p_avs_enc);
//int_32_t mpeg2_decoder_create( mpeg2_dec_create_t * create);
//int_32_t mpeg2_decoder_destroy();
//int_32_t mpeg2_decoder_decode( mpeg2_dec_frame_t * frame, mpeg2_dec_stats_t * stats);

#endif