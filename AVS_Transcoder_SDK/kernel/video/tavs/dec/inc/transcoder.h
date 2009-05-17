#ifndef _TRANSCODER_H_
#define _TRANSCODER_H_

//#define DLL_EXPORT __declspec(dllexport) //标志着这个这个函数将成为对外的接口

#include "transcoding_type.h"
#include "global.h"

void avs_encoder_create  (c_avs_enc* p_avs_enc);
void avs_encoder_encode  (c_avs_enc* p_avs_enc);
void avs_encoder_destroy(c_avs_enc* p_avs_enc);

#endif
