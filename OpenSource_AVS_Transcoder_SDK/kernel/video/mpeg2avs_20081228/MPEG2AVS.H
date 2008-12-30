#ifndef MPEG2AVS_H_
#define MPEG2AVS_H_
#include "typedef.h"
#include "global.h"
#define DLL_EXPORT __declspec(dllexport) //标志着这个这个函数将成为对外的接口

#ifdef WIN32
#pragma warning (disable:4996)
#endif

#ifdef _cplusplus
extern "C"{
#endif
	


DLL_EXPORT int_32_t trans_create   (trans_create_t * create, char* encoder_cfg_file_name, c_avs_enc* p_avs_enc);
DLL_EXPORT int_32_t trans_transcode(trans_create_t * create, int_32_t transcoding_type,   c_avs_enc* p_avs_enc);
DLL_EXPORT int_32_t trans_destroy  (trans_create_t *create,  c_avs_enc* p_avs_enc);

#ifdef _cplusplus
}
#endif  // _cplusplus
#endif