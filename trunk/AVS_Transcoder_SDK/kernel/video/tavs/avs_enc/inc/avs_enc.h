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
#ifndef _AVS_ENC_H_
#define _AVS_ENC_H_
#include "global.h"
#define DLL_EXPORT __declspec(dllexport) //标志着这个这个函数将成为对外的接口
DLL_EXPORT void avs_encoder_create  (c_avs_enc* p_avs_enc);
DLL_EXPORT void avs_encoder_encode  (c_avs_enc* p_avs_enc);
DLL_EXPORT void avs_encoder_destroy(c_avs_enc* p_avs_enc);
#endif