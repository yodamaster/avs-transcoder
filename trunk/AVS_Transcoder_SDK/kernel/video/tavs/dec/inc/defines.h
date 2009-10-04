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


#ifndef _DEFINES_H_
#define _DEFINES_H_

#define snprintf  _snprintf
/*
 * define FastME ;
 * #define FIELDINTE
 */
#define FastME
#define _THREE_STEP_MOTION_SEARCH_
#define AVS_OUT_BUFFER_SIZE  (1024 * 1024 * 4)

#define ET_SIZE        300  /* !< size of error text buffer */
#define NUM_PIC_TYPE     5

#define max(a, b)      (((a) > (b)) ? (a) : (b))
#define min(a, b)      (((a) < (b)) ? (a) : (b))
#define clamp(a, b, c)    ((a) < (b) ? (b) : ((a) > (c) ? (c) : (a)))  /* !< clamp a to the range of [b; c] */
#define Clip3(min, max, val)  (((val) < (min)) ? (min) : (((val) > (max)) ? (max) : (val)))

#define MAXMODE            13
#define INIT_FRAME_RATE    30
#define EOSEOS      1    /* !< End Of Sequence */

#define MAX_SYMBOLS_PER_MB  1200    /* !< Maximum number of different syntax elements for one MB */
#define STREAM_BUF_SIZE 1024 
#endif
