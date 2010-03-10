/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
* edited by
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the users without any
* license fee or royalty on an "as is" basis. The AVS disclaims
* any and all warranties, whether express, implied, or statutory,
* including any implied warranties of merchantability or of fitness
* for a particular purpose. In no event shall the contributors or 
* the AVS be liable for any incidental, punitive, or consequential
* damages of any kind whatsoever arising from the use of this program.
*
* This disclaimer of warranty extends to the user of this program
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The AVS does not represent or warrant that the program furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy is available from the AVS Web site at
* http://www.avs.org.cn
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

#ifndef POLAR_STEREO
#define POLAR_STEREO

#ifdef __cplusplus
extern "C" {
#endif 

#include "av3enc.h"

#define LIMIT 1.9999
#define COSANGLE .9

void polar_stereo(CoderInfo*    coderInfo,
                  ChannelInfo*  channelInfo,
                  double*       spectrum[MAX_CHANNELS],
                  int           max_channel,
                  int           allowSquarePolar);

void polar_reconstruct(CoderInfo*   coderInfo,
                       ChannelInfo* channelInfo,
					   int          max_channel);

void PQSPSC_stereo(int		*quantFreqL, 
				   int		*quantFreqR, 
				   PSInfo	*psInfoM,
				   PSInfo	*psInfoA,
				   int		*sfb_offset, 
				   int		nr_of_sfb);
                
#ifdef __cplusplus
}
#endif 

#endif 