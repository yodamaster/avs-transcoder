/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
*
* Lei Miao (win.miaolei@samsung.com), Samsung AIT
* Lei Miao, CBC Multi-channel extension, 2005-09-19
*
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

#ifndef _SAM_ENCODE_H_
#define _SAM_ENCODE_H_

#include <stdio.h>
#include "av3enc.h"

void sam_cbc_init(int fsidx, int bitrate, FILE *outFile);

int sam_cbc(
			ChannelInfo channelInfo[],
			int	nch,
			int sampleRateIdx,	
			int outputformat,
			int signal_type,
			int	windowShape,
			int num_window_groups,
			int	window_group_length[],
			short *swb_offset[],
			int	*scalefac[][8],
			int	*quant[][8],
			int	maxSfb,
			int stereo_mode,
			int	*stereo_info[],
			int frame_length,
			int wflag,
			int isExtended,
			int channelIdx,
			int numChannels,
			int mc_present,
			int extended_bytes,
			AV3EncFramePtr hEncoder);

#endif