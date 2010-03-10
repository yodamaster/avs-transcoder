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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "psych.h"
#include "av3enc.h"
#include "block_switch.h"

void PsyInit()
{
	check_short_init();
}

 /* Do psychoacoustical analysis */
void PsyCalculate(double *samplebuff[],ChannelInfo * channelInfo,
			 int *psy_signal_type,unsigned int numChannels)
{
  unsigned int channel;

  for (channel = 0; channel < numChannels; channel++)
  {
    if (channelInfo[channel].present)
    {
		if (channelInfo[channel].paired_ch == -1)
			check_short(channel,samplebuff,(WINDOW_TYPE *)&psy_signal_type[channel]);
		else if(channelInfo[channel].ch_is_left)
		{
			int leftChan = channel;
			int rightChan = channelInfo[channel].paired_ch;
			check_short(leftChan,samplebuff,(WINDOW_TYPE *)&psy_signal_type[leftChan]);
			check_short(rightChan,samplebuff,(WINDOW_TYPE *)&psy_signal_type[rightChan]);
		}
	}
  }
}


void BlockSwitch(CoderInfo * coderInfo, int *psy_block_type, unsigned int numChannels)
{
  unsigned int channel;
  int desire = ONLY_LONG_WINDOW;

  /* Use the same block type for all channels
     If there is 1 channel that wants a short block,
     use a short block on all channels.
   */
  for (channel = 0; channel < numChannels; channel++)
  {
    if (psy_block_type[channel] == ONLY_SHORT_WINDOW)
      desire = ONLY_SHORT_WINDOW;
  }

  for (channel = 0; channel < numChannels; channel++)
  {
    int lasttype = coderInfo[channel].block_type;

    if (desire == ONLY_SHORT_WINDOW	|| coderInfo[channel].desired_block_type == ONLY_SHORT_WINDOW)
    {
      if (lasttype == ONLY_LONG_WINDOW || lasttype == SHORT_LONG_WINDOW)
		coderInfo[channel].block_type = LONG_SHORT_WINDOW;
      else
		coderInfo[channel].block_type = ONLY_SHORT_WINDOW;
    }
    else
    {
      if (lasttype == ONLY_SHORT_WINDOW || lasttype == LONG_SHORT_WINDOW)
		coderInfo[channel].block_type = SHORT_LONG_WINDOW;
      else
		coderInfo[channel].block_type = ONLY_LONG_WINDOW;
    }
    coderInfo[channel].desired_block_type = desire;
  }
}


