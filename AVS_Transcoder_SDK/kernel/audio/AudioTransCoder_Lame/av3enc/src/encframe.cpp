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
#include <assert.h>
#include <math.h>
#include <memory.h>

#include "av3enc.h"
#include "polar.h"
#include "av3enc.h"
#include "quantization.h"
#include "psych.h"
#include "sam_encode.h"
#include "intMDCT.h"
#include "block_switch.h"
#include "mdct.h"
#include "bew_flpvq.h"
#include "bew_mr_analysis.h"

#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) > (y) ? (y) : (x))

extern void FlushBufferInit();

/*base bandwidth for q=100*/
static const int bwbase = 16000;
/*bandwidth multiplier (for quantiser)*/
static const int bwmult = 120;
/*max bandwidth/samplerate ratio*/
static const double bwfac = 0.45;

static struct {
	 int rate; 
	 int cutoff;
} rates[] = {
	{29500, 5000},
	{37500, 7000},
	{47000, 10000},
	{64000, 16000},
	{76000, 20000},
	{0, 0}
};


/* Scalefactorband data table */
 
static SR_INFO srInfo[12+1] =
{
    { 96000, 41, 12,
        {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 16, 16, 24, 28,
            36, 44, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
        },{
            4, 4, 4, 4, 4, 4, 8, 8, 8, 16, 28, 36
        }
    }, { 88200, 41, 12,
        {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 16, 16, 24, 28,
            36, 44, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
        },{
            4, 4, 4, 4, 4, 4, 8, 8, 8, 16, 28, 36
        }
    }, { 64000, 47, 12,
        {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            8, 8, 8, 8, 12, 12, 12, 16, 16, 16, 20, 24, 24, 28,
            36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40
        },{
            4, 4, 4, 4, 4, 4, 8, 8, 8, 16, 28, 32
        }
    }, { 48000, 49, 14,
        {
            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8,  8,  8,  8,
            12, 12, 12, 12, 16, 16, 20, 20, 24, 24, 28, 28, 32, 32, 32, 32, 32, 32,
            32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 96
        }, {
            4,  4,  4,  4,  4,  8,  8,  8, 12, 12, 12, 16, 16, 16
        }
    }, { 44100, 49, 14,
        {
            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8,  8,  8,  8,
            12, 12, 12, 12, 16, 16, 20, 20, 24, 24, 28, 28, 32, 32, 32, 32, 32, 32,
            32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 96
        }, {
            4,  4,  4,  4,  4,  8,  8,  8, 12, 12, 12, 16, 16, 16
        }
    }, { 32000, 51, 14,
        {
            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8,
            8,  8,  8,  12, 12, 12, 12, 16, 16, 20, 20, 24, 24, 28,
            28, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
            32, 32, 32, 32, 32, 32, 32, 32, 32
        },{
            4,  4,  4,  4,  4,  8,  8,  8,  12, 12, 12, 16, 16, 16
        }
    }, { 24000, 47, 15,
        {
            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8,  8,  8,  8,
            8,  8,  8,  12, 12, 12, 12, 16, 16, 16, 20, 20, 24, 24, 28, 28, 32,
            36, 36, 40, 44, 48, 52, 52, 64, 64, 64, 64, 64
        }, {
            4,  4,  4,  4,  4,  4,  4,  8,  8,  8, 12, 12, 16, 16, 20
        }
    }, { 22050, 47, 15,
        {
            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  8,  8,  8,  8,  8,  8,  8,
            8,  8,  8,  12, 12, 12, 12, 16, 16, 16, 20, 20, 24, 24, 28, 28, 32,
            36, 36, 40, 44, 48, 52, 52, 64, 64, 64, 64, 64
        }, {
            4,  4,  4,  4,  4,  4,  4,  8,  8,  8, 12, 12, 16, 16, 20
        }
    }, { 16000, 43, 15,
        {
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 20, 20, 20, 24,
            24, 28, 28, 32, 36, 40, 40, 44, 48, 52, 56, 60, 64, 64, 64
        }, {
            4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 12, 12, 16, 20, 20
        }
    }, { 12000, 43, 15,
        {
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 20, 20, 20, 24,
            24, 28, 28, 32, 36, 40, 40, 44, 48, 52, 56, 60, 64, 64, 64
        }, {
            4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 12, 12, 16, 20, 20
        }
    }, { 11025, 43, 15,
        {
            8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 20, 20, 20, 24,
            24, 28, 28, 32, 36, 40, 40, 44, 48, 52, 56, 60, 64, 64, 64
        }, {
            4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 12, 12, 16, 20, 20
        }
    }, { 8000, 40, 15,
        {
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 16,
            16, 16, 16, 16, 16, 16, 20, 20, 20, 20, 24, 24, 24, 28,
            28, 32, 36, 36, 40, 44, 48, 52, 56, 60, 64, 80
        }, {
            4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 12, 16, 20, 20
        }
    },
    { -1 }
};

/* Returns the sample rate index */
int GetSRIndex(unsigned int sampleRate)
{
    if (92017 <= sampleRate) return 0;	//96k
    if (75132 <= sampleRate) return 1;	//88.2k
    if (55426 <= sampleRate) return 2;	//64k
    if (46009 <= sampleRate) return 3;	//48k
    if (37566 <= sampleRate) return 4;	//44.1k
    if (27713 <= sampleRate) return 5;	//32k
    if (23004 <= sampleRate) return 6;	//24k
    if (18783 <= sampleRate) return 7;	//22.05k
    if (13856 <= sampleRate) return 8;	//16k
    if (11502 <= sampleRate) return 9;	//12k
    if (9391 <= sampleRate) return 10;	//11.025k

    return 11;		//8k
}

int AV3EncSetConfiguration(AV3EncFramePtr hEncoder,
                                    AV3EncCfgPtr config)
{
	int i;
	
    if (config->bitRate > (6144.0 * (double)config->sampleRate/1024.0 + .5))
	  return 0;

	/* set quantization quality */
	hEncoder->quality = 100;

    if (config->bitRate && !config->bandWidth)
    {	
        int f0, f1;
		int r0, r1;
		// bit rate is not changed while over 79kbps/ch, YAN  2006-04-27
		int bitRateTemp;   

        bitRateTemp = (double)config->bitRate * 44100 / config->sampleRate;

		f0 = f1 = rates[0].cutoff;
		r0 = r1 = rates[0].rate;
		
		for (i = 0; rates[i].rate; i++)
		{
			f0 = f1;
			f1 = rates[i].cutoff;
			r0 = r1;
			r1 = rates[i].rate;
			if (rates[i].rate >= bitRateTemp)
				break;
		}
/*
		if (config->bitRate > r1)
			config->bitRate = r1;
		if (config->bitRate < r0)
			config->bitRate = r0;
*/
		if (f1 > f0)
			config->bandWidth =
					pow((double)bitRateTemp / r1,
					log((double)f1 / f0) / log ((double)r1 / r0)) * (double)f1;
		else
			config->bandWidth = f1;

        config->bandWidth =
				(double)config->bandWidth * config->sampleRate / 44100;
//				config->bitRate = (double)config->bitRate * config->sampleRate / 44100;

		if (config->bandWidth > bwbase)
		  config->bandWidth = bwbase;
	}

    if (!config->bandWidth)
    {
		config->bandWidth = (hEncoder->quality - 100) * bwmult + bwbase;
    }

	if(hEncoder->config.isLFE && hEncoder->config.LFE_3SFB)
		hEncoder->config.desbits = (hEncoder->config.numChannels-1) * (hEncoder->config.bitRate * 1024)
				/ hEncoder->config.sampleRate;
	else
		hEncoder->config.desbits = hEncoder->config.numChannels * (hEncoder->config.bitRate * 1024)
				/ hEncoder->config.sampleRate;

	return 1;
}

AV3EncFramePtr AV3EncOpen(unsigned long sampleRate,
                                  unsigned int numChannels,
                                  unsigned long *inputSamples,
                                  unsigned long *maxOutputBytes)
{
    unsigned int channel;
    AV3EncFramePtr hEncoder;

    *inputSamples = 1024*numChannels;
    *maxOutputBytes = (6144/8)*numChannels;

    hEncoder = (AV3EncFrame*)malloc(sizeof(AV3EncFrame));
    memset(hEncoder, 0, sizeof(AV3EncFrame));

    hEncoder->config.numChannels = numChannels;
    hEncoder->config.sampleRate = sampleRate;
    hEncoder->config.sampleRateIdx = GetSRIndex(sampleRate);

    /* Initialize variables to default values */
    hEncoder->frameNum = 0;
    hEncoder->flushFrame = 0;

    /* Default configuration */
	hEncoder->config.allowFLPVQ = 1;
    hEncoder->config.allowSPSC = 1;
    hEncoder->config.bitRate = 0; /* default bitrate / channel */
    hEncoder->config.bandWidth = bwfac * hEncoder->config.sampleRate;
    if (hEncoder->config.bandWidth > bwbase)
		hEncoder->config.bandWidth = bwbase;
    hEncoder->config.shortctl = SHORTCTL_NORMAL;

	/* default channel map is straight-through */
	for( channel = 0; channel <8; channel++ )
		hEncoder->config.channel_map[channel] = channel;
	
    /*
        by default we have to be compatible with all previous software
        which assumes that we will generate ADTS
        /AV
    */
    hEncoder->config.outputFormat = 0;

    /* find correct sampling rate depending parameters */
    hEncoder->srInfo = &srInfo[hEncoder->config.sampleRateIdx];

    for (channel = 0; channel < numChannels; channel++) 
	{
        hEncoder->coderInfo[channel].window_shape = SINE_WINDOW;
        hEncoder->coderInfo[channel].block_type = ONLY_LONG_WINDOW;
        hEncoder->coderInfo[channel].num_window_groups = 1;
        hEncoder->coderInfo[channel].window_group_length[0] = 1;

        hEncoder->freqBuff[channel] = (double*)malloc(2*FRAME_LEN*sizeof(double));
        hEncoder->overlapBuff[channel] = (double*)malloc(FRAME_LEN*sizeof(double));
        memset(hEncoder->overlapBuff[channel], 0, FRAME_LEN*sizeof(double));
   
        hEncoder->sampleBuff[channel] = NULL;
        hEncoder->nextSampleBuff[channel] = NULL;
    }

    /* Initialize coder functions */
	PsyInit();
	check_short_init();
	InitFLPVQStructure( hEncoder->config.allowFLPVQ,
		hEncoder->channelInfo,
		hEncoder->config.sampleRateIdx,
		hEncoder->config.numChannels );

	InitMultiResolution( hEncoder->channelInfo,
		hEncoder->config.numChannels,
		hEncoder->config.sampleRate,
        hEncoder->srInfo->num_cb_short,
		hEncoder->srInfo->cb_width_short );
	
	

    quantizeInit(hEncoder->coderInfo, hEncoder->config.numChannels);

	/* Return handle */
    return hEncoder;
}

int AV3EncClose(AV3EncFramePtr hEncoder)
{
    unsigned int channel;

 	quantizeEnd(hEncoder->coderInfo, hEncoder->config.numChannels);

    /* Free remaining buffer memory */
    for (channel = 0; channel < hEncoder->config.numChannels; channel++) 
	{
	    if (hEncoder->freqBuff[channel]) 
	        free(hEncoder->freqBuff[channel]);
	    if (hEncoder->overlapBuff[channel]) 
            free(hEncoder->overlapBuff[channel]); 
		if (hEncoder->sampleBuff[channel])
			free(hEncoder->sampleBuff[channel]);
		if (hEncoder->nextSampleBuff[channel])
			free(hEncoder->nextSampleBuff[channel]);
    }

    /* Free handle */
    if (hEncoder) 
		free(hEncoder);

    return 0;
}

void GetChannelInfo(ChannelInfo *channelInfo, int numChannels)
{
    int numChannelsLeft = numChannels;

	/* mono */
    if(numChannelsLeft == 1) {
        channelInfo[numChannels-numChannelsLeft].present = 1;
		/* no paired channel */
		channelInfo[numChannels-numChannelsLeft].paired_ch = -1;
        numChannelsLeft--;
    }
    /* Stereo, multi-channel extension */
    else if(numChannelsLeft >= 2) {
        /* Left Front channel info */
        channelInfo[numChannels-numChannelsLeft].present = 1;
        channelInfo[numChannels-numChannelsLeft].common_window = 0;
        channelInfo[numChannels-numChannelsLeft].ch_is_left = 1;
        channelInfo[numChannels-numChannelsLeft].paired_ch = numChannels-numChannelsLeft+1;
        numChannelsLeft--;

        /* Right Front channel info */
        channelInfo[numChannels-numChannelsLeft].present = 1;
        channelInfo[numChannels-numChannelsLeft].common_window = 0;
        channelInfo[numChannels-numChannelsLeft].ch_is_left = 0;
        channelInfo[numChannels-numChannelsLeft].paired_ch = numChannels-numChannelsLeft-1;
        numChannelsLeft--;

		/* multi-channel extension */
		if(numChannelsLeft) {
			/* Center channel info */
			channelInfo[numChannels-numChannelsLeft].present = 1;
			/* no paired channel */
			channelInfo[numChannels-numChannelsLeft].paired_ch = -1;
			numChannelsLeft--;

			while(numChannelsLeft > 1) {
				/* Left Surround channel info */
				channelInfo[numChannels-numChannelsLeft].present = 1;
				channelInfo[numChannels-numChannelsLeft].common_window = 0;
				channelInfo[numChannels-numChannelsLeft].ch_is_left = 1;
				channelInfo[numChannels-numChannelsLeft].paired_ch = numChannels-numChannelsLeft+1;
				numChannelsLeft--;

				/* Right Surround channel info */
				channelInfo[numChannels-numChannelsLeft].present = 1;
				channelInfo[numChannels-numChannelsLeft].common_window = 0;
				channelInfo[numChannels-numChannelsLeft].ch_is_left = 0;
				channelInfo[numChannels-numChannelsLeft].paired_ch = numChannels-numChannelsLeft-1;
				numChannelsLeft--;
			}
		}

		/* LFE channel */
		if(numChannelsLeft)
		{
			channelInfo[numChannels-1].present = 1;
			/* no paired channel */
			channelInfo[numChannels-1].paired_ch = -1;
		}
    }
}

int AV3EncEncode(AV3EncFramePtr hEncoder,
                          int *inputBuffer,
                          unsigned int samplesInput,
                          unsigned char *outputBuffer)
{
    unsigned int channel, i;
    int sb, frameBytes, temp;	
    unsigned int offset;

    /* local copy's of parameters */
    ChannelInfo *channelInfo = hEncoder->channelInfo;
    CoderInfo *coderInfo = hEncoder->coderInfo;
    unsigned int numChannels = hEncoder->config.numChannels;
    unsigned int sampleRate = hEncoder->config.sampleRate;
    unsigned int AVSVersion = hEncoder->config.AVSVersion;
    unsigned int allowSPSC= hEncoder->config.allowSPSC;
    unsigned int bandWidth = hEncoder->config.bandWidth;
    unsigned int shortctl = hEncoder->config.shortctl;
	unsigned int outputformat = hEncoder->config.outputFormat;

    /* Increase frame number */
    hEncoder->frameNum++;

    if (samplesInput == 0)
        hEncoder->flushFrame++;

    /* After 2 flush frames all samples have been encoded,
       return 0 bytes written */
    if (hEncoder->flushFrame > 2) 
        return 0;

    /* Determine the channel configuration */
    GetChannelInfo(channelInfo, numChannels);

    /* Update current sample buffers & remap channels*/
    for (channel = 0; channel < numChannels; channel++) 
	{
		double *tmp;

		if (!hEncoder->sampleBuff[channel])
			hEncoder->sampleBuff[channel] = (double*)malloc(FRAME_LEN*sizeof(double));
		
		tmp = hEncoder->sampleBuff[channel];
        hEncoder->sampleBuff[channel] = hEncoder->nextSampleBuff[channel];
		hEncoder->nextSampleBuff[channel] = tmp;

        if (samplesInput == 0)
        {
            /* start flushing*/
            for (i = 0; i < FRAME_LEN; i++)
                hEncoder->nextSampleBuff[channel][i] = 0.0;
        }
        else
        {
			int samples_per_channel = samplesInput/numChannels;

			{
			  float *input_channel;
			  
			  if(!hEncoder->config.isLFE && channel>2)
				  input_channel = (float*)inputBuffer + hEncoder->config.channel_map[channel]-1;
			  else
				  input_channel = (float*)inputBuffer + hEncoder->config.channel_map[channel];

		      for (i = 0; i < samples_per_channel; i++)
				{
					hEncoder->nextSampleBuff[channel][i] = (double)*input_channel;
					input_channel += numChannels;
				}

			}

            for (i = (int)(samplesInput/numChannels); i < FRAME_LEN; i++)
                hEncoder->nextSampleBuff[channel][i] = 0.0;
		}
    }

    if (hEncoder->frameNum <= 1) /* Still filling up the buffers */
        return 0;

    /* Psychoacoustics */
	PsyCalculate(hEncoder->nextSampleBuff,channelInfo,hEncoder->psy_block_type,numChannels);

	for (channel = 0; channel < numChannels; channel++)
	{
		if (hEncoder->psy_block_type[channel] == ONLY_LONG_WINDOW)
			hEncoder->coderInfo[channel].signal_type = STATIONARY_TYPE;
		else
			hEncoder->coderInfo[channel].signal_type = TRANSIENT_TYPE;
	}

	/* set common type for channel pair */
	if ((numChannels == 2 || numChannels == 6)&&(hEncoder->coderInfo[0].signal_type == TRANSIENT_TYPE || 
		hEncoder->coderInfo[1].signal_type == TRANSIENT_TYPE))
	{
		hEncoder->coderInfo[0].signal_type = TRANSIENT_TYPE;
		hEncoder->coderInfo[1].signal_type = TRANSIENT_TYPE;
	}
	if ((numChannels == 6) && (hEncoder->coderInfo[3].signal_type == TRANSIENT_TYPE ||
		hEncoder->coderInfo[4].signal_type == TRANSIENT_TYPE))
	{
		hEncoder->coderInfo[3].signal_type = TRANSIENT_TYPE;
		hEncoder->coderInfo[4].signal_type = TRANSIENT_TYPE;
	}
//	/* 2005-12-8 xuhengu */
//	hEncoder->coderInfo[0].signal_type = STATIONARY_TYPE;
//	hEncoder->coderInfo[1].signal_type = STATIONARY_TYPE;
	
	
	/* force signal type */
    if (shortctl == SHORTCTL_NOSHORT)
    {
		for (channel = 0; channel < numChannels; channel++)
		{
			coderInfo[channel].signal_type = STATIONARY_TYPE;
		}
    }
    if (shortctl == SHORTCTL_NOLONG)
    {
		for (channel = 0; channel < numChannels; channel++)
		{
			coderInfo[channel].signal_type = TRANSIENT_TYPE;
		}
    }

	/* LFE: force to stationary type */
	if(hEncoder->config.isLFE==1 && hEncoder->config.LFE_3SFB) {
		coderInfo[numChannels-1].signal_type = STATIONARY_TYPE;
	}
    /* AV3 Filterbank, MDCT with overlap and add */
    for (channel = 0; channel < numChannels; channel++) {
        int lowpass, i;
		
		FilterBank(hEncoder->sampleBuff[channel],
				hEncoder->overlapBuff[channel],
				hEncoder->freqBuff[channel]);
		/* calculate the last line which is not zero */
		lowpass = (bandWidth * BLOCK_LEN_LONG) / (sampleRate>>1) + 1;
		
		for (i = lowpass; i < BLOCK_LEN_LONG; i++)
			hEncoder->freqBuff[channel][i] = 0.0;		
    }

    /* TMP: Build sfb offset table and other stuff */
    for (channel = 0; channel < numChannels; channel++) {

        channelInfo[channel].psInfo.present = 0;

		if (coderInfo[channel].signal_type == TRANSIENT_TYPE) {
			coderInfo[channel].max_sfb = hEncoder->srInfo->num_cb_short;
            coderInfo[channel].nr_of_sfb = hEncoder->srInfo->num_cb_short;
			
            coderInfo[channel].swb = hEncoder->srInfo->num_cb_short;

            coderInfo[channel].num_window_groups = 3;
            coderInfo[channel].window_group_length[0] = 3;
            coderInfo[channel].window_group_length[1] = 3;
            coderInfo[channel].window_group_length[2] = 2;
            coderInfo[channel].window_group_length[3] = 0;
            coderInfo[channel].window_group_length[4] = 0;
            coderInfo[channel].window_group_length[5] = 0;
            coderInfo[channel].window_group_length[6] = 0;
            coderInfo[channel].window_group_length[7] = 0;

            offset = 0;
            for (sb = 0; sb < coderInfo[channel].nr_of_sfb; sb++) {
                coderInfo[channel].sfb_offset[sb] = offset;
                offset += hEncoder->srInfo->cb_width_short[sb];
            }
            coderInfo[channel].sfb_offset[coderInfo[channel].nr_of_sfb] = offset;
        } else {
            coderInfo[channel].max_sfb = hEncoder->srInfo->num_cb_long;
            coderInfo[channel].nr_of_sfb = hEncoder->srInfo->num_cb_long;

			coderInfo[channel].swb = hEncoder->srInfo->num_cb_long;

            coderInfo[channel].num_window_groups = 1;
            coderInfo[channel].window_group_length[0] = 1;

            offset = 0;
            for (sb = 0; sb < coderInfo[channel].nr_of_sfb; sb++) {
                coderInfo[channel].sfb_offset[sb] = offset;
                offset += hEncoder->srInfo->cb_width_long[sb];
            }
            coderInfo[channel].sfb_offset[coderInfo[channel].nr_of_sfb] = offset;
        }
    }
	/* FLPVQ encode */

	for(channel=0; channel< numChannels; channel++) {
		ChannelInfo* pChInfo = &(hEncoder->channelInfo[channel]);
		if (hEncoder->config.allowFLPVQ /*&& coderInfo[channel].signal_type==STATIONARY_TYPE*/) {
			FlpVQEncode(&(pChInfo->flpInfo), 
				coderInfo[channel].nr_of_sfb,
				coderInfo[channel].sfb_offset,
				hEncoder->freqBuff[channel]);
		}
	}	
	
	/* Multi-Resolution */
	for(channel=0; channel< numChannels; channel++) {
		ChannelInfo* pChInfo = &(hEncoder->channelInfo[channel]);
		MultiResolutionEncode(&(pChInfo->mrInfo),
								hEncoder->coderInfo[channel].signal_type,
								hEncoder->freqBuff[channel]);
	}			
	
    for (channel = 0; channel < numChannels; channel++) {
		if (coderInfo[channel].signal_type == TRANSIENT_TYPE) {
			SortForGrouping(&coderInfo[channel],
					&channelInfo[channel],
					hEncoder->srInfo->cb_width_short,
					hEncoder->freqBuff[channel]);
		}
		CalcAvgEnrg(&coderInfo[channel], hEncoder->freqBuff[channel]);
	}

	/* LFE limitation */
	if(hEncoder->config.isLFE==1 && hEncoder->config.LFE_3SFB)
		coderInfo[numChannels-1].nr_of_sfb = coderInfo[numChannels-1].max_sfb = 3;

LOOP_BUFFER_CONTROL:
   /* Quantize and code the signal */
    for (channel = 0; channel < numChannels; channel++) {
		quantize(&coderInfo[channel], 
			&channelInfo[channel],
			hEncoder->freqBuff[channel],
			hEncoder->quality);
    }
	
	// move this part below in order to enable PQSPSC in multichannel case
	// xun 2006-4-16
//	if (allowSPSC && numChannels == 2)
//		PQSPSC_stereo( coderInfo[0].quantFreq, 
//					   coderInfo[1].quantFreq, 
//					   &channelInfo[0].psInfo,
//					   &channelInfo[1].psInfo,
//					   coderInfo->sfb_offset,
//					   coderInfo->nr_of_sfb);
//	else
//		for (channel = 0; channel < numChannels; channel++)
//			channelInfo[channel].psInfo.present = 0;

	{
		int	*sample[2][8];
		int	sf_buffer[2][8][100];
		int	*scalefactors[2][8];
		short swb_offset_buf[8][64];
		short *swb_offset[8];
		int	stereo_mode;
		int	stereo_info_buf[8*16];
		int	*stereo_info[8];
		int maxSfb;
		int numChannelsLeft;
		int numChannelsCoded;
		int isExtended;
		int mc_present;

		ChannelInfo *localChannelInfo;
		CoderInfo *localCoderInfo;

		numChannelsLeft = numChannels;
		if(numChannels > 2)
			mc_present = 1;
		else
			mc_present = 0;
		frameBytes = 0;
		isExtended = 0;

		FlushBufferInit();

		if(numChannels > 2)
		{
			/* extended CBC raw data */
			numChannelsLeft -= 2;
			isExtended = 1;
		}
		else
			isExtended = 0;
		while(numChannelsLeft && isExtended) {
			localChannelInfo = &channelInfo[numChannels-numChannelsLeft];
			localCoderInfo = &coderInfo[numChannels-numChannelsLeft];

			/* current encoded channel number */
			if(localChannelInfo->paired_ch == -1)
				numChannelsCoded = 1;               /* no pair channel */
			else
				numChannelsCoded = 2;

			/* the same window used or not */
			if(numChannelsCoded == 2) {
				if(localCoderInfo->block_type != (localCoderInfo+1)->block_type)  {
					printf("Error in common window\n");
					exit(0);
				}
			}
			
			// move this part here to enable PQSPSC in multichannel case
			// xun 2006-4-16
			if (allowSPSC && numChannelsCoded == 2)
				PQSPSC_stereo( localCoderInfo[0].quantFreq, 
							   localCoderInfo[1].quantFreq, 
							   &localChannelInfo[0].psInfo,
							   &localChannelInfo[1].psInfo,
							   localCoderInfo[0].sfb_offset,
							   localCoderInfo[0].nr_of_sfb);
			else
				for (channel = 0; channel < numChannelsCoded; channel++)
					localChannelInfo[channel].psInfo.present = 0;

			for (channel = 0; channel < numChannelsCoded; channel++) {
				if ((localCoderInfo+channel)->signal_type == TRANSIENT_TYPE) {
					ReSortForGrouping(localCoderInfo+channel,
							localChannelInfo+channel,
							hEncoder->srInfo->cb_width_short,
							(localCoderInfo+channel)->quantFreq);

				}
			}

			for (channel = 0; channel < numChannelsCoded; channel++) {
				int length = 0;
				/* 2005-11-28 xuhengyu bug*/
//				for(i = 0; i < coderInfo->num_window_groups; i++)
				for(i = 0; i < localCoderInfo->num_window_groups; i++)
				{
					sample[channel][i] = (localCoderInfo+channel)->quantFreq+length*128;
					length += (localCoderInfo+channel)->window_group_length[i];
				}
			}

			for (channel = 0; channel < numChannelsCoded; channel++) {
				int length = 0;
				for(i = 0; i < (localCoderInfo+channel)->num_window_groups; i++)
				{
					int sfb_temp;
					for(sfb_temp=0;sfb_temp<(localCoderInfo+channel)->nr_of_sfb;sfb_temp++)
					{
						/* 2005-11-26 xuhengyu */
//						sf_buffer[channel][i][sfb_temp] = (localCoderInfo+channel)->scale_factor[(localCoderInfo+channel)->nr_of_sfb*length+sfb_temp];
						sf_buffer[channel][i][sfb_temp] = (localCoderInfo+channel)->scale_factor[(localCoderInfo+channel)->nr_of_sfb*i+sfb_temp];
					}
					scalefactors[channel][i] = sf_buffer[channel][i];
					length += (localCoderInfo+channel)->window_group_length[i];
				}
			}
			{
				int length = 0;
				for(i = 0; i < localCoderInfo->num_window_groups; i++)
				{
					int sfb_temp;
					for(sfb_temp=0;sfb_temp<localCoderInfo->nr_of_sfb;sfb_temp++)
					{
						swb_offset_buf[i][sfb_temp] = localCoderInfo->sfb_offset[sfb_temp]*localCoderInfo->window_group_length[i];

						/* 2005-11-29 xuhengyu */
						if(localChannelInfo->psInfo.present)   
							stereo_info_buf[16*i+sfb_temp] = localChannelInfo->psInfo.ps_used[i*localCoderInfo->nr_of_sfb+sfb_temp];
//							stereo_info_buf[16*i+sfb_temp] = localChannelInfo->psInfo.ps_used[length*localCoderInfo->nr_of_sfb+sfb_temp];
					}
					/* the last scalefactor band */
					swb_offset_buf[i][sfb_temp] = localCoderInfo->sfb_offset[sfb_temp]*localCoderInfo->window_group_length[i];
					length += localCoderInfo->window_group_length[i];
				}
			}

			stereo_mode = localChannelInfo->psInfo.present;
			
			for(i = 0; i < localCoderInfo->num_window_groups; i++)
			{
				swb_offset[i] = swb_offset_buf[i];
				stereo_info[i] = stereo_info_buf+16*i;
			}

			if(hEncoder->config.channel_map[numChannels-numChannelsLeft]!=3) {
				for(maxSfb=localCoderInfo->nr_of_sfb-1; ;maxSfb--)
				{
					if(scalefactors[0][0][maxSfb] != 90)
						break;
				}
				maxSfb += 1;
				if(localCoderInfo->signal_type == TRANSIENT_TYPE)
				{
					if(maxSfb < 10)
						maxSfb = 10;
				}
				else
				{
					if(maxSfb < 38)
						maxSfb = 38;
				}
			}
			else {
				maxSfb = localCoderInfo->nr_of_sfb;
			}


			/* extendec CBC encoding */
			temp = 
			sam_cbc(hEncoder->channelInfo, numChannelsCoded, hEncoder->config.sampleRateIdx, 
					outputformat, localCoderInfo->signal_type, 0, localCoderInfo->num_window_groups,
					localCoderInfo->window_group_length, swb_offset, scalefactors, sample,
					maxSfb, stereo_mode, stereo_info, 20000, 1, isExtended,
					hEncoder->config.channel_map[numChannels-numChannelsLeft], 
					numChannels, hEncoder->config.isLFE, 0,					
					hEncoder				    
					);

			if (temp == -1)
			{
				hEncoder->quality *= 0.707;
				for (channel = 0; channel < numChannels; channel++) {
					if (coderInfo[channel].signal_type == TRANSIENT_TYPE) {
						SortForGrouping(&coderInfo[channel],
										&channelInfo[channel],
										hEncoder->srInfo->cb_width_short,
										NULL);	
					}
				}
				goto LOOP_BUFFER_CONTROL;
			}

			frameBytes += temp /8;
			numChannelsLeft -= numChannelsCoded;
		}

		isExtended = 0;

		/* stereo/mono coding */
		localChannelInfo = &channelInfo[0];
		localCoderInfo = &coderInfo[0];

		/* current encoded channel number */
		if(localChannelInfo->paired_ch == -1)
			numChannelsCoded = 1;               /* no pair channel */
		else
			numChannelsCoded = 2;

		/* the same window used or not */
		if(numChannelsCoded == 2) {
			if (localCoderInfo->signal_type != (localCoderInfo + 1)->signal_type) {
				printf("Error in common window\n");
				exit(0);
			}
		}

		// move this part here to enable PQSPSC in multichannel case
		// xun 2006-4-16
		if (allowSPSC && numChannelsCoded == 2)
			PQSPSC_stereo( localCoderInfo[0].quantFreq, 
						   localCoderInfo[1].quantFreq, 
						   &localChannelInfo[0].psInfo,
						   &localChannelInfo[1].psInfo,
						   localCoderInfo[0].sfb_offset,
						   localCoderInfo[0].nr_of_sfb);
		else
			for (channel = 0; channel < numChannelsCoded; channel++)
				localChannelInfo[channel].psInfo.present = 0;
		
		for (channel = 0; channel < numChannelsCoded; channel++) {
			if ((localCoderInfo + channel)->signal_type == TRANSIENT_TYPE) {
				ReSortForGrouping(localCoderInfo+channel,
						localChannelInfo+channel,
						hEncoder->srInfo->cb_width_short,
						(localCoderInfo+channel)->quantFreq);
			}
		}

		for (channel = 0; channel < numChannelsCoded; channel++) {
			int length = 0;
			/* 2005-11-28 xuhengyu bug*/
//			for(i = 0; i < coderInfo->num_window_groups; i++)
			for(i = 0; i < localCoderInfo->num_window_groups; i++)
			{
				sample[channel][i] = (localCoderInfo+channel)->quantFreq+length*128;
				length += (localCoderInfo+channel)->window_group_length[i];
			}
		}

		for (channel = 0; channel < numChannelsCoded; channel++) {
			int length = 0;
			for(i = 0; i < (localCoderInfo+channel)->num_window_groups; i++)
			{
				int sfb_temp;
				/* 2005-11-26 xuhengyu */
//				for(sfb_temp=0;sfb_temp<(localCoderInfo+channel)->nr_of_sfb;sfb_temp++)
//				{
//					sf_buffer[channel][i][sfb_temp] = 
//						(localCoderInfo+channel)->scale_factor[(localCoderInfo+channel)->nr_of_sfb*length+sfb_temp];
//				}
				for(sfb_temp=0;sfb_temp<(localCoderInfo+channel)->nr_of_sfb;sfb_temp++)
				{
					sf_buffer[channel][i][sfb_temp] = 
						(localCoderInfo+channel)->scale_factor[(localCoderInfo+channel)->nr_of_sfb*i+sfb_temp];
				}
				scalefactors[channel][i] = sf_buffer[channel][i];
				length += (localCoderInfo+channel)->window_group_length[i];
			}
		}
		
		{			
			for(i = 0; i < localCoderInfo->num_window_groups; i++)
			{
				int sfb_temp;
				for(sfb_temp=0;sfb_temp<localCoderInfo->nr_of_sfb;sfb_temp++)
				{
					swb_offset_buf[i][sfb_temp] = 
						localCoderInfo->sfb_offset[sfb_temp]*localCoderInfo->window_group_length[i];

					/* 2005-11-29 xuhengyu */
					if(localChannelInfo->psInfo.present)   
						stereo_info_buf[16*i+sfb_temp] = 
						localChannelInfo->psInfo.ps_used[i*localCoderInfo->nr_of_sfb+sfb_temp];
				}
				/* the last scalefactor band */
				swb_offset_buf[i][sfb_temp] = 
					localCoderInfo->sfb_offset[sfb_temp]*localCoderInfo->window_group_length[i];				
			}	
		}


		stereo_mode = localChannelInfo->psInfo.present;
		
		for(i = 0; i < localCoderInfo->num_window_groups; i++)
		{
			swb_offset[i] = swb_offset_buf[i];
			stereo_info[i] = stereo_info_buf+16*i;
		}

		for(maxSfb=localCoderInfo->nr_of_sfb-1; ;maxSfb--)
		{
			if(scalefactors[0][0][maxSfb] != 90)
				break;
		}
		maxSfb += 1;
		if (localCoderInfo->signal_type == TRANSIENT_TYPE) {
			if(maxSfb < 10)
				maxSfb = 10;
		}
		else {
			if(maxSfb < 38)
				maxSfb = 38;
		}

		/* stereo/mono encoding */
		temp = 
		sam_cbc(hEncoder->channelInfo, numChannelsCoded, hEncoder->config.sampleRateIdx, 
				outputformat, localCoderInfo->signal_type, 0, localCoderInfo->num_window_groups,
				localCoderInfo->window_group_length, swb_offset, scalefactors, sample,
				maxSfb, stereo_mode, stereo_info, 20000, 1, isExtended, 0, numChannels,
				mc_present, frameBytes,				
				hEncoder
				);

		if (temp == -1)
		{
			hEncoder->quality *= 0.707;
			for (channel = 0; channel < numChannels; channel++) {
				if (coderInfo[channel].signal_type == TRANSIENT_TYPE) {
					SortForGrouping(&coderInfo[channel],
									&channelInfo[channel],
									hEncoder->srInfo->cb_width_short,
									NULL);	
				}
			}
			goto LOOP_BUFFER_CONTROL;
		}
		frameBytes += temp / 8;
	}
    /* Adjust quality to get correct average bitrate */
    if (hEncoder->config.bitRate)
	{
		double fix;
        int diff = (frameBytes * 8) - hEncoder->config.desbits;
		
		hEncoder->bitDiff += diff;
      
		fix = (double)hEncoder->bitDiff / hEncoder->config.desbits;
		fix *= 0.01;
		fix = max(fix, -0.2);
		fix = min(fix, 0.2);

		if ( diff*fix > 0.0 )
		{
           hEncoder->quality *= (1.0 - fix);
            hEncoder->quality = min(hEncoder->quality, 300);
			hEncoder->quality = max(hEncoder->quality, 50);
		}
    }

	// buffer status updates - Xun 2006-4-16
	hEncoder->fltBufferFullness -= (float)(frameBytes * 8);
	// simulate sending
	hEncoder->fltBufferFullness += hEncoder->fltAvgBitPerFrm;

    return frameBytes;
}


