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

#include <math.h>
#include "polar.h"
#include "av3enc.h"

/* caculate the average energy of a given vector upto its last nonzero item */
static double averEnrg(double* spectrum, int length)
{
    double sumEnrg = 0;
    int i;
    int lastx;
    
    lastx = 0;
    for(i = 0; i < length; i++)
    {
        if(spectrum[i]){
            sumEnrg += spectrum[i]*spectrum[i];
            /* find the last nonzero item of the vector */
			lastx = i;
        }
    }
    
    return sumEnrg/(++lastx);
}
 
/* caculate the cos(angle between two vectors) */
static double cosAngle(double* spectrum_pair0, double* spectrum_pair1, int length)
{
    double crossProduct;
    double norm_pair0;
    double norm_pair1;
    int i;

    crossProduct = 0;
    norm_pair0 = 0;
    norm_pair1 = 0;
    
    for(i = 0; i < length; i++)
    {
        crossProduct += spectrum_pair0[i] * spectrum_pair1[i];
        norm_pair0 += spectrum_pair0[i] * spectrum_pair0[i];
        norm_pair1 += spectrum_pair1[i] * spectrum_pair1[i];
    }
    
    if(norm_pair0 != 0 && norm_pair1 != 0)
        return crossProduct/sqrt(norm_pair0*norm_pair1);
    else
        return -2; /*signal an error*/    
    
}

/* check for the consistence of scale_factor_grouping of a channel pair */
static int common_grouping(CoderInfo* coderInfoL, CoderInfo* coderInfoR)
{
    if(coderInfoL->window_shape != coderInfoR->window_shape)
        return 0;
    if(coderInfoL->block_type != coderInfoR->block_type)
        return 0;
    if(coderInfoL->block_type == ONLY_SHORT_WINDOW){
        if(coderInfoL->num_window_groups != coderInfoR->num_window_groups)
            return 0;
        if(coderInfoL->window_group_length[0] != coderInfoR->window_group_length[0])
            return 0;
        if(coderInfoL->window_group_length[1] != coderInfoR->window_group_length[1])
            return 0;
        if(coderInfoL->window_group_length[2] != coderInfoR->window_group_length[2])
            return 0;
        if(coderInfoL->window_group_length[3] != coderInfoR->window_group_length[3])
            return 0;
        if(coderInfoL->window_group_length[4] != coderInfoR->window_group_length[4])
            return 0;
        if(coderInfoL->window_group_length[5] != coderInfoR->window_group_length[5])
            return 0;
        if(coderInfoL->window_group_length[6] != coderInfoR->window_group_length[6])
            return 0;
        /* the following check is redudent
        if(coderInfoL->window_group_length[7] != coderInfoR->window_group_length[7])
            return 0;
        */
    }
    return 1;
}    
         
/* square_polar coupling */
static void couple_lossless(double* spectrum_pair0, double* spectrum_pair1, int length)
{
    int i;

	for(i = 0; i < length; i++)
    {
        int     cmp;
        double  diff;
        
        cmp = fabs(spectrum_pair0[i])>fabs(spectrum_pair1[i]);
        diff = spectrum_pair0[i] - spectrum_pair1[i];       
        
        if(cmp == 1){
            spectrum_pair1[i] 
                = (spectrum_pair0[i]>0)? diff : -diff;
        }else{
            spectrum_pair0[i] = spectrum_pair1[i];
            spectrum_pair1[i]
                = (spectrum_pair1[i]>0)? diff : -diff;
        }
   }
} 

static void decouple_lossless(double* spectrum_pair0, double* spectrum_pair1, int length)
{
    int i;

	for(i = 0; i < length; i++)
    {
        double temp;
        if(spectrum_pair0[i] > 0){
            if(spectrum_pair1[i] > 0){
                spectrum_pair1[i] = spectrum_pair0[i] - spectrum_pair1[i];
            }else{
                temp = spectrum_pair1[i];
                spectrum_pair1[i] = spectrum_pair0[i];
                spectrum_pair0[i] = spectrum_pair0[i] + temp;
            }
        }else{
            if(spectrum_pair1[i] > 0){
                spectrum_pair1[i] = spectrum_pair0[i] + spectrum_pair1[i];
            }else{
                temp = spectrum_pair1[i];
                spectrum_pair1[i] = spectrum_pair0[i];
                spectrum_pair0[i] = spectrum_pair0[i] - temp;
            }
        }
    }
} 

void polar_stereo(CoderInfo*    coderInfo,
                  ChannelInfo*  channelInfo,
                  double*       spectrum[MAX_CHANNELS],
                  int           max_channel,
                  int           allowSquarePolar)
{
    int chn;

	for(chn = 0; chn < max_channel; chn++)
    {
        if(channelInfo[chn].present){
            if(/*channelInfo[chn].cpe &&*/ channelInfo[chn].ch_is_left){
                int rch = channelInfo[chn].paired_ch;
                int nr_of_sfb = coderInfo[chn].nr_of_sfb;

                channelInfo[chn].psInfo.present = 0;
                channelInfo[rch].psInfo.present = 0;
                /* check if common setting exits between the paired channels */                
                channelInfo[chn].common_window = allowSquarePolar && common_grouping(&coderInfo[chn], &coderInfo[rch]);
                channelInfo[rch].common_window = channelInfo[chn].common_window;
                
                if(channelInfo[chn].common_window){
                    int ps_used = 0;
                    int always_used = 1;
                    int sfb;

					PSInfo* psInfoM = &(channelInfo[chn].psInfo);
                    PSInfo* psInfoA = &(channelInfo[rch].psInfo);
                                        
                    for(sfb = 0; sfb < nr_of_sfb; sfb++)
                    {
                        int ps;
                        int start;
                        int length;
						double cor;
                        
                        start = coderInfo[chn].sfb_offset[sfb];
                        length = coderInfo[chn].sfb_offset[sfb+1] - start;
                        /* check if strong correlation exists between two channels 
                           within the current sfb
                        */ 
						cor = cosAngle(spectrum[chn]+start,spectrum[rch]+start,length);
                        ps = ( cor > COSANGLE)? 1 : 0;                        
                        
                        /* register whether  polar stereo is used for the current sfb */
                        psInfoM->ps_used[sfb] = ps;
                        psInfoA->ps_used[sfb] = ps; 
                   
                        if(ps){
                            /* perform channel pair lossless coupling */
                            couple_lossless(spectrum[chn]+start, spectrum[rch]+start, length);
                            ps_used = 1;
                        }else{
                            always_used = 0;
                        }
                                                
                    } 
                                       
                    /* psInfo.present = 0 : polar stereo not used;
                                        1 : polar stereo used for some (not all) sfbs;
                                        2 : polar stereo used for all sfbs;
                                        3 : reserved
                    */
                    channelInfo[chn].psInfo.present = ps_used + always_used;
                    channelInfo[rch].psInfo.present = ps_used + always_used;
                    
                    /* adjust avarage energy */
                    if(ps_used){
                        double avgenrg0;
                        double avgenrg1;

                        avgenrg0 = averEnrg(spectrum[chn], 
                                            coderInfo[chn].sfb_offset[nr_of_sfb]);
                        avgenrg1 = averEnrg(spectrum[rch],
                                            coderInfo[rch].sfb_offset[nr_of_sfb]);
                        
                        if(avgenrg0 != 0 && avgenrg1 !=0){
                            double scale;
                            scale = 2 * avgenrg0/(avgenrg0 + avgenrg1);
                            coderInfo[chn].avgenrg *= 2 - scale;
                            coderInfo[rch].avgenrg *= scale;													 
                        }						
                    } 
                }                    
            }
        }
    }
}

void polar_reconstruct(CoderInfo*   coderInfo,
                       ChannelInfo* channelInfo,
					   int          max_channel)
{
    int chn;
	
	for(chn = 0; chn < max_channel; chn++)
    {
        if(channelInfo[chn].ch_is_left){
            int rch = channelInfo[chn].paired_ch;
            PSInfo* psInfoL = &(channelInfo[chn].psInfo);
            if(psInfoL->present){
                int sfb;

				for(sfb = 0; sfb < coderInfo[chn].nr_of_sfb; sfb++)
                {
                    if(psInfoL->ps_used[sfb]){
                        int start; 
                        int length;
                        double* requantL;
                        double* requantR;
                         
                        start = coderInfo[chn].sfb_offset[sfb];
                        length = coderInfo[chn].sfb_offset[sfb+1] - start;
                        requantL = coderInfo[chn].requantFreq + start;
                        requantR = coderInfo[rch].requantFreq + start;
                        /* perform decoupling for channel pair */
                        decouple_lossless(requantL, requantR, length);
                    }
                }
            }
        }
    }
}


// This function use a simple but efficient criteria to make decision use PQSPC or not. 
// More precise criteria is to compare the final encoding bits as that described in standard, 
// but the later is more complex than the former.
void PQSPSC_stereo(int		*quantFreqL, 
				   int		*quantFreqR, 
				   PSInfo	*psInfoM,
				   PSInfo	*psInfoA,
				   int		*sfb_offset, 
				   int		nr_of_sfb)
{
	int sfb, i;
	int coupled_lineL[FRAME_LEN], coupled_lineR[FRAME_LEN];
	int sum0, sum1;
	int iCoupleCompare;
	int iCmpResult = 0;
	int iAlwaysUsed = 1;

	for (sfb = 0; sfb < nr_of_sfb; sfb++)
	{
		sum0 = 0;
		sum1 = 0;		
		for (i = sfb_offset[sfb]; i < sfb_offset[sfb+1]; i++)
		{
			int     cmp;
			double  diff;

			sum0 += abs(quantFreqL[i]) +abs(quantFreqR[i]);

			cmp = abs(quantFreqL[i])>abs(quantFreqR[i]);
			diff = quantFreqL[i] - quantFreqR[i];			

			if(cmp == 1){
				coupled_lineL[i] = quantFreqL[i];
				coupled_lineR[i] = (quantFreqL[i]>0) ? diff : -diff;
			}else{
				coupled_lineL[i] = quantFreqR[i];
				coupled_lineR[i] = quantFreqR[i] > 0 ? diff : -diff;
			}

			sum1 += abs(coupled_lineL[i]) + abs(coupled_lineR[i]);
		}

		if (sum0 != 0)
		{		
			iCoupleCompare = sum1 < sum0 ? 1 : 0;
			iCmpResult += iCoupleCompare;
			iAlwaysUsed &= iCoupleCompare;
			
			if (iCoupleCompare)
			{
				psInfoM->ps_used[sfb] = 1;
				psInfoA->ps_used[sfb] = 1;
			}else
			{
				psInfoM->ps_used[sfb] = 0;
				psInfoA->ps_used[sfb] = 0;
			}
		}
		else
		{
			psInfoM->ps_used[sfb] = psInfoA->ps_used[sfb] = 0;
		}
	}

	if (iCmpResult <= 5)
	{
		psInfoM->present = psInfoA->present = 0;
	}
	else if (iAlwaysUsed)
	{
		psInfoM->present = psInfoA->present = 2;
		for (i = 0; i < FRAME_LEN; i++)
		{
			quantFreqL[i] = coupled_lineL[i];
			quantFreqR[i] = coupled_lineR[i];
		}
	}
	else
	{
		psInfoM->present = psInfoA->present = 1;
		for (sfb = 0; sfb < nr_of_sfb; sfb++)
		{
			if (psInfoM->ps_used[sfb])
			{
				for (i = sfb_offset[sfb]; i < sfb_offset[sfb+1]; i++)
				{
					quantFreqL[i] = coupled_lineL[i];
					quantFreqR[i] = coupled_lineR[i];
				}
			}
		}
	}

	return;
}