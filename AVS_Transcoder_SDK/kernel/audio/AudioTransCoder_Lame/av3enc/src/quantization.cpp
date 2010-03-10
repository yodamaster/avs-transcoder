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
#include <stdlib.h>

#include "av3enc.h"
#include "quantization.h"
#include "psych.h"
#include <memory.h>

static double adj43[PRECALC_SIZE];

void  quantizeInit(CoderInfo *coderInfo, unsigned int numChannels)
{
    unsigned int channel; 
    unsigned int i;
    double pi0;
    double pi1;

    adj43[0] = 0.0;
    pi1 = 0;         
    
	for (i = 1; i < PRECALC_SIZE; i++)
    {
        pi0 = pi1;
        pi1 = pow(i, 4.0/3.0);
        adj43[i] =  i - .5 - pow(.5 * (pi0 + pi1), .75);
	}
    
    for (channel = 0; channel < numChannels; channel++) 
		coderInfo[channel].quantFreq = (int*)malloc(BLOCK_LEN_LONG*sizeof(int));

}

void quantizeEnd(CoderInfo *coderInfo, unsigned int numChannels)
{
    unsigned int channel;

   for (channel = 0; channel < numChannels; channel++) 
		    if (coderInfo[channel].quantFreq) 
		        free(coderInfo[channel].quantFreq);    
}

int quantize(CoderInfo *coderInfo, ChannelInfo *channelInfo, double *xr,double quality)
{
    int sb, i, do_q = 0;
    int bits = 0, sign;
    double xr_pow[FRAME_LEN];
    double xmin[MAX_SCFAC_BANDS];
    int xi[FRAME_LEN];

    /* Use local copy's */
    int *scale_factor = coderInfo->scale_factor;
	
	memset(xi, 0, FRAME_LEN*sizeof(int));

    /* Set all scalefactors to 0 */
    coderInfo->global_gain = 0;
    for (sb = 0; sb < coderInfo->nr_of_sfb; sb++)
        scale_factor[sb] = 0;

    /* Compute xr_pow */
    for (i = 0; i < FRAME_LEN; i++) {
        double temp = fabs(xr[i]);
        xr_pow[i] = sqrt(temp * sqrt(temp));
        do_q += (temp > 1E-20);
    }

    if (do_q) {
        CalcAllowedDist(coderInfo, xr, xmin,quality);
	      
		coderInfo->global_gain = 0;
		FixNoise(coderInfo, xr_pow, xi, xmin);
		
		
		  for ( i = 0; i < FRAME_LEN; i++ )  {
            sign = (xr[i] < 0) ? -1 : 1;
            xi[i] *= sign;
        }
    } else {
        coderInfo->global_gain = 0;
        memset(xi, 0, FRAME_LEN*sizeof(int));
    }

	  for (i = 0; i < FRAME_LEN; i++) {
        coderInfo->quantFreq[i] = xi[i];
	  }

	  for (i = 0; i < coderInfo->nr_of_sfb; i++) {
		    scale_factor[i] = coderInfo->global_gain - scale_factor[i] + SF_OFFSET;
    }
    return bits;
}

static void QuantizeBand(const double *xp, int *pi, double istep, int offset, int end)
{
    int j;
    fi_union *fi;

    fi = (fi_union *)pi;
    for (j = offset; j < end; j++)
    {
        double x0 = istep * xp[j];

        x0 += MAGIC_FLOAT; fi[j].f = (float)x0;
        fi[j].f = x0 + (adj43 - MAGIC_INT)[fi[j].i];
        fi[j].i -= MAGIC_INT;
    }
}

static void CalcAllowedDist(CoderInfo *coderInfo, double *xr, double *xmin, int quality)
{
    int sfb, start, end, l;
    const double globalthr = 196.0 / (double)quality;
    int last = coderInfo->lastx;
    int lastsb = 0;
    int *cb_offset = coderInfo->sfb_offset;
    int num_cb = coderInfo->nr_of_sfb;
    double avgenrg = coderInfo->avgenrg;

    for (sfb = 0; sfb < num_cb; sfb++)
    if (last > cb_offset[sfb])
        lastsb = sfb;
 
    for (sfb = 0; sfb < num_cb; sfb++)
    {
        double thr, tmp;
        double enrg = 0.0;

        start = cb_offset[sfb];
        end = cb_offset[sfb + 1];

        if (sfb > lastsb){
            xmin[sfb] = 0;
            continue;
        }

        if (coderInfo->signal_type != TRANSIENT_TYPE) {
            double enmax = -1.0;
            double lmax;

            lmax = start;
            for (l = start; l < end; l++)
            {
	              if (enmax < (xr[l] * xr[l])){
	                  enmax = xr[l] * xr[l];
	                  lmax = l;
	              }
            }

            start = lmax - 2;
            end = lmax + 3;
      
            if (start < 0)
	              start = 0;
            if (end > last)
	              end = last;
        }

        for (l = start; l < end; l++)
            enrg += xr[l]*xr[l];

        thr = enrg/((double)(end-start)*avgenrg);
        thr = pow(thr, 0.1*(lastsb-sfb)/lastsb + 0.3);

        tmp = 1.0 - ((double)start / (double)last);
        tmp = tmp * tmp * tmp + 0.075;

        thr = 1.0 / (1.4*thr + tmp);

        xmin[sfb] = ((coderInfo->signal_type == TRANSIENT_TYPE) ? 0.65 : 1.12)
                    * globalthr * thr;
    }
}

static int FixNoise(CoderInfo* coderInfo, double* xr_pow, int* xi, double* xmin)
{
    int   sb;
    int   start;
    int   end;
        
    end = coderInfo->sfb_offset[0];
    for (sb = 0; sb < coderInfo->nr_of_sfb; sb++)
    {
        double  maxx;
        double  Xsum0;
		    double  Xsum1;
		    double  fac;
        double  dfac;
		    int     sfac;
        int     delta0;
		    int     delta1;
		    int     k;
		    int     i;

		start = end;
        end = coderInfo->sfb_offset[sb+1];

        if (!xmin[sb]){
	          for (i = start; i < end; i++)
	              xi[i] = 0;
	          coderInfo->scale_factor[sb] = 10;
	          continue;
	      }    
      
        maxx = 0.0;
        for (i = start; i < end; i++)
        {
	          if (xr_pow[i] > maxx)
	              maxx = xr_pow[i];
        }

        if (maxx < 10.0){
	          for (i = start; i < end; i++)
	              xi[i] = 0;
	          coderInfo->scale_factor[sb] = 10;
	          continue;
        }
        
		    Xsum0 = 0.0;
		    Xsum1 = 0.0;
		    for(i = start; i < end; i++)
			      Xsum0 +=xr_pow[i]; 
		    		    
		    fac = -1.5*log(xmin[sb]) + 2*log(double(end - start)) - 2*log(Xsum0);
		    fac *= LOG2_8_3;
		    sfac = (int)(fac + .5);
		
		    k = 0;
		    delta0 = 0;
		    delta1 = 0;
		    dfac = 1;
		    while(k++ < 4){
			      coderInfo->scale_factor[sb] = sfac;
			      if(1 == k){	
				        fac = pow(2, .1875*sfac);
				        dfac = fac;
			      }else{
				        dfac = pow(2, .1875*delta0);
				        fac *= dfac;
			      }			
			      for(i = start; i < end; i++)
				        xr_pow[i] *= dfac;
			      QuantizeBand(xr_pow, xi, 1, start, end);
			      for(i = start; i < end; i++)
				        Xsum1 += xi[i];
			      Xsum1 = (Xsum1+.01)/fac;
			      delta1 = (int)(2*log(Xsum0/Xsum1)*LOG2_8_3+.5);
			      delta1 = (delta1>4)?  4 : delta1;
			      delta1 = (delta1<-4)? -4 : delta1;
			      if(!delta1)
				        break;
			      if(delta1 * delta0 < 0 && -delta1 >= abs(delta0))
				        break;
			      sfac += delta1;
			      Xsum0 = Xsum1;
			      Xsum1 = 0;
			      delta0 = delta1;
			      continue;
	      }
    } 
    return 0;
}

int SortForGrouping(CoderInfo* coderInfo,
                    ChannelInfo *channelInfo,
                    int *sfb_width_table,
                    double *xr)
{
    int i,j,ii;
    int index = 0;
    double xr_tmp[1024];
    int group_offset=0;
    int k=0;
    int windowOffset = 0;


    /* set up local variables for used quantInfo elements */
    int* sfb_offset = coderInfo->sfb_offset;
    int* nr_of_sfb = &(coderInfo->nr_of_sfb);
    int* window_group_length;
    int num_window_groups;
    *nr_of_sfb = coderInfo->max_sfb;              /* Init to max_sfb */
    window_group_length = coderInfo->window_group_length;
    num_window_groups = coderInfo->num_window_groups;

    /* calc org sfb_offset just for shortblock */
    sfb_offset[k]=0;
    for (k=1 ; k <*nr_of_sfb+1; k++) {
        sfb_offset[k] = sfb_offset[k-1] + sfb_width_table[k-1];
    }

    /* sort the input spectral coefficients */
    index = 0;
    group_offset=0;
	if (xr != NULL)
	{
		for (i=0; i< num_window_groups; i++) {
			for (k=0; k<*nr_of_sfb; k++) {
				for (j=0; j < window_group_length[i]; j++) {
					for (ii=0;ii< sfb_width_table[k];ii++)
						xr_tmp[index++] = xr[ii+ sfb_offset[k] + 128*j +group_offset];
				}
			}
			group_offset +=  128*window_group_length[i];
		}

		for (k=0; k<1024; k++){
			xr[k] = xr_tmp[k];
		}
	}


    /* now calc the new sfb_offset table for the whole p_spectrum vector*/
    index = 0;
    sfb_offset[index++] = 0;
    windowOffset = 0;
    for (i=0; i < num_window_groups; i++) {
        for (k=0 ; k <*nr_of_sfb; k++) {
            sfb_offset[index] = sfb_offset[index-1] + sfb_width_table[k]*window_group_length[i] ;
            index++;
        }
        windowOffset += window_group_length[i];
    }

    *nr_of_sfb = *nr_of_sfb * num_window_groups;  /* Number interleaved bands. */

    return 0;
}

int ReSortForGrouping(CoderInfo* coderInfo,
					  ChannelInfo *channelInfo,
					  int *sfb_width_table,
					  int *xr)
{
    int i,j,ii;
    int index = 0;
    int xr_tmp[1024];
	int group_offset=0;
    int k=0;
    int windowOffset = 0;
	
    /* set up local variables for used quantInfo elements */
    int* sfb_offset = coderInfo->sfb_offset;
    int* nr_of_sfb = &(coderInfo->nr_of_sfb);
    int* window_group_length;
    int num_window_groups;
    *nr_of_sfb = coderInfo->max_sfb;              /* Init to max_sfb */
    window_group_length = coderInfo->window_group_length;
    num_window_groups = coderInfo->num_window_groups;
	
	/* 2005-11-26 xuhengyu */
//	*nr_of_sfb = *nr_of_sfb/num_window_groups;  /* Number interleaved bands. */
	
	index = 0;
    sfb_offset[index++] = 0;
    for (k=0 ; k <*nr_of_sfb; k++) {
		sfb_offset[index] = sfb_offset[index-1] + sfb_width_table[k];
		index++;
    }
	
    /* sort the input spectral coefficients */
    index = 0;
    group_offset=0;
    for (i=0; i< num_window_groups; i++) {
        for (k=0; k<*nr_of_sfb; k++) {
            for (j=0; j < window_group_length[i]; j++) {
                for (ii=0;ii< sfb_width_table[k];ii++)
                    xr_tmp[ii+ sfb_offset[k] + 128*j +group_offset] = xr[index++];
            }
        }
        group_offset +=  128*window_group_length[i];
    }
	
	/* resort for CBC */
	ii = 0;
    for(index = 0; index<num_window_groups; index++) {
		for (i = 0; i < 128; i+=4) {
			group_offset = (128 * ii) + (i * window_group_length[index]);
			for (k = 0; k < 4; k++) 
				for(j = 0; j < window_group_length[index]; j++)  
					xr[group_offset+4*j+k] = xr_tmp[128*(j+ii)+i+k];
		}
		ii += window_group_length[index];
	}
	
	return 0;
}


void resortSpectrum(int *pFreqCoef, int startLine,
                           int nSfb, int *sfbw)
{
	double tmpSpec[FRAME_LEN];
	int sb, i;
	int sfbOffs[MAX_BANDS];
	int sf, offsSrc, offsDest;
	
	/* sort scalefactor bands from 01234567 01234567 ... 
	to 00... 11... 22... 33... 44... 55... 66... 77...*/
	sfbOffs[0] = 0;
	for (sf=0; sf<nSfb; sf++)
		sfbOffs[sf+1] = sfbOffs[sf] + sfbw[sf];
	
	offsSrc = 0;
	for (sf=0; sf<nSfb; sf++) {
		for (sb=0; sb<MAX_SHORT_WINDOWS; sb++) {
			offsDest = BLOCK_LEN_SHORT*sb + sfbOffs[sf];
			for (i=0; i<sfbw[sf]; i++)
				tmpSpec[offsDest+i] = pFreqCoef[offsSrc+i];
			offsSrc += sfbw[sf];
		}
	}
	
	/* sort coefficients from 00... 11... 22... 33... 44... 55... 66... 77... 
	to 7654321076543210 ... */
	for (sb=0; sb<MAX_SHORT_WINDOWS; sb++) {
		int k=0;
		while (startLine+sb+MAX_SHORT_WINDOWS*k < FRAME_LEN) {
			pFreqCoef[startLine+sb+MAX_SHORT_WINDOWS*k] =
				tmpSpec[startLine + (MAX_SHORT_WINDOWS-1-sb) *
                (FRAME_LEN-startLine)/MAX_SHORT_WINDOWS + k];
			k++;
		}
	}
}


void CalcAvgEnrg(CoderInfo *coderInfo, const double *xr)
{
    int end, l;
    int last = 0;
    double totenrg = 0.0;

    end = coderInfo->sfb_offset[coderInfo->nr_of_sfb];
    for (l = 0; l < end; l++)
    {
        if (xr[l]){
            last = l;
            totenrg += xr[l] * xr[l];
        }
    }
    last++;

    coderInfo->lastx = last;
    coderInfo->avgenrg = totenrg / last;
}
