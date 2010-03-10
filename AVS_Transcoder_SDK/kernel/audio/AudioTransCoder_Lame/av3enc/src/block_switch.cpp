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
#include <malloc.h>
#include <stdio.h>
#include "av3enc.h"
#include "block_switch.h"
#include "psychfft.h"
#include "block_switch.h"

static double  sample[PSY_CHANS][3*SUBBLOCK_LEN+BLOCK_LEN_LONG];
static Complex cmpw[129],data[129];
static int nexp;
static VAR_SUBBLOCK    var_subblock;

void check_short(int ch,
				  double *p_time_signal[],
				  enum	WINDOW_TYPE *block_type)
{
	double maxe,maxp=0;
	
	/*caculate in time field*/
    step1(p_time_signal,sample,ch);
	       /* input sample value */
	step2(sample,&var_subblock,&maxe,ch);
	       /* calculate subblock enery and enery diff rate,find max diff rate*/

	/*caculate in freqency field*/
	if(maxe>=E_SWITCH&&(&var_subblock)->last_blocktype[ch]==0)
	{
		step3(sample,&var_subblock,ch);
		  /* subblock FFT*/
		step4(&var_subblock);
		  /* calculate predicted values*/
	    step5(&var_subblock,&maxp);
		  /* calculate each subblock unpredictability */
		
	}

   /* making decision*/
	if((maxp>=P_SWITCH)||(maxe>=(2*E_SWITCH/3)&&(&var_subblock)->last_blocktype[ch]==2))
		*block_type=ONLY_SHORT_WINDOW;
	else
		*block_type=ONLY_LONG_WINDOW;

	(&var_subblock)->last_blocktype[ch]=*block_type;
	
}

void check_short_init(void)
{
	prepsychfft(128,0,&nexp,cmpw);
    subblock_calc_init(sample,&var_subblock);
}
 

void subblock_calc_init( double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN],
		                 VAR_SUBBLOCK *var_subblock)
{ 
  int ch;
  int i;

       /* init var_subblock */
    for(i = 0; i < 2*SUBBLOCK_LEN; ++i)
	var_subblock->hw[i] = (0.5 - 0.5 * cos( M_PI * (double)(i + 0.5)/(double)(SUBBLOCK_LEN) ) );
    for(ch = 0; ch < PSY_CHANS; ++ch)
	{
	  var_subblock->last_subblock_enery[ch]=100;
	  var_subblock->last_blocktype[ch]=0;
	}
    for(ch = 0; ch < PSY_CHANS; ++ch){
		for(i=0;i<BLOCK_LEN_LONG+3*SUBBLOCK_LEN;++i)
			sample[ch][i]=0.0;
		}
}


void step1(double* p_time_signal[],
	       double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN], 
	       int ch) 
{
	int i;

    for(i = 0; i < 3*SUBBLOCK_LEN; ++i){
	sample[ch][i] = sample[ch][i+BLOCK_LEN_LONG];
	}
	for(i=0;i<BLOCK_LEN_LONG;i++)
	{
	sample[ch][i+3*SUBBLOCK_LEN] = (*p_time_signal)[i]/(double)32768.0;
    }
}


void step2(double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN],
					  VAR_SUBBLOCK *var_subblock,double *maxe_r,int ch)
{
	int j,i;
	double temprate,temp;

	for(j = 0; j < NUM_SUBBLOCK; ++j){
		var_subblock->subblock_enery[j]=0.0;

        /* calculate subblock enery */        
        for(i = 0; i < SUBBLOCK_LEN*2; ++i){
        var_subblock->subblock_enery[j]+=sample[ch][2*SUBBLOCK_LEN+SUBBLOCK_LEN * j + i]*sample[ch][2*SUBBLOCK_LEN+SUBBLOCK_LEN * j + i];
	    }
	}
	   /* calculate enery diff rate*/
	temprate=(var_subblock->subblock_enery[0]-var_subblock->last_subblock_enery[ch])/var_subblock->last_subblock_enery[ch];
	if(temprate<0)
		temprate=-temprate;
	temp=temprate;
	for(j=1;j<NUM_SUBBLOCK;++j)
	{
		temprate=(var_subblock->subblock_enery[j]-var_subblock->subblock_enery[j-1])/var_subblock->subblock_enery[j-1];
		if(temprate<0)
			temprate=-temprate;
		if(temprate>temp)
			temp=temprate;
		if(temp>=E_SWITCH)
			break;
	}
	*maxe_r=temp;
	var_subblock->last_subblock_enery[ch]=var_subblock->subblock_enery[NUM_SUBBLOCK-1];
}



void step3(double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN], 
               VAR_SUBBLOCK *var_subblock,
	           int ch)
{
    int w,i,j;
    double *xl;
	
	        
    xl = (double *)malloc( sizeof(double) * 2*SUBBLOCK_LEN);
    for(j = 0; j < NUM_SUBBLOCK+2; ++j){
	    /* windowing */        
        for(i = 0; i < 2*SUBBLOCK_LEN; ++i){
        xl[i] =var_subblock->hw[i]* sample[ch][SUBBLOCK_LEN *j + i];
	    data[i+1]=set(xl[i],0.);
	}
		
	

        /* FFT for subblock */
	psychfft(2*SUBBLOCK_LEN,0,1.0,nexp,cmpw,data);

	    /*cal fft_r and fft_f*/
    for(w = 0; w < SUBBLOCK_LEN; w++){
	    var_subblock->fft_r[j][w] 
	        = sqrt(data[w+1].real*data[w+1].real + data[w+1].imag*data[w+1].imag);

	    if( data[w].real  > 0.0 ){
	        if(data[w].imag  >= 0.0 )
		    var_subblock->fft_f[j][w] = atan( data[w+1].imag  /data[w+1].real  );
		else
		    var_subblock->fft_f[j][w] = atan(data[w+1].imag /data[w+1].real  )+ M_PI * 2.0;
	    } else if(data[w].real < 0.0 ) {
	        var_subblock->fft_f[j][w] = atan(data[w+1].imag  /data[w+1].real  ) + M_PI;
	    } else {
	        if(data[w].imag > 0.0 )
		    var_subblock->fft_f[j][w] = M_PI * 0.5;
		else if(data[w].imag < 0.0 )
		    var_subblock->fft_f[j][w] = M_PI * 1.5;
		else
		    var_subblock->fft_f[j][w] = 0.0; 
	    }
	}
	    
    }

    free(xl);
}


void step4(VAR_SUBBLOCK *var_subblock)
{
	int w,j;
	double r,f,rp,fp;

	for(j = 0; j < NUM_SUBBLOCK; ++j)
	{
		for(w = 0; w < SUBBLOCK_LEN; ++w)
			{
				var_subblock->r_pred[j][w] = 2.0 * var_subblock->fft_r[j+1][w] - var_subblock->fft_r[j][w];
	            var_subblock->f_pred[j][w] = 2.0 * var_subblock->fft_f[j+1][w] - var_subblock->fft_f[j][w];
			}
	}
	for(j = 0; j < NUM_SUBBLOCK; ++j){
		for(w = 0; w < SUBBLOCK_LEN; ++w){
				r = var_subblock->fft_r[j+2][w];
	            f = var_subblock->fft_f[j+2][w];
	            rp =var_subblock->r_pred[j][w];
	            fp =var_subblock->f_pred[j][w];
				if( r + fabs(rp) != 0.0 )
					var_subblock->c[j][w] = sqrt( psy_sqr(r*cos(f) - rp*cos(fp))
				      +psy_sqr(r*sin(f) - rp*sin(fp)) )/ ( r + fabs(rp) ) ;
				else
					var_subblock->c[j][w] = 0.0; 
			}
	}
}

void step5(VAR_SUBBLOCK *var_subblock,double *max_p)
{
	int w,j;
	double tmp_cb,tempp;

	for(j = 0; j < NUM_SUBBLOCK; ++j){
		tmp_cb = 0.0;
		for(w =0; w <SUBBLOCK_LEN ; ++w){
			tmp_cb += psy_sqr(var_subblock->fft_r[j+2][w]) * var_subblock->c[j][w]; 
		}
		var_subblock->cw[j]= tmp_cb;
	}
    tempp=var_subblock->cw[0];
	for(j=1;j<NUM_SUBBLOCK;++j)
	{
		if(var_subblock->cw[j]>tempp)
			tempp=var_subblock->cw[j];
		if(tempp>=P_SWITCH)
			break;
	}
	*max_p=tempp;
}