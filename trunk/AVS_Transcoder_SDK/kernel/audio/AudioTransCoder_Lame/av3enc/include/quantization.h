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

#ifndef QUANT_H
#define QUANT_H

#ifdef __cplusplus
extern "C" {
#endif 

#include "av3enc.h"

#define IXMAX_VAL 8191
#define PRECALC_SIZE (IXMAX_VAL+2)
#define LARGE_BITS 100000
#define SF_OFFSET 100
#define MAGIC_FLOAT (65536*(128))
#define MAGIC_INT 0x4b000000
#define LOG2_8_3 3.847186775703902419626465816005

typedef union {
    float f;
    int i;
} fi_union;

void quantizeInit(CoderInfo *coderInfo, 
                  unsigned int numChannels
                 );

void quantizeEnd(CoderInfo *coderInfo, 
                 unsigned int numChannels
		         );

int quantize(CoderInfo *coderInfo,
             ChannelInfo *channelInfo,
             double *xr,
			 double quantlity
            );

int SortForGrouping(CoderInfo* coderInfo,
		            ChannelInfo *channelInfo,
		            int *sfb_width_table,
		            double *xr);
  
int ReSortForGrouping(CoderInfo* coderInfo,
                      ChannelInfo *channelInfo,
                      int *sfb_width_table,
                      int *xr);
void resortSpectrum(int *pFreqCoef, int startLine,
                           int nSfb, int *sfbw);

void CalcAvgEnrg(CoderInfo *coderInfo,
		         const double *xr);

static int FixNoise(CoderInfo *coderInfo, double *xr_pow, int *xi, double *xmin);
static void CalcAllowedDist(CoderInfo *coderInfo, double *xr, double *xmin, int quality);

#ifdef __cplusplus
}
#endif 

#endif 
