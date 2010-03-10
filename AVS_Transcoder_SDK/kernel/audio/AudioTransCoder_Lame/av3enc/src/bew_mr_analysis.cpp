/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by Beijing E-world and modified 
* for AVS standardization
* 
* Contact Information: Lei Wang (lwang@davworld.net), E-World, Beijing
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
#include "av3enc.h"
#include "bew_mr_analysis.h"
#include <memory.h>
 
#ifndef M_PI /* PI */
#define M_PI            3.14159265358979323846
#endif

#define QUICK_MDCT8

void InitMultiResolution(ChannelInfo* pChInfo,
						 int nChanNums,
						 int nSamplingRate,
                         int nCBNum,
                         int *anCBWidth)
{
	int ch;

	MRInfo* pMrInfo ;

	for (ch = 0; ch < nChanNums; ch++) {
		pMrInfo = &pChInfo[ch].mrInfo;
		pMrInfo->nCBNum = nCBNum;
		pMrInfo->anCBWidth = anCBWidth;
		pMrInfo->anLayoutOfSubband = 1024;
		pMrInfo->nMRSwitchOn = 0;
    }

}

const double COS_T[4] = 
{
	0.99879545620517241,
	0.90398929312344334,
	0.67155895484701833,
	0.33688985339222005
};

const double SIN_T[4] =
{
	0.04906767432741802,
	0.42755509343028208,
	0.74095112535495911,
	0.94154406518302081
};

// cos(a) * 2
const double COS_TN[4] = 
{
	1.99759091241034482,
	1.80797858624688668,
	1.34311790969403666,
	0.6737797067844401
};

// sin(a) * 2
const double SIN_TN[4] = 
{
	0.09813534865483604,
	0.85511018686056416,
	1.48190225070991822,
	1.88308813036604162
};

const double win[16] = 
{
	0.09801714032956,  0.29028467725446,  0.47139673682600,  0.63439328416365,
	0.77301045336274,  0.88192126434835,  0.95694033573221,  0.99518472667220,
	0.99518472667220,  0.95694033573221,  0.88192126434836,  0.77301045336274,
	0.63439328416365,  0.47139673682600,  0.29028467725446,  0.09801714032956
};


static void mdct8(double* pFreqCoef,
		   float* pTimeCoef)
{
	int N = 16;
	double tempr, tempi; 
	double t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i;
	double pr[4], pi[4];
	
	//i = 0
	tempr = pTimeCoef[11] + pTimeCoef[12];
	tempi = pTimeCoef[4] - pTimeCoef[3];
	
	pr [0] = tempr * COS_T[0] +  tempi * SIN_T[0];
	pi [0] = tempi * COS_T[0] -  tempr * SIN_T[0];
	
	//i = 1
	tempr = pTimeCoef[9] + pTimeCoef[14];
	tempi = pTimeCoef[6] - pTimeCoef[1];

	pr[1] = tempr * COS_T[1] + tempi * SIN_T[1];
	pi[1] = tempi * COS_T[1] - tempr * SIN_T[1];

	//i = 2
	tempr = pTimeCoef[7] - pTimeCoef[0];
	tempi = pTimeCoef[8] + pTimeCoef[15];

	pr[2] = tempr * COS_T[2] + tempi * SIN_T[2];
	pi[2] = tempi * COS_T[2] - tempr * SIN_T[2];

	//i = 3
	tempr = pTimeCoef[5] - pTimeCoef[2];
	tempi = pTimeCoef[10] + pTimeCoef[13];

	pr[3] = tempr * COS_T[3] + tempi * SIN_T[3];
	pi[3] = tempi * COS_T[3] - tempr * SIN_T[3];
	
	
	t0r = pr[0] + pr[2];
	t0i = pi[0] + pi[2];
	
	t1r = pr[1] + pr[3];
	t1i = pi[1] + pi[3];
	
	t2r = pr[0] - pr[2];
	t2i = pi[0] - pi[2];
	
	t3r = pr[1] - pr[3];
	t3i = pi[1] - pi[3];
	
	pr[0] = t0r + t1r;
	pi[0] = t0i + t1i;
	
	pr[1] = t2r + t3i;
	pi[1] = t2i - t3r;
	
	pr[2] = t0r - t1r;
	pi[2] = t0i - t1i;
	
	pr[3] = t2r - t3i;
	pi[3] = t2i + t3r;
	
	//=========================================================================
	// i == 0
	tempr = pr[0] * COS_TN[0] + pi[0] * SIN_TN[0];
	tempi = pi[0] * COS_TN[0] - pr[0] * SIN_TN[0];
	
	pFreqCoef[0] = -tempr;
	pFreqCoef[7] = tempi;
	
	// i == 1
	tempr = pr[1] * COS_TN[1] + pi[1] * SIN_TN[1];
	tempi = pi[1] * COS_TN[1] - pr[1] * SIN_TN[1];
	
	pFreqCoef[2] = -tempr;
	pFreqCoef[5] = tempi;
	
	
	// i == 2
	tempr = pr[2] * COS_TN[2] + pi[2] * SIN_TN[2];
	tempi = pi[2] * COS_TN[2] - pr[2] * SIN_TN[2];
	
	pFreqCoef[4] = -tempr;
	pFreqCoef[3] = tempi;
	
	// i == 3
	tempr = pr[3] * COS_TN[3] + pi[3] * SIN_TN[3];
	tempi = pi[3] * COS_TN[3] - pr[3] * SIN_TN[3];
	
	pFreqCoef[6] = -tempr;
	pFreqCoef[1] = tempi;
}


static void sortSpectrum(double *pFreqCoef, int startLine,
                         int nSfb, int *sfbw)
{
  double tmpSpec[FRAME_LEN];
  int sb;

  for (sb=0; sb<MAX_SHORT_WINDOWS; sb++) {
    int k=0;
    while (startLine+sb+MAX_SHORT_WINDOWS*k < FRAME_LEN) {
      tmpSpec[startLine + (MAX_SHORT_WINDOWS-1-sb) *
              (FRAME_LEN-startLine)/MAX_SHORT_WINDOWS + k] =
        pFreqCoef[startLine+sb+MAX_SHORT_WINDOWS*k];
      k++;
    }
  }

  for(sb=0; sb<FRAME_LEN; sb++)
	  pFreqCoef[sb] = tmpSpec[sb];
}

void MultiResolutionEncode(MRInfo* pMrInfo,
                             enum SIGNAL_TYPE nSignalType,
                             double* pFreqCoef)
{
  int i;
  double ddddTemp[16];

  if (nSignalType == TRANSIENT_TYPE)
  {
    pMrInfo->nMRSwitchOn = 1;
  }
  else
  {
    pMrInfo->nMRSwitchOn = 0;
  }
  
  if (pMrInfo->nMRSwitchOn)
  {
    int Ns = 8; 
    int startLine = 0;
    int l;
    float windowedSamples[16], lastSamples[4];


    for (l=0; l<Ns/2; l++)
      lastSamples[l] = 0.0f;

    for (l=startLine; l<=pMrInfo->anLayoutOfSubband-Ns; l+=Ns) {
      if (l==startLine) {
        for (i=0; i<Ns/2; i++)
          windowedSamples[i] = 0.0f;
        for (i=Ns/2; i<Ns; i++)
          windowedSamples[i] = pFreqCoef[l+i-Ns/2];
        for (i=Ns; i<2*Ns; i++)
          windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
        for (i=0; i<Ns/2; i++)
          lastSamples[i] = pFreqCoef[l+i+Ns/2];
      }
      else 
	  {
		  if (l==pMrInfo->anLayoutOfSubband-Ns) 
		  {
          for (i=0; i<Ns/2; i++)
            windowedSamples[i] = lastSamples[i] * win[i];
          for (i=Ns/2; i<Ns; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
          for (i=Ns; i<3*Ns/2; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2];
          for (i=3*Ns/2; i<Ns; i++)
            windowedSamples[i] = 0.0f;
          for (i=0; i<Ns/2; i++)
            lastSamples[i] = 0.0f;
		  }
          else 
		  {
          for (i=0; i<Ns/2; i++)
            windowedSamples[i] = lastSamples[i] * win[i];
          for (i=Ns/2; i<2*Ns; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
          for (i=0; i<Ns/2; i++)
            lastSamples[i] = pFreqCoef[l+i+Ns/2];
		  }
      }
      mdct8(pFreqCoef+l, windowedSamples);
	  memcpy(ddddTemp, pFreqCoef+l, 16 * sizeof(double));
    }

    sortSpectrum(pFreqCoef, startLine, pMrInfo->nCBNum, pMrInfo->anCBWidth);
  }
}

