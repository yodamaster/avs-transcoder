/************************* EAC  Audio Decoder ************************************************
	This software module was originally developed by Beijing E-world technology Co., Ltd. 
This work(and including software and documentation) is provided by the copyright hoder
under the following license:By obtaining, using and/or copying this work, you (the licensee) 
agree that you have read, understood, and will comply with the following terms and conditions.
without permission from Beijing E-world technology Co.,Ltd, any forms of copy,modification 
and distribution are forbidden. The name and trademarks of copyright holders may NOT be used 
in advertising or publicity pertaining to the software without specific, written prior 
permission. Title to copyright in this software and any associated documentation will at all 
times remain with copyright holders, and all right reserved.
***********************************************************************************************/

#include "av3enc.h"
#include "multiresolutionanalysis.h"
 
#ifndef M_PI /* PI */
#define M_PI            3.14159265358979323846
#endif
#define QUICK_MDCT8

void InitMultiResolution(ChannelInfo* pChInfo,
						   int nChanNums,
						   int nSamplingRate,
                           int nCBNum,
                           int *anCBWidth
                           )
{
	int ch;

	pMrInfo->nMRSwitchOn = 0;
	pMrInfo->anLayoutOfSubband = 0;

	for (ch = 0; ch < nChanNums; ch++) {
		MRInfo* pMrInfo = &pChInfo[ch].mrInfo;
		pMrInfo->nCBNum = nCBNum;
		pMrInfo->anCBWidth = anCBWidth;
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
	
	//    fft_qh(pr, pi, 4, 1);
	//=============  4 point FFT ==============================================
	
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
//	pFreqCoef[8] = -tempi;
//	pFreqCoef[15] = tempr;
	
	// i == 1
	tempr = pr[1] * COS_TN[1] + pi[1] * SIN_TN[1];
	tempi = pi[1] * COS_TN[1] - pr[1] * SIN_TN[1];
	
	pFreqCoef[2] = -tempr;
	pFreqCoef[5] = tempi;
//	pFreqCoef[10] = -tempi;
//	pFreqCoef[13] = tempr;
	
	
	// i == 2
	tempr = pr[2] * COS_TN[2] + pi[2] * SIN_TN[2];
	tempi = pi[2] * COS_TN[2] - pr[2] * SIN_TN[2];
	
	pFreqCoef[4] = -tempr;
	pFreqCoef[3] = tempi;
//	pFreqCoef[12] = -tempi;
//	pFreqCoef[11] = tempr;
	
	// i == 3
	tempr = pr[3] * COS_TN[3] + pi[3] * SIN_TN[3];
	tempi = pi[3] * COS_TN[3] - pr[3] * SIN_TN[3];
	
	pFreqCoef[6] = -tempr;
	pFreqCoef[1] = tempi;
//	pFreqCoef[14] = -tempi;
//	pFreqCoef[9] = tempr;	
}


static void sortSpectrum(double *pFreqCoef, int startLine,
                         int nSfb, int *sfbw)
{
  double tmpSpec[FRAME_LEN];
  int sb, i;
  int sfbOffs[MAX_BAND_NUM/MAX_SHORT_WINDOWS+1];
  int sf, offsSrc, offsDest;

  /* sort coefficients from 7654321076543210 ...
                         to 00... 11... 22... 33... 44... 55... 66... 77... */
  for (sb=0; sb<MAX_SHORT_WINDOWS; sb++) {
    int k=0;
    while (startLine+sb+MAX_SHORT_WINDOWS*k < FRAME_LEN) {
      tmpSpec[startLine + (MAX_SHORT_WINDOWS-1-sb) *
              (FRAME_LEN-startLine)/MAX_SHORT_WINDOWS + k] =
        pFreqCoef[startLine+sb+MAX_SHORT_WINDOWS*k];
      k++;
    }
  }

  /* sort scalefactor bands from 00... 11... 22... 33... 44... 55... 66... 77...
                              to 01234567 01234567 ... */
  sfbOffs[0] = 0;
  for (sf=0; sf<nSfb/MAX_SHORT_WINDOWS; sf++)
    sfbOffs[sf+1] = sfbOffs[sf] + sfbw[MAX_SHORT_WINDOWS*sf];

  offsDest = 0;
  for (sf=0; sf<nSfb/MAX_SHORT_WINDOWS; sf++) {
    for (sb=0; sb<MAX_SHORT_WINDOWS; sb++) {
      offsSrc = BLOCK_LEN_SHORT*sb + sfbOffs[sf];
      for (i=0; i<sfbw[MAX_SHORT_WINDOWS*sf]; i++)
        pFreqCoef[offsDest+i] = tmpSpec[offsSrc+i];
      offsDest += sfbw[MAX_SHORT_WINDOWS*sf];
    }
  }


}

void MultiResolutionEncode(MRVQInfo* pMrInfo,
                             enum SIGNAL_TYPE nSignalType,
                             double* pFreqCoef)
{
  int i;
  double temp1, temp2;	
  double ddddTemp[16];

  if (nSignalType == TRANSIENT_TYPE)
  {
    pMrInfo->nMRSwitchOn = 1;
  }
  else
  {
    pMrInfo->nMRSwitchOn = 0;
  }

#ifdef TEST_MR_ALWAYS_ON
    pMrInfo->nMRSwitchOn = 1;
#endif

#ifdef TEST_MR_ALWAYS_OFF
    pMrInfo->nMRSwitchOn = 0;
#endif

  if (pMrInfo->nMRSwitchOn)
  {
    int Ns = 8; /* length of short mdct */
    /*int startLine = (int)floor(pMrInfo->anLayoutOfSubband[1] / Ns) * Ns;*/
    int startLine = 0;
    int l;
    float windowedSamples[16], lastSamples[4];

#ifndef QUICK_MDCT8
	float win[16];
    for (l=0; l<2*Ns; l++) /* sine window */
      win[l] = (float)sin(M_PI/(2*Ns) * (l+0.5f));
#endif

    for (l=0; l<Ns/2; l++)
      lastSamples[l] = 0.0f;

    for (l=startLine; l<=pMrInfo->anLayoutOfSubband-Ns; l+=Ns) {
      if (l==startLine) {
        /* first block is special */
        for (i=0; i<Ns/2; i++)
          windowedSamples[i] = 0.0f;
        for (i=Ns/2; i<Ns; i++)
          windowedSamples[i] = pFreqCoef[l+i-Ns/2];
        for (i=Ns; i<2*Ns; i++)
          windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
        /* save Ns/2 lastSamples for next transform */
        for (i=0; i<Ns/2; i++)
          lastSamples[i] = pFreqCoef[l+i+Ns/2];
      }
      else {
        if (l==pMrInfo->anLayoutOfSubband-Ns) {
          /* last block is special, too */
          for (i=0; i<Ns/2; i++)
            windowedSamples[i] = lastSamples[i] * win[i];
          for (i=Ns/2; i<Ns; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
          for (i=Ns; i<3*Ns/2; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2];
          for (i=3*Ns/2; i<Ns; i++)
            windowedSamples[i] = 0.0f;
          /* no next transform */
          for (i=0; i<Ns/2; i++)
            lastSamples[i] = 0.0f;
        }
        else {
          /* normal block */
          for (i=0; i<Ns/2; i++)
            windowedSamples[i] = lastSamples[i] * win[i];
          for (i=Ns/2; i<2*Ns; i++)
            windowedSamples[i] = pFreqCoef[l+i-Ns/2] * win[i];
          /* save Ns/2 lastSamples for next transform */
          for (i=0; i<Ns/2; i++)
            lastSamples[i] = pFreqCoef[l+i+Ns/2];
        }
      }
      /* transform */
      mdct8(pFreqCoef+l, windowedSamples);
	  memcpy(ddddTemp, pFreqCoef+l, 16 * sizeof(double));
    }

    /* sort spectrum to 8 subblocks */
    sortSpectrum(pFreqCoef, startLine, pMrInfo->nCBNum, pMrInfo->anCBWidth);
  }
}
