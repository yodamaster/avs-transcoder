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
// flpvq.h
#ifndef __FLP_H__
#define __FLP_H__

#include "av3enc.h"

#define GRID_POINTS				120
#define FIRST_FIRSTSUBVECTOR	1024
#define FIRST_SECONDSUBVECTOR   1024
#define SECOND_FIRSTSUBVECTOR   512
#define SECOND_SECONDSUBVECTOR  256
#define NUMBER_SUBVECTOR        6


void InitFLPVQStructure(int nUseFlp,
					  ChannelInfo* pChInfo,
					  int nSampleIndex,
					  int nChannelNum);

void FlpVQEncode(FLPVQInfo* pFlpInfo,
			   int nNumOfSfb,
			   int anSfbOffset[],
			   double adPWaveletCoef[]);

double LevinsonDurbin(int nPredOrder,   
                      int nPredLength,       
                      double* pPredCoef,        
                      double* pReflectCoef);

void AutoCorrelation(int nPredOrder, 
					 int nPredLength,
					 double* pPredCoef,
					 double* pAutoCorCoef);

void Levinson(double* pAutoCorCoef, 
			  int nPredOrder,
			  double* pReflectCoef,
			  double* pErrorPower);

void CalculatePredCoef(int nCuttedOrder,
					   double* pReflectCoef,
					   double* pPredCoef);

void FlpVQFilter(int nLength,
			   double* pPWaveletCoef,
			   FLPVQInfo* pFlpInfo);

void FlpVQDecode(FLPVQInfo* pFlpInfo,
               int numberOfBands,         
               int anSfbOffset[],
			   double adPWaveletCoef[]);

void IFlpVQFilter(int nLength, 
				double* pPWaveletCoef, 
				FLPVQInfo* pFlpInfo);

void PredCoefToLineSpectrumFreq(int nPredOrder, 
								double* pPredCoef,
								double* pLineSpectrumFreq);

double BET_Chebps(double x, double f[], int n);

void VectorQuantLineSpectrumFreq(FLPVQInfo* pFlpInfo);

void LineSpectrumPairToPredCoef(double* pLineSpectrumPair,
								double* pPredCoef);

void GetLspPol(double *lsp, double *f);

#endif //__FLP_H__

