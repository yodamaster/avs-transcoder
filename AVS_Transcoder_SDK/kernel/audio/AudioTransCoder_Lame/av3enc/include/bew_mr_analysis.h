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
// MultiResolutionAnalysis.h

#ifndef __MULTI_RESOLUTION_ANALYSIS_H__
#define __MULTI_RESOLUTION_ANALYSIS_H__


void InitMultiResolution(ChannelInfo* pChInfo,
						   int nChanNums,
						   int nSamplingRate,
                           int nCBNum,
                           int *anCBWidth
                           );

void MultiResolutionEncode(MRInfo* pMrInfo,
						   enum SIGNAL_TYPE nSignalType,
						   double* pFreqCoef);

void MultiResolutionDecode(MRInfo* pMrInfo,
						   double* pFreqCoef);

#endif //__MULTI_RESOLUTION_ANALYSIS_H__

