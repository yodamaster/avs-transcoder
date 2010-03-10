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

#ifndef _block_switch_h
#define _block_switch_h
#include "psychfft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PSY_CHANS 2
#define BLOCK_LEN_LONG  1024
#define NUM_SUBBLOCK 16
#define SUBBLOCK_LEN  64
#define E_SWITCH 2.5
#define P_SWITCH 20

#define psy_sqr(x) ((x)*(x)) 
 
typedef struct {
  double fft_r[NUM_SUBBLOCK+2][SUBBLOCK_LEN];
  double fft_f[NUM_SUBBLOCK+2][SUBBLOCK_LEN];
  double subblock_enery[NUM_SUBBLOCK];
  double r_pred[NUM_SUBBLOCK][SUBBLOCK_LEN];
  double f_pred[NUM_SUBBLOCK][SUBBLOCK_LEN];
  double c[NUM_SUBBLOCK][SUBBLOCK_LEN];
  double cw[NUM_SUBBLOCK];
  double last_subblock_enery[PSY_CHANS];
  double hw[2*SUBBLOCK_LEN];
  int    last_blocktype[PSY_CHANS];
  } VAR_SUBBLOCK;


void subblock_calc_init( double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN],
	                     VAR_SUBBLOCK *var_subblock);
void check_short_init(void);
void step1(double* p_time_signal[],
	       double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN], 
	       int ch);
void step2(double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN],
		   VAR_SUBBLOCK *var_subblock,
		   double *maxe_r,int ch);

void step3(double sample[][BLOCK_LEN_LONG+3*SUBBLOCK_LEN], 
           VAR_SUBBLOCK *var_subblock,
		   int ch);
void step4(VAR_SUBBLOCK *var_subblock);
void step5(VAR_SUBBLOCK *var_subblock,
		   double *maxp);
void check_short(int ch,
				 double *p_time_signal[],
				 enum	WINDOW_TYPE *block_type);
#endif
