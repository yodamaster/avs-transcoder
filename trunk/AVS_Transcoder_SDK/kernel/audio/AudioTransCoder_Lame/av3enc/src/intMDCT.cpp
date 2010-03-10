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

#include "intMDCT.h"
#include "av3enc.h"
#include <stdio.h>
#include <stdlib.h>


int log2int(int x)
{
	int m;
  
	if (x<0)
		x = -x;
	for (m=0; x>1; x>>=1)
		m++;
	return (m);
}


int msbHeadroomINT32(int *x, int len)
{
	int i,max;

	max = 0;
	for (i=0; i<len; i++)
		max |= ABS(x[i]);
	return (30-log2int(max));
}


void shiftLeftINT32(int* x, int len, int shift)
{
	int i;

	if (shift>0)
		for (i=0; i<len; i++)
			x[i] <<= shift;
}


void shiftRightINT32(int* x, int len, int shift)
{
	int i;

	if (shift>0)
		for (i=0; i<len; i++)
			x[i] = (((x[i]>>(shift-1))+1)>>1);
}


int calPreShiftINT32(int *x, int len)
{
	int m;
	
	m = msbHeadroomINT32(x,len);
	if (m>8)
		return (m-6);
	else
		return (2);
}


int multiINT32 (int x, int y)
{
	return ((int)((((__int64)x*y)>>(SHIFT-1))+1)>>1);
}


void hardamard(int* u0, int* u1, int* u2, int* u3)
{
	*u2 = *u2 + *u0;
	*u0 = -*u0;
	*u3 = *u3 - *u1;
	*u1 = -*u1;
	*u0 = *u0 + ((*u2+*u3+1)>>1);
	*u1 = *u1 + ((*u2-*u3+1)>>1);
	*u2 = *u2 - *u1;
	*u3 = *u3 - *u0;
	*u0 = -*u0;
}


void dctII_4p(int *x)
{
	int tmp;

	hardamard(&x[0],&x[1],&x[2],&x[3]);
	
	/* Lifting rotation
	   angle = pi/8
	   | c -s | = | 1  0 | | 1 -s | | 1  0 |
	   | s  c |   | t  1 | | 0  1 | | t  1 | */
	x[0] = x[0] + multiINT32(TAN_PI_16,x[1]);
	x[1] = x[1] - multiINT32(SIN_PI_8,x[0]);
	x[0] = x[0] + multiINT32(TAN_PI_16,x[1]);
	
	tmp = x[0]; x[0] = x[2]; x[2] = x[3];
	x[3] = x[1]; x[1] = tmp;
}

	
void dctII_8p(int *x)
{
	int u[8];
	int tmp,tmp1;

	u[0] = x[0]; u[1] = x[7]; u[2] = x[3]; u[3] = x[4];
	u[4] = x[1]; u[5] = x[6]; u[6] = x[2]; u[7] = x[5];

	hardamard(&u[0],&u[1],&u[2],&u[3]);
	hardamard(&u[4],&u[5],&u[6],&u[7]);
	
	/*  Rotation A1
	    sin(pi/4) * | 1 -1 |
	                | 1  1 |*/
	tmp = multiINT32((u[1]-u[3]),SIN_PI_4);
	tmp1 = multiINT32((u[1]+u[3]),SIN_PI_4);
	u[1] = tmp; u[3] = tmp1;
	
	hardamard(&u[1],&u[7],&u[5],&u[3]);
	
	/* Rotation A2
	   sin(pi/4) * | 1  1 |
	               | 1 -1 |*/
	x[0] = multiINT32((u[2]+u[6]),SIN_PI_4);
	x[4] = multiINT32((u[2]-u[6]),SIN_PI_4);
	
	/* Lifting rotation B
	   angle = pi/8
	   | s  c | = |-t  1 | |-1  s | | 1  0 |
	   |-c  s |   | 1  0 | | 0  1 | | t  1 | */
	u[0] = u[0] + multiINT32(TAN_PI_16,u[4]);
	x[6] = multiINT32(SIN_PI_8,u[0]) - u[4];
	x[2] = u[0] - multiINT32(TAN_PI_16,x[6]);
		
	/* Lifting rotation C
	   angle = 3pi/16
	   | c -s | = | t  1 | |-1 -s | | 0  1 |
	   |-s -c |   | 1  0 | | 0  1 | | 1 -t | */
	u[5] = u[5] - multiINT32(TAN_3PI_32,u[1]);
	x[7] = -u[1] - multiINT32(SIN_3PI_16,u[5]);
	x[1] = u[5] + multiINT32(TAN_3PI_32,x[7]);

	/* Lifting rotation D
	   angle = pi/16
	   | c  s | = |-t  1 | |-1  s | | 0  1 |
	   | s -c |   | 1  0 | | 0  1 | | 1  t |*/
	u[3] = u[3] + multiINT32(TAN_PI_32,u[7]);
	x[3] = multiINT32(SIN_PI_16,u[3]) - u[7];
	x[5] = u[3] - multiINT32(TAN_PI_32,x[3]);
}


void flipud(int *u, int len)
{
	int i,len1,tmp;
	
	len1 = len>>1;
	for (i=0; i<len1; i++)
	{
		tmp = u[i]; u[i] = u[len-1-i]; u[len-1-i] = tmp;
	}
}


void flipudMinus(int *u, int len)
{
	int i,len1,tmp;
	
	len1 = len>>1;
	for (i=0; i<len1; i++)
	{
		tmp = u[i]; u[i] = -u[len-1-i]; u[len-1-i] = -tmp;
	}
}


void swapBuf(int *u, int len)
{
	int i,len1,len2,tmp;
	len1 = len>>1; len2 = len>>2;
	for (i=0; i<len2; i++)
	{
		tmp = u[len2+i]; u[len2+i] = u[len1+i]; u[len1+i] = tmp;
	}
}


void permute(int *u, int logm, int dir)
{		
	int i,j,k,len;
	int *p;
	
	if (dir==1)
		for (i=0; i<logm-1; i++)
		{
			len = 1<<(logm-i); k = 1<<i;
			for (p=u,j=0; j<k; j++)
			{			
				swapBuf(p,len); p += len;
			}
		}
	else
		for (i=0; i<logm-1; i++)
		{
			k = 1<<(logm-i-2); len = 1<<(i+2);
			for (p=u,j=0; j<k; j++)
			{			
				swapBuf(p,len); p += len;
			}
		}
}


void rec_dctII(int *u, int logm);


void rec_dctIV(int *u, int logm)
{
	static int i,j,m,m2;
	static int tmp,tmp1;
	static int s,t,step;
	
	m = 1<<logm; m2 = m>>1;
	
	/* Matrix T - Rotation*/
	step = 1<<(LOG2_N_MDCT-logm);
	for (i=0; i<m2; i++)
	{
		s = sin_tab[(2*i+1)*step]; t = tan_tab[(2*i+1)*step];
		
		/* Lifting rotation
		   angle = (2i+1)pi/4m
		   | c  s | = |-t  1 | |-1  s | | 0  1 |
		   | s -c |   | 1  0 | | 0  1 | | 1  t |*/
		u[i] = u[i] + multiINT32(t,u[m-1-i]);
		u[m-1-i] = multiINT32(s,u[i]) - u[m-1-i];
		u[i] = u[i] - multiINT32(t,u[m-1-i]);
	}
	
	/* Matrix J	- reverse order*/
	flipud(&u[m2],m2);
	
	/* Matrix D -  change sign alternatively*/
	for (i=m2+1; i<m; i+=2)
		u[i] = -u[i];
	
	/* Matrices C2 - DCT-II*/
	rec_dctII(&u[m2],logm-1);
	rec_dctII(u,logm-1);

	m = 1<<logm; m2 = m>>1;

	/* Matrix P - alternative implementation*/
	flipud(&u[m2],m2);
	permute(u,logm,1);
			
	/* Transpose of matrix B */
	for (i=1; i<m-1; i+=2)
	{
		/* Rotation
		   sin(pi/4) * | 1  1 |
		               |-1  1 | */
		tmp = multiINT32((u[i]+u[i+1]),SIN_PI_4);
		tmp1 = multiINT32((u[i+1]-u[i]),SIN_PI_4);
		u[i] = tmp;
		u[i+1] = tmp1;
	}
}
		

void rec_dctII(int *u, int logm)
{
	static int i,j,m,m2;
	static int tmp,tmp1;

	if (logm == 3)
	{
		dctII_8p(u);
		return;
	}

	if (logm == 2)
	{
		dctII_4p(u);
		return;
	}

	m = 1<<logm; m2 = m>>1;

	/* Matrix A2*/
	for (i=0; i<m2; i++)
	{
		/* Rotation
		   sin(pi/4) * | 1  1 |
		               | 1 -1 |*/
		tmp = multiINT32((u[i]+u[m-1-i]),SIN_PI_4);
		tmp1 = multiINT32((u[i]-u[m-1-i]),SIN_PI_4);
		u[i] = tmp;
		u[m-1-i] = tmp1;
	}
	
	/* Matrix J*/
	flipud(&u[m2],m2);
	
	/* Matrix C4 - DCT-IV*/
	rec_dctIV(&u[m2],logm-1);

	/* Matrix C2 - DCT-II*/
	rec_dctII(u,logm-1);
	
	m = 1<<logm; m2 = m>>1;
	
	/* Matrix J*/
	flipud(&u[m2],m2);

	/* Matrix P - alternative implementation*/
	flipud(&u[m2],m2);
	permute(u,logm,1);
}


void MDCT_window(int *signal, int len, int dir)
{
	int i,len1,step;
	
	len1 = (len>>1);
	
	if (len==N_MDCT)
		step = 1;			/*Long block*/
	else
		step = 8;			/*Short block*/
	
	/* Sine window */
	for (i=0; i<len1; i++)
	{
		signal[len1-1-i] += dir*multiINT32(tan_tab[N_MDCT-step*(2*i+1)],signal[len1+i]);
		signal[len1+i] -= dir*multiINT32(sin_tab[N_MDCT-step*(2*i+1)],signal[len1-1-i]);
		signal[len1-1-i] += dir*multiINT32(tan_tab[N_MDCT-step*(2*i+1)],signal[len1+i]);
	}
}


void liftingStepT2(int *signal0, int *signal1, int *liftBuffer, int len)
{
	int i,len1,len2;
	int pre_shift;
	
	len1 = len>>1; len2 = len>>2;
	
	for (i=0; i<len2; i++)				/* D: altenatively change sign*/
	{
		liftBuffer[2*i] = signal0[2*i];
		liftBuffer[2*i+1] = -signal0[2*i+1];
	}
	
	pre_shift = calPreShiftINT32(liftBuffer,len1);
	shiftLeftINT32(liftBuffer,len1,pre_shift);

	rec_dctIV(liftBuffer,log2int(len1));			/* DCT-IV */

	for (i=0; i<len1; i++)
		liftBuffer[i] = multiINT32(SQRT2,liftBuffer[i]);
	
	shiftRightINT32(liftBuffer,len1,pre_shift);
}


void liftingStepT1(int *signal0, int *signal1, int *liftBuffer, int len)
{
	int i,len1,len2;
	int pre_shift;
	
	len1 = len>>1; len2 = len>>2;
	
	for (i=0; i<len1; i++)
		liftBuffer[i] = signal1[i];

	pre_shift = calPreShiftINT32(liftBuffer,len1);
	shiftLeftINT32(liftBuffer,len1,pre_shift);

	rec_dctIV(liftBuffer,log2int(len1));			/* DCT-IV */ 

	for (i=0; i<len1; i++)
		liftBuffer[i] = multiINT32(liftBuffer[i],SIN_PI_4);

	shiftRightINT32(liftBuffer,len1,pre_shift);
}


void liftingStepS(int *signal0, int *signal1, int *liftBuffer, int len)
{
	int i,len1,len2;
	int pre_shift;
	int *v1,*v2;
	int step;
	
	len1 = len>>1; len2 = len>>2;
	
	if (len==N_MDCT)
		step = 1;			/*Long block*/
	else
		step = 8;			/*Short block*/
		
	v1 = &liftBuffer[len1]; v2 = &liftBuffer[len];

	for (i=0; i<len1; i++)
		liftBuffer[i] = signal0[i];

	pre_shift = calPreShiftINT32(liftBuffer,len1);
	shiftLeftINT32(liftBuffer,len1,pre_shift);

	for (i=0; i<len1; i++)
	{
		v1[i] = - multiINT32(tan_tab[N_MDCT-step*(2*i+1)],liftBuffer[len1-1-i]);
		v2[i] = liftBuffer[i];
	}

	rec_dctIV(v2,log2int(len1));			/*DCT-IV*/

	for (i=0; i<len1; i++)
	{
		liftBuffer[i] = v2[i];
		v2[i] = multiINT32(SQRT2,v2[i]);
	}

	for (i=0; i<len2; i++)					/*D: altenatively change sign*/
		liftBuffer[2*i+1] = -liftBuffer[2*i+1];

	rec_dctIV(liftBuffer,log2int(len1));	/* DCT-IV*/ 

	for (i=0; i<len1; i++)
	{
		v2[i] = -v2[i] - liftBuffer[i];
		liftBuffer[i] = v1[i] + v2[i];
	}

	shiftRightINT32(liftBuffer,len1,pre_shift);
}


void intDCTIV(int *signal0, int *signal1, int len)
{
	int liftBuffer[N_MDCT+N_MDCT_2];
	int i,len1,len2;
	int step;
	
	len1 = len>>1; len2 = len>>2;
	
	if (len==N_MDCT)
		step = 1;			/*Long block*/
	else
		step = 8;			/*Short block*/
	
	/* Forward integer DCT-IV*/
	/* C4 = R1*R2*S*T1*T2*Peo*/

	/* Step(1) Even-odd permutation matrix Peo*/ 
	permute(signal0,log2int(len),-1);
	
	/* Step(2) Lifting matrix T2*/
	liftingStepT2(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal1[i] += liftBuffer[i] + signal0[i];

	/* Step(3) Lifting matrix T1*/
	for (i=0; i<len2; i++)				/* -D: altenatively change sign */
		signal0[2*i] = -signal0[2*i];
	liftingStepT1(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal0[i] += liftBuffer[i];

	/* Step(4) Lifting matrix S*/
	liftingStepS(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal1[i] += liftBuffer[i];

	/* Step(5) Lifing matrix R2*/
	for (i=0; i<len1; i++)
		signal0[i] += multiINT32(sin_tab[step*(2*i+1)],signal1[len1-1-i]);

	/* Step(6) Lifting matrix R1*/
	for (i=0; i<len1; i++)
		signal1[i] -= multiINT32(tan_tab[N_MDCT-step*(2*i+1)],signal0[len1-1-i]);
	/* End of forward integer DCT-IV*/
}


void intIDCTIV(int *signal0, int *signal1, int len)
{
	int liftBuffer[N_MDCT+N_MDCT_2];
	int i,len1,len2;
	int step;
	
	len1 = len>>1; len2 = len>>2;
	
	if (len==N_MDCT)
		step = 1;			/*Long block*/
	else
		step = 8;			/*Short block*/
	
	/* Integer inverse DCT-IV */
	/* Inv(C4) = Inv(R1*R2*S*T1*T2*Peo) */
	
	/* Inverse step(6) Lifting matrix R1 */
	for (i=0; i<len1; i++)
		signal1[i] += multiINT32(tan_tab[N_MDCT-step*(2*i+1)],signal0[len1-1-i]);
		
	/* Inverse step(5) Lifing matrix R2 */
	for (i=0; i<len1; i++)
		signal0[i] -= multiINT32(sin_tab[step*(2*i+1)],signal1[len1-1-i]);

	/* Inverse step(4) Lifting matrix S */
	liftingStepS(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal1[i] -= liftBuffer[i];	

	/* Inverse step(3) Lifting matrix T1 */
	liftingStepT1(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal0[i] -= liftBuffer[i];
	for (i=0; i<len2; i++)				/* -D: altenatively change sign */
		signal0[2*i] = -signal0[2*i];
	
	/* Inverse step(2) Lifting matrix T2 */
	liftingStepT2(signal0,signal1,liftBuffer,len);
	for (i=0; i<len1; i++)
		signal1[i] -= liftBuffer[i] + signal0[i];

	/* Inverse step(1) Even-odd permutation matrix Peo */
	permute(signal0,log2int(len),1);
	/* End of Integer inverse DCT-IV */
}


void intMDCT(double* p_in, double* p_overlap, double* p_out, int blockType)
{
	int i,j;
	int transFac = 8;
	int startOnesLength = N_MDCT_2-N_MDCT_SHORT_2;  /* number of ones in start window */
	int signal[N_MDCT+N_MDCT_2];

	/* Copy buffers */
	for (i=0; i<N_MDCT_2; i++)
		signal[i] = p_overlap[i];

	for (i=0; i<N_MDCT; i++)
    	signal[N_MDCT_2+i] = p_in[i];
  
	/* Window switching */
	switch(blockType)
	{
		case ONLY_LONG_WINDOW:
		case SHORT_LONG_WINDOW:
    		MDCT_window(&signal[N_MDCT_2],N_MDCT,1);
    		break;
  		case LONG_SHORT_WINDOW:
    		MDCT_window(&signal[N_MDCT_2+startOnesLength],N_MDCT_SHORT,1);
    		break;
  		case ONLY_SHORT_WINDOW:
    		for (j=0; j<transFac; j++)
 				MDCT_window(&signal[N_MDCT_SHORT_2+j*N_MDCT_SHORT],N_MDCT_SHORT,1);
    		break;
  		default:
  			printf("blockType not implemented!\n");
	    	break;
  }
    
	/* DCT-IV */
	switch(blockType)
	{
		case ONLY_LONG_WINDOW:
		case LONG_SHORT_WINDOW:
		case SHORT_LONG_WINDOW:
			flipudMinus(signal,N_MDCT);				/* reverse order and change sign */
    		intDCTIV(signal,&signal[N_MDCT_2],N_MDCT);	 /* Integer DCT-IV */
     		break;
		case ONLY_SHORT_WINDOW:
    		for (j=0; j<transFac; j++)
    		{
    			flipudMinus(&signal[j*N_MDCT_SHORT],N_MDCT_SHORT);
      			intDCTIV(&signal[j*N_MDCT_SHORT],&signal[j*N_MDCT_SHORT+N_MDCT_SHORT_2],N_MDCT_SHORT);
			}
    		break;  
		default:
			printf("blockType not implemented!\n");
			break;
	}

	/* copy buffers */
	for (i=0; i<N_MDCT; i++)
		p_out[i] = signal[i];
  
	for (i=0; i<N_MDCT_2; i++)
		p_overlap[i] = signal[N_MDCT+i];
}


void intIMDCT(double* p_in, double* p_overlap, double* p_out, int blockType)
{
	int i,j;
	int transFac = 8;
	int startOnesLength = N_MDCT_2-N_MDCT_SHORT_2;  /* number of ones in start window*/
	int signal[N_MDCT+N_MDCT_2];

	/* Copy buffers */
	for (i=0; i<N_MDCT_2; i++)
		signal[i] = p_overlap[i];

	for (i=0; i<N_MDCT; i++)
    	signal[N_MDCT_2+i] = p_in[i];
  
  	/* InvDCT-IV */
	switch(blockType)
	{
		case ONLY_LONG_WINDOW:
		case LONG_SHORT_WINDOW:
		case SHORT_LONG_WINDOW:
			intIDCTIV(&signal[N_MDCT_2],&signal[N_MDCT],N_MDCT);	/* Integer inverse DCT-IV */
			flipudMinus(&signal[N_MDCT_2],N_MDCT);			/* reverse order and change sign */
     		break;
		case ONLY_SHORT_WINDOW:
    		for (j=0; j<transFac; j++)
    		{
    			intIDCTIV(&signal[j*N_MDCT_SHORT+N_MDCT_2],&signal[j*N_MDCT_SHORT+N_MDCT_2+N_MDCT_SHORT_2],N_MDCT_SHORT);
    			flipudMinus(&signal[j*N_MDCT_SHORT+N_MDCT_2],N_MDCT_SHORT);
			}
    		break;  
		default:
			printf("blockType not implemented!\n");
			break;
	}
	
	/* Window switching */
	switch(blockType)
	{
		case ONLY_LONG_WINDOW:
		case SHORT_LONG_WINDOW:
    		MDCT_window(signal,N_MDCT,-1);
    		break;
  		case LONG_SHORT_WINDOW:
    		MDCT_window(&signal[startOnesLength],N_MDCT_SHORT,-1);
    		break;
  		case ONLY_SHORT_WINDOW:
    		for (j=0; j<transFac; j++)
 				MDCT_window(&signal[j*N_MDCT_SHORT+startOnesLength],N_MDCT_SHORT,-1);
    		break;
  		default:
  			printf("blockType not implemented!\n");
	    	break;
	}
    
	/* copy buffers */
	for (i=0; i<N_MDCT; i++)
		p_out[i] = signal[i];
  
	for (i=0; i<N_MDCT_2; i++)
		p_overlap[i] = signal[N_MDCT+i];
}