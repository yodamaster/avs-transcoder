/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
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
/*
 * $Id: filtbank.c,v 1.9 2001/09/04 18:39:35 menno Exp $
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mdct.h"


#ifndef M_PI
#define  M_PI 3.14159265358979323846
#endif

#define		TWOPI			2*M_PI
#define		BLOCK_LEN_LONG	1024


//extern double max_e, min_e;

static void exac_fft(double *pr,double *pi,int n,int l)
{ 
    int it,m,is,i,j,nv,l0,nt,k;
    double p,q,s,vr,vi,poddr,poddi;
	double *fr,*fi;
	
    k = 0;
	
	nt = n;
    for (;;)
	{
		m = nt >> 1;  
		if ((m << 1) == nt) 
		{
			k++;
			nt = m;
		} 
		else
		{
			break;
		}
	}
	
	fr = (double *)malloc( sizeof(double) * n );
	fi = (double *)malloc( sizeof(double) * n );
	
    for (it=0; it<=n-1; it++)
	{ 
		m=it;
		is=0;
		
		for (i=0; i<=k-1; i++)
		{
			j=m/2;
			is=2*is+(m-2*j); 
			m=j;
		}
		fr[it]=pr[is];
		fi[it]=pi[is];
	}

    pr[0]=1.0;
	pi[0]=0.0;
    p=6.2831853071795865/(1.0*n);
    pr[1]=cos(p); 
	pi[1]=-sin(p);

    if (l!=1) pi[1]=-pi[1];

    for (i=2; i<=n-1; i++)
	{
		p=pr[i-1]*pr[1];
		q=pi[i-1]*pi[1];
		s=(pr[i-1]+pi[i-1])*(pr[1]+pi[1]);
		pr[i]=p-q; pi[i]=s-p-q;
	}

    for (it=0; it<=n-2; it=it+2)
	{ 
		vr=fr[it];
		vi=fi[it];
		fr[it]=vr+fr[it+1];
		fi[it]=vi+fi[it+1];
		fr[it+1]=vr-fr[it+1]; 
		fi[it+1]=vi-fi[it+1];
	}

    m=n/2; 
	nv=2;
    for (l0=k-2; l0>=0; l0--)
	{
		m=m/2; 
		nv=2*nv;
		for (it=0; it<=(m-1)*nv; it=it+nv)
		{
			for (j=0; j<=(nv/2)-1; j++)
			{ 
				p=pr[m*j]*fr[it+j+nv/2];
				q=pi[m*j]*fi[it+j+nv/2];
				s=pr[m*j]+pi[m*j];
				s=s*(fr[it+j+nv/2]+fi[it+j+nv/2]);
				poddr=p-q; poddi=s-p-q;
				fr[it+j+nv/2]=fr[it+j]-poddr;
				fi[it+j+nv/2]=fi[it+j]-poddi;
				fr[it+j]=fr[it+j]+poddr;
				fi[it+j]=fi[it+j]+poddi;


//            if(q>=max_e) max_e=q;
//    	       if(q<min_e)  min_e=q;

			}
		}
	}

	for (i=0; i<=n-1; i++)
	{
		pr[i]=fr[i];
		pi[i]=fi[i];
	}
	
	free(fr);
	free(fi);
	
	return;
}


static void MDCT(double *data, int N)
{
    double *xi, *xr;
    double tempr, tempi, c, s, cold, cfreq, sfreq; /* temps for pre and post twiddle */
    double freq = TWOPI / N;
    double cosfreq8, sinfreq8;
    int i, n;

    xi = (double*)malloc((N >> 2)*sizeof(double));
    xr = (double*)malloc((N >> 2)*sizeof(double));

    /* prepare for recurrence relation in pre-twiddle */
    cfreq = cos (freq);
    sfreq = sin (freq);
    cosfreq8 = cos (freq * 0.125);
    sinfreq8 = sin (freq * 0.125);
    c = cosfreq8;
    s = sinfreq8;

    for (i = 0; i < (N >> 2); i++) {
        /* calculate real and imaginary parts of g(n) or G(p) */
        n = (N >> 1) - 1 - 2 * i;

        if (i < (N >> 3))
            tempr = data [(N >> 2) + n] + data [N + (N >> 2) - 1 - n]; /* use second form of e(n) for n = N / 2 - 1 - 2i */
        else
            tempr = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n = N / 2 - 1 - 2i */

        n = 2 * i;
        if (i < (N >> 3))
            tempi = data [(N >> 2) + n] - data [(N >> 2) - 1 - n]; /* use first form of e(n) for n=2i */
        else
            tempi = data [(N >> 2) + n] + data [N + (N >> 2) - 1 - n]; /* use second form of e(n) for n=2i*/

        /* calculate pre-twiddled FFT input */
        xr[i] = tempr * c + tempi * s;
        xi[i] = tempi * c - tempr * s;

        /* use recurrence to prepare cosine and sine for next value of i */
        cold = c;
        c = c * cfreq - s * sfreq;
        s = s * cfreq + cold * sfreq;
    }

    /* Perform in-place complex FFT of length N/4 */
	exac_fft(xr, xi, N / 4, 1);

    /* prepare for recurrence relations in post-twiddle */
    c = cosfreq8;
    s = sinfreq8;

    /* post-twiddle FFT output and then get output data */
    for (i = 0; i < (N >> 2); i++) {
        /* get post-twiddled FFT output  */
        tempr = 2. * (xr[i] * c + xi[i] * s);
        tempi = 2. * (xi[i] * c - xr[i] * s);

        /* fill in output values */
        data [2 * i] = -tempr;   /* first half even */
        data [(N >> 1) - 1 - 2 * i] = tempi;  /* first half odd */
        data [(N >> 1) + 2 * i] = -tempi;  /* second half even */
        data [N - 1 - 2 * i] = tempr;  /* second half odd */

        /* use recurrence to prepare cosine and sine for next value of i */
        cold = c;
        c = c * cfreq - s * sfreq;
        s = s * cfreq + cold * sfreq;
    }

    if (xr) free(xr);
    if (xi) free(xi);
}

static void IMDCT(double *data, int N)
{
    double *xi, *xr;
    double tempr, tempi, c, s, cold, cfreq, sfreq; /* temps for pre and post twiddle */
    double freq = 2.0 * M_PI / N;
    double fac, cosfreq8, sinfreq8;
    int i;

    xi = (double*)malloc((N >> 2)*sizeof(double));
    xr = (double*)malloc((N >> 2)*sizeof(double));

    /* Choosing to allocate 2/N factor to Inverse Xform! */
	
    fac = 2. / N; /* remaining 2/N from 4/N IFFT factor */

    /* prepare for recurrence relation in pre-twiddle */
    cfreq = cos (freq);
    sfreq = sin (freq);
    cosfreq8 = cos (freq * 0.125);
    sinfreq8 = sin (freq * 0.125);
    c = cosfreq8;
    s = sinfreq8;

    for (i = 0; i < (N >> 2); i++) {
        /* calculate real and imaginary parts of g(n) or G(p) */
        tempr = -data[2 * i];
        tempi = data[(N >> 1) - 1 - 2 * i];

        /* calculate pre-twiddled FFT input */
        xr[i] = tempr * c - tempi * s;
        xi[i] = tempi * c + tempr * s;

        /* use recurrence to prepare cosine and sine for next value of i */
        cold = c;
        c = c * cfreq - s * sfreq;
        s = s * cfreq + cold * sfreq;
    }

    /* Perform in-place complex IFFT of length N/4 */
	exac_fft(xr, xi, N / 4, -1);

    /* prepare for recurrence relations in post-twiddle */
    c = cosfreq8;
    s = sinfreq8;

    /* post-twiddle FFT output and then get output data */
    for (i = 0; i < (N >> 2); i++) {

        /* get post-twiddled FFT output  */
        tempr = fac * (xr[i] * c - xi[i] * s);
        tempi = fac * (xi[i] * c + xr[i] * s);

        /* fill in output values */
        data [(N >> 1) + (N >> 2) - 1 - 2 * i] = tempr;
        if (i < (N >> 3))
            data [(N >> 1) + (N >> 2) + 2 * i] = tempr;
        else
            data [2 * i - (N >> 2)] = -tempr;

        data [(N >> 2) + 2 * i] = tempi;
        if (i < (N >> 3))
            data [(N >> 2) - 1 - 2 * i] = -tempi;
        else
            data [(N >> 2) + N - 1 - 2*i] = tempi;

        /* use recurrence to prepare cosine and sine for next value of i */
        cold = c;
        c = c * cfreq - s * sfreq;
        s = s * cfreq + cold * sfreq;
    }

    if (xr) free(xr);
    if (xi) free(xi);
}

void FilterBank(double *p_in_data,
                double *p_overlap,
                double *p_out_mdct)
{
    double *transf_buf;
	int i;

	transf_buf = (double*)malloc(2 * BLOCK_LEN_LONG * sizeof(double));

	memcpy(transf_buf, p_overlap, BLOCK_LEN_LONG * sizeof(double));
	memcpy(transf_buf + BLOCK_LEN_LONG, p_in_data, BLOCK_LEN_LONG * sizeof(double));
	memcpy(p_overlap, p_in_data, BLOCK_LEN_LONG * sizeof(double));
	
	for ( i = 0 ; i < BLOCK_LEN_LONG * 2; i++){
		transf_buf[i] = transf_buf[i] * adAnaAndSynWindow[i];
	}
	
	MDCT(transf_buf, 2 * BLOCK_LEN_LONG);
	memcpy(p_out_mdct, transf_buf, BLOCK_LEN_LONG * sizeof(double));

	if (transf_buf)
		free(transf_buf);
}

