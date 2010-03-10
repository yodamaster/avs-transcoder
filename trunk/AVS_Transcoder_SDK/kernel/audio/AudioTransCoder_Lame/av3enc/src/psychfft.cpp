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
#include "psychfft.h"

Complex set(double r,double i)
{
	Complex temp;
	temp.real=r;
	temp.imag=i;
	return temp;
}

Complex conjg(Complex x)
{
	Complex temp;
	temp.real=x.real;
	temp.imag=-x.imag;
	return temp;
}

Complex  cmulc(Complex c1,Complex c2)
{
	Complex temp;
	temp.real=c1.real*c2.real -c1.imag *c2.imag;
	temp.imag =c1.imag*c2.real +c1.real *c2.imag ;
	return temp;
}

Complex cmulf(Complex c1,float c2)
{
	Complex temp;
	temp.real =c1.real *c2;
	temp.imag =c1.imag *c2;
	return temp;
}

Complex caddc(Complex c1,Complex c2)
{
	Complex temp;
	temp.real =c1.real +c2.real ;
	temp.imag =c1.imag +c2.imag ;
	return temp;
}
  
Complex csubc(Complex c1,Complex c2)
{
	Complex temp;
	temp.real =c1.real -c2.real ;
	temp.imag =c1.imag -c2.imag ;
	return temp;
}

void prepsychfft(int n,int mode,int *nexp,Complex *w)
{
	int i,k,nt,nexp1;
	float s;
	Complex c1,c2;
	nexp1=1;nt=1;
	do{
		nt=1;
		for(i=1;i<=nexp1;i++)
			nt=2*nt;
		if(nt>=n)
			break;
		nexp1=nexp1+1;
	}while(1);
	
	if(nt==n)
	{
		s=8*atan(1.0)/(float)(nt);
		c1=set(cos(s),-sin(s));
		if(mode!=0) c1=conjg(c1);
		c2=set(1.,0);
		for(k=1;k<=nt;k++)
		{
			w[k]=c2;
			c2=cmulc(c2,c1);
		}
	}
	else
	{
		nexp1=-1;
	}
	*nexp=nexp1;

}

void psychfft(int n,int mode,float t,int nexp,Complex *w,Complex *x)
{
	Complex c1,c2;
	int k,mm,ll,j,jj,kk,i,nn,nv2,nm1;
	float s;
	mm=1;
	ll=n;
	for(k=1;k<=nexp;k++)
	{
		nn=ll/2;
		jj=mm+1;
		for(i=1;i<=n;i=i+ll)
		{
			kk=i+nn;
			c1=caddc(x[i],x[kk]);
			x[kk]=csubc(x[i],x[kk]);
			x[i]=c1;
		}
		if(nn==1)
			continue;
		else
		{
			for(j=2;j<=nn;j++)
			{
				c2=w[jj];
				for(i=j;i<=n;i=i+ll)
				{
					kk=i+nn;
					c1=caddc(x[i],x[kk]);
					x[kk]=cmulc(csubc(x[i],x[kk]),c2);
					x[i]=c1;
				}
				jj=jj+mm;
			}
			ll=nn;
			mm=mm*2;
		}
	}
	nv2=n/2;
	nm1=n-1;
	j=1;
	for(i=1;i<nm1;i++)
	{
		if(i>=j);
		else
		{
			c1=x[j];
			x[j]=x[i];
			x[i]=c1;
		}
		k=nv2;
		do
		{
			if(k>=j)
				break;
			else
			{
				j=j-k;
				k=k/2;
			}
		}while(1);
		j=j+k;
	}
	if(mode==0)
		s=t;
	else
		s=1./(t*(float)n);
	for(i=1;i<=n;i++)
		x[i]=cmulf(x[i],s);
}


