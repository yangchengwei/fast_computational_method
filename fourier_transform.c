#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int main(void)
{
	int k, n, N=1 << 15;
	double *X_re, *X_im, *x_re, *x_im;
	double a, c, s, ca, sa;
	double tempc;
	
	X_re = (double *) malloc(N*sizeof(double));
	X_im = (double *) malloc(N*sizeof(double));
	x_re = (double *) malloc(N*sizeof(double));
	x_im = (double *) malloc(N*sizeof(double));
	
	for(n=0;n<N;n++)
	{
		x_re[n] = n+1;
		x_im[n] = 0.0;
	}
	
	for(k=0;k<N;k++)
	{
		X_re[k] = 0.0;
		X_im[k] = 0.0;
		
		a = 2*M_PI*k/N;
		ca = cos(a);
		sa = sin(a);
		c = 1.0;
		s = 0.0;
		for(n=0;n<N;n++)
		{
			//X_re[k] += x_re[n]*cos(2*M_PI*k*n/N)+x_im[n]*sin(2*M_PI*k*n/N);
			//X_im[k] += x_im[n]*cos(2*M_PI*k*n/N)-x_re[n]*sin(2*M_PI*k*n/N);
			
			X_re[k] += x_re[n]*c+x_im[n]*s;
			X_im[k] += x_im[n]*c-x_re[n]*s;
			tempc = c;
			c = c*ca - s*sa;
			s = s*ca + tempc*sa;
		}
	}
	
	for(k=0;k<N;k++)
	{
		if (X_im[k]>=0)
		{
			printf("%.9f+%.9fi\n", X_re[k], X_im[k]);
		}
		else
		{
			printf("%.9f%.9fi\n", X_re[k], X_im[k]);
		}
	}
}
