#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
	int m, k, N=1 << 4;
	double *x, *y_DCT, *y_DST;
	double c_DCT, s_DCT, c_DST, s_DST;
	double ca_DCT, sa_DCT, ca_DST, sa_DST, temp;
	
	x     = (double *) malloc(N*sizeof(double));
	y_DCT = (double *) malloc(N*sizeof(double));
	y_DST = (double *) malloc(N*sizeof(double));
	
	for(k=0;k<N;k++)
	{
		x[k] = k+1;
	}
	
	for(m=0;m<N;m++)
	{
		y_DST[m] = 0.0;
		y_DCT[m] = 0.0;
		
		temp = M_PI*m/N;
		ca_DCT = cos(temp);
		sa_DCT = sin(temp);
		temp = M_PI*m*0.5/N;
		c_DCT = cos(temp);
		s_DCT = sin(temp);
		temp = M_PI*(m+1)/(N+1);
		ca_DST = cos(temp);
		sa_DST = sin(temp);
		c_DST = ca_DST;
		s_DST = sa_DST;
		for(k=0;k<N;k++)
		{
			y_DCT[m] += x[k]*c_DCT;
			y_DST[m] += x[k]*s_DST;
			temp = c_DCT;
			c_DCT = c_DCT*ca_DCT - s_DCT*sa_DCT;
			s_DCT = s_DCT*ca_DCT + temp*sa_DCT;
			temp = c_DST;
			c_DST = c_DST*ca_DST - s_DST*sa_DST;
			s_DST = s_DST*ca_DST + temp*sa_DST;
		}
	}
	printf("DCT-II result:\n");
	for(k=0;k<N;k++)
		printf("DCT[%d] = %.9f\n", k, y_DCT[k]);
	
	printf("\nDST-I result:\n");
	for(k=0;k<N;k++)
		printf("DST[%d] = %.9f\n", k, y_DST[k]);
	
	return 0;
}
