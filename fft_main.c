#include <time.h>
#include "fast_fourier_transform.c"

#define DEBUG_OTHER 1
#define PRINT_RESULT 10

int main()
{
	int i, N = (int)(pow(2,6)*pow(3,5)*pow(5,5));
	clock_t t1, t2;
	double *x_re, *x_im;
	
	x_re = (double *) malloc(N*sizeof(double));
	x_im = (double *) malloc(N*sizeof(double));
	
	
	/* main */
	
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	
	FFT_general(x_re, x_im, N);
	
	t2 = clock();
	
	printf("time = %f\n",1.0*(t2-t1)/(double) CLOCKS_PER_SEC);
	
	#if PRINT_RESULT
	for(i=0;i<N;i+=N/PRINT_RESULT)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif
	
	
	/* other */
	
	#if DEBUG_OTHER 
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	
	FFT_general_np_separated(x_re, x_im, N);
	
	t2 = clock();
	
	printf("time = %f\n",1.0*(t2-t1)/(double) CLOCKS_PER_SEC);
	
	#if PRINT_RESULT
	for(i=0;i<N;i+=N/10)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif
	#endif
	return 0;
}

