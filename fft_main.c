
#include <time.h>
#include "fast_fourier_transform.c"

#define DEBUG_OTHER 1
#define FFT_TIMES 0
#define PRINT_RESULT 0

int main()
{
	int i;
	//int N = (int)(pow(2,2)*pow(3,1)*pow(5,0));
	//int N = (int)(pow(2,5)*pow(3,5)*pow(5,5));
	//int N = 134217728; //70.78
	//int N = 33554432; //10....
	//int N = 43046721 //7.765
	int N = 14348907; //2.375
	double *x_re, *x_im;
	clock_t t1, t2, t3, t4;
	x_re = (double *) malloc(N*sizeof(double));
	x_im = (double *) malloc(N*sizeof(double));
	
	/* main */
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	t1 = clock();
	#if FFT_TIMES 
	for(i=0;i<FFT_TIMES;++i)
	#endif
	{
		FFT_general(x_re, x_im, N);
		//rearrange(x_re, x_im, N);
		//butterfly0(x_re, x_im, N);
	}
	t2 = clock();
	printf("time = %f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
	#if PRINT_RESULT
	for(i=0;i<N;++i)
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
	#if FFT_TIMES 
	for(i=0;i<FFT_TIMES;++i)
	#endif
	{
		//bit_reverse0(x_re, x_im, N);
		//bit_reverse2(x_re, x_im, N);
		rearrange(x_re, x_im, N);
		t3 = clock();
		t4 = clock();
		butterfly0(x_re, x_im, N);
	}
	t2 = clock();
	printf("time2= %f\n",(t2-t4)/(double) CLOCKS_PER_SEC);
	printf("time3= %f\n",(t3-t1)/(double) CLOCKS_PER_SEC);
	#if PRINT_RESULT
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif
	#endif
	
	return;
}

