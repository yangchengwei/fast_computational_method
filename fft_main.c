#include "fast_fourier_transform.c"
#define DEBUG_OTHER 1
#define FFT_TIMES 2500
#define PRINT_RESULT 0

int main()
{
	int i;
	int N = (int)(pow(2,0)*pow(3,9)*pow(5,0));
	double x_re[N], x_im[N];
	clock_t t1, t2;
	
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
		bit_reverse(x_re, x_im, N);
		butterfly(x_re, x_im, N);
	}
	t2 = clock();
	printf("time=%f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
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
		bit_reverse(x_re, x_im, N);
		//bit_reverse2(x_re, x_im, N);
		//rearrange(x_re, x_im, N);
		butterfly0(x_re, x_im, N);
	}
	t2 = clock();
	printf("time=%f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
	#if PRINT_RESULT
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif
	#endif
	
	return;
}

