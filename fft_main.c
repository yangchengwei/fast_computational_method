#include <time.h>
#include "fast_fourier_transform.c"

#define DEBUG_OTHER 0
#define PRINT_RESULT 1

int main()
{
	int i, j, N;
	double *x_re, *x_im;
	clock_t t1, t2;
	
	printf("N=");
	scanf("%d",&N);
	printf("N=%d\n",N);
	
	if (N==-1){
		int *N_array;
		
		N_array = (int *) malloc(7*sizeof(double));
		N_array[0]=16777216;
		N_array[1]=14348907;
		N_array[2]=48828125;
		N_array[3]=10077696;
		N_array[4]=11390625;
		N_array[5]=10000000;
		N_array[6]=10935000;
		
		for(j=0;j<7;j++)
		{
			N=N_array[j];
			printf("N=%d\n",N);
			
			x_re = (double *) malloc(N*sizeof(double));
			x_im = (double *) malloc(N*sizeof(double));
			
			for(i=0;i<N;++i)
			{
				x_re[i] = i;
				x_im[i] = 0.0;
			}
			
			t1 = clock();
			
			FFT_general(x_re, x_im, N);
			
			t2 = clock();
			
			printf("time = %f\n",1.0*(t2-t1)/(double) CLOCKS_PER_SEC);
			
			free(x_re);
			free(x_im);
		}
	}
	else
	{
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
		system("pause");
		for(i=0;i<N;i++)
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
		system("pause");
		for(i=0;i<N;i++)
		{
			printf("%f + %f i\n", x_re[i], x_im[i]);
		}
		#endif
		
		#endif
	}
	
	return 0;
}

