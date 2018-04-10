#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h> 

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double **A, *x, *b, T1;
	double first_time[4];
	int i, j, k=0, N;

	for(N=2000;N<=20000;N*=2)	// without parallel computing
	{
		A = (double **) malloc( N * sizeof(double*) );
		A[0] = (double *) malloc( N*N*sizeof(double));
		
		//for(i=1;i<N;++i) A[i] = A[i-1] + N; 	// cannot be parallel computed 
		#pragma omp parallel for
		for(i=1;i<N;++i) A[i] = A[0] + i*N;   	// can be parallel computed 
		
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		
		#pragma omp parallel for private(j)
		for (i=0;i<N;i++)
		{
			srand(time(NULL)+j);
			#pragma omp parallel for
			for(j=0;j<N;++j)
			{
				A[i][j] = rand();
			}
			x[i] = rand();
		}
				
		t1 = clock();
		for(i=0;i<N;++i) 
		{
			b[i] = 0.0;
			for(j=0;j<N;++j)
			{
				b[i] += A[i][j]*x[j];
			}
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("N= %d , time: %f \n", N, T1);
		first_time[k++]=T1;
		
		free(b);
		free(x);
		free(A[0]);
		free(A);	
	} 
	k=0;
	for(N=2000;N<=20000;N*=2)	// with parallel computing
	{
		A = (double **) malloc( N * sizeof(double*) );
		A[0] = (double *) malloc( N*N*sizeof(double));
		
		//for(i=1;i<N;++i) A[i] = A[i-1] + N; 	// cannot be parallel computed 
		#pragma omp parallel for
		for(i=1;i<N;++i) A[i] = A[0] + i*N;   	// can be parallel computed 
		
		x = (double *) malloc( N * sizeof(double) );
		b = (double *) malloc( N * sizeof(double) );
		
		#pragma omp parallel for private(j)
		for (i=0;i<N;i++)
		{
			srand(time(NULL)+j);
			#pragma omp parallel for
			for(j=0;j<N;++j)
			{
				A[i][j] = rand();
			}
			x[i] = rand();
		}
				
		t1 = clock();
		
		#pragma omp parallel for private(j)
		for(i=0;i<N;++i) 
		{
			b[i] = 0.0;
			for(j=0;j<N;++j)
			{
				b[i] += A[i][j]*x[j];
			}
		}
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("N= %d , time: %f (parallel computing: %f seconds faster.)\n",
				N, T1, first_time[k++]-T1);
		
		free(b);
		free(x);
		free(A[0]);
		free(A);	
	} 
	
	return 0;
} 
