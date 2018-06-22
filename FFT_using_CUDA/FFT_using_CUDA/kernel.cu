
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

#define THREADS_PER_BLOCK 512

#define PRINT_RESULT 1
#define DEBUG_OTHER 0
#define DEBUG 0

/* constant */

__constant__ double w3_im2;
__constant__ double w5_re1;
__constant__ double w5_re2;
__constant__ double w5_im3;
__constant__ double w5_im4;
__constant__ double w7_re1;
__constant__ double w7_im1;
__constant__ double w7_re2;
__constant__ double w7_im2;
__constant__ double w7_re3;
__constant__ double w7_im3;

/* functions declaration */

cudaError_t fftCuda(double *x_re, double *x_im, int N);

/* kernel */

__global__ void butterflyKernel_7(double *x_re, double *x_im, double *w_re, double *w_im, const int m, const int s, const int M)
{
	int k, A, B, C, D, E, F, G;
	double tA_re, tA_im, tB_re, tB_im, tw_re, tw_im, tC_re, tC_im;
	double tD_re, tD_im, tE_re, tE_im, tF_re, tF_im, tG_re, tG_im, t;
	
	if ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) >= M) return;
	
	k = (blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) % m;
	A = ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) / m) * s + k;
	B = A + m;
	C = B + m;
	D = C + m;
	E = D + m;
	F = E + m;
	G = F + m;

	tA_re = x_re[A];
	tA_im = x_im[A];
	tB_re = w_re[k] * x_re[B] - w_im[k] * x_im[B];
	tB_im = w_re[k] * x_im[B] + w_im[k] * x_re[B];
	tw_re = w_re[k] * w_re[k] - w_im[k] * w_im[k];
	tw_im = 2 * w_re[k] * w_im[k];
	tC_re = tw_re*x_re[C] - tw_im*x_im[C];
	tC_im = tw_re*x_im[C] + tw_im*x_re[C];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tD_re = tw_re*x_re[D] - tw_im*x_im[D];
	tD_im = tw_re*x_im[D] + tw_im*x_re[D];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tE_re = tw_re*x_re[E] - tw_im*x_im[E];
	tE_im = tw_re*x_im[E] + tw_im*x_re[E];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tF_re = tw_re*x_re[F] - tw_im*x_im[F];
	tF_im = tw_re*x_im[F] + tw_im*x_re[F];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tG_re = tw_re*x_re[G] - tw_im*x_im[G];
	tG_im = tw_re*x_im[G] + tw_im*x_re[G];
	
	x_re[A] = tA_re + tB_re + tC_re + tD_re + tE_re + tF_re + tG_re;
	x_re[B] = tA_re + (tB_re + tG_re)*w7_re1 + (tG_im - tB_im)*w7_im1 + (tC_re + tF_re)*w7_re2 + (tF_im - tC_im)*w7_im2 + (tD_re + tE_re)*w7_re3 + (tE_im - tD_im)*w7_im3;
	x_re[C] = tA_re + (tB_re + tG_re)*w7_re2 + (tG_im - tB_im)*w7_im2 + (tC_re + tF_re)*w7_re3 + (tC_im - tF_im)*w7_im3 + (tD_re + tE_re)*w7_re1 + (tD_im - tE_im)*w7_im1;
	x_re[D] = tA_re + (tB_re + tG_re)*w7_re3 + (tG_im - tB_im)*w7_im3 + (tC_re + tF_re)*w7_re1 + (tC_im - tF_im)*w7_im1 + (tD_re + tE_re)*w7_re2 + (tE_im - tD_im)*w7_im2;
	x_re[E] = tA_re + (tB_re + tG_re)*w7_re3 + (tB_im - tG_im)*w7_im3 + (tC_re + tF_re)*w7_re1 + (tF_im - tC_im)*w7_im1 + (tD_re + tE_re)*w7_re2 + (tD_im - tE_im)*w7_im2;
	x_re[F] = tA_re + (tB_re + tG_re)*w7_re2 + (tB_im - tG_im)*w7_im2 + (tC_re + tF_re)*w7_re3 + (tF_im - tC_im)*w7_im3 + (tD_re + tE_re)*w7_re1 + (tE_im - tD_im)*w7_im1;
	x_re[G] = tA_re + (tB_re + tG_re)*w7_re1 + (tB_im - tG_im)*w7_im1 + (tC_re + tF_re)*w7_re2 + (tC_im - tF_im)*w7_im2 + (tD_re + tE_re)*w7_re3 + (tD_im - tE_im)*w7_im3;
	x_im[A] = tA_im + tB_im + tC_im + tD_im + tE_im + tF_im + tG_im;
	x_im[B] = tA_im + (tB_im + tG_im)*w7_re1 + (tB_re - tG_re)*w7_im1 + (tC_im + tF_im)*w7_re2 + (tC_re - tF_re)*w7_im2 + (tD_im + tE_im)*w7_re3 + (tD_re - tE_re)*w7_im3;
	x_im[C] = tA_im + (tB_im + tG_im)*w7_re2 + (tB_re - tG_re)*w7_im2 + (tC_im + tF_im)*w7_re3 + (tF_re - tC_re)*w7_im3 + (tD_im + tE_im)*w7_re1 + (tE_re - tD_re)*w7_im1;
	x_im[D] = tA_im + (tB_im + tG_im)*w7_re3 + (tB_re - tG_re)*w7_im3 + (tC_im + tF_im)*w7_re1 + (tF_re - tC_re)*w7_im1 + (tD_im + tE_im)*w7_re2 + (tD_re - tE_re)*w7_im2;
	x_im[E] = tA_im + (tB_im + tG_im)*w7_re3 + (tG_re - tB_re)*w7_im3 + (tC_im + tF_im)*w7_re1 + (tC_re - tF_re)*w7_im1 + (tD_im + tE_im)*w7_re2 + (tE_re - tD_re)*w7_im2;
	x_im[F] = tA_im + (tB_im + tG_im)*w7_re2 + (tG_re - tB_re)*w7_im2 + (tC_im + tF_im)*w7_re3 + (tC_re - tF_re)*w7_im3 + (tD_im + tE_im)*w7_re1 + (tD_re - tE_re)*w7_im1;
	x_im[G] = tA_im + (tB_im + tG_im)*w7_re1 + (tG_re - tB_re)*w7_im1 + (tC_im + tF_im)*w7_re2 + (tF_re - tC_re)*w7_im2 + (tD_im + tE_im)*w7_re3 + (tE_re - tD_re)*w7_im3;
}
__global__ void butterflyKernel_5(double *x_re, double *x_im, double *w_re, double *w_im, const int m, const int s, const int M)
{
	int k, A, B, C, D, E;
	double tA_re, tA_im, tB_re, tB_im, tw_re, tw_im, tC_re, tC_im;
	double tD_re, tD_im, tE_re, tE_im, t;

	if ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) >= M) return;

	k = (blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) % m;
	A = ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) / m) * s + k;
	B = A + m;
	C = B + m;
	D = C + m;
	E = D + m;

	tA_re = x_re[A];
	tA_im = x_im[A];
	tB_re = w_re[k] * x_re[B] - w_im[k] * x_im[B];
	tB_im = w_re[k] * x_im[B] + w_im[k] * x_re[B];
	tw_re = w_re[k] * w_re[k] - w_im[k] * w_im[k];
	tw_im = 2 * w_re[k] * w_im[k];
	tC_re = tw_re*x_re[C] - tw_im*x_im[C];
	tC_im = tw_re*x_im[C] + tw_im*x_re[C];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tD_re = tw_re*x_re[D] - tw_im*x_im[D];
	tD_im = tw_re*x_im[D] + tw_im*x_re[D];
	t = tw_re;
	tw_re = tw_re*w_re[k] - tw_im*w_im[k];
	tw_im = t    *w_im[k] + tw_im*w_re[k];
	tE_re = tw_re*x_re[E] - tw_im*x_im[E];
	tE_im = tw_re*x_im[E] + tw_im*x_re[E];

	x_re[A] = tA_re + tB_re + tC_re + tD_re + tE_re;
	x_re[B] = tA_re + w5_re1*(tB_re + tE_re) + w5_im4*(tB_im - tE_im) + w5_re2*(tC_re + tD_re) + w5_im3*(tC_im - tD_im);
	x_re[C] = tA_re + w5_re2*(tB_re + tE_re) + w5_im3*(tB_im - tE_im) + w5_re1*(tC_re + tD_re) + w5_im4*(tD_im - tC_im);
	x_re[D] = tA_re + w5_re2*(tB_re + tE_re) + w5_im3*(tE_im - tB_im) + w5_re1*(tC_re + tD_re) + w5_im4*(tC_im - tD_im);
	x_re[E] = tA_re + w5_re1*(tB_re + tE_re) + w5_im4*(tE_im - tB_im) + w5_re2*(tC_re + tD_re) + w5_im3*(tD_im - tC_im);
	x_im[A] = tA_im + tB_im + tC_im + tD_im + tE_im;
	x_im[B] = tA_im + w5_re1*(tB_im + tE_im) + w5_im4*(tE_re - tB_re) + w5_re2*(tC_im + tD_im) + w5_im3*(tD_re - tC_re);
	x_im[C] = tA_im + w5_re2*(tB_im + tE_im) + w5_im3*(tE_re - tB_re) + w5_re1*(tC_im + tD_im) + w5_im4*(tC_re - tD_re);
	x_im[D] = tA_im + w5_re2*(tB_im + tE_im) + w5_im3*(tB_re - tE_re) + w5_re1*(tC_im + tD_im) + w5_im4*(tD_re - tC_re);
	x_im[E] = tA_im + w5_re1*(tB_im + tE_im) + w5_im4*(tB_re - tE_re) + w5_re2*(tC_im + tD_im) + w5_im3*(tC_re - tD_re);
}
__global__ void butterflyKernel_3(double *x_re, double *x_im, double *w_re, double *w_im, const int m, const int s, const int M)
{
	int k, A, B, C;
	double tA_re, tA_im, tB_re, tB_im, tw_re, tw_im, tC_re, tC_im;

	if ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) >= M) return;

	k = (blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) % m;
	A = ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) / m) * s + k;
	B = A + m;
	C = B + m;

	tA_re = x_re[A];
	tA_im = x_im[A];
	tB_re = w_re[k] * x_re[B] - w_im[k] * x_im[B];
	tB_im = w_re[k] * x_im[B] + w_im[k] * x_re[B];
	tw_re = w_re[k] * w_re[k] - w_im[k] * w_im[k];
	tw_im = 2 * w_re[k] * w_im[k];
	tC_re = tw_re*x_re[C] - tw_im*x_im[C];
	tC_im = tw_re*x_im[C] + tw_im*x_re[C];

	x_re[A] = tA_re + tB_re + tC_re;
	x_im[A] = tA_im + tB_im + tC_im;
	x_re[B] = tA_re - 0.5*(tB_re + tC_re) + w3_im2*(tB_im - tC_im);
	x_im[B] = tA_im + w3_im2*(tC_re - tB_re) - 0.5*(tB_im + tC_im);
	x_re[C] = tA_re - 0.5*(tB_re + tC_re) + w3_im2*(tC_im - tB_im);
	x_im[C] = tA_im + w3_im2*(tB_re - tC_re) - 0.5*(tB_im + tC_im);
}
__global__ void butterflyKernel_2(double *x_re, double *x_im, double *w_re, double *w_im, const int m, const int s, const int M)
{
	int k, A, B;
	double tA_re, tA_im, tB_re, tB_im;

	if ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) >= M) return;

	k = (blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) % m;
	A = ((blockIdx.x * THREADS_PER_BLOCK + threadIdx.x) / m) * s + k;
	B = A + m;

	tA_re = x_re[A];
	tA_im = x_im[A];
	tB_re = w_re[k] * x_re[B] - w_im[k] * x_im[B];
	tB_im = w_re[k] * x_im[B] + w_im[k] * x_re[B];

	x_re[A] = tA_re + tB_re;
	x_re[B] = tA_re - tB_re;
	x_im[A] = tA_im + tB_im;
	x_im[B] = tA_im - tB_im;
}

/* main */

int main()
{
	int i, N;
	double *x_re, *x_im, t;
	clock_t t1, t2;
	cudaError_t cudaStatus;

	printf("N=");
	scanf("%d", &N);
	printf("N=%d\n", N);

	x_re = (double *)malloc(N * sizeof(double));
	x_im = (double *)malloc(N * sizeof(double));

	/* initial CUDA */

	// Choose which GPU to run on.
	cudaFree(0);
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }

	// constant 
	t = sqrt(3.0) / 2.0;	cudaStatus = cudaMemcpyToSymbol(w3_im2, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w3_im2\n"); system("pause"); exit(cudaStatus); }
	t = cos(2.0*M_PI / 5.0);	cudaStatus = cudaMemcpyToSymbol(w5_re1, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = cos(4.0*M_PI / 5.0);	cudaStatus = cudaMemcpyToSymbol(w5_re2, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = sin(4.0*M_PI / 5.0);	cudaStatus = cudaMemcpyToSymbol(w5_im3, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = sin(2.0*M_PI / 5.0);	cudaStatus = cudaMemcpyToSymbol(w5_im4, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = cos(2.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_re1, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = -sin(2.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_im1, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = cos(4.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_re2, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = -sin(4.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_im2, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = cos(6.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_re3, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }
	t = -sin(6.0*M_PI / 7.0);	cudaStatus = cudaMemcpyToSymbol(w7_im3, &t, 8);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess w5_re1\n"); system("pause"); exit(cudaStatus); }

	/* main */

	for (i = 0; i<N; ++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}

	t1 = clock();

	fftCuda(x_re, x_im, N);

	t2 = clock();

	printf("time = %f\n", 1.0*(t2 - t1) / (double)CLOCKS_PER_SEC);

#if PRINT_RESULT
	system("pause");
	for (i = 0; i<N; i++)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
#endif

	/* other */

#if DEBUG_OTHER 

	for (i = 0; i<N; ++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}

	t1 = clock();

	//fftHost(x_re, x_im, N);

	t2 = clock();

	printf("time = %f\n", 1.0*(t2 - t1) / (double)CLOCKS_PER_SEC);

#if PRINT_RESULT
	system("pause");
	for (i = 0; i<N; i++)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
#endif

#endif

	system("pause");
	return 0;
}

/* functions definition */

cudaError_t fftCuda(double *x_re, double *x_im, int N)
{
	/* initial */

	// Device variable
	int memorySize = N * sizeof(double);
	double *dev_x_re, *dev_x_im;
	double *dev_w_re, *dev_w_im;
	cudaError_t cudaStatus;

	// Device memory
	cudaStatus = cudaMalloc((void**)&dev_x_re, memorySize);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }
	cudaStatus = cudaMalloc((void**)&dev_x_im, memorySize);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }
	cudaStatus = cudaMalloc((void**)&dev_w_re, memorySize);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega allocate %s\n", cudaGetErrorString(cudaStatus)); system("pause"); exit(cudaStatus); }
	cudaStatus = cudaMalloc((void**)&dev_w_im, memorySize);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega allocate %s\n", cudaGetErrorString(cudaStatus)); system("pause"); exit(cudaStatus); }

	// Host variable
	int copy_N = N, power_2 = 0, power_3 = 0, power_5 = 0, power_7 = 0, power_75 = 0, power_753 = 0, power_sum = 0;
	int p, q, i, k, m, s, M, step, gate, add, *order;
	double w_N_re, w_N_im, t, *temp_re, *temp_im;

	// Power Computation
	while (copy_N % 7 == 0) { power_7++; copy_N /= 7; }
	while (copy_N % 5 == 0) { power_5++; copy_N /= 5; }
	while (copy_N % 3 == 0) { power_3++; copy_N /= 3; }
	while (copy_N % 2 == 0) { power_2++; copy_N /= 2; }
	if (copy_N != 1) { printf("ERROR: N is not radix-2,3,5,7 !\n"); system("pause"); exit(EXIT_FAILURE);}
	power_sum = power_7 + power_5 + power_3 + power_2;
	power_753 = power_7 + power_5 + power_3;
	power_75 = power_7 + power_5;

	// Host memory
	temp_re = (double *)malloc(N * sizeof(double));
	if (temp_re == NULL) { printf("Failed to allocate host memory temp_re!\n"); system("pause"); exit(EXIT_FAILURE); }
	temp_im = (double *)malloc(N * sizeof(double));
	if (temp_im == NULL) { printf("Failed to allocate host memory temp_im!\n"); system("pause"); exit(EXIT_FAILURE); }
	order = (int *)malloc(power_sum * sizeof(int));
	if (order == NULL) { printf("Failed to allocate host memory order!\n"); system("pause"); exit(EXIT_FAILURE); }

	// order
	for (i = 0; i<power_7; i++)				order[i] = 7;
	for (i = power_7; i<power_75; i++)		order[i] = 5;
	for (i = power_75; i<power_753; i++)	order[i] = 3;
	for (i = power_753; i<power_sum; i++)	order[i] = 2;



	/* FFT */

	/* bit reverse */

	// copy x
	for (i = 0; i<N; i++)
	{
		temp_re[i] = x_re[i];
		temp_im[i] = x_im[i];
	}

#if DEBUG
	clock_t T1, T2;
	T1 = clock();
#endif

	// bit reverse main
	step = N / order[0];
	q = step;			// first change
	for (p = 1; p<N - 1; p++)
	{
		// change value
		x_re[p] = temp_re[q];
		x_im[p] = temp_im[q];

		// compute next place
		i = 0;
		add = step;
		gate = (order[i++] - 1)*add;
		while (q >= gate && gate > 0)
		{
			q = q - gate;
			add = add / order[i];
			gate = (order[i++] - 1)*add;
		}
		q = q + add;
	}

#if DEBUG
	T2 = clock();
	printf("bit reverse time = %f\n", 1.0*(T2 - T1) / (double)CLOCKS_PER_SEC);
#endif

	// length N array copy
	cudaStatus = cudaMemcpy(dev_x_re, x_re, memorySize, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }
	cudaStatus = cudaMemcpy(dev_x_im, x_im, memorySize, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }



	/* butterfly */

	m = 1;
	s = 1;
	// parallel
	for (i = 0; i<power_7; i++)
	{
		s *= order[i];

		// omega computation
		temp_re[0] = 1.0;
		temp_im[0] = 0.0;
		w_N_re = cos(2.0*M_PI / s);
		w_N_im = -sin(2.0*M_PI / s);
		for (k = 1; k<m; ++k)
		{
			temp_re[k] = w_N_re*temp_re[k - 1] - w_N_im*temp_im[k - 1];
			temp_im[k] = w_N_re*temp_im[k - 1] + w_N_im*temp_re[k - 1];
		}

		// omega copy
		cudaStatus = cudaMemcpy(dev_w_re, temp_re, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }
		cudaStatus = cudaMemcpy(dev_w_im, temp_im, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }

		// kernel
		M = N / order[i];
#if DEBUG
		printf("blockSize, threadSize: %d, %d\n", M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK);
#endif
		butterflyKernel_7 << < M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK >> >
			(dev_x_re, dev_x_im, dev_w_re, dev_w_im, m, s, M);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: kernel\n"); system("pause"); exit(cudaStatus); }

		m *= order[i];
	}
	for (i = power_7; i<power_75; i++)
	{
		s *= order[i];

		// omega computation
		temp_re[0] = 1.0;
		temp_im[0] = 0.0;
		w_N_re = cos(2.0*M_PI / s);
		w_N_im = -sin(2.0*M_PI / s);
		for (k = 1; k<m; ++k)
		{
			temp_re[k] = w_N_re*temp_re[k - 1] - w_N_im*temp_im[k - 1];
			temp_im[k] = w_N_re*temp_im[k - 1] + w_N_im*temp_re[k - 1];
		}

		// omega copy
		cudaStatus = cudaMemcpy(dev_w_re, temp_re, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }
		cudaStatus = cudaMemcpy(dev_w_im, temp_im, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }

		// kernel
		M = N / order[i];
#if DEBUG
		printf("blockSize, threadSize: %d, %d\n", M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK);
#endif
		butterflyKernel_5 << < M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK >> >
			(dev_x_re, dev_x_im, dev_w_re, dev_w_im, m, s, M);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: kernel\n"); system("pause"); exit(cudaStatus); }

		m *= order[i];
	}
	for (i = power_75; i<power_753; i++)
	{
		s *= order[i];

		// omega computation
		temp_re[0] = 1.0;
		temp_im[0] = 0.0;
		w_N_re = cos(2.0*M_PI / s);
		w_N_im = -sin(2.0*M_PI / s);
		for (k = 1; k<m; ++k)
		{
			temp_re[k] = w_N_re*temp_re[k - 1] - w_N_im*temp_im[k - 1];
			temp_im[k] = w_N_re*temp_im[k - 1] + w_N_im*temp_re[k - 1];
		}

		// omega copy
		cudaStatus = cudaMemcpy(dev_w_re, temp_re, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }
		cudaStatus = cudaMemcpy(dev_w_im, temp_im, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }

		// kernel
		M = N / order[i];
#if DEBUG
		printf("blockSize, threadSize: %d, %d\n", M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK);
#endif
		butterflyKernel_3 << < M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK >> >
			(dev_x_re, dev_x_im, dev_w_re, dev_w_im, m, s, M);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: kernel\n"); system("pause"); exit(cudaStatus); }

		m *= order[i];
	}
	for (i = power_753; i<power_sum; i++)
	{
		s *= order[i];

		// omega computation
		temp_re[0] = 1.0;
		temp_im[0] = 0.0;
		w_N_re = cos(2.0*M_PI / s);
		w_N_im = -sin(2.0*M_PI / s);
		for (k = 1; k<m; ++k)
		{
			temp_re[k] = w_N_re*temp_re[k - 1] - w_N_im*temp_im[k - 1];
			temp_im[k] = w_N_re*temp_im[k - 1] + w_N_im*temp_re[k - 1];
		}

		// omega copy
		cudaStatus = cudaMemcpy(dev_w_re, temp_re, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }
		cudaStatus = cudaMemcpy(dev_w_im, temp_im, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: omega copy\n"); system("pause"); exit(cudaStatus); }

		// kernel
		M = N / order[i];
#if DEBUG
		printf("blockSize, threadSize: %d, %d\n", M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK);
#endif
		butterflyKernel_2 << < M / THREADS_PER_BLOCK + 1, THREADS_PER_BLOCK >> >
			(dev_x_re, dev_x_im, dev_w_re, dev_w_im, m, s, M);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess: kernel\n"); system("pause"); exit(cudaStatus); }

		m *= order[i];
	}

	/* FINISH */

	// Memory
	cudaStatus = cudaMemcpy(x_re, dev_x_re, memorySize, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }
	cudaStatus = cudaMemcpy(x_im, dev_x_im, memorySize, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) { printf("cudaNotSuccess\n"); system("pause"); exit(cudaStatus); }

	return cudaStatus;
}