#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

#define DEBUG 0


/* Fast Fourier Transform (parallel) (all-in-one) */

int FFT_general(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2 + power_3 + power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d, power_sum=%d\n", power_2, power_3, power_5, power_sum);
    #endif
    
    int i, *order;
    int power_53 = power_5+power_3;
    order = (int *) malloc(power_sum*sizeof(int));
    for (i=0;i<power_5;i++)				order[i]=5;
    for (i=power_5;i<power_53;i++)		order[i]=3;
    for (i=power_53;i<power_sum;i++)	order[i]=2;
    
    /* FFT */
    
    /* bit reverse */
    
    int p,q;
    int step,gate,add;
	double t, *copy_re, *copy_im;
	copy_re = (double *) malloc(N*sizeof(double));
	copy_im = (double *) malloc(N*sizeof(double));
	#pragma omp parallel for
	for(i=0;i<N;i++)
	{
		copy_re[i] = x_re[i];
		copy_im[i] = x_im[i];
	}
	
    step = N/order[0];
    q = step;			// first change
    for(p=1;p<N-1;p++)
    {
		// change value
	    x_re[p] = copy_re[q];
	    x_im[p] = copy_im[q];
        
        // compute next place
        i = 0;
        add = step;
        gate = (order[i++]-1)*add;
        while(q >= gate && gate > 0)
        {
            #if DEBUG
            printf("q=%d,gate=%d,add=%d,i=%d\n",q,gate,add,i);
			printf("order[i]=%d\n",order[i]);
            #endif
            q = q-gate;
            add = add/order[i];
            gate = (order[i++]-1)*add;
        }
        q = q+add;
    }
    
	/* butterfly */
	
    int k, m, s;
	int A, B, C, D, E, M;
	double w_re, w_im, w_N_re, w_N_im;
	double tA_re, tB_re, tC_re, tD_re, tE_re;
	double tA_im, tB_im, tC_im, tD_im, tE_im;
	double tw_re, tw_im;
	
	double sqrt3_2 = sqrt(3.0)/2.0;
	double cos2pi_5 = cos(2.0*M_PI/5.0);
	double sin2pi_5 = sin(2.0*M_PI/5.0);
	double cos4pi_5 = cos(4.0*M_PI/5.0);
	double sin4pi_5 = sin(4.0*M_PI/5.0);
	
	m = 1;
	s = 1;
	for(i=0;i<power_5;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/(double)s);
		w_N_im = -sin(2.0*M_PI/(double)s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,C,D,E,tA_re,tB_re,tC_re,tD_re,tE_re,tA_im,tB_im,tC_im,tD_im,tE_im,tw_re,tw_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			C = B + m;
			D = C + m;
			E = D + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
			
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			tw_re = copy_re[k]*copy_re[k] - copy_im[k]*copy_im[k];
			tw_im = 2*copy_re[k]*copy_im[k];
			
			tC_re = tw_re*x_re[C] - tw_im*x_im[C];
			tC_im = tw_re*x_im[C] + tw_im*x_re[C];
			
			t = tw_re;
			tw_re = tw_re*copy_re[k] - tw_im*copy_im[k];
			tw_im = t    *copy_im[k] + tw_im*copy_re[k];
			
			tD_re = tw_re*x_re[D] - tw_im*x_im[D];
			tD_im = tw_re*x_im[D] + tw_im*x_re[D];
			
			t = tw_re;
			tw_re = tw_re*copy_re[k] - tw_im*copy_im[k];
			tw_im = t    *copy_im[k] + tw_im*copy_re[k];
			
			tE_re = tw_re*x_re[E] - tw_im*x_im[E];
			tE_im = tw_re*x_im[E] + tw_im*x_re[E];

			x_re[A] = tA_re + tB_re + tC_re + tD_re + tE_re;
			x_re[B] = tA_re + cos2pi_5*(tB_re+tE_re) + sin2pi_5*(tB_im-tE_im) + cos4pi_5*(tC_re+tD_re) + sin4pi_5*(tC_im-tD_im);
			x_re[C] = tA_re + cos4pi_5*(tB_re+tE_re) + sin4pi_5*(tB_im-tE_im) + cos2pi_5*(tC_re+tD_re) + sin2pi_5*(tD_im-tC_im); 
			x_re[D] = tA_re + cos4pi_5*(tB_re+tE_re) + sin4pi_5*(tE_im-tB_im) + cos2pi_5*(tC_re+tD_re) + sin2pi_5*(tC_im-tD_im);
			x_re[E] = tA_re + cos2pi_5*(tB_re+tE_re) + sin2pi_5*(tE_im-tB_im) + cos4pi_5*(tC_re+tD_re) + sin4pi_5*(tD_im-tC_im);
			
			x_im[A] = tA_im + tB_im + tC_im + tD_im + tE_im;
			x_im[B] = tA_im + cos2pi_5*(tB_im+tE_im) + sin2pi_5*(tE_re-tB_re) + cos4pi_5*(tC_im+tD_im) + sin4pi_5*(tD_re-tC_re);
			x_im[C] = tA_im + cos4pi_5*(tB_im+tE_im) + sin4pi_5*(tE_re-tB_re) + cos2pi_5*(tC_im+tD_im) + sin2pi_5*(tC_re-tD_re);
			x_im[D] = tA_im + cos4pi_5*(tB_im+tE_im) + sin4pi_5*(tB_re-tE_re) + cos2pi_5*(tC_im+tD_im) + sin2pi_5*(tD_re-tC_re);
			x_im[E] = tA_im + cos2pi_5*(tB_im+tE_im) + sin2pi_5*(tB_re-tE_re) + cos4pi_5*(tC_im+tD_im) + sin4pi_5*(tC_re-tD_re);
		}
		m *= order[i];
	}
	for(i=power_5;i<power_53;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,C,tA_re,tB_re,tC_re,tA_im,tB_im,tC_im,tw_re,tw_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			C = B + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
					
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			tw_re = copy_re[k]*copy_re[k] - copy_im[k]*copy_im[k];
			tw_im = 2*copy_re[k]*copy_im[k];
						
			tC_re = tw_re*x_re[C] - tw_im*x_im[C];
			tC_im = tw_re*x_im[C] + tw_im*x_re[C];
					
			x_re[A] = tA_re + tB_re + tC_re;
			x_re[B] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tB_im-tC_im);
			x_re[C] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tC_im-tB_im);
					
			x_im[A] = tA_im + tB_im + tC_im;
			x_im[B] = tA_im + sqrt3_2*(tC_re-tB_re) - 0.5*(tB_im+tC_im);
			x_im[C] = tA_im + sqrt3_2*(tB_re-tC_re) - 0.5*(tB_im+tC_im);
			
			/* teacher version (correctness test) */
			//double a = -0.5, b = -sqrt(3)/2;
			//x_re[A] = tA_re + tB_re + tC_re;
			//x_re[B] = tA_re + (a*tB_re - b*tB_im) + (a*tC_re + b*tC_im);
			//x_re[C] = tA_re + (a*tC_re - b*tC_im) + (a*tB_re + b*tB_im);
			//x_im[A] = tA_im + tB_im + tC_im;
			//x_im[B] = tA_im + (a*tB_im + b*tB_re) + (a*tC_im - b*tC_re);
			//x_im[C] = tA_im + (a*tB_im - b*tB_re) + (a*tC_im + b*tC_re);
		}
		m *= order[i];
	}
	for(i=power_53;i<power_sum;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,tA_re,tB_re,tA_im,tB_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
					
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			x_re[A] = tA_re + tB_re;
			x_re[B] = tA_re - tB_re;
			
			x_im[A] = tA_im + tB_im;
			x_im[B] = tA_im - tB_im;
		}
		m *= order[i];
	}
	
    /* FINISH */
    
    return 0;
}


/* Fast Fourier Transform (parallel) (separated) */

int bit_reverse(double *x_re, double *x_im, int N, int *order)
{
    int i,p,q;
    int step,gate,add;
	double t, *copy_re, *copy_im;
	copy_re = (double *) malloc(N*sizeof(double));
	copy_im = (double *) malloc(N*sizeof(double));
	#pragma omp parallel for
	for(i=0;i<N;i++)
	{
		copy_re[i] = x_re[i];
		copy_im[i] = x_im[i];
	}
	
    step = N/order[0];
    q = step;			// first change
    for(p=1;p<N-1;p++)
    {
		// change value
	    x_re[p] = copy_re[q];
	    x_im[p] = copy_im[q];
        
        // compute next place
        i = 0;
        add = step;
        gate = (order[i++]-1)*add;
        while(q >= gate && gate > 0)
        {
            #if DEBUG
            printf("q=%d,gate=%d,add=%d,i=%d\n",q,gate,add,i);
			printf("order[i]=%d\n",order[i]);
            #endif
            q = q-gate;
            add = add/order[i];
            gate = (order[i++]-1)*add;
        }
        q = q+add;
    }
	return 0;
}

int butterfly(double *x_re, double *x_im, int N, int power_5, int power_53, int power_sum, int *order)
{
    int i, k, m, s;
	int A, B, C, D, E, M;
	double w_re, w_im, w_N_re, w_N_im;
	double tA_re, tB_re, tC_re, tD_re, tE_re;
	double tA_im, tB_im, tC_im, tD_im, tE_im;
	double tw_re, tw_im;
	
	double sqrt3_2 = sqrt(3)/2.0;;
	double cos2pi_5 = cos(2.0*M_PI/5.0);
	double sin2pi_5 = sin(2.0*M_PI/5.0);
	double cos4pi_5 = cos(4.0*M_PI/5.0);
	double sin4pi_5 = sin(4.0*M_PI/5.0);
	
	double t, *copy_re, *copy_im;
	copy_re = (double *) malloc(N*sizeof(double));
	copy_im = (double *) malloc(N*sizeof(double));
	
	m = 1;
	s = 1;
	for(i=0;i<power_5;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,C,D,E,tA_re,tB_re,tC_re,tD_re,tE_re,tA_im,tB_im,tC_im,tD_im,tE_im,tw_re,tw_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			C = B + m;
			D = C + m;
			E = D + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
			
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			tw_re = copy_re[k]*copy_re[k] - copy_im[k]*copy_im[k];
			tw_im = 2*copy_re[k]*copy_im[k];
			
			tC_re = tw_re*x_re[C] - tw_im*x_im[C];
			tC_im = tw_re*x_im[C] + tw_im*x_re[C];
			
			t = tw_re;
			tw_re = tw_re*copy_re[k] - tw_im*copy_im[k];
			tw_im = t    *copy_im[k] + tw_im*copy_re[k];
			
			tD_re = tw_re*x_re[D] - tw_im*x_im[D];
			tD_im = tw_re*x_im[D] + tw_im*x_re[D];
			
			t = tw_re;
			tw_re = tw_re*copy_re[k] - tw_im*copy_im[k];
			tw_im = t    *copy_im[k] + tw_im*copy_re[k];
			
			tE_re = tw_re*x_re[E] - tw_im*x_im[E];
			tE_im = tw_re*x_im[E] + tw_im*x_re[E];

			x_re[A] = tA_re + tB_re + tC_re + tD_re + tE_re;
			x_re[B] = tA_re + cos2pi_5*(tB_re+tE_re) + sin2pi_5*(tB_im-tE_im) + cos4pi_5*(tC_re+tD_re) + sin4pi_5*(tC_im-tD_im);
			x_re[C] = tA_re + cos4pi_5*(tB_re+tE_re) + sin4pi_5*(tB_im-tE_im) + cos2pi_5*(tC_re+tD_re) + sin2pi_5*(tD_im-tC_im); 
			x_re[D] = tA_re + cos4pi_5*(tB_re+tE_re) + sin4pi_5*(tE_im-tB_im) + cos2pi_5*(tC_re+tD_re) + sin2pi_5*(tC_im-tD_im);
			x_re[E] = tA_re + cos2pi_5*(tB_re+tE_re) + sin2pi_5*(tE_im-tB_im) + cos4pi_5*(tC_re+tD_re) + sin4pi_5*(tD_im-tC_im);
			
			x_im[A] = tA_im + tB_im + tC_im + tD_im + tE_im;
			x_im[B] = tA_im + cos2pi_5*(tB_im+tE_im) + sin2pi_5*(tE_re-tB_re) + cos4pi_5*(tC_im+tD_im) + sin4pi_5*(tD_re-tC_re);
			x_im[C] = tA_im + cos4pi_5*(tB_im+tE_im) + sin4pi_5*(tE_re-tB_re) + cos2pi_5*(tC_im+tD_im) + sin2pi_5*(tC_re-tD_re);
			x_im[D] = tA_im + cos4pi_5*(tB_im+tE_im) + sin4pi_5*(tB_re-tE_re) + cos2pi_5*(tC_im+tD_im) + sin2pi_5*(tD_re-tC_re);
			x_im[E] = tA_im + cos2pi_5*(tB_im+tE_im) + sin2pi_5*(tB_re-tE_re) + cos4pi_5*(tC_im+tD_im) + sin4pi_5*(tC_re-tD_re);
		}
		m *= order[i];
	}
	for(i=power_5;i<power_53;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,C,tA_re,tB_re,tC_re,tA_im,tB_im,tC_im,tw_re,tw_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			C = B + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
					
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			tw_re = copy_re[k]*copy_re[k] - copy_im[k]*copy_im[k];
			tw_im = 2*copy_re[k]*copy_im[k];
						
			tC_re = tw_re*x_re[C] - tw_im*x_im[C];
			tC_im = tw_re*x_im[C] + tw_im*x_re[C];
					
			x_re[A] = tA_re + tB_re + tC_re;
			x_re[B] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tB_im-tC_im);
			x_re[C] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tC_im-tB_im);
					
			x_im[A] = tA_im + tB_im + tC_im;
			x_im[B] = tA_im + sqrt3_2*(tC_re-tB_re) - 0.5*(tB_im+tC_im);
			x_im[C] = tA_im + sqrt3_2*(tB_re-tC_re) - 0.5*(tB_im+tC_im);
		}
		m *= order[i];
	}
	for(i=power_53;i<power_sum;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			copy_re[k] = w_re;
			copy_im[k] = w_im;
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		#pragma omp parallel for private(k,A,B,tA_re,tB_re,tA_im,tB_im)
		for(M=0;M<N/order[i];++M)
		{
			k = M % m;
			A = (M / m)*s + k;
			B = A + m;
			
			tA_re = x_re[A];
			tA_im = x_im[A];
					
			tB_re = copy_re[k]*x_re[B] - copy_im[k]*x_im[B];
			tB_im = copy_re[k]*x_im[B] + copy_im[k]*x_re[B];
			
			x_re[A] = tA_re + tB_re;
			x_re[B] = tA_re - tB_re;
			
			x_im[A] = tA_im + tB_im;
			x_im[B] = tA_im - tB_im;
		}
		m *= order[i];
	}
	return 0;
}

int FFT_general_separated(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2 + power_3 + power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d, power_sum=%d\n", power_2, power_3, power_5, power_sum);
    #endif
    
    int i, *order;
    int power_53 = power_5+power_3;
    order = (int *) malloc(power_sum*sizeof(int));
    for (i=0;i<power_5;i++)				order[i]=5;
    for (i=power_5;i<power_53;i++)		order[i]=3;
    for (i=power_53;i<power_sum;i++)	order[i]=2;
    
    /* FFT */
    
    bit_reverse(x_re, x_im, N, order);
    
	butterfly  (x_re, x_im, N, power_5, power_53, power_sum, order);
	
    /* FINISH */
    
    return 0;
}


/* Fast Fourier Transform (not parallel) (separated) */

int bit_reverse_np(double *x_re, double *x_im, int N, int *order)
{
    int i,p,q;
    int step,gate,add;
	double t, *copy_re, *copy_im;
	copy_re = (double *) malloc(N*sizeof(double));
	copy_im = (double *) malloc(N*sizeof(double));
	clock_t t1, t2;
	t1 = clock();
	for(i=0;i<N;i++)
	{
		copy_re[i] = x_re[i];
		copy_im[i] = x_im[i];
	}
	t2 = clock();
	printf("bit time = %f\n",1.0*(t2-t1)/CLOCKS_PER_SEC);
	
    step = N/order[0];
    q = step;			// first change
    for(p=1;p<N-1;p++)
    {
		// change value
    	//printf("%d -> %d\n", p,q);
	    x_re[p] = copy_re[q];
	    x_im[p] = copy_im[q];
        
        // compute next place
        i = 0;
        add = step;
        gate = (order[i++]-1)*add;
        while(q >= gate && gate > 0)
        {
            #if DEBUG
            printf("q=%d,gate=%d,add=%d,i=%d\n",q,gate,add,i);
			printf("order[i]=%d\n",order[i]);
            #endif
            q = q-gate;
            add = add/order[i];
            gate = (order[i++]-1)*add;
        }
        q = q+add;
    }
	return 0;
}

int butterfly_np(double *x_re, double *x_im, int N, int power_5, int power_53, int power_sum, int *order)
{
    int i, k, m, s;
	int A, B, C, D, E;
	double w_re, w_im, w_N_re, w_N_im, t;
	double tA_re, tB_re, tC_re, tD_re, tE_re;
	double tA_im, tB_im, tC_im, tD_im, tE_im;
	double tw_re, tw_im;
	
	double sqrt3_2 = sqrt(3)/2.0;
	double cos2pi_5 = cos(2.0*M_PI/5.0);
	double sin2pi_5 = sin(2.0*M_PI/5.0);
	double cos4pi_5 = cos(4.0*M_PI/5.0);
	double sin4pi_5 = sin(4.0*M_PI/5.0);
	
	m = 1;
	s = 1;
	for(i=0;i<power_5;i++)
	{
		s = s*order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			for(A=k;A<N;A+=order[i]*m)
			{
				B = A + m;
				C = A + 2*m;
				D = A + 3*m;
				E = A + 4*m;
				
				t = x_re[B];
				x_re[B] = w_re*x_re[B] - w_im*x_im[B];
				x_im[B] = w_re*x_im[B] + w_im*t;
				
				tw_re = w_re*w_re - w_im*w_im;
				tw_im = w_re*w_im + w_im*w_re;
				
				t = x_re[C];
				x_re[C] = tw_re*x_re[C] - tw_im*x_im[C];
				x_im[C] = tw_re*x_im[C] + tw_im*t;
				
				t = tw_re;
				tw_re = w_re*tw_re - w_im*tw_im;
				tw_im = w_re*tw_im + w_im*t;
				
				t = x_re[D];
				x_re[D] = tw_re*x_re[D] - tw_im*x_im[D];
				x_im[D] = tw_re*x_im[D] + tw_im*t;
				
				t = tw_re;
				tw_re = w_re*tw_re - w_im*tw_im;
				tw_im = w_re*tw_im + w_im*t;
				
				t = x_re[E];
				x_re[E] = tw_re*x_re[E] - tw_im*x_im[E];
				x_im[E] = tw_re*x_im[E] + tw_im*t;
				
				tA_re = x_re[A] + x_re[B] + x_re[C] + x_re[D] + x_re[E];
				tB_re = x_re[A] + cos2pi_5*x_re[B] + sin2pi_5*x_im[B]
								+ cos4pi_5*x_re[C] + sin4pi_5*x_im[C]
								+ cos4pi_5*x_re[D] - sin4pi_5*x_im[D]
								+ cos2pi_5*x_re[E] - sin2pi_5*x_im[E];
				tC_re = x_re[A] + cos4pi_5*x_re[B] + sin4pi_5*x_im[B]
								+ cos2pi_5*x_re[C] - sin2pi_5*x_im[C]
								+ cos2pi_5*x_re[D] + sin2pi_5*x_im[D]
								+ cos4pi_5*x_re[E] - sin4pi_5*x_im[E];
				tD_re = x_re[A] + cos4pi_5*x_re[B] - sin4pi_5*x_im[B]
								+ cos2pi_5*x_re[C] + sin2pi_5*x_im[C]
								+ cos2pi_5*x_re[D] - sin2pi_5*x_im[D]
								+ cos4pi_5*x_re[E] + sin4pi_5*x_im[E];
				tE_re = x_re[A] + cos2pi_5*x_re[B] - sin2pi_5*x_im[B]
								+ cos4pi_5*x_re[C] - sin4pi_5*x_im[C]
								+ cos4pi_5*x_re[D] + sin4pi_5*x_im[D]
								+ cos2pi_5*x_re[E] + sin2pi_5*x_im[E];
				
				tA_im = x_im[A] + x_im[B] + x_im[C] + x_im[D] + x_im[E];
				tB_im = x_im[A] + cos2pi_5*x_im[B] - sin2pi_5*x_re[B]
								+ cos4pi_5*x_im[C] - sin4pi_5*x_re[C]
								+ cos4pi_5*x_im[D] + sin4pi_5*x_re[D]
								+ cos2pi_5*x_im[E] + sin2pi_5*x_re[E];
				tC_im = x_im[A] + cos4pi_5*x_im[B] - sin4pi_5*x_re[B]
								+ cos2pi_5*x_im[C] + sin2pi_5*x_re[C]
								+ cos2pi_5*x_im[D] - sin2pi_5*x_re[D]
								+ cos4pi_5*x_im[E] + sin4pi_5*x_re[E];
				tD_im = x_im[A] + cos4pi_5*x_im[B] + sin4pi_5*x_re[B]
								+ cos2pi_5*x_im[C] - sin2pi_5*x_re[C]
								+ cos2pi_5*x_im[D] + sin2pi_5*x_re[D]
								+ cos4pi_5*x_im[E] - sin4pi_5*x_re[E];
				x_im[E]=x_im[A] + cos2pi_5*x_im[B] + sin2pi_5*x_re[B]
								+ cos4pi_5*x_im[C] + sin4pi_5*x_re[C]
								+ cos4pi_5*x_im[D] - sin4pi_5*x_re[D]
								+ cos2pi_5*x_im[E] - sin2pi_5*x_re[E];
				#if DEBUG
				printf("SS. %f,%f,%f,%f,%f\n",x_re[A],x_re[B],x_re[C],x_re[D],x_re[E]);
				printf("SS. %f,%f,%f,%f,%f\n",x_im[A],x_im[B],x_im[C],x_im[D],x_im[E]);
				#endif
				x_re[A] = tA_re; x_re[B] = tB_re; x_re[C] = tC_re; x_re[D] = tD_re; x_re[E] = tE_re;
				x_im[A] = tA_im; x_im[B] = tB_im; x_im[C] = tC_im; x_im[D] = tD_im;
				#if DEBUG
				printf("EE. %f,%f,%f,%f,%f\n",x_re[A],x_re[B],x_re[C],x_re[D],x_re[E]);
				printf("EE. %f,%f,%f,%f,%f\n\n",x_im[A],x_im[B],x_im[C],x_im[D],x_im[E]);
				#endif
			}
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		m = m * order[i];
	}
	for(i=power_5;i<power_53;i++)
	{
		s *= order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			tw_re = w_re*w_re - w_im*w_im;
			tw_im = w_re*w_im + w_im*w_re;
			for(A=k;A<N;A+=order[i]*m)
			{
				B = A + m;
				C = A + 2*m;
				
				tA_re = x_re[A];
				tA_im = x_im[A];
				
				tB_re = w_re*x_re[B] - w_im*x_im[B];
				tB_im = w_re*x_im[B] + w_im*x_re[B];
				
				tC_re = tw_re*x_re[C] - tw_im*x_im[C];
				tC_im = tw_re*x_im[C] + tw_im*x_re[C];
				
				x_re[A] = tA_re + tB_re + tC_re;
				x_re[B] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tB_im-tC_im);
				x_re[C] = tA_re - 0.5*(tB_re+tC_re) + sqrt3_2*(tC_im-tB_im);
				
				x_im[A] = tA_im + tB_im + tC_im;
				x_im[B] = tA_im + sqrt3_2*(tC_re-tB_re) - 0.5*(tB_im+tC_im);
				x_im[C] = tA_im + sqrt3_2*(tB_re-tC_re) - 0.5*(tB_im+tC_im);
			}
			t    = w_re;
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		m *= order[i];
	}
	for(i=power_53;i<power_sum;i++)
	{
		s = s*order[i];
		w_re = 1.0;
		w_im = 0.0;
		w_N_re = cos(2.0*M_PI/s);
		w_N_im = -sin(2.0*M_PI/s);
		for(k=0;k<m;++k)
		{
			for(A=k;A<N;A+=order[i]*m)
			{
				B = A + m;
				
				tA_re = x_re[A];
				tA_im = x_im[A];
				
				tB_re = w_re*x_re[B] - w_im*x_im[B];
				tB_im = w_re*x_im[B] + w_im*x_re[B];
				
				x_re[A] = tA_re + tB_re;
				x_re[B] = tA_re - tB_re;
				
				x_im[A] = tA_im + tB_im;
				x_im[B] = tA_im - tB_im;
			}
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		m = m * order[i];
	}
	return 0;
}

int FFT_general_np_separated(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2 + power_3 + power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d, power_sum=%d\n", power_2, power_3, power_5, power_sum);
    #endif
    
    int i, *order;
    int power_53 = power_5+power_3;
    order = (int *) malloc(power_sum*sizeof(int));
    for (i=0;i<power_5;i++)				order[i]=5;
    for (i=power_5;i<power_53;i++)		order[i]=3;
    for (i=power_53;i<power_sum;i++)	order[i]=2;
    
    /* FFT */
	clock_t t1, t2, t3;
	t1 = clock();
    bit_reverse_np(x_re, x_im, N, order);
	t2 = clock();
	butterfly_np  (x_re, x_im, N, power_5, power_53, power_sum, order);
	t3 = clock();
	printf("time1 = %f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
	printf("time2 = %f\n",(t3-t2)/(double) CLOCKS_PER_SEC);
    /* FINISH */
    
    return 0;
}


/* other bitreverse */

int rearrange(double *x_re, double *x_im, int N)
{
    int copy_N = N;
    int power_2=0, power_3=0, power_5=0;
    int radix_2=1, radix_3=1, radix_5=1;
    while(copy_N%2==0){ power_2++; copy_N/=2; radix_2*=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; radix_3*=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; radix_5*=5; }
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d\n", power_2, power_3, power_5);
    printf("radix_2=%d, radix_3=%d, radix_5=%d\n", radix_2, radix_3, radix_5);
    #endif
    
    int i, j;
	int m,p,q,k;
	int radix_3_5 = radix_3*radix_5;
	double *out_re, *out_im;
	out_re = (double *) malloc(N*sizeof(double));
	out_im = (double *) malloc(N*sizeof(double));
    
    if (power_2>0 && N!=2)
    {
		m = radix_2/2;
	    q = 0;
	    for(p=0;p<radix_2;++p)
		{
			for(i=0;i<radix_3_5;i++)
			{
				#if DEBUG
		    	printf("%d -> %d\n", p*radix_3_5+i, q+radix_2*i);	// FIXME: do not compute everytime
			    #endif
		    	out_re[p*radix_3_5+i] = x_re[q+radix_2*i];
				out_im[p*radix_3_5+i] = x_im[q+radix_2*i];
		    }
	        k = m;
	        while(q >= k & k > 0)
	        {
	            q = q-k;				// 1->0
	            k = k/2;				// next
		        #if DEBUG
		        //printf("q=%d, k=%d\n",q ,k);
		        #endif
	        }
	        q = q+k;
	        #if DEBUG
	        //printf("q=%d\n",q);
	        #endif
		}
		for(i=0;i<N;i++)
		{
			x_re[i] = out_re[i];
			x_im[i] = out_im[i];
		}
	}
	
	
	if (power_3>0 && !(power_3==1 && power_5==0) )
    {
		m = radix_3/3;
	    q = 0;
	    for(j=0;j<N;j+=radix_3_5)
		{
			for(p=0;p<radix_3;++p)
			{
				for(i=0;i<radix_5;i++)
				{
					#if DEBUG
			    	printf("%d -> %d\n", p*radix_5+i+j, q+radix_3*i+j);
			    	#endif
			    	out_re[p*radix_5+i+j] = x_re[q+radix_3*i+j];
					out_im[p*radix_5+i+j] = x_im[q+radix_3*i+j];
			    }
		        k = m;
		        while(q >= 2*k & k > 0)
		        {
		            q = q-2*k;				// 2->0
		            k = k/3;				// next
			        #if DEBUG
			        //printf("q=%d, k=%d\n",q ,k);
			        #endif
		        }
		        q = q+k;
		        #if DEBUG
		        //printf("q=%d\n",q);
		        #endif
			}
		}
		for(i=0;i<N;i++)
		{
			x_re[i] = out_re[i];
			x_im[i] = out_im[i];
		}
	}
	
	if (power_5>1)
    {
		m = radix_5/5;
	    q = 0;
	    for(j=0;j<N;j+=radix_5)
		{
			for(p=0;p<radix_5;++p)
			{
				#if DEBUG
				printf("%d -> %d\n", p+j, q+j);
				#endif
		    	out_re[p+j] = x_re[q+j];
				out_im[p+j] = x_im[q+j];
				
		        k = m;
		        while(q >= 4*k & k > 0)
		        {
		            q = q-4*k;				// 2->0
		            k = k/5;				// next
			        #if DEBUG
			        //printf("q=%d, k=%d\n",q ,k);
			        #endif
		        }
		        q = q+k;
		        #if DEBUG
		        //printf("q=%d\n",q);
		        #endif
			}
		}
		for(i=0;i<N;i++)
		{
			x_re[i] = out_re[i];
			x_im[i] = out_im[i];
		}
	}
    return 0;
}

