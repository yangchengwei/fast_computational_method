#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#define DEBUG 1

int main()
{
	int i;
	int N = (int)(pow(2,2)*pow(3,1)*pow(5,0));
	double x_re[N], x_im[N];
	
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	rearrange(x_re, x_im, N);
	/*butterfly(x_re, x_im, N);*/
	
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	return;
}

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
	double out_re[N], out_im[N];
    
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

int butterfly(double *x_re, double *x_im, int N)
{
	int k, p, q, m;
	double w_re, w_im, w_N_re, w_N_im, t; 
	m = 1;
	while(m<N)
	{
		w_re = 1.0;
		w_im = 0.0;
		w_N_re =  cos(M_PI/m);
		w_N_im = -sin(M_PI/m);
		for(k=0;k<m;++k)
		{
			for(p=k;p<N;p+=2*m)
			{
				q = p + m;
				// printf("(%d,%d) (%f,%f) \n", p,q, w_re, w_im);
				// multiply (w_re + w_im * i) on x[q]
				t = x_re[q]; 
				x_re[q] = w_re*x_re[q] - w_im*x_im[q];
				x_im[q] = w_re*x_im[q] + w_im*t; 
				
				t = x_re[p];
				x_re[p] = x_re[p] + x_re[q];
				x_re[q] = t       - x_re[q]; 
				t = x_im[p];
				x_im[p] = x_im[p] + x_im[q];
				x_im[q] = t       - x_im[q];
			}
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		m = m * 2;
	}
	
	return;
}
