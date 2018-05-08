#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int main()
{
	int i;
	int N = 8;
	double y_re[N], y_im[N], x_re[N], x_im[N];
	
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	bit_reverse(x_re, x_im, N);	
	butterfly(x_re, x_im, N);
	
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	
	return;
}

int bit_reverse(double *x_re, double *x_im, int N)
{
    int m,p,q,k;
    int N235=N;
    int N2=0, N3=0, N5=0;
    double t;
    
    while(N235%2==0){ N235=N235/2; N2++;}
	while(N235%3==0){ N235=N235/3; N3++;}
	while(N235%5==0){ N235=N235/5; N5++;}
	printf("N2=%d,N3=%d,N5=%d \n",N2,N3,N5);
    
    m = pow(2,N2)/2;				// first bit after bit reverse
    q = m;							// first change
    for(p=1;p<N-1;++p)				// from 1 to N-2 (0 and N-1 don't need to change) 
    {
        printf("%d <-> %d\n", p,q);
        if(p < q)
        {
            t = x_re[p];
            x_re[p] = x_re[q];
			x_re[q] = t;
            t = x_im[p];
            x_im[p] = x_im[q];
			x_im[q] = t;			 
        }
        
        k = m;						// first bit after bit reverse
        
        while(q >= k & k > 0)		// if == 1
        {
            q = q-k;				// 1->0
            k = k/2;				// next bit
        }
        q = q+k;					// 0->1
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
				printf("(%d,%d) (%f,%f) \n", p,q, w_re, w_im);
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
