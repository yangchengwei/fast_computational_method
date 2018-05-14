#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#define DEBUG 0

int main()
{
	int i;
	int N = (int)(pow(2,0)*pow(3,3)*pow(5,0));
	int times=1;
	double x_re[N], x_im[N];
	clock_t t1, t2;
	
	///////////////////////// bit_reverse(x_re, x_im, N);
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	for(i=0;i<times;++i)
	bit_reverse(x_re, x_im, N);
	t2 = clock();
	butterfly(x_re, x_im, N);
	printf("bit_reverse \t time=%f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
	#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif
	
	///////////////////////// bit_reverse2(x_re, x_im, N);
	/*for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	t1 = clock();
	for(i=0;i<times;++i)
	bit_reverse2(x_re, x_im, N);
	t2 = clock();
	// butterfly(x_re, x_im, N);
	printf("bit_reverse2 \t time=%f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
	#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif*/
	
	///////////////////////// rearrange(x_re, x_im, N);
	/*for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	t1 = clock();
	for(i=0;i<times;++i)
	rearrange(x_re, x_im, N);
	t2 = clock();
	printf("rearrange \t time=%f\n",(t2-t1)/(double) CLOCKS_PER_SEC);
    #if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	#endif*/
	
	return;
}

int bit_reverse(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2+power_3+power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d\n", power_2, power_3, power_5);
    #endif
    
    int i;
    int p,q,temp;
    int step,gate,add;
	double t;
    double copy_re[N], copy_im[N];
	for(i=0;i<N;i++)
	{
		copy_re[i] = x_re[i];
		copy_im[i] = x_im[i];
	}
	
    int value[N];
    for(i=0;i<N;i++) value[i]=i;
    
    int order[power_sum];
    i=0;
    while(i<power_5) order[i++]=5;
    while(i<power_5+power_3) order[i++]=3;
    while(i<power_sum) order[i++]=2;
    
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
        while(q >= gate & gate > 0)
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

int bit_reverse2(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2+power_3+power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d\n", power_2, power_3, power_5);
    #endif
    
    int i;
    int p,q,temp;
    int step,gate,add;
	double t;
    
    int value[N];
    for(i=0;i<N;i++) value[i]=i;
    
    int order[power_sum];
    i=0;
    while(i<power_5) order[i++]=5;
    while(i<power_5+power_3) order[i++]=3;
    while(i<power_sum) order[i++]=2;
    
    step = N/order[0];
    q = step;			// first change
    for(p=1;p<N-1;p++)
    {
		// change value
    	//printf("%d -> %d\n", p,q);
    	if (value[q]!=p)
    	{
    		i=q;
	    	if (value[q]!=q)
			{
				i=value[q];
			}
			// printf("(p,q,i)=(%d,%d,%d)\n",p,q,i);
			// printf("value[p,q,i]=(%d,%d,%d)\n",value[p],value[q],value[i]);
			t = x_re[p];
	        x_re[p] = x_re[i];
			x_re[i] = t;
	        t = x_im[p];
	        x_im[p] = x_im[i];
			x_im[i] = t;
			temp = value[p];
			value[p] = value[q];
			value[q] = temp;
			// printf("(p,q,i)=(%d,%d,%d)\n",p,q,i);
			// printf("value[p,q,i]=(%d,%d,%d)\n",value[p],value[q],value[i]);
			// for (i=0;i<N;i++) printf("%d ",value[i]);
			// printf("\n\n");
		}
        
        // compute next place
        i = 0;
        add = step;
        gate = (order[i++]-1)*add;
        while(q >= gate & gate > 0)
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
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
    while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }
    power_sum = power_2+power_3+power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d\n", power_2, power_3, power_5);
    #endif
    
    int i, j, order[power_sum];
    i=0;
    while(i<power_5) order[i++]=5;
    while(i<power_5+power_3) order[i++]=3;
    while(i<power_sum) order[i++]=2;
    
	int k, m, s;
	int A, B, C, D, E;
	double w_re, w_im, w_N_re, w_N_im, t;
	double tA_re, tA_im, tB_re, tB_im, tC_re, tC_im;
	double tA, tB, tC, tD, tE;
	m = 1;
	s = 1;
	for(i=0;i<power_sum;i++)
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
				switch (order[i]){
					case 2:
						B = A + m;
						// printf("(%d,%d) (%f,%f) %d \n", p, q, w_re, w_im, k);
						// multiply (w_re + w_im * i) on x[q]
						t = x_re[B]; // outside loop
						x_re[B] = w_re*x_re[B] - w_im*x_im[B];
						x_im[B] = w_re*x_im[B] + w_im*t;
						
						t = x_re[A]; // butterfly
						x_re[A] = x_re[A] + x_re[B];
						x_re[B] = t       - x_re[B]; 
						t = x_im[A];
						x_im[A] = x_im[A] + x_im[B];
						x_im[B] = t       - x_im[B];
						break;
			    	case 3:
						B = A + m;
						C = A + 2*m;
						tB_re = w_re*x_re[B] - w_im*x_im[B];
						tB_im = w_re*x_im[B] + w_im*x_re[B];
						x_re[B] = tB_re; x_im[B] = tB_im;
						
						tC_re = w_re*(w_re*x_re[C] - w_im*x_im[C]) - w_im*(w_re*x_im[C] + w_im*x_re[C]);
						tC_im = w_re*(w_re*x_im[C] + w_im*x_re[C]) + w_im*(w_re*x_re[C] - w_im*x_im[C]);
						x_re[C] = tC_re; x_im[C] = tC_im;
						
						tA_re = x_re[A] + x_re[B] + x_re[C];
						tB_re = x_re[A] + ((-1.0/2.0)*x_re[B]-(-sqrt(3)/2.0)*x_im[B]) + ((-1.0/2.0)*x_re[C]-(sqrt(3)/2.0)*x_im[C]);
						tC_re = x_re[A] + ((-1.0/2.0)*x_re[B]-(sqrt(3)/2.0)*x_im[B]) + ((-1.0/2.0)*x_re[C]-(-sqrt(3)/2.0)*x_im[C]);
						
						tA_im = x_im[A] + x_im[B] + x_im[C];
						tB_im = x_im[A] + ((-sqrt(3)/2.0)*x_re[B]+(-1.0/2.0)*x_im[B]) + ((sqrt(3)/2.0)*x_re[C]+(-1.0/2.0)*x_im[C]);
						tC_im = x_im[A] + ((sqrt(3)/2.0)*x_re[B]+(-1.0/2.0)*x_im[B]) + ((-sqrt(3)/2.0)*x_re[C]+(-1.0/2.0)*x_im[C]);
						
						#if DEBUG
						printf("SS. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("SS. %f,%f,%f\n",x_im[A],x_im[B],x_im[C]);
						#endif
						x_re[A] = tA_re; x_re[B] = tB_re; x_re[C] = tC_re;
						x_im[A] = tA_im; x_im[B] = tB_im; x_im[C] = tC_im;
						#if DEBUG
						printf("EE. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("EE. %f,%f,%f\n\n",x_im[A],x_im[B],x_im[C]);
						#endif
						break;
				    case 5:
						break;
				    default:
						printf("error\n");
				}
			}
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		m = m * order[i];
	}
	
	return;
}
