
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define DEBUG 0


/* bit_reverse */

int bit_reverse(double *x_re, double *x_im, int N, int power_2, int power_3, int power_5, int power_sum, int *order)
{
    int i,p,q;
    int step,gate,add;
	double t, *copy_re, *copy_im;
	copy_re = (double *) malloc(N*sizeof(double));
	copy_im = (double *) malloc(N*sizeof(double));
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

int bit_reverse0(double *x_re, double *x_im, int N)
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
    
	double copy_re[N], copy_im[N];
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


/* butterfly */

int butterfly(double *x_re, double *x_im, int N, int power_2, int power_3, int power_5, int power_sum, int *order)
{
    int i, j, k, m, s;
	int A, B, C, D, E;
	double w_re, w_im, w_N_re, w_N_im, t;
	double tA_re, tB_re, tC_re, tD_re, tE_re;
	double tA_im, tB_im, tC_im, tD_im, tE_im;
	double tA, tB, tC, tD, tE;
	double tw_re, tw_im;
	
	double sqrt3_2 = sqrt(3)/2.0;
	double cos2pi_5 = cos(2.0*M_PI/5.0);
	double sin2pi_5 = sin(2.0*M_PI/5.0);
	double cos4pi_5 = cos(4.0*M_PI/5.0);
	double sin4pi_5 = sin(4.0*M_PI/5.0);
	
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
						
						t = x_re[B];
						x_re[B] = w_re*x_re[B] - w_im*x_im[B];
						x_im[B] = w_re*x_im[B] + w_im*t;
						
						t = x_re[A];
						x_re[A] = x_re[A] + x_re[B];
						x_re[B] = t       - x_re[B]; 
						t = x_im[A];
						x_im[A] = x_im[A] + x_im[B];
						x_im[B] = t       - x_im[B];
						break;
						
			    	case 3:
						B = A + m;
						C = A + 2*m;
						
						t = x_re[B];
						x_re[B] = w_re*x_re[B] - w_im*x_im[B];
						x_im[B] = w_re*x_im[B] + w_im*t;
						
						t = x_re[C];
						x_re[C] = (w_re*w_re - w_im*w_im)*x_re[C] - (w_re*w_im + w_im*w_re)*x_im[C];
						x_im[C] = (w_re*w_re - w_im*w_im)*x_im[C] + (w_re*w_im + w_im*w_re)*t;
						
						tA_re = x_re[A] + x_re[B] + x_re[C];
						tB_re = x_re[A] - 0.5*x_re[B] + sqrt3_2*x_im[B]
										- 0.5*x_re[C] - sqrt3_2*x_im[C];
						tC_re = x_re[A] - 0.5*x_re[B] - sqrt3_2*x_im[B]
										- 0.5*x_re[C] + sqrt3_2*x_im[C];
						
						tA_im = x_im[A] + x_im[B] + x_im[C];
						tB_im = x_im[A] - sqrt3_2*x_re[B] - 0.5*x_im[B]
										+ sqrt3_2*x_re[C] - 0.5*x_im[C];
						x_im[C]=x_im[A] + sqrt3_2*x_re[B] - 0.5*x_im[B]
										- sqrt3_2*x_re[C] - 0.5*x_im[C];
						
						#if DEBUG
						printf("SS. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("SS. %f,%f,%f\n",x_im[A],x_im[B],x_im[C]);
						#endif
						x_re[A] = tA_re; x_re[B] = tB_re; x_re[C] = tC_re;
						x_im[A] = tA_im; x_im[B] = tB_im;
						#if DEBUG
						printf("EE. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("EE. %f,%f,%f\n\n",x_im[A],x_im[B],x_im[C]);
						#endif
						break;
						
				    case 5:
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
	
	return 0;
}

int butterfly0(double *x_re, double *x_im, int N)
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
	double tA_re, tB_re, tC_re, tD_re, tE_re;
	double tA_im, tB_im, tC_im, tD_im, tE_im;
	double tA, tB, tC, tD, tE;
	double tw_re, tw_im;
	
	double sqrt3_2 = sqrt(3)/2.0;
	double cos2pi_5 = cos(2.0*M_PI/5.0);
	double sin2pi_5 = sin(2.0*M_PI/5.0);
	double cos4pi_5 = cos(4.0*M_PI/5.0);
	double sin4pi_5 = sin(4.0*M_PI/5.0);
	
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
						
						t = x_re[B];
						x_re[B] = w_re*x_re[B] - w_im*x_im[B];
						x_im[B] = w_re*x_im[B] + w_im*t;
						
						t = x_re[A];
						x_re[A] = x_re[A] + x_re[B];
						x_re[B] = t       - x_re[B]; 
						t = x_im[A];
						x_im[A] = x_im[A] + x_im[B];
						x_im[B] = t       - x_im[B];
						break;
						
			    	case 3:
						B = A + m;
						C = A + 2*m;
						
						t = x_re[B];
						x_re[B] = w_re*x_re[B] - w_im*x_im[B];
						x_im[B] = w_re*x_im[B] + w_im*t;
						
						t = x_re[C];
						x_re[C] = (w_re*w_re - w_im*w_im)*x_re[C] - (w_re*w_im + w_im*w_re)*x_im[C];
						x_im[C] = (w_re*w_re - w_im*w_im)*x_im[C] + (w_re*w_im + w_im*w_re)*t;
						
						tA_re = x_re[A] + x_re[B] + x_re[C];
						tB_re = x_re[A] - 0.5*x_re[B] + sqrt3_2*x_im[B]
										- 0.5*x_re[C] - sqrt3_2*x_im[C];
						tC_re = x_re[A] - 0.5*x_re[B] - sqrt3_2*x_im[B]
										- 0.5*x_re[C] + sqrt3_2*x_im[C];
						
						tA_im = x_im[A] + x_im[B] + x_im[C];
						tB_im = x_im[A] - sqrt3_2*x_re[B] - 0.5*x_im[B]
										+ sqrt3_2*x_re[C] - 0.5*x_im[C];
						x_im[C]=x_im[A] + sqrt3_2*x_re[B] - 0.5*x_im[B]
										- sqrt3_2*x_re[C] - 0.5*x_im[C];
						
						#if DEBUG
						printf("SS. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("SS. %f,%f,%f\n",x_im[A],x_im[B],x_im[C]);
						#endif
						x_re[A] = tA_re; x_re[B] = tB_re; x_re[C] = tC_re;
						x_im[A] = tA_im; x_im[B] = tB_im;
						#if DEBUG
						printf("EE. %f,%f,%f\n",x_re[A],x_re[B],x_re[C]);
						printf("EE. %f,%f,%f\n\n",x_im[A],x_im[B],x_im[C]);
						#endif
						break;
						
				    case 5:
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
	
	return 0;
}

int butterfly2(double *x_re, double *x_im, int N)
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
	double tD_re, tD_im, tE_re, tE_im;
	double tA, tB, tC, tD, tE;
	double tw_re, tw_im;
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
						B = A + m;
						C = A + 2*m;
						D = A + 3*m;
						E = A + 4*m;
						
						tB_re = w_re*x_re[B] - w_im*x_im[B];
						tB_im = w_re*x_im[B] + w_im*x_re[B];
						x_re[B] = tB_re; x_im[B] = tB_im;
						
						tw_re = w_re*w_re - w_im*w_im;
						tw_im = w_re*w_im + w_im*w_re;
						
						tC_re = tw_re*x_re[C] - tw_im*x_im[C];
						tC_im = tw_re*x_im[C] + tw_im*x_re[C];
						x_re[C] = tC_re; x_im[C] = tC_im;
						
						t = tw_re;
						tw_re = w_re*tw_re - w_im*tw_im;
						tw_im = w_re*tw_im + w_im*t;
						
						tD_re = tw_re*x_re[D] - tw_im*x_im[D];
						tD_im = tw_re*x_im[D] + tw_im*x_re[D];
						x_re[D] = tD_re; x_im[D] = tD_im;
						
						t = tw_re;
						tw_re = w_re*tw_re - w_im*tw_im;
						tw_im = w_re*tw_im + w_im*t;
						
						tE_re = tw_re*x_re[E] - tw_im*x_im[E];
						tE_im = tw_re*x_im[E] + tw_im*x_re[E];
						x_re[E] = tE_re; x_im[E] = tE_im;
						
						tA_re = x_re[A] + x_re[B] + x_re[C] + x_re[D] + x_re[E];
						tB_re = x_re[A] + (cos(2.0*M_PI/5)*x_re[B])-(-sin(2.0*M_PI/5)*x_im[B])
										+ (cos(4.0*M_PI/5)*x_re[C])-(-sin(4.0*M_PI/5)*x_im[C])
										+ (cos(6.0*M_PI/5)*x_re[D])-(-sin(6.0*M_PI/5)*x_im[D])
										+ (cos(8.0*M_PI/5)*x_re[E])-(-sin(8.0*M_PI/5)*x_im[E]);
						tC_re = x_re[A] + (cos(4.0*M_PI/5)*x_re[B])-(-sin(4.0*M_PI/5)*x_im[B])
										+ (cos(8.0*M_PI/5)*x_re[C])-(-sin(8.0*M_PI/5)*x_im[C])
										+ (cos(2.0*M_PI/5)*x_re[D])-(-sin(2.0*M_PI/5)*x_im[D])
										+ (cos(6.0*M_PI/5)*x_re[E])-(-sin(6.0*M_PI/5)*x_im[E]);
						tD_re = x_re[A] + (cos(6.0*M_PI/5)*x_re[B])-(-sin(6.0*M_PI/5)*x_im[B])
										+ (cos(2.0*M_PI/5)*x_re[C])-(-sin(2.0*M_PI/5)*x_im[C])
										+ (cos(8.0*M_PI/5)*x_re[D])-(-sin(8.0*M_PI/5)*x_im[D])
										+ (cos(4.0*M_PI/5)*x_re[E])-(-sin(4.0*M_PI/5)*x_im[E]);
						tE_re = x_re[A] + (cos(8.0*M_PI/5)*x_re[B])-(-sin(8.0*M_PI/5)*x_im[B])
										+ (cos(6.0*M_PI/5)*x_re[C])-(-sin(6.0*M_PI/5)*x_im[C])
										+ (cos(4.0*M_PI/5)*x_re[D])-(-sin(4.0*M_PI/5)*x_im[D])
										+ (cos(2.0*M_PI/5)*x_re[E])-(-sin(2.0*M_PI/5)*x_im[E]);
						
						tA_im = x_im[A] + x_im[B] + x_im[C] + x_im[D] + x_im[E];
						tB_im = x_im[A] + (cos(2.0*M_PI/5)*x_im[B])+(-sin(2.0*M_PI/5)*x_re[B])
										+ (cos(4.0*M_PI/5)*x_im[C])+(-sin(4.0*M_PI/5)*x_re[C])
										+ (cos(6.0*M_PI/5)*x_im[D])+(-sin(6.0*M_PI/5)*x_re[D])
										+ (cos(8.0*M_PI/5)*x_im[E])+(-sin(8.0*M_PI/5)*x_re[E]);
						tC_im = x_im[A] + (cos(4.0*M_PI/5)*x_im[B])+(-sin(4.0*M_PI/5)*x_re[B])
										+ (cos(8.0*M_PI/5)*x_im[C])+(-sin(8.0*M_PI/5)*x_re[C])
										+ (cos(2.0*M_PI/5)*x_im[D])+(-sin(2.0*M_PI/5)*x_re[D])
										+ (cos(6.0*M_PI/5)*x_im[E])+(-sin(6.0*M_PI/5)*x_re[E]);
						tD_im = x_im[A] + (cos(6.0*M_PI/5)*x_im[B])+(-sin(6.0*M_PI/5)*x_re[B])
										+ (cos(2.0*M_PI/5)*x_im[C])+(-sin(2.0*M_PI/5)*x_re[C])
										+ (cos(8.0*M_PI/5)*x_im[D])+(-sin(8.0*M_PI/5)*x_re[D])
										+ (cos(4.0*M_PI/5)*x_im[E])+(-sin(4.0*M_PI/5)*x_re[E]);
						tE_im = x_im[A] + (cos(8.0*M_PI/5)*x_im[B])+(-sin(8.0*M_PI/5)*x_re[B])
										+ (cos(6.0*M_PI/5)*x_im[C])+(-sin(6.0*M_PI/5)*x_re[C])
										+ (cos(4.0*M_PI/5)*x_im[D])+(-sin(4.0*M_PI/5)*x_re[D])
										+ (cos(2.0*M_PI/5)*x_im[E])+(-sin(2.0*M_PI/5)*x_re[E]);
						#if DEBUG
						printf("SS. %f,%f,%f,%f,%f\n",x_re[A],x_re[B],x_re[C],x_re[D],x_re[E]);
						printf("SS. %f,%f,%f,%f,%f\n",x_im[A],x_im[B],x_im[C],x_im[D],x_im[E]);
						#endif
						x_re[A] = tA_re; x_re[B] = tB_re; x_re[C] = tC_re; x_re[D] = tD_re; x_re[E] = tE_re;
						x_im[A] = tA_im; x_im[B] = tB_im; x_im[C] = tC_im; x_im[D] = tD_im; x_im[E] = tE_im;
						#if DEBUG
						printf("EE. %f,%f,%f,%f,%f\n",x_re[A],x_re[B],x_re[C],x_re[D],x_re[E]);
						printf("EE. %f,%f,%f,%f,%f\n\n",x_im[A],x_im[B],x_im[C],x_im[D],x_im[E]);
						#endif
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
	
	return 0;
}


/* Fast Fourier Transform */

int FFT_general(double *x_re, double *x_im, int N)
{
    int copy_N = N, power_sum = 0;
    int power_2=0, power_3=0, power_5=0;
	while(copy_N>1) /* slower */
	{
    	if		(copy_N%2==0) { power_2++; copy_N/=2; }
    	else if (copy_N%3==0) { power_3++; copy_N/=3; }
    	else   /*copy_N%5==0*/{ power_5++; copy_N/=5; }
	}
    /*while(copy_N%2==0){ power_2++; copy_N/=2; }
    while(copy_N%3==0){ power_3++; copy_N/=3; }
    while(copy_N%5==0){ power_5++; copy_N/=5; }*/ 
    power_sum = power_2 + power_3 + power_5;
    #if DEBUG
    printf("power_2=%d, power_3=%d, power_5=%d, power_sum=%d\n", power_2, power_3, power_5, power_sum);
    #endif
    
    int i, order[power_sum];
    int t = power_5+power_3;
    for (i=0;i<power_5;i++)		order[i]=5;
    for (i=power_5;i<t;i++)		order[i]=3;
    for (i=t;i<power_sum;i++)	order[i]=2;
    
    /* FFT */
    bit_reverse(x_re, x_im, N, power_2, power_3, power_5, power_sum, order);
    butterfly  (x_re, x_im, N, power_2, power_3, power_5, power_sum, order);
    /* FINISH */
    
    return 0;
}
