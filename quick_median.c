#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 0

int quick_find(int *x, int left, int right, int goal);
int quicksort(int *x, int left, int right);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, median;				// array for sorting & median
	double T1;					// final computing time
	int i, N;

	srand( time(NULL) );

	for(N=160000;N<=10240000;N*=2)
	{
		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		//#pragma omp parallel
		{
			//#pragma omp parallel for
			for(i=0;i<N;++i)
			{
				y[i] = x[i] = rand() % N;
			}
		}
		
		#if DEBUG        // if DEBUG == 1, then compile the following codes 
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif			// end of if block
		
		t1 = clock();
		median = quick_find(x, 0, N, N/2);
		
		#if DEBUG
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif
		
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("median= %d ; Quick Median  %d elements: %f\n",median, N, T1);

		// by quick sorting
		for(i=0;i<N;++i) y[i] = x[i];
		
		t1 = clock();
		quicksort(x, 0, N);
		median = x[N/2];
		
		#if DEBUG
		for(i=0;i<N;++i)
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif
		
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("median= %d ; Quick Sorting %d elements: %f\n\n",median, N, T1);

		free(x);
		free(y);
	} 
	return 0;
}

int quick_find(int *x, int left, int right, int goal)
{
	int i, j, k;
	int pivot, t;
	
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
    	while(1)
		{
      		while(i < right && pivot >= x[i]) i++;
      		while(j >  left && pivot <  x[j]) j--;
      		if(i>=j) break;
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
        }
        x[left] = x[j];
        x[j] = pivot;
        #if DEBUG
        printf("i= %d ,j= %d ,goal= %d ,pivot= %d\n",i,j,goal,pivot);
		for(k=left;k<right;++k)
		{
			printf("x[%d]=%d\n",k,x[k]);
		}
		system("pause");
        #endif
        
        if (j==goal)		return x[goal];
		else if (j>goal)	return quick_find(x, left, j, goal);
		else				return quick_find(x, j+1, right, goal);
    }
    else
    {
    	return x[goal];
	}
}

int quicksort(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t;
	
	//  
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
    	// 
    	while(1)
		{	// x: 5(pivot) 4 8 3 2 7 3 2 10  -> i = 2, j = 7 (ユ传 x[i], x[j]) 
			//    5(pivot) 4 2 3 2 7 3 8 10  -> i = 5, j = 6 (ユ传 x[i], x[j])
			//    5(pivot) 4 2 3 2 3 7 8 10  -> i = 6, j = 5 (ぃユ传!!!) 
			//    3                5 (location: j) 
			// x: 4(pivot) 4 8 3 2 7 3 4 10  -> i = 2, j = 7 (ユ传 x[i], x[j]) 
			//    4(pivot) 4 4 3 2 7 3 8 10  -> i = 5, j = 6 (ユ传 x[i], x[j])
			//    4(pivot) 4 4 3 2 3 7 8 10  -> i = 6, j = 5 (ぃユ传!!!) 
			//    3                4 (location: j) 			
			// x: 8 2 1 8 7 8 9 4 5 8      i = 3, j = 9 (ユ传 x[i], x[j])
			// x: 8 2 1 8 7 8 9 4 5 8
      		while(i < right && pivot >= x[i]) i++; // ┕娩т材  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // ┕オ娩т材  pivot >= x[j] 
      		#if DEBUG
			printf("%d %d %d\n", i,j,pivot);
			#endif
      		if(i>=j) break;
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		#if DEBUG
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
        }
        //t = x[left];
        x[left] = x[j];
        x[j] = pivot;
        #if DEBUG
        printf("i=%d,j=%d\n",i,j);
		for(k=left;k<right;++k)
		{
			printf("x[%d]=%d\n",k,x[k]);
		}
		system("pause");
        #endif
		quicksort(x, left, j);
		quicksort(x, j+1, right);
    }
    else 
    {
    	return 1;
	}
}
