#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(void)
{
	clock_t time_start, time_end;
	clock_t empty_loop;
	int i;
	int A = 9;
	double B = 9.0;

	// Empty
	time_start = clock();
	for (i = 0; i < 1000000000; i++){}
	time_end = clock();
	printf("Empty_loop:\n\t%.4f (seconds)\n", (time_end - time_start) / (double)(CLOCKS_PER_SEC));
	printf("\n");
	empty_loop = time_end - time_start;

	// Addition
	time_start = clock();
	for (i = 0; i < 1000000000; i++) { A = A + 9; }
	time_end = clock();
	printf("Addition_loop:\n\t%.4f (seconds)\n", (time_end - time_start) / (double)(CLOCKS_PER_SEC));
	printf("One_addition_time:\n\t%e (seconds)\n", (time_end - time_start - empty_loop) / (double)(CLOCKS_PER_SEC) / 1000000000);
	printf("\n");

	// Multiplication
	time_start = clock();
	for (i = 0; i < 1000000000; i++) { A = A * 9; }
	time_end = clock();
	printf("Multiplication_loop:\n\t%.4f (seconds)\n", (time_end - time_start) / (double)(CLOCKS_PER_SEC));
	printf("One_multiplication_time:\n\t%e (seconds)\n", (time_end - time_start - empty_loop) / (double)(CLOCKS_PER_SEC) / 1000000000);
	printf("\n");

	// Sine
	time_start = clock();
	for (i = 0; i < 1000000000; i++) { B = sin(B); }
	time_end = clock();
	printf("Sine_loop:\n\t%.4f (seconds)\n", (time_end - time_start) / (double)(CLOCKS_PER_SEC));
	printf("One_sine_time:\n\t%e (seconds)\n", (time_end - time_start - empty_loop) / (double)(CLOCKS_PER_SEC) / 1000000000);
	printf("\n");

	printf("A= %d \n", A);
	printf("B= %f \n", B);
	system("pause");
	return 0;
}