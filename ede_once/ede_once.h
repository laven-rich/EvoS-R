#pragma once
#include <corecrt_math.h>
#include <time.h>
#include <ctime>
#include <malloc.h>
#include <float.h>
#include <queue>  
#include <iostream>

#define URAND  ((double)rand()/((double)RAND_MAX+1.0))  
#define PI 3.141592654
#define LOW_VAL -10000

const double e = exp(1.0);

double **data0;
using namespace std;

//to load and arrange the data
double gaussrand(double U, double V);
double generateGaussianNoise(double mu, double sigma);
double **loadData(int choose, int *d, int *n);
void malloc1D(double *&a, int D);
void malloc1E(int *&a, int D);
void malloc2D(int **&a, int xDim, int yDim);
void malloc2E(double **&a, int xDim, int yDim);
double min_fun(double *a, int b);
double max_fun(double *a, int b);

//functions to implement the EvoS-R algorithm
double getDistance(double *avector, double *bvector, int n);
double getDistance2(double *avector, double *bvector, int n);
double **lhs(int choose, int n);
int de(FILE *fp, FILE *fp2, double *p1, double *p2, double *p3, int *p4, int *valids, int N, int D, int Gmax);
void qsort(double *a, int L, int R);
void quicksort(double *v, int N);

double generateGaussianNoise(double mu, double sigma)
{
	const double epsilon = DBL_MIN;
	const double two_pi = 2.0*3.14159265358979323846;
	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	}while (u1 <= epsilon);
	
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

/*****************************************************************
*min
*****************************************************************/
double min_fun(double *a, int b) {
	double min = a[0];
	for (int i = 1; i < b; i++) {
		if (a[i] < min)
			min = a[i];
	}

	return min;
}

/*****************************************************************
*max
*****************************************************************/
double max_fun(double *a, int b) {
	double max = a[0];
	for (int i = 1; i < b; i++) {
		if (a[i] > max)
			max = a[i];
	}

	return max;
}

/*****************************************************************
*sum
*****************************************************************/
double sum_fun(double *a, int b) {
	double sum = 0;
	for (int i = 0; i < b; i++) {
		sum += a[i];
	}

	return sum;
}

/*****************************************************************
*mean
*****************************************************************/
double sum_mean(double *a, int b) {
	double sum = 0;
	for (int i = 0; i < b; i++) {
		sum += a[i];
	}
	sum = sum / (b * 1.0);

	return sum;
}

/*****************************************************************
*stand deviation
*****************************************************************/
double sum_std(double *a, int b) {
	double avg = 0, sum = 0, standard;
	for (int i = 0; i < b; i++) {
		avg += a[i];
	}
	avg = avg / (b * 1.0);  

	for (int i = 0; i < b; i++) {
		sum += pow(a[i]-avg, 2);
	}
	sum = sum / (b*1.0);      
	standard = pow(sum, 0.5); 

	return standard;
}


/*****************************************************************
*nk2_func
*input£º cluster£¬data p_pars£¬population individual p_index£¬division in_clu3£¬num. of data N£¬dimension D
*output£ºfitness value
*****************************************************************/
double func_nk2(double k_num, int *in_cluster, double *a1, double *a2, double *a3, int *relation, int N, int D, int nm) {
	int i, j;
	double sum = 0;

	for (i = 0; i < N; i++) {
		//01.noise
		if (*(in_cluster + i) == 0) {
			sum = sum + *(a1 + i * N + *(relation + i * 2 + 0)) + *(a1 + i * N + *(relation + i * 2 + 1));
			continue;
		}
		//02.first 
		if (*(in_cluster + i) == *(in_cluster + *(relation + i * 2 + 0))) {
			sum = sum+ *(a2 + i * N + *(relation + i * 2 + 0));
		}
		else {
			sum = sum + *(a3 + i * N + *(relation + i * 2 + 0));
		}
		//03.second
		if (*(in_cluster + i) == *(in_cluster + *(relation + i * 2 + 1))) {
			sum = sum + *(a2 + i * N + *(relation + i * 2 + 1));
		}
		else {
			sum = sum + *(a3 + i * N + *(relation + i * 2 + 1));
		}
	}

	sum = sum/(N*2.0);

	return sum;
}

/*****************************************************************
*function of allocation
*****************************************************************/
void malloc1D(double *&a, int D) {
	a = (double *)malloc(D * sizeof(double));
	if (a == NULL)
		perror("malloc");
}

void malloc1E(int *&a, int D) {
	a = (int *)malloc(D * sizeof(int));
	if (a == NULL)
		perror("malloc");
}

void malloc2D(int **&a, int xDim, int yDim)
{
	a = (int **)malloc(xDim * sizeof(int *));
	a[0] = (int *)malloc(xDim * yDim * sizeof(int));
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;
	}
	if (a == NULL)
		perror("malloc");
}

void malloc2E(double **&a, int xDim, int yDim)
{
	a = (double **)malloc(xDim * sizeof(double *));
	a[0] = (double *)malloc(xDim * yDim * sizeof(double));
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;
	}
	if (a == NULL)
		perror("malloc");
}

/*****************************************************************
*function of distance calculation
*****************************************************************/
double getDistance(double *avector, double *bvector, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i<n; i++)
		sum = sum + pow(*(avector + i) - *(bvector + i), 2);

	return sqrt(sum);
}

double getDistance2(double *avector, double *bvector, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i<n; i++)
		sum = sum + pow(*(avector + i) - *(bvector + i), 2);

	return sum;
}

/*****************************************************************
*permutation function
*****************************************************************/
void permutate(int *sequence, int N)
{
	int i, itr = N;       //N denote the numbers of the data
	for (i = 0; i<N; i++) {
		sequence[i] = i;
	}

	for (i = 0; i<itr; i++) {
		int a = rand() % N, b = rand() % N;
		int tmp = sequence[a];
		sequence[a] = sequence[b];
		sequence[b] = tmp;
	}
}

/******************************************************************************\
*								 Quicksort: qsort							 *
\******************************************************************************/
void qsort(double *a, int L, int R)
{
	int i, j;
	double x, w;

	i = L;
	j = R;
	x = a[(L + R) / 2];
	do {
		while (a[i] < x) {
			i = i + 1;
		}
		while (x < a[j]) {
			j = j - 1;
		}
		if (i <= j) {
			w = a[i];
			a[i] = a[j];
			a[j] = w;
			i = i + 1;
			j = j - 1;
		}
	} while (i <= j);
	if (L < j) {
		qsort(a, L, j);
	}
	if (i < R) {
		qsort(a, i, R);
	}
	return;
}

/******************************************************************************\
*								 Quicksort: qsort							 *
\******************************************************************************/
void quicksort(double *v, int N) {

	qsort(v, 0, N - 1);

}

/*****************************************************************
*latin hyper_lhs
*****************************************************************/
double **lhs(int choose, int n)
{
	int i, j;
	double **arraydata2;  
	FILE *fp=NULL;
	//adjust the path according to your pc
	if (choose == 1) {
		if ((fp = fopen("C:/../sample/d1.txt", "r")) == NULL)
			fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 2) {
		if ((fp = fopen("C:/../sample/d2.txt", "r")) == NULL)
			fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 3) {
		if ((fp = fopen("C:/../sample/d3.txt", "r")) == NULL)
			fprintf(stderr, "cannot open data.txt!\n");
	}

	malloc2E(arraydata2, n, choose);  
	for (i = 0; i<n; i++)
		for (j = 0; j<choose; j++)
			fscanf(fp, "%lf", &arraydata2[i][j]);  
	return arraydata2;
}

