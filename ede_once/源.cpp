#include <stdio.h>
#include <stdlib.h>
//#include "boost/random.hpp"
#include "ede_once.h"


/*****************************************************************
random gaussian number 
*****************************************************************/
double gaussrand(double U, double V)
{
	//static double U, V;
	static int phase = 0;
	double Z;

	if (phase == 0)
	{
		U = (rand() + 1.1) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-1 * log(U))* sin(2 * PI * V);
	}
	else
	{
		Z = sqrt(-2 * long(U)) * cos(2 * PI * V);
	}

	phase = 1 - phase;
	Z = Z * 0.1 + 0;
	return Z;
}

/*****************************************************************
*load function
*input：data.txt，
D-data dimension, N-num. of all data
*****************************************************************/
double **loadData(int choose, int *d, int *n)
{
	int i, j;
	double **arraydata;  
	FILE *fp;
	if (choose == 1) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/twospiral.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 2) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/spiral.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 3) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/corner.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 4) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/cluster2.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 5) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/half.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	else if (choose == 6) {
		if ((fp = fopen("C:/Users/richard/source/repos/ede_once/moon.txt", "r")) == NULL)    fprintf(stderr, "cannot open data.txt!\n");
	}
	if (fscanf(fp, "D=%d,N=%d\n", d, n) != 2)        fprintf(stderr, "load error!\n");

	malloc2E(arraydata, *n, *d);  
	for (i = 0; i<*n; i++)
		for (j = 0; j<*d; j++)
			fscanf(fp, "%lf", &arraydata[i][j]);  
	return arraydata;
}

/*****************************************************************
*evos-r function
*input： population p，dimension D，iter time Gmax
*output：best results 
*****************************************************************/
int de(FILE *fp, FILE *fp2, double *p1, double *p2, double *p3, int *p4, int *valids, double *alpha1,
	double *alpha2, double *alpha3, int *relat1, int N, int D, int Gmax) {
	int i, j, k, l, m, r1, r2, r3, r4, is_exist, clu_num2, top_ele2;
	int NP = 100, numofE = 0, noise_num2, count_arr, inclu_num2, index;
	int *used1, *used2, **clu_div2, **inval2, **pop_inclu, **reach2; 

	double F = 0.5, CR = 0.5, best_val=0.0, best_val2=DBL_MAX, wheel=0.0;
	double **next_index, **next_param;
	queue<int> q2;

	malloc1E(used1, N);
	malloc1E(used2, N);
	malloc2D(clu_div2, NP, N);      //label record
	malloc2D(inval2, NP, 20);       //invalid location
	malloc2D(pop_inclu, NP, 20);    //inclu num.
	malloc2D(reach2, N, N);         //reachable
	malloc2E(next_index, NP, 40);   //mutant encoding
	malloc2E(next_param, NP, 2);    //results of mutants

	double **alpha_no2, **alpha_in2, **alpha_out2;
	int **con_fun2;

	malloc2E(alpha_no2, N, N);
	malloc2E(alpha_in2, N, N);
	malloc2E(alpha_out2, N, N);
	malloc2D(con_fun2, N, 2);

	//self_incre
	double r_mean2, r_std2, *group_radius2;
	malloc1D(group_radius2, 20);

	//arrays to save the better cr and f values
	int better_count = -1;
	double mu_cr = 0.5, mu_f = 0.5;
	double *s_cr, *s_f;
	malloc1D(s_cr, 100);
	malloc1D(s_f, 100);

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			alpha_no2[i][j] = *(alpha1 + i * N + j);
			alpha_in2[i][j] = *(alpha2 + i * N + j);
			alpha_out2[i][j] = *(alpha3 + i * N + j);
		}
		for (j = 0; j < 2; j++) {
			con_fun2[i][j] = *(relat1 + i * 2 + j);
		}
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			reach2[i][j] = 0;
		}
	}
	for (i = 0; i < NP; i++) {
		for (j = 0; j < 20; j++) {
			inval2[i][j] = 0;
			pop_inclu[i][j] = 0;
		}
	}

	for (i = 0; i < 100; i++) {
		s_cr[i] = 0;
		s_f[i] = 0;
	}

	printf("begin the DE algorithm...\n");
	for (i = 0; i < Gmax; i++) {
		for (j = 0; j < NP; j++) {
			//mutation_two parts
			do {
				r1 = (int)(NP*URAND);
			} while (r1 == j);
			do {
				r2 = (int)(NP*URAND);
			} while (r2 == j || r2 == r1);
			do {
				r3 = (int)(NP*URAND);
			} while (r3 == j || r3 == r1 || r3 == r2);

			//operations based on sets
			for (k = 0; k < 20; k++) {
				is_exist = 1;
				for (l = 0; l < 20; l++) {
					if ((int)*(p1 + r2 * 40 + k) == (int)*(p1 + r3 * 40 + l)) {
						next_index[j][k] = -1;
						is_exist = 0; 
					}
				}
				if (is_exist == 1) {
					next_index[j][k] = *(p1 + r2 * 40 + k);
				}
			}

			//SD
			double F_i;
			do {
				F_i = generateGaussianNoise(mu_f, 0.1);
			} while (F_i <= 0 || F_i >= 1);

			for (k = 0; k < 20; k++) {
				if ((int)next_index[j][k] == -1)
					continue;
				if (URAND >= F_i) {
					next_index[j][k] = -1;
				}
				else {
					used1[(int)next_index[j][k]] = 1;
				}	
			}

			//V=r1+SD
			for (k = 0; k < 20; k++) {
				if ((int)next_index[j][k] == -1) {   //add from r1
					if (used1[(int)*(p1 + r1 * 40 + k)] == 1) {  //rerandom
						do {
							r4 = (int)(N*URAND);    
						} while (used1[r4] == 1);
						next_index[j][k] = r4;
						used1[r4] = 1;  
					}
					else {
						next_index[j][k] = *(p1 + r1 * 40 + k);  
						used1[(int) next_index[j][k]] = 1;
					}
				}
			}

			//upper and lower bounds
			for (k = 0; k < 20; k++) {
				next_index[j][k + 20] = *(p1 + r1 * 40 + k + 20) + F * (*(p1 + r2 * 40 + k + 20) - *(p1 + r3 * 40 + k + 20));
				if (next_index[j][k + 20] > sqrt(D) / 2) {
					next_index[j][k + 20] = sqrt(D) / 2;
				}
				else if(next_index[j][k + 20] <= 0){
					next_index[j][k + 20] = *(p1 + r1 * 40 + k + 20);   
				}
			}

			//crossover
			double CR_i; 
			do {
				CR_i = generateGaussianNoise(mu_cr, 0.1);
			} while (CR_i<=0 || CR_i>=1);

			for (k = 0; k < N; k++)
				used1[k] = 0;
			for (k = 0; k < 20; k++) {
				if (URAND >= CR_i) {  
					next_index[j][k] = *(p1 + j * 40 + k);	
					next_index[j][k+20] = *(p1 + j * 40 + k + 20);
				}
				if (used1[(int) next_index[j][k]] == 1) {
					do {
						r4 = (int)(N*URAND);    
					} while (used1[r4] == 1);
					next_index[j][k] = r4;
					used1[r4] = 1;  
				}
				else {
					used1[(int)next_index[j][k]] = 1;
				}
			}
		
			/* data arrangement */
			for (k = 0; k < N; k++)
				used2[k] = 0;   
			clu_num2 = 0;
			for (k = 0; k < N; k++)
				clu_div2[j][k] = 0;   

			for (k = 0; k < 20; k++) {
				inval2[j][k] = 0;  
				if (used2[(int)next_index[j][k]] == 1) {
					inval2[j][k] = 1;   
					continue;
				}

				clu_num2++;
				clu_div2[j][(int) *(p1 + j * 40 + k)] = clu_num2;   
				inclu_num2 = 1;
				q2.push((int) *(p1 + j * 40 + k));
				used2[(int) *(p1 + j * 40 + k)] = 1;

				for (l = 0; l < N; l++) {
					for (m = 0; m < N; m++) {
						if (*(p3 + l * N + m) <= *(p1 + j * 40 + k + 20)) {  
							reach2[l][m] = 1;
						}
						else {
							reach2[l][m] = 0;
						}
					}
				}

				while (!q2.empty()) {
					top_ele2 = q2.front();
					q2.pop();
					for (l = 0; l < N; l++) {
						if (l == top_ele2 or used2[l] == 1)
							continue;
						if (reach2[top_ele2][l] == 1) {
							clu_div2[j][l] = clu_num2;
							q2.push(l);
							used2[l] = 1;
							inclu_num2++;
						}
					}
				}

				pop_inclu[j][k] = inclu_num2;
				if (inclu_num2 < 3) {
					for (l = 0; l < N; l++) {
						if (clu_div2[j][l] == clu_num2) {
							clu_div2[j][l] = 0;
							used2[l] = 0;
						}
					}
					clu_num2--;
					inval2[j][k] = 1;   
				}
			}

			next_param[j][0] = clu_num2;
			noise_num2 = 0;
			for (k = 0; k < N; k++) {
				if (clu_div2[j][k] == 0)
					noise_num2++;
			}

			count_arr = 0;
			while (noise_num2 > (N / 10)) {
				count_arr++;
				if (count_arr > 3) {
					break;
				}

				int count_inval = 0;
				for (k = 0; k < 20; k++) {
					if (inval2[j][k] == 1) { 
						group_radius2[k] = 0;
						count_inval++;
					}
					else {
						group_radius2[k] = next_index[j][k + 20];
					}
				}
				
				r_mean2 = sum_mean(group_radius2, 20);
				r_std2 = sum_std(group_radius2, 20);

				for (k = 0; k < 20; k++) {
					if (inval2[j][k] == 1)
						continue;
					wheel = (N - pop_inclu[j][k]) *1.0 / (N*1.0);
					if (URAND < wheel) {
						next_index[j][k + 20] = next_index[j][k + 20] + r_std2;
					}
						
				}

				for (k = 0; k < N; k++)
					used2[k] = 0;  
				clu_num2 = 0;
				for (k = 0; k < N; k++)
					clu_div2[j][k] = 0;   

				for (k = 0; k < 20; k++) {
					inval2[j][k] = 0;  
					if (used2[(int)next_index[j][k]] == 1) {
						inval2[j][k] = 1;   
						continue;
					}

					clu_num2++;
					clu_div2[j][(int) *(p1 + j * 40 + k)] = clu_num2;   
					inclu_num2 = 1;
					q2.push((int) *(p1 + j * 40 + k));
					used2[(int) *(p1 + j * 40 + k)] = 1;

					for (l = 0; l < N; l++) {
						for (m = 0; m < N; m++) {
							if (*(p3 + l * N + m) <= *(p1 + j * 40 + k + 20)) {  
								reach2[l][m] = 1;
							}
							else {
								reach2[l][m] = 0;
							}
						}
					}

					while (!q2.empty()) {
						top_ele2 = q2.front();
						q2.pop();
						for (l = 0; l < N; l++) {
							if (l == top_ele2 or used2[l] == 1)
								continue;
							if (reach2[top_ele2][l] == 1) {
								clu_div2[j][l] = clu_num2;
								q2.push(l);
								used2[l] = 1;
								inclu_num2++;
							}
						}
					}

					pop_inclu[j][k] = inclu_num2;
					if (inclu_num2 < 3) {
						for (l = 0; l < N; l++) {
							if (clu_div2[j][l] == clu_num2) {
								clu_div2[j][l] = 0;
								used2[l] = 0;
							}
						}
						clu_num2--;
						inval2[j][k] = 1;   
					}
				}
				next_param[j][0] = clu_num2;
				noise_num2 = 0;
				for (k = 0; k < N; k++) {
					if (clu_div2[j][k] == 0)
						noise_num2++;
				}
			}


			//fitness
			if ((int) next_param[j][0] == 1) {
				next_param[j][1] = DBL_MAX;
			}
			else {
				next_param[j][1] = func_nk2(next_param[j][0], clu_div2[j], *alpha_no2, *alpha_in2, *alpha_out2, *con_fun2, N, D, noise_num2);
			}

			numofE++;
			if (numofE % 500 == 0)
				printf("the %d th fitness evaluations...\n", numofE);
			
			//selection(greedy)
			if (next_param[j][1] < *(p2 + j * 2 + 1)) {
				for (k = 0; k < 40; k++) {
					*(p1 + j * 40 + k) = next_index[j][k];
				}
				*(p2 + j * 2 + 0) = next_param[j][0];
				*(p2 + j * 2 + 1) = next_param[j][1];

				for (k = 0; k < N; k++) {
					*(p4 + j * N + k) = clu_div2[j][k];
				}

				for (k = 0; k < 20; k++) {
					*(valids + j * 20 + k) = inval2[j][k];
				}

				//fitness better than parents need to be updated(the self adaptive strategy)
				better_count++;
				if (better_count < 100) {
					s_cr[better_count] = CR_i;
					s_f[better_count] = F_i;
					mu_cr = 0.9*mu_cr + 0.1*sum_mean(s_cr, better_count + 1);
					mu_f = 0.9*mu_f + 0.1*sum_mean(s_f, better_count + 1);
				}
				else {
					s_cr[better_count%100] = CR_i;
					s_f[better_count%100] = F_i;
					mu_cr = 0.9*mu_cr + 0.1*sum_mean(s_cr, 100);
					mu_f = 0.9*mu_f + 0.1*sum_mean(s_f, 100);
				}
			}

			if (*(p2 + j * 2 + 1) < best_val2) {
				best_val = *(p2 + j * 2 + 0);
				best_val2 = *(p2 + j * 2 + 1);
				printf("The %d th evaluation: %d, NK2_index: %f\n", numofE, (int)best_val, best_val2);
				index = j;
			}
			fprintf(fp, "%f\n", best_val2);
			fprintf(fp2, "%d\n", (int) best_val);
		}
	}
	printf("end the DE algorithm...\n");

	return index;
}

int  main()
{
	srand((unsigned int)(time(NULL)));   //seed

	//file record, adjust the path according to your pc setting
	FILE *fshow = NULL, *fshow2 = NULL, *fshow3 = NULL, * fshow4 = NULL, *fshow5 = NULL, *fshow6 = NULL;
	fshow = fopen("C:/../results/index_results.txt", "w");
	if (!fshow) {
		perror("cannot open file!\n");
		return -1;
	}
	fshow2 = fopen("C:/../results/partition_results.txt", "w");
	if (!fshow2) {
		perror("cannot open file!\n");
		return -1;
	}
	fshow3 = fopen("C:/../results/count_results.txt", "w");
	if (!fshow3) {
		perror("cannot open file!\n");
		return -1;
	}
	fshow4 = fopen("C:/../results/radius_results.txt", "w");
	if (!fshow4) {
		perror("cannot open file!\n");
		return -1;
	}
	fshow5 = fopen("C:/../results/fitness_iteration.txt", "w");
	if (!fshow5) {
		perror("cannot open file!\n");
		return -1;
	}
	fshow6 = fopen("C:/../results/cluster_iteration.txt", "w");
	if (!fshow6) {
		perror("cannot open file!\n");
		return -1;
	}

	//parameters and data structure
	int i, j, k, l, D, N, Gmax, NP, best = 0, tmp = 0, clu_num, top_ele1, noise_num, count_arr, label, inclu_num;
	int choose, k_correct;

	double data_min, data_max, dist, dist_min = DBL_MAX, tmp2=0.0, wheel0=0.0, param1 = 0, param2 = 0, param3 = 0;
	double *uk, *lk;
	double intra_dist = 0, inter_dist = 0, intra_sum = 0, totaltime = 0.0;

	double **pop_cod1, **pop_cod2;   //encoding
	double **all_dist, **all_dist3, **reach;  //distance between data
	double **data_smp;  
	int *used, *used2, *lhs_label, *disturb0, *valids0, **clu_div, **inval, **org_inclu; 
	queue<int> q1;

	//nk2
	int **den_cmp, **con_fun;
	double den_mean, den_std, den_sum, dist_mean, dist_std, c1, c2, c3, c_den;
	double *den_data, *all_dist2, **alpha_no, **alpha_in, **alpha_out, **group_p, **group_c;

	//self_increase
	double r_mean, r_std;
	double *group_radius;

	printf("Choose the dataset:1.twospiral 2.spiral 3.corner 4.cluster2 5.half 6.moon\n");
	scanf("%d", &choose);

	//load	
	data0 = loadData(choose, &D, &N);      

	//definition
	NP = 100;  
	Gmax = 1000/ NP; //set the iteration time
	double theta; 
	int percent = 2;
	int position = (((N * N) / 2)*percent) / 100;
	printf("The times of iteration(Gmax):%d\n", Gmax);

	
	malloc1D(uk, D);
	malloc1D(lk, D);
	malloc1E(used, N);
	malloc1E(used2, N);
	malloc1E(lhs_label, 20);      //num. of sample
	malloc1E(disturb0, 20);
	malloc1E(valids0, 20);
	malloc2E(pop_cod1, NP, 40);   //population encoding
	malloc2E(pop_cod2, NP, 2);    //result
	malloc2D(clu_div, NP, N);     //acutual label
	malloc2D(inval, NP, 20);      //location of invalid 
	malloc2D(org_inclu, NP, 20);  //inclu num. of points
	malloc2E(all_dist, N, N);     
	malloc2E(all_dist3, N, N);     
	malloc2E(reach, N, N);		  //is_reachable
	malloc2E(data_smp, 20, D);    

	//new allocation
	malloc1D(den_data, N);
	malloc1D(all_dist2, N*N);
	malloc2D(den_cmp, N, N);
	malloc2D(con_fun, N, 2);
	malloc2E(alpha_no, N, N);
	malloc2E(alpha_in, N, N);
	malloc2E(alpha_out, N, N);
	malloc2E(group_p, N, 3);   
	malloc2E(group_c, N, 3);   //c1-c3 of each group
	malloc1D(group_radius, 20);

	clock_t starttime, endtime;
	data_smp = lhs(D, 20);

	//initialization
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			reach[i][j] = 0;
		}
	}
	for (i = 0; i < NP; i++) {
		for (j = 0; j < 20; j++) {
			inval[i][j] = 0;
			org_inclu[i][j] = 0;
		}
	}

	//prepare
	for (i = 0; i<D; i++) {
		data_min = DBL_MAX, data_max = LOW_VAL;
		for (j = 0; j<N; j++) {
			if (data0[j][i] > data_max)
				data_max = data0[j][i];
			if (data0[j][i] < data_min)
				data_min = data0[j][i];
		}
		uk[i] = data_max;
		lk[i] = data_min;
		for (j = 0; j < N; j++)
			data0[j][i] = (data0[j][i] - lk[i]) / (uk[i] - lk[i]);
	}

	for (i = 0; i < D; i++) {
		uk[i] = 1;
		lk[i] = 0;
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			all_dist2[i*N+j]=getDistance2(data0[i], data0[j], D);
		}
	}
	quicksort(all_dist2, N*N);
	theta = all_dist2[position];

	double dt_min, dt_mid;
	dt_min = all_dist2[0];
	dt_mid = all_dist2[N*N / 2];

	//individual density, all_dist and all_dist2 arrays
	for (i = 0; i < N; i++) {
		den_sum = 0.0;
		for (j = 0; j < N; j++) {
			if (j == i) {
				all_dist[i][j] = 0;
				den_sum = den_sum + 1;
				continue;
			}
			all_dist[i][j] = getDistance2(data0[i], data0[j], D);
			all_dist3[i][j]= getDistance(data0[i], data0[j], D);
			den_sum = den_sum + pow(e, ((-1)*all_dist[i][j]* all_dist[i][j])/ (theta*theta));
		}
		den_data[i] = den_sum;
	}
	double den_max = max_fun(den_data, N);
	for (i = 0; i < N; i++) {
		den_data[i] = den_data[i] / den_max;
	}

	//density, std
	den_mean = sum_mean(den_data, N);
	den_std = sum_std(den_data, N);

	if (den_mean > den_std) {
		c_den = den_mean - den_std;
	}
	else {
		c_den = den_mean / 2.0;
	}

	//density compare
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (den_data[i] > den_data[j])
				den_cmp[i][j] = 1;
			else
				den_cmp[i][j] = 0;
		}
	}

	//01.twp related objects(nk2)
	int sel_index, pair;
	double sel_min;
	for (i = 0; i < N; i++) {
		pair = 0;
		sel_min = DBL_MAX;
		for (j = 0; j < N; j++) {
			if (j == i || den_cmp[i][j] == 1) {
				continue;
			}
			if (all_dist[i][j] < sel_min) {
				pair = 1;
				sel_min = all_dist[i][j];
				sel_index = j;
			}
		}
		if (pair == 1)
			con_fun[i][0] = sel_index;
		else
			con_fun[i][0] = -1;   
	}
	//02.find the value of the second location(nk2)
	for (i = 0; i < N; i++) {
		if (con_fun[i][0] == -1) {   
			sel_min = DBL_MAX;
			for (j = 0; j < N; j++) {
				if (j == i)
					continue;
				if (all_dist[i][j] < sel_min) {
					sel_min = all_dist[i][j];
					sel_index = j;
				}
			}
			con_fun[i][0] = sel_index;
		}
		sel_min = DBL_MAX;
		for (j = 0; j < N; j++) {
			if ((j == i) || (j == con_fun[i][0]))
				continue;
			if (all_dist[i][j] < sel_min) {
				sel_min = all_dist[i][j];
				sel_index = j;
			}
		}
		con_fun[i][1] = sel_index;
	}
	//03.calculation of the c1-c3 and group_c(nk2)
	for (i = 0; i < N; i++) {
		for (j = 0; j < 3; j++) {
			group_p[i][0] = all_dist[i][con_fun[i][0]];
			group_p[i][1] = all_dist[i][con_fun[i][1]];
			group_p[i][2] = all_dist[con_fun[i][0]][con_fun[i][1]];
		}
		dist_mean = sum_mean(group_p[i], 3);
		dist_std=sum_std(group_p[i], 3);
		group_c[i][0] = dist_mean;
		group_c[i][1] = dist_mean+dist_std;
		group_c[i][2] = dist_mean+2.0*dist_std;
	}

	//the three alpha arrays
	//01.alpha_noise
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (j == i) {
				alpha_no[i][j] = 0;  
				continue;
			}
			if ((all_dist[i][j] > group_c[i][1]) && (den_data[i] <= c_den)) {
				alpha_no[i][j] = 0;
			}
			else {
				alpha_no[i][j] = den_data[j];
			}
		}
	}
	//02.alpha_in
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (j == i) {
				alpha_in[i][j] = 0;  
				continue;
			}
			if (all_dist[i][j] < group_c[i][0]) {
				alpha_in[i][j] = 0;
			}
			else if (all_dist[i][j] >= group_c[i][0] && all_dist[i][j] <= group_c[i][2]) {
				alpha_in[i][j] = (all_dist[i][j] - group_c[i][0]) / (group_c[i][2] - group_c[i][0]) * den_data[j];
			}
			else {
				alpha_in[i][j] = den_data[j];
			}
		}
	}
	//03.alpha_out
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (j == i) {
				alpha_out[i][j] = 0; 
				continue;
			}
			if (all_dist[i][j] < group_c[i][0]) {
				alpha_out[i][j] = den_data[j];
			}
			else if (all_dist[i][j] >= group_c[i][0] && all_dist[i][j] <= group_c[i][2]) {
				alpha_out[i][j] = (group_c[i][2] - all_dist[i][j]) / (group_c[i][2] - group_c[i][0]) * den_data[j];
			}
			else {
				alpha_out[i][j] = 0;
			}
		}
	}

	for (i = 0; i < N; i++) {
		used[i] = 0;  
	}
	//
	for (i = 0; i < 20; i++) {
		dist_min = DBL_MAX;
		for (j = 0; j < N; j++) {
			if (used[j] == 1)
				continue;
			dist = getDistance(data_smp[i], data0[j], D);
			if (dist < dist_min) {
				dist_min = dist;
				label = j;
			}
		}
		lhs_label[i] = label;
		used[label] = 1;
	}

	//iteration settings
	int iter,iter_times=1;
	double *iter_fitness, *iter_runtime;
	malloc1D(iter_fitness, iter_times);
	malloc1D(iter_runtime, iter_times);

	fprintf(fshow3, "the cluster number in %d times runs:\n", iter_times);
	for (iter = 0; iter < iter_times; iter++) {
		starttime = clock();   //start
		printf("the %d th run:\n", iter + 1);
		//population initialization
		for (i = 0; i < NP; i++) {
			//selection of radii
			for (j = 0; j < N; j++) {
				used[j] = 0;   
			}

			for (j = 0; j < 20; j++) {
				disturb0[j] = j;
			}
			permutate(disturb0, 20);

			for (j = 0; j < 20; j++) {
				pop_cod1[i][j] = lhs_label[disturb0[j]];
				//small random value 
				do {
					tmp2 = URAND * (sqrt(D)) / (30 * 1.0);
				} while (tmp2 == 0.0);
				pop_cod1[i][j + 20] = tmp2;
			}

			/* data arrangement */

			for (j = 0; j < N; j++)
				used2[j] = 0;   
			clu_num = 0;        
			for (j = 0; j < N; j++)
				clu_div[i][j] = 0;   //noise at first

			for (j = 0; j < 20; j++) {
				//start from the first point
				inval[i][j] = 0;
				if (used2[(int)pop_cod1[i][j]] == 1) {
					inval[i][j] = 1;   //invalid record
					continue;
				}

				clu_num++;
				clu_div[i][(int)pop_cod1[i][j]] = clu_num;   //label 
				inclu_num = 1;
				q1.push((int)pop_cod1[i][j]);
				used2[(int)pop_cod1[i][j]] = 1;

				for (k = 0; k < N; k++) {
					for (l = 0; l < N; l++) {
						if (all_dist3[k][l] <= pop_cod1[i][j + 20]) {
							reach[k][l] = 1;
						}
						else {
							reach[k][l] = 0;
						}
					}
				}

				while (!q1.empty()) {
					top_ele1 = q1.front();
					q1.pop();
					for (k = 0; k < N; k++) {
						if (k == top_ele1 or used2[k] == 1)
							continue;
						if (reach[top_ele1][k] == 1) {
							clu_div[i][k] = clu_num;
							q1.push(k);
							used2[k] = 1;
							inclu_num++;
						}
					}
				}

				//cluster judgement
				org_inclu[i][j] = inclu_num;
				if (inclu_num < 3) {
					for (k = 0; k < N; k++) {
						if (clu_div[i][k] == clu_num) {
							clu_div[i][k] = 0;
							used2[k] = 0;
						}
					}
					clu_num--;
					inval[i][j] = 1;
				}
			}
			pop_cod2[i][0] = clu_num;

			//sum of the num. of noise
			noise_num = 0;
			for (j = 0; j < N; j++) {
				if (clu_div[i][j] == 0)
					noise_num++;
			}

			//统计对应的平均半径和标准差，作为自增的参考
			//只统计有效的半径
			int count_inval = 0;
			for (j = 0; j < 20; j++) {
				if (inval[i][j] == 1) {  
					group_radius[j] = 0;
					count_inval++;
				}
				else {
					group_radius[j] = pop_cod1[i][j + 20];
				}
			}

			r_mean = sum_mean(group_radius, 20);
			r_std = sum_std(group_radius, 20);

			//most points are labelled as noises, rearrangement
			count_arr = 0;
			while (noise_num > (N / 10)) {
				count_arr++;
				if (count_arr > 3) {
					break;
				}

				for (j = 0; j < 20; j++) {
					if (inval[i][j] == 1)
						continue;
					wheel0 = (N - org_inclu[i][j])*1.0 / (N*1.0);
					if (URAND < wheel0) {
						pop_cod1[i][j + 20] = pop_cod1[i][j + 20] + r_std;  //self-incre
					}
				}

				for (j = 0; j < N; j++)
					used2[j] = 0;   
				clu_num = 0;        
				for (j = 0; j < N; j++)
					clu_div[i][j] = 0;   

				for (j = 0; j < 20; j++) {
					inval[i][j] = 0;  
					if (used2[(int)pop_cod1[i][j]] == 1) {
						inval[i][j] = 1;   
						continue;
					}

					clu_num++;
					clu_div[i][(int)pop_cod1[i][j]] = clu_num;   
					inclu_num = 1;
					q1.push((int)pop_cod1[i][j]);
					used2[(int)pop_cod1[i][j]] = 1;

					for (k = 0; k < N; k++) {
						for (l = 0; l < N; l++) {
							if (all_dist3[k][l] <= pop_cod1[i][j + 20]) {
								reach[k][l] = 1;
							}
							else {
								reach[k][l] = 0;
							}
						}
					}

					while (!q1.empty()) {
						top_ele1 = q1.front();
						q1.pop();
						for (k = 0; k < N; k++) {
							if (k == top_ele1 or used2[k] == 1)
								continue;
							if (reach[top_ele1][k] == 1.0) {
								clu_div[i][k] = clu_num;
								q1.push(k);
								used2[k] = 1;
								inclu_num++;
							}
						}
					}

					org_inclu[i][j] = inclu_num;
					if (inclu_num < 3) {
						for (k = 0; k < N; k++) {
							if (clu_div[i][k] == clu_num) {
								clu_div[i][k] = 0;
								used2[k] = 0;
							}
						}
						clu_num--;
						inval[i][j] = 1;
					}
				}
				pop_cod2[i][0] = clu_num;

				noise_num = 0;
				for (j = 0; j < N; j++) {
					if (clu_div[i][j] == 0)
						noise_num++;
				}
			}

			//fitness
			if ((int)pop_cod2[i][0] == 1) {
				pop_cod2[i][1] = DBL_MAX;
			}
			else {
				pop_cod2[i][1] = func_nk2(pop_cod2[i][0], clu_div[i], *alpha_no, *alpha_in, *alpha_out, *con_fun, N, D, noise_num);
			}
		}

		printf("initial...\n");

		fprintf(fshow5, "the %d th run:\n", iter + 1);
		fprintf(fshow6, "the %d th run:\n", iter + 1);
		best = de(fshow5, fshow6, *pop_cod1, *pop_cod2, *all_dist3, *clu_div, *inval, *alpha_no, *alpha_in, *alpha_out, *con_fun, N, D, Gmax);

		printf("done...\n");
		endtime = clock();  	//end
		totaltime = (double)(endtime - starttime);

		//results output
		int is_noise = 0, is_true = 0;
		iter_fitness[iter] = pop_cod2[best][1];
		iter_runtime[iter] = totaltime / (double)CLOCKS_PER_SEC;
		
		fprintf(fshow, "the %d th run:\n",iter+1);
		printf("Execution time: %.3f s\n", totaltime / (double)CLOCKS_PER_SEC);
		printf("The best numbers of the clusters are:%d\n", (int)pop_cod2[best][0]);
		printf("\nThe best(min) NK2 index value is: ");
		printf("%f\n", pop_cod2[best][1]);

		fprintf(fshow, "// Counts // Fitness // Execution time :\n");
		fprintf(fshow, "%d\n", (int)pop_cod2[best][0]);
		fprintf(fshow, "%f\n", pop_cod2[best][1]);
		fprintf(fshow, "%f\n", totaltime / (double)CLOCKS_PER_SEC);
		fprintf(fshow3, "%d ", (int)pop_cod2[best][0]);
	}

	double mean_fitnesss, std_fitness, mean_runtime, std_runtime;
	mean_fitnesss = sum_mean(iter_fitness, iter_times);
	std_fitness = sum_std(iter_fitness, iter_times);
	mean_runtime = sum_mean(iter_runtime, iter_times);
	std_runtime = sum_std(iter_runtime, iter_times);

	fprintf(fshow3, "\n");
	fprintf(fshow3, "the average fitness is:\n");
	fprintf(fshow3, "%f %f\n", mean_fitnesss, std_fitness);
	fprintf(fshow3, "the average run time is:\n");
	fprintf(fshow3, "%f %f\n", mean_runtime, std_runtime);

	system("pause");
	return 0;
}


