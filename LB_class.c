/* ================================================================= */
/*   Subgradient method for Goods Reccomendation Problem             */
/*                                                                   */
/*   Programmed by Takahiro Kan <takahiro.kan@ist.osaka-u.ac.jp>     */
/*   Date: 2016/11/15                                                */
/* ================================================================= */

/* including files ------------------------------------------------- */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "mt19937ar.h"

/* definitions(constant) ------------------------------------------- */

#define TRUE 1
#define FALSE 0
#define INVALID -1
#define NUM_EPSILON 1e-10
#define NUM_INF DBL_MAX 
#define NUM_MIN DBL_MIN 
#define N 5 /* numbers of reccomendation Goods 	*/
#define J 310
#define SAMPLE_for_av 10000

#define LB_TIME 300
#define ITR_NUM 0
     
/* definitions(type) ----------------------------------------------- */

/* data of grp problem */
typedef struct{
	int I;				/*	numbers of the customers */
	double Cmax;
	double* P;
	double** cost;
	double* res; 
} Grpdata;

/* various data often necessary during the search */
typedef struct{
  	double start_time;  	/* the time the search started */
  	double finish_time;  	/* the time the search finished */

  	int** sol;  			/* solution */
	double ub;
	double lb;
} Vdata;

/* local file */
void* malloc_e(size_t size);

void Construction(Grpdata* grpdata_ptr,Vdata* vdata_ptr); 
void quicksort_mirror(int i,int j,double* A_original,int* idx);
void quicksort(int i,int j,double* A,int* idx);
int partiton(int i,int j,double a,double* A,int* idx);
int pivot(int i,int j,double* A);
void swap_d(int i,int j,double* A);
void swap_i(int i,int j,int* A);
void shuffle_for_av(int I,int* A);
void de_select(double* cf_original,int num,int* idx);
double a_select(double* cf_original,int num);
void swap(int i,int j,double* A,int* B);
int check_feasible(Grpdata* grpdata_ptr,int** x);
void LH(Grpdata* grpdata_ptr,Vdata* vdata_ptr);


/* cpu_time.c */
double cpu_time();

/* functions ------------------------------------------------------- */

void Construction(Grpdata* grpdata_ptr,Vdata* vdata_ptr)
{
	vdata_ptr->lb = 0.0;
	printf("LB_TIME:%d ITR_NUM:%d\n",LB_TIME,ITR_NUM);
	double c_av = 0.0;
	double Cr;
	double P_zero = 0.0;
	double C_rest = grpdata_ptr->Cmax;

	int* I_idx;
	I_idx = (int*)malloc_e(sizeof(int)*grpdata_ptr->I);

	int* J_idx;
	J_idx = (int*)malloc_e(sizeof(int)*J);

	int* J_idx_2;
	J_idx_2 = (int*)malloc_e(sizeof(int)*J);

	double* P_check;
	P_check = (double*)malloc_e(sizeof(double)*J);

/*	make averafe of C */
	for (int j = 0; j < J; ++j)
	{
		P_zero += grpdata_ptr->P[j];
	}

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		I_idx[i] = i;
	}
	shuffle_for_av(grpdata_ptr->I,I_idx);

	for (int i = 0; i < SAMPLE_for_av; ++i)
	{
		c_av += a_select(grpdata_ptr->cost[I_idx[i]],J);	
	}

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		I_idx[i] = i;
	}



	quicksort_mirror(0,grpdata_ptr->I-1,grpdata_ptr->res,I_idx);



	Cr = P_zero/grpdata_ptr->res[I_idx[0]] + c_av*(double)grpdata_ptr->I/(double)SAMPLE_for_av;
	//C_rest = -C_rest*0.6;



	for (int j = 0; j < J; ++j)
	{
		P_check[j] = grpdata_ptr->P[j];
	}

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{	


		if (Cr <= C_rest)
		{
/*	step 1	*/


/*	calculation of P_pnealty		*/
			for (int j = 0; j < J; ++j)
			{
				J_idx[j] = j;
			}

			de_select(grpdata_ptr->cost[I_idx[i]],J,J_idx);

			for (int n = 0; n < N; ++n)
			{

				vdata_ptr->sol[I_idx[i]][n] = J_idx[n];
				Cr += grpdata_ptr->cost[I_idx[i]][J_idx[n]];
				C_rest = C_rest - grpdata_ptr->cost[I_idx[i]][J_idx[n]];
			}
	

			double P_temp1 = 0.0;

			for (int j = 0; j < J; ++j)
			{
				for (int n = 0; n < N; ++n)
				{
					if (j == vdata_ptr->sol[I_idx[i]][n])
					{
						P_check[j] -= grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][j];
					}
				}

				if ( P_check[j] > 0.0)
				{ 
					P_temp1 += P_check[j];
				}		
			}

			Cr = P_temp1/grpdata_ptr->res[I_idx[i]] + c_av*(double)(grpdata_ptr->I - i)/(double)SAMPLE_for_av;

		}else
		{
/*	step 2-A	*/

			int k = 0;

			for (int j = 0; j < J; ++j)
			{
				J_idx_2[j] = j;
			}
			quicksort_mirror(0,J-1,grpdata_ptr->cost[I_idx[i]],J_idx_2);

			for (int j = 0; j < J; ++j)
			{
				if (P_check[J_idx_2[j]] > 0.0)
				{
					vdata_ptr->sol[I_idx[i]][k] = J_idx_2[j];
					P_check[J_idx_2[j]] -= grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][J_idx_2[j]];
					k++;
				}

				if (k == N)
				{
					break;
				}
			}

/*	step 2-B	*/

			if (k < N)
			{
				for (int j = J-1; j >= 0; --j)
				{
					if (P_check[J_idx_2[j]] <= 0.0)
					{
						vdata_ptr->sol[I_idx[i]][k] = J_idx_2[j];
						P_check[J_idx_2[j]] -= grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][J_idx_2[j]];
						k++;
					}
					if (k == N)
					{
						break;
					}
				}
			}
		}

	}
	if(check_feasible(grpdata_ptr,vdata_ptr->sol) == TRUE)
	{
		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			for (int n = 0; n < N; ++n)
			{
				vdata_ptr->lb += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
			}
		}
		printf("first lb:%lf\n",vdata_ptr->lb);
	}

	free(I_idx);
	free(J_idx);
	free(J_idx_2);
	free(P_check);
	
	return;

}

/* ------------------------------------------------------- */
/*   			descending_select		   				   */
/*                                                         */
/* ------------------------------------------------------- */

void de_select(double* cf_original,int num,int* idx)
{

	double* cf;
	cf = (double*)malloc_e(sizeof(double)*num);

	memcpy(cf,cf_original,sizeof(double)*num);

	int i,j,left,right;
	int k = N - 1;
	double x;
	left = 0;
	right = J - 1;
	while(left < right){
		x = cf[k];

		i = left;
		j = right;
		while(1){
			while(cf[i] > x) i++;
			while(cf[j] < x) j--;
			if(i >= j) {break;}
			else{
				swap(i,j,cf,idx);
				i++; j--;
			}
		}
		if(i <= k) left = j + 1;
		if(k <= j) right = i - 1;
	}
	free(cf);
	return;
}

/* ------------------------------------------------------- */
/*   			ascending_select		   				   */
/*                                                         */
/* ------------------------------------------------------- */

double a_select(double* cf_original,int num){

	double c_sum = 0.0;

	double* cf;
	cf = (double*)malloc_e(sizeof(double)*num);

	memcpy(cf,cf_original,sizeof(double)*num);


	int i,j,left,right;
	int k = N - 1;

	double x;
	left = 0;
	right = J - 1;
	while(left < right){
		x = cf[k];

		i = left;
		j = right;
		while(1){
			while(cf[i] < x) i++;
			while(cf[j] > x) j--;
			if(i >= j) {break;}
			else{
				swap_d(i,j,cf);
				i++; j--;
			}
		}
		if(i <= k) left = j + 1;
		if(k <= j) right = i - 1;
	}

	for (int n = 0; n < N; ++n)
	{
		c_sum += cf[n];
	}
	free(cf);
	return(c_sum);
}


/* ------------------------------------------------------- */
/*   		quicksort(desending)				      	   */
/*                                                         */
/* ------------------------------------------------------- */

void quicksort(int i,int j,double* A,int* idx)
{
	int pv,k;
	double a;

	pv = pivot(i,j,A);
	if(pv != -1)
	{
		a = A[pv];
		k = partiton(i,j,a,A,idx);
		quicksort(i,k-1,A,idx);
		quicksort(k,j,A,idx);
	}
	return;
}


void quicksort_mirror(int i,int j,double* A_original,int* idx)
{
	double* A;
	A = (double*)malloc_e(sizeof(double)*(j+1));
	memcpy(A,A_original,sizeof(double)*(j+1));

	int pv,k;
	double a;

	pv = pivot(i,j,A);
	if(pv != -1)
	{
		a = A[pv];
		k = partiton(i,j,a,A,idx);
		quicksort(i,k-1,A,idx);
		quicksort(k,j,A,idx);
	}
	free(A);
	return;
	
}

int partiton(int i,int j,double a,double* A,int* idx)
{
	int l,r,k;

	l = i;
	r = j;
	while(1)
	{
		while(A[l] >= a)  l = l + 1;
		while(A[r] < a) r = r - 1;
		if(l <= r )
		{
			swap_d(l,r,A);
			swap_i(l,r,idx);
			l++;
			r--;
		}
		else break;
	}
	k = l;
	return(k);
}

int pivot(int i,int j,double* A)
{
	int pv,k;

	k= i + 1;
	while(k <= j &&  A[i]== A[k]) k = k + 1;
	if(k > j) pv = -1;
	else if(A[i] >= A[k]) pv = i;	
	else pv = k;
	return(pv);
}

void swap_d(int i,int j,double* A)
{
	double temp;
	temp = A[i];
	A[i] = A[j];
	A[j] = temp;
	return;
}

void swap_i(int i,int j,int* idx)
{
	int temp;
	temp = idx[i];
	idx[i] = idx[j];
	idx[j] = temp;
	return;
}


/* ------------------------------------------------------- */
/*  		 Check feasible solutions	                   */
/*                                                         */
/*           grpdata_ptr(I):data of the GRP                */
/*		     x(I):array of solutions					   */
/* ------------------------------------------------------- */

int check_feasible(Grpdata* grpdata_ptr,int** x){

	double SUM_C = 0.0;
	double* SUM_G;
	int P_sum = 0;
	SUM_G = (double*)malloc_e(sizeof(double)*(J));
	for (int j = 0; j < J; ++j)
	{
		SUM_G[j] = 0.0;
	}

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		for (int n = 0; n < N; ++n)
		{
			SUM_C += grpdata_ptr->cost[i][x[i][n]];
			SUM_G[x[i][n]] += grpdata_ptr->res[i]*grpdata_ptr->cost[i][x[i][n]];
		}

	}

	for (int j = 0; j < J; ++j)
	{
		if (SUM_G[j] >= grpdata_ptr->P[j])
		{
			P_sum++;
		}
	}

	free(SUM_G);
	//printf("P_sum:%d\n",P_sum);
	//printf("delta:%lf\n",SUM_C - grpdata_ptr->Cmax);

	if (SUM_C <= grpdata_ptr->Cmax && P_sum == J)
	{
		//printf("feasilbe!\n");
		return TRUE;
	}else
	{
		//printf("violation!\n");
		return FALSE;
	}
}

/* ------------------------------------------------------- */
/*  		 Lagragian Heuristics 		                   */
/*					1& 2								   */
/*           grpdata_ptr(I):data of the GRP                */
/*			 vdata_ptr(I/O):various data				   */
/*			 exchange in two customers					   */
/* ------------------------------------------------------- */

void LH(Grpdata* grpdata_ptr,Vdata* vdata_ptr){
	int dumy;
	double fpenalty = 0.0;
	double penalty;
	double P_rpre;
	double P_pre;
	double P_pre1;
	double P_pre2;
	double C_pre;
	double P_sec = 0.0;
	double P_sec_pre;
	double temp_lb;
	double temp;

/*	make sol_tld	*/
	int**  sol_tld;
	sol_tld = (int**)malloc_e(sizeof(int*)*grpdata_ptr->I);

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		sol_tld[i] = (int*)malloc_e(sizeof(int)*N);
		for (int n = 0; n < N; ++n)
		{
			sol_tld[i][n] = vdata_ptr->sol[i][n];
		}
	}

	double* P_cp;
	P_cp = (double*)malloc_e(sizeof(double)*(J + 1));


	double sum_C;
	double sum_G;

	double* sum_Gi;
	sum_Gi = (double*)malloc_e(sizeof(double)*J);

	for (int j = 0; j < J; ++j)
	{
		sum_Gi[j] = 0.0;
	}
	sum_C = 0.0;
	sum_G = 0.0;

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		for (int n = 0;  n < N; ++n)
		{
			sum_C += grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
			sum_G += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
			sum_Gi[vdata_ptr->sol[i][n]] += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
		}
	}

	P_cp[0] = sum_C - grpdata_ptr->Cmax;

/*	make first penalty(fpenalty)	*/
	if (P_cp[0] > 0.0)
	{
		fpenalty += (double)J*P_cp[0]/grpdata_ptr->Cmax;
	}
	
	for (int j = 1; j < J + 1; ++j)
	{
		P_cp[j] = grpdata_ptr->P[j-1] - sum_Gi[j-1];
		if (P_cp[j] > 0.0)
		{
			fpenalty += P_cp[j]/grpdata_ptr->P[j-1];
			P_sec += P_cp[j]/grpdata_ptr->P[j-1];
		}

	}

	//printf("P_sec:%.30lf\n",P_sec);



/*	make neighbor of 2 customers' sol	*/

	double s_time = cpu_time();
	int r,r_n,r_pre,r_sum,i1,i2,r1,r2,r_check;
	int iteration = 0;
	while(vdata_ptr->lb == 0.0)
	{
		temp_lb = 0.0;
		penalty = 0.0;

		i1 = (int)(genrand_real2() * grpdata_ptr->I);

/*	single swap 	*/


			while(1){			
				r_sum = 0;
				r = (int)(genrand_real2() * J);

				for (int n = 0; n < N; ++n)
				{ 
					if (r != sol_tld[i1][n])
					{
						r_sum++;
					}else
					{ 
						r_sum = 0;
						break;
					}
				}

				if (r_sum == N)
				{
					r_n = (int)(genrand_real2() * N);
					r_pre = sol_tld[i1][r_n];
					sol_tld[i1][r_n] = r;
					break;
				}

			}

			C_pre = P_cp[0];
			P_cp[0] = P_cp[0] - grpdata_ptr->cost[i1][r_pre] + grpdata_ptr->cost[i1][r];

			P_rpre = P_cp[r_pre+1];
			P_cp[r_pre+1] = P_cp[r_pre+1] + grpdata_ptr->res[i1]*grpdata_ptr->cost[i1][r_pre];


			P_pre = P_cp[r+1];
			P_cp[r+1] = P_cp[r+1]  - grpdata_ptr->res[i1]*grpdata_ptr->cost[i1][r];


			if (P_cp[0] > 0.0)
			{
				penalty += (double)J*P_cp[0]/grpdata_ptr->Cmax;
			}

			P_sec_pre = P_sec;

			if (P_rpre > 0.0)
			{
				P_sec -= P_rpre/grpdata_ptr->P[r_pre];

			}

			if (P_pre > 0.0)
			{
				P_sec -= P_pre/grpdata_ptr->P[r];
			}

			if (P_cp[r_pre+1] > 0.0)
			{
				P_sec += P_cp[r_pre+1]/grpdata_ptr->P[r_pre];
			}

			if (P_cp[r+1] > 0.0)
			{
				P_sec += P_cp[r+1]/grpdata_ptr->P[r];
			}

			if (P_sec > 0.0)
			{
				penalty += P_sec;
			}

				//printf("penalty:%.30lf\n",penalty);

			if (penalty < NUM_EPSILON)
			{

/*	calculation of lb	*/	
				for (int i = 0; i < grpdata_ptr->I; ++i)
				{
					for (int n = 0; n < N; ++n)
					{
						temp_lb += grpdata_ptr->res[i]*grpdata_ptr->cost[i][sol_tld[i][n]];
					}
				}

					//dumy = check_feasible(grpdata_ptr,sol_tld);
				if (temp_lb > vdata_ptr->lb && check_feasible(grpdata_ptr,sol_tld) == TRUE)
				{
					//printf("1111111111\n");
					//printf("penalty:%.30lf P_cp[0]:%lf P_sec:%.30lf\n",penalty,P_cp[0],P_sec);
					vdata_ptr->lb = temp_lb;
					//printf("lb:%lf\n",temp_lb);
					//printf("ITR_NUM:%d\n",iteration);
					iteration = 0;

				}else
				{
					P_cp[r+1] = P_pre;
					P_cp[r_pre+1] = P_rpre;
					sol_tld[i1][r_n] = r_pre;
					P_cp[0] = C_pre;
					P_sec = P_sec_pre;
					++iteration;
				}


			}else if(penalty < fpenalty)
			{
				fpenalty = penalty;
				//printf("111 penalty:%.30lf\n",penalty);
			}else
			{
				P_cp[r+1] = P_pre;
				P_cp[r_pre+1] = P_rpre;
				sol_tld[i1][r_n] = r_pre;
				P_cp[0] = C_pre;
				P_sec = P_sec_pre;
				++iteration;
			}

/*	double swap		*/
		temp_lb = 0.0;
		penalty = 0.0;

		while(1){

			r_check = 0;
			i2 = (int)(genrand_real2() * grpdata_ptr->I);
			
			while(i1 == i2)
			{
				i2 = (int)(genrand_real2() * grpdata_ptr->I);
			}

			r1 = (int)(genrand_real2() * N);
			r2 = (int)(genrand_real2() * N);

			for (int n = 0; n < N; ++n)
			{
				if (sol_tld[i1][r1] != sol_tld[i2][n] && sol_tld[i2][r2] != sol_tld[i1][n])
				{
					++r_check;
				}else
				{
					break;
				}
			}

			if (r_check == N)
			{
				break;
			}
		}

		temp = sol_tld[i1][r1];
		sol_tld[i1][r1] = sol_tld[i2][r2];
		sol_tld[i2][r2] = temp;

		C_pre = P_cp[0];
		P_cp[0] = P_cp[0] - grpdata_ptr->cost[i1][sol_tld[i2][r2]] + grpdata_ptr->cost[i1][sol_tld[i1][r1]]- grpdata_ptr->cost[i2][sol_tld[i1][r1]] + grpdata_ptr->cost[i2][sol_tld[i2][r2]];
		
		P_pre1 = P_cp[sol_tld[i1][r1]+1];
		P_cp[sol_tld[i1][r1]+1] = P_cp[sol_tld[i1][r1]+1] + grpdata_ptr->res[i2]*grpdata_ptr->cost[i2][sol_tld[i1][r1]] - grpdata_ptr->res[i1]*grpdata_ptr->cost[i1][sol_tld[i1][r1]];

		P_pre2 = P_cp[sol_tld[i2][r2]+1];
		P_cp[sol_tld[i2][r2]+1] = P_cp[sol_tld[i2][r2]+1] + grpdata_ptr->res[i1]*grpdata_ptr->cost[i1][sol_tld[i2][r2]] - grpdata_ptr->res[i2]*grpdata_ptr->cost[i2][sol_tld[i2][r2]];

		if (P_cp[0] > 0.0)
		{
			penalty += (double)J*P_cp[0]/grpdata_ptr->Cmax;
		}




		P_sec_pre = P_sec;

		if (P_pre1 > 0.0)
		{
			P_sec +=  - P_pre1/grpdata_ptr->P[sol_tld[i1][r1]];
		}

		if (P_pre2 > 0.0)
		{
			P_sec += - P_pre2/grpdata_ptr->P[sol_tld[i2][r2]];
		}
		
		if (P_cp[sol_tld[i1][r1]+1] > 0.0)
		{
			P_sec += P_cp[sol_tld[i1][r1]+1]/grpdata_ptr->P[sol_tld[i1][r1]];
		}

		if (P_cp[sol_tld[i2][r2]+1] > 0.0)
		{
			P_sec += P_cp[sol_tld[i2][r2]+1]/grpdata_ptr->P[sol_tld[i2][r2]];
		}


		if (P_sec > 0.0)
		{
			penalty += P_sec;
		}
		//printf("P_sec:%.30lf\n",P_sec);
		

		//printf("penalty:%.30lf\n",penalty);


		if (penalty < NUM_EPSILON)
		{

/*	calculation of lb	*/	
			for (int i = 0; i < grpdata_ptr->I; ++i)
			{
				for (int n = 0; n < N; ++n)
				{
					temp_lb += grpdata_ptr->res[i]*grpdata_ptr->cost[i][sol_tld[i][n]];
				}
			}
			//printf("delta lb:%lf\n",vdata_ptr->lb - temp_lb);
			//dumy = check_feasible(grpdata_ptr,sol_tld);
			if (temp_lb > vdata_ptr->lb && check_feasible(grpdata_ptr,sol_tld) == TRUE)
			{
				//printf("penalty:%.30lf P_cp[0]:%lf P_sec:%.30lf\n",penalty,P_cp[0],P_sec);
				//printf("22222222222\n");
				vdata_ptr->lb = temp_lb;
				//printf("lb:%lf\n",temp_lb);
				//printf("ITR_NUM:%d\n",iteration);
				iteration = 0;

			}else
			{
				P_cp[sol_tld[i1][r1]+1] = P_pre1;
				P_cp[sol_tld[i2][r2]+1] = P_pre2;
				sol_tld[i2][r2] = sol_tld[i1][r1];
				sol_tld[i1][r1] = temp;
				P_cp[0] = C_pre;
				P_sec = P_sec_pre;
				++iteration;
			}

		}else if(penalty < fpenalty)
		{
			fpenalty = penalty;
			//printf("222 penalty:%.30lf\n",penalty);
		}else
		{
			P_cp[sol_tld[i1][r1]+1] = P_pre1;
			P_cp[sol_tld[i2][r2]+1] = P_pre2;
			sol_tld[i2][r2] = sol_tld[i1][r1];
			sol_tld[i1][r1] = temp;
			P_cp[0] = C_pre;
			P_sec = P_sec_pre;
			++iteration;
		}

	}
	

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		free(sol_tld[i]);
	}
	free(sol_tld);
	free(P_cp);
	return;	
}

/* ------------------------------------------------------- */
/*   		shuffle of I 						      	   */
/*                                                         */
/* ------------------------------------------------------- */

void shuffle_for_av(int I,int* A)
{
	int temp;
	for (int i = 0; i < SAMPLE_for_av; ++i)
	{
		int p = (int)(genrand_real2() * I);
		temp = A[i];
		A[i] = A[p];
		A[p] = temp;

	}
	return;
}
