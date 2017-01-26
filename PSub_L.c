/* ================================================================= */
/*   Subgradient method for Goods Reccomendation Problem             */
/*                                                                   */
/*   Programmed by Takahiro Kan <takahiro.kan@ist.osaka-u.ac.jp>     */
/*   Date: 2016/11/15                                                */
/* ================================================================= */

/* including files ------------------------------------------------- */

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
#define NUM_EPSILON 0.001
#define NUM_INF DBL_MAX 
#define NUM_MIN DBL_MIN 
#define N 5 /* numbers of reccomendation Goods 	*/
#define J 310
#define SAMPLE 5000
#define RAN_NUM 0
     
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
void PSubgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr); 
void k_select(double* cf,int* idx);
void swap(int i,int j,double* A,int* B);
void LH(Grpdata* grpdata_ptr,Vdata* vdata_ptr,double* d_ptr);
int check_feasible(Grpdata* grpdata_ptr,int** x);
void shuffle(int I,int* A);



/* cpu_time.c */
double cpu_time();

/* functions ------------------------------------------------------- */


int main(int argc, char *argv[])
{ 	
	init_genrand(RAN_NUM);
	Grpdata grpdata;
	Vdata vdata;
	vdata.ub = NUM_INF; 
	vdata.lb = 0.0;

	FILE* input_file;

	if(argc != 2){
	  fprintf(stderr,"Please input the name of data file!\n");
	  		exit(1);
	    }

/* scan I,res and cost */
	input_file = fopen(argv[1], "r");
	fscanf(input_file,"%d",&grpdata.I);
	grpdata.res = (double*)malloc_e(sizeof(double)*grpdata.I);
	grpdata.cost = (double**)malloc_e(sizeof(double*)*grpdata.I);
    for (int i = 0; i < grpdata.I; ++i)
    {
    	fscanf(input_file,"%lf",&grpdata.res[i]); 
    	grpdata.cost[i] = (double*)malloc_e(sizeof(double)*J);
    	for (int j = 0; j < J; ++j)
    	{
    		fscanf(input_file,"%lf",&grpdata.cost[i][j]); 
    	}
    }
	fclose(input_file); 

/*	make Cmax & P[j]	*/
    int* r_check;
    r_check =(int*)malloc_e(sizeof(int)*N);
    int r_sum;
    int ch_sum;
   	int a;

/*	make Cmax	*/	
	grpdata.Cmax = 0.0;
    for (int i = 0; i < grpdata.I; ++i)
    {
	    r_sum = 0;
	    for (int j = 0; j < N; ++j)
	    {
	    	r_check[j] = -1;
	    }
	    
		while(r_sum < N)
	    {
	    	int ch_sum = 0;
	    	a = (int)(genrand_real1() * J);
	    	for (int k = 0; k < N; ++k)
	    	{
	    		if (a == r_check[k])
	    		{
			    	break;
	    		}
	    		else
	    		{
	    			++ch_sum;
	    		}
	    	}	

		   	if (ch_sum==N)
		   	{
			    ++r_sum;
			    grpdata.Cmax += grpdata.cost[i][a];
			    r_check[r_sum-1] = a;
	    	}
	   	}
	}
    //grpdata.Cmax = grpdata.Cmax*1.2;
    

/*	make P[J]*/
	double P_sum;
	grpdata.P = (double*)malloc_e(sizeof(double)*J);
	for (int j = 0; j < J; ++j)
	{
	    r_sum = 0;
	    P_sum = 0.0;
	    for (int i = 0; i < N; ++i)
	    {
	    	r_check[i] = -1;
	    }

		while(r_sum < N)
	    {
	    	ch_sum = 0;
	    	a = (int)(genrand_real1() * grpdata.I);
	    	for (int k = 0; k < N; ++k)
	    	{
	    		if (a == r_check[k])
	    		{
			    	break;
	    		}
	    		else
	    		{
	    			++ch_sum;
	    		}
	    	}	
		   	if (ch_sum==N)
		   	{
			    ++r_sum;
			    r_check[r_sum-1] = a;
			    P_sum += grpdata.res[a]*grpdata.cost[a][j];
	    	}

	   	}
	   	grpdata.P[j] = P_sum;
	   	//grpdata.P[j] = P_sum*1.1;
	}

/*	make sol	*/
	vdata.sol = (int**)malloc_e(sizeof(int*)*grpdata.I);
	for (int i = 0; i < grpdata.I; ++i)
	{
		vdata.sol[i] = (int*)malloc_e(sizeof(int)*N);
	}

	printf("I = %d\n",grpdata.I);
	printf("Cmax = %lf\n",grpdata.Cmax);
/*	
	printf("res[0] = %lf\n",grpdata.res[0]);
	printf("cost[0][0] = %lf\n",grpdata.cost[0][0]);
	printf("res[I-1] = %lf\n",grpdata.res[grpdata.I-1]);
	printf("cost[I-1][J-1] = %lf\n",grpdata.cost[grpdata.I-1][J-1]);
*/

/*
	for (int j = 0; j < J; ++j)
	{
		printf("%lf ",grpdata.P[j]);
		if (j%10==9)
		{
			printf("\n");
		}
	}
*/

	//printf("parameter:alpha\n");
	printf("parameter:beta\n");

	printf("RAN_NUM:%d\n",RAN_NUM);

/* subgrident method */
	vdata.start_time = cpu_time();
	PSubgrad(&grpdata,&vdata);
	vdata.finish_time = cpu_time();
	
	printf("CPU time: %.3f seconds\n",vdata.finish_time - vdata.start_time); 

/* free memory	*/
	for (int i = 0; i < grpdata.I; ++i)
	{
		free(grpdata.cost[i]);
	}
	free(grpdata.cost);

	for (int i = 0; i < grpdata.I; ++i)
	{
		free(vdata.sol[i]);
	}
	free(vdata.sol);

	free(grpdata.res);
	free(grpdata.P);
	return 0;
}


/* ------------------------------------------------------- */
/*  		 PSubgradient Methods 		                   */
/*														   */
/*           grpdata_ptr(I):data of the GRP                */
/*			 vdata_ptr(I/O):various data				   */
/* ------------------------------------------------------- */

void PSubgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr){
	double best_ub = NUM_INF;
	double GAP = NUM_INF;
	int k = 1;
	int t = 0;
	double lmd = 0.0;

	double* mu;
	mu = (double*)malloc_e(sizeof(double)*J);
	for (int j = 0; j < J; ++j)
	{
		mu[j] = 0.0;
	}

	double s;
	double** cf;
	cf = (double**)malloc_e(sizeof(double*)*grpdata_ptr->I);

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		cf[i] = (double*)malloc_e(sizeof(double)*J);
	}

	int* idx;
	idx = (int*)malloc_e(sizeof(int)*J);

	double* sum_Gi;
	sum_Gi = (double*)malloc_e(sizeof(double)*J);

	double* d;
	d = (double*)malloc_e(sizeof(double)*(J + 1));

	double* d_all;
	d_all = (double*)malloc_e(sizeof(double)*(J + 1));

	int* I_idx;
	I_idx = (int*)malloc_e(sizeof(int)*grpdata_ptr->I);

	double sum_C;
	double sum_G;
	double sum_muG;
	double d_nol;
	double pre_ub = NUM_INF;
	double gg; 
	
	do{
		if (t > 0)
		{
			pre_ub = vdata_ptr->ub;
		}

		//printf("t=%d\n",t);

		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			I_idx[i] = i;
		}

		for (int j = 0; j < J; ++j)
		{
			sum_Gi[j] = 0.0;
		}
		sum_C = 0.0;
		sum_G = 0.0;
		sum_muG = 0.0;
		d_nol = 0.0;


/*	sampling of I  */
		shuffle(grpdata_ptr->I,I_idx);

/*	sort cf 	*/
		for (int i = 0; i < SAMPLE; ++i)
		{ 

			for (int j = 0; j < J; ++j)
			{
				idx[j] = j;
				cf[i][j] = (1 + mu[j])*grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][j] -lmd*grpdata_ptr->cost[I_idx[i]][j];
			}
	

			k_select(cf[i],idx);
			for (int n = 0;n < N; ++n)
			{
				vdata_ptr->sol[i][n] = idx[n];
			}
		}


/*	calculation of d	 */
		for (int i = 0; i < SAMPLE; ++i)
		{
			for (int n = 0;  n < N; ++n)
			{
				sum_C += grpdata_ptr->cost[I_idx[i]][vdata_ptr->sol[i][n]];
				sum_G += grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][vdata_ptr->sol[i][n]];
				sum_Gi[vdata_ptr->sol[i][n]] += grpdata_ptr->res[I_idx[i]]*grpdata_ptr->cost[I_idx[i]][vdata_ptr->sol[i][n]];
			}
		}
		d[0] = grpdata_ptr->Cmax - grpdata_ptr->I*sum_C/(double)SAMPLE;
		d_nol += d[0]*d[0];
		for (int j = 1; j < J+1; ++j)
		{
			d[j] = grpdata_ptr->I*sum_Gi[j-1]/(double)SAMPLE - grpdata_ptr->P[j-1];
			sum_muG += mu[j-1]*d[j];
			d_nol += d[j]*d[j];
		}

/*	recalculation of lmd,mu 	*/
		d_nol = (double)sqrt(d_nol);

/* ------------------------------------------------------- */
/*   			parameter settings		        	       */
/*                                                         */
/* ------------------------------------------------------- */

		//s = pow(0.9,t)/d_nol;		//alpha
/**/
		if (t==0)					//beta
		{
			s = (double)1/(double)d_nol;
		}else
		{
			s =(double)1/(double)((t)*d_nol);
		}
/**/	

		if (lmd-s*d[0] > 0.0)
		{
			lmd = lmd-s*d[0];
		}else{
			lmd = 0.0;
		}

		for (int j = 1; j < J+1; ++j)
		{
			if(mu[j-1]-s*d[j] > 0.0)
			{
				mu[j-1] = mu[j-1]-s*d[j];
			}else{
				mu[j-1]= 0.0;
			}
		}
		t++;
		s =(double)1/(double)(t*d_nol);

		if (t*SAMPLE >= k*grpdata_ptr->I)
		{ 

			for (int j = 0; j < J; ++j)
			{
				sum_Gi[j] = 0.0;
			}
			sum_C = 0.0;
			sum_G = 0.0;
			sum_muG = 0.0;
/*	sort cf[] 	*/
			for (int i = 0; i < grpdata_ptr->I; ++i)
			{
				for (int j = 0; j < J; ++j)
				{
					idx[j] = j;
					cf[i][j] = (1 + mu[j])*grpdata_ptr->res[i]*grpdata_ptr->cost[i][j] - lmd*grpdata_ptr->cost[i][j];
				}

				k_select(cf[i],idx);
				for (int n = 0;n < N; ++n)
				{
					vdata_ptr->sol[i][n] = idx[n];
				}
			}

/*	calculation of ub	 */
			for (int i = 0; i < grpdata_ptr->I; ++i)
			{
				for (int n = 0;  n < N; ++n)
				{
					sum_C += grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
					sum_G += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
					sum_Gi[vdata_ptr->sol[i][n]] += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
				}
			}

			d_all[0] = grpdata_ptr->Cmax - sum_C;
			for (int j = 1; j < J+1; ++j)
			{
				d_all[j] = sum_Gi[j-1] - grpdata_ptr->P[j-1];
				sum_muG += mu[j-1]*d[j];
			}
			vdata_ptr->ub = sum_G + lmd*d_all[0] + sum_muG;
			if (vdata_ptr->ub < best_ub)
			{
				best_ub = vdata_ptr->ub;
			}
			//printf("ub:%lf\n",best_ub);

			LH(grpdata_ptr,vdata_ptr,d_all);
			gg = GAP;
			GAP = 1 - vdata_ptr->lb/best_ub;
			if (gg != GAP)
			{
				printf("GAP:%lf \n",GAP);
			}
			k++;
		}

	}while(GAP >= NUM_EPSILON);
		printf("number of occurence:%d\n",t);
		printf("best_ub:%lf\n",best_ub);
		printf("best_lb:%lf\n",vdata_ptr->lb);
/* free memory	*/
		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			free(cf[i]);
		}
		free(mu);
		free(cf);
		free(idx);
		free(sum_Gi);
		free(d);
		free(I_idx);
		return;
	}


/* ------------------------------------------------------- */
/*  		 Lagragian Heuristics 		                   */
/*														   */
/*           grpdata_ptr(I):data of the GRP                */
/*			 vdata_ptr(I/O):various data				   */
/*			 d_ptr(I):array of Subgrad 					   */
/* ------------------------------------------------------- */

void LH(Grpdata* grpdata_ptr,Vdata* vdata_ptr,double* d_ptr){
	double fpenalty = 0.0;
	double penalty;
	int r;
	int r_pre;
	int r_n;
	int r_sum;
	int test;
	double P_pre;
	double C_pre;
	double P_sec = 0.0;
	double P_sec_pre;
	double temp_lb = 0.0;

/*	make Lsol_tld	*/
	int**  Lsol_tld;
	Lsol_tld = (int**)malloc_e(sizeof(int*)*grpdata_ptr->I);
	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		Lsol_tld[i] = (int*)malloc_e(sizeof(int)*N);
		for (int n = 0; n < N; ++n)
		{
			Lsol_tld[i][n] = vdata_ptr->sol[i][n];
		}
	}

	double* P_cp;
	P_cp = (double*)malloc_e(sizeof(double)*(J + 1));
	P_cp[0] = -d_ptr[0];

/*	make first penalty(fpenalty)	*/
	if (P_cp[0] >= 0.0)
	{
		fpenalty += (double)J*P_cp[0]/grpdata_ptr->Cmax;
	}
	
	for (int j = 1; j < J + 1; ++j)
	{
		P_cp[j] = -d_ptr[j];
		if (-d_ptr[j] >= 0.0)
		{
			fpenalty += P_cp[j]/grpdata_ptr->P[j-1];
			P_sec += P_cp[j]/grpdata_ptr->P[j-1];
		}

	}

	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		penalty = 0.0;
/*	make neighbor of sol	*/
		r_sum = 0;
		while(1){
			r = (int)(genrand_real1() * J);
			for (int n = 0; n < N; ++n)
			{ 
				if (r != Lsol_tld[i][n])
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
				r_n = (int)(genrand_real1() * N);
				r_pre = Lsol_tld[i][r_n];
				Lsol_tld[i][r_n] = r;
				break;
			}
		}

		C_pre = P_cp[0];
		P_cp[0] = P_cp[0] - grpdata_ptr->cost[i][r_pre] + grpdata_ptr->cost[i][r];
		P_pre = P_cp[r+1];
		P_cp[r+1] = P_cp[r+1] + grpdata_ptr->res[i]*grpdata_ptr->cost[i][r_pre] - grpdata_ptr->res[i]*grpdata_ptr->cost[i][r];


		if (P_cp[0] >= 0.0)
		{
			penalty += (double)J*P_cp[0]/grpdata_ptr->Cmax;
		}

		P_sec_pre = P_sec;
		P_sec += - P_pre/grpdata_ptr->P[r];

		if (P_cp[r+1] >= 0.0)
		{
			P_sec += P_cp[r+1]/grpdata_ptr->P[r];
		}
		penalty  += P_sec;


		if (penalty < fpenalty)
		{  
			fpenalty = penalty;
		}else
		{
			Lsol_tld[i][r_n] = r_pre;
			P_cp[0] = C_pre;
			P_cp[r+1] = P_pre;
			P_sec = P_sec_pre;

		}


		if (fpenalty < NUM_EPSILON && fpenalty > -NUM_EPSILON)
		{
/*	calculation of lb	*/		
			if (check_feasible(grpdata_ptr,Lsol_tld) == TRUE)
			{

				for (int i = 0; i < grpdata_ptr->I; ++i)
				{
					for (int n = 0; n < N; ++n)
					{
						temp_lb += grpdata_ptr->res[i]*grpdata_ptr->cost[i][Lsol_tld[i][n]];
					}
				}
				if (temp_lb > vdata_ptr->lb)
				{
					vdata_ptr->lb = temp_lb;
				}

				break;
			}			
		}
	}
	//test = check_feasible(grpdata_ptr,Lsol_tld);
	printf("fpenalty:%lf\n",fpenalty);	
	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		free(Lsol_tld[i]);
	}
	free(Lsol_tld);
	free(P_cp);
	return;
	
}



/* ------------------------------------------------------- */
/*   			rearrange of array		 				   */
/*                                                         */
/* ------------------------------------------------------- */

void k_select(double* cf,int* idx){
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

	return;
}


/* ------------------------------------------------------- */
/*   		swap of array(both double & int)	      	   */
/*                                                         */
/* ------------------------------------------------------- */

void swap(int i,int j,double* A,int* B)
{
	double temp1;
	int temp2;

	temp1 = A[i];
	A[i] = A[j];
	A[j] = temp1;

	temp2 = B[i];
	B[i] = B[j];
	B[j] = temp2;
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
	//printf("P_sum%d\n",P_sum);
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
/*   Memory allocation with error check                    */
/*                                                         */
/*   size(I): size of allocating memory                    */
/*   return_value: pointer to allocated memory             */
/* ------------------------------------------------------- */
void *malloc_e(size_t size){
  
  void *s;

  if((s = malloc(size)) == NULL){
    fprintf(stderr,"malloc_e: Not enough memory.\n");
    exit(EXIT_FAILURE);
  }
  return(s);
}

void shuffle(int I,int* A)
{
	int temp;
	for (int i = 0; i < SAMPLE; ++i)
	{
		int p = (int)(genrand_real1() * I);
		temp = A[i];
		A[i] = A[p];
		A[p] = temp;

	}

	return;
}
