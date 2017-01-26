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
#define N 5 			/* numbers of reccomendation Goods for one customer	*/
#define J 310			/* numbers of reccomendation Goods */
#define RAN_NUM1 0
#define Rc 0.5
#define Rp 1.5

#define UB_TIME 1800
     
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
void Subgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr); 
void k_select(double* cf,int* idx);
void swap(int i,int j,double* A,int* B);

/* cpu_time.c */
double cpu_time();

/*	LB_class.c 	*/
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
int check_feasible(Grpdata* grpdata_ptr,int** x);
void LH(Grpdata* grpdata_ptr,Vdata* vdata_ptr);



/* functions ------------------------------------------------------- */


int main(int argc, char *argv[])
{ 	
	double input_start_time = cpu_time();		
	init_genrand(RAN_NUM1);
	Grpdata grpdata;
	Vdata vdata;
	vdata.ub = NUM_INF; 

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

/*	make Cmax	*/	
	
	grpdata.Cmax = 0.0;
	for (int i = 0; i < grpdata.I; ++i)
	{
		for (int j = 0; j < J; ++j)
		{
			grpdata.Cmax += grpdata.cost[i][j];
		}
	}
	grpdata.Cmax = (double)N*grpdata.Cmax/(double)J;
	grpdata.Cmax = grpdata.Cmax*Rc;


/*	make P[J]*/


	double P_sum;
	grpdata.P = (double*)malloc_e(sizeof(double)*J);

	for (int j = 0; j < J; ++j)
	{
		P_sum =0.0;
		for (int i = 0; i < grpdata.I; ++i)
		{
			P_sum += grpdata.res[i]*grpdata.cost[i][j];
		}
		grpdata.P[j] = (double)N*P_sum/(double)J;
		grpdata.P[j] = grpdata.P[j]*Rp;

	}



/*	make sol	*/
	vdata.sol = (int**)malloc_e(sizeof(int*)*grpdata.I);
	for (int i = 0; i < grpdata.I; ++i)
	{
		vdata.sol[i] = (int*)malloc_e(sizeof(int)*N);
	}

	printf("SM---------\n");
	printf("I = %d\n",grpdata.I);
	//printf("Cmax = %lf\n",grpdata.Cmax);

/* ------------------------------------------------------- */
/*   			parameter settings		        	       */
/*                                                         */
/* ------------------------------------------------------- */

	printf("LB&UB TIME:%d seconds\n",UB_TIME);
	printf("RAN_NUM1:%d\n",RAN_NUM1);
	printf("Rc:%.1lf Rp:%.1lf\n",Rc,Rp);
	printf("input time:%.3lf seconds\n",cpu_time() - input_start_time);

/* subgrident method */
	vdata.start_time = cpu_time();
	Construction(&grpdata,&vdata);

	if (vdata.lb == 0.0)
	{
		LH(&grpdata,&vdata);
	}

	printf("LB time: %.3f seconds\n",cpu_time() - vdata.start_time);
	printf("UB time: %.3f seconds\n",UB_TIME - (cpu_time() - vdata.start_time));
	//printf("best lb:%lf\n",vdata.lb);


	Subgrad(&grpdata,&vdata);
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
/*  		 Subgradient Methods 		                   */
/*														   */
/*           grpdata_ptr(I):data of the GRP                */
/*			 vdata_ptr(I/O):various data				   */
/* ------------------------------------------------------- */

void Subgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr){

/* ------------------------------------------------------- */
/*   			output file  			        	       */
/*                                                         */
/* ------------------------------------------------------- */
	FILE* output_file; 
	output_file = fopen("Sub.csv","w");
	double best_ub = NUM_INF;
	int t = 0;
	double lmd = 0.0;

	double* mu;
	mu = (double*)malloc_e(sizeof(double)*J);
	for (int j = 0; j < J; ++j)
	{
		mu[j] = 0.0;
	}

	double s;
	double* cf;
	cf = (double*)malloc_e(sizeof(double)*J);

	int* idx;
	idx = (int*)malloc_e(sizeof(int)*J);

	double* sum_Gi;
	sum_Gi = (double*)malloc_e(sizeof(double)*J);

	double* d;
	d = (double*)malloc_e(sizeof(double)*(J + 1));

	double sum_C;
	double sum_G;
	double sum_muG;
	double d_nol;
	double sub_start_time = cpu_time();
	
	do{
		for (int j = 0; j < J; ++j)
		{
			sum_Gi[j] = 0.0;
		}
		sum_C = 0.0;
		sum_G = 0.0;
		sum_muG = 0.0;
		d_nol = 0.0;

		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			for (int j = 0; j < J; ++j)
			{
				idx[j] = j;
				cf[j] = (1 + mu[j])*grpdata_ptr->res[i]*grpdata_ptr->cost[i][j] - lmd*grpdata_ptr->cost[i][j];				
			}

			k_select(cf,idx);

			
			
			for (int n = 0;n < N; ++n)
			{  
				vdata_ptr->sol[i][n] = idx[n];
			}
		}



/*	calculation of ub & d	 */

		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			for (int n = 0;  n < N; ++n)
			{
				sum_C += grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
				sum_G += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
				sum_Gi[vdata_ptr->sol[i][n]] += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
			}
		}
		d[0] = grpdata_ptr->Cmax - sum_C;
		d_nol += d[0]*d[0];
		
		for (int j = 1; j < J+1; ++j)
		{
			d[j] = sum_Gi[j-1] - grpdata_ptr->P[j-1];
			sum_muG += mu[j-1]*d[j];
			d_nol += d[j]*d[j];
		}
		vdata_ptr->ub = sum_G + lmd*d[0] + sum_muG;

		if (vdata_ptr->ub < best_ub)
		{ 

			if (cpu_time() - vdata_ptr->start_time < UB_TIME)
			{
				fprintf(output_file,"%lf,%lf\n",cpu_time() - vdata_ptr->start_time,(vdata_ptr->ub - vdata_ptr->lb)*100.0/vdata_ptr->ub);
				//printf("gap:%lf\n",(vdata_ptr->ub - vdata_ptr->lb)*100.0/vdata_ptr->ub);
			}
			best_ub = vdata_ptr->ub;

		}
		

/*	recalculation of lmd,mu 	*/
		d_nol = (double)sqrt(d_nol);

/* ------------------------------------------------------- */
/*   			parameter settings		        	       */
/*                                                         */
/* ------------------------------------------------------- */
		//s = pow(0.9,t)/d_nol;		//alpha


		s = (double)1/(double)((t+1)*d_nol); //beta
		

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


	}while(cpu_time() -  vdata_ptr->start_time < UB_TIME);
	printf("number of occurence:%d\n",t);
	printf("iteration per second:%.3lf\n",(double)t/(cpu_time() - sub_start_time));
	printf("best ub:%.3lf\n",best_ub);
	printf("best lb:%.3lf\n",vdata_ptr->lb);
	printf("best gap:%.3lf %%\n",(best_ub - vdata_ptr->lb)*100.0/best_ub);
	fclose(output_file);



/* free memory	*/
	free(mu);
	free(cf);
	free(idx);
	free(sum_Gi);
	free(d);
	return;

}


/* ------------------------------------------------------- */
/*   			rearrange of array		   				   */
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


