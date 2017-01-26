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

/* definitions(constant) ------------------------------------------- */

#define TRUE 1
#define FALSE 0
#define INVALID -1
#define NUM_EPSILON 0.0001
#define NUM_INF DBL_MAX 
     
/* definitions(type) ----------------------------------------------- */

/* data of grp problem */
typedef struct{
	int I;	/*	numbers of the customers */
	int J;	/*	numbers of the goods */
	int N;	/* numbers of reccomendation Goods 	*/
	double Cmax;
	double* P;
	double** cost;
	double* res;
} Grpdata;

/* various data often necessary during the search */
typedef struct{
  double start_time;  /* the time the search started */
  double finish_time;  /* the time the search finished */
  int **sol;  /* solution */
  double ub;
  double lb;
} Vdata;

/* local file */
void* malloc_e(size_t size);
void Subgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr); 
void k_select(double* cf,int* idx,int N,int J);
void swap_d(int i,int j,double* A);
void swap_i(int i,int j,int* A);

/* cpu_time.c */
double cpu_time();

/* functions ------------------------------------------------------- */


int main(int argc, char *argv[])
{
	Grpdata grpdata;
	Vdata vdata;
	grpdata.N = 5;
	grpdata.J = 310;
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
    	grpdata.cost[i] = (double*)malloc_e(sizeof(double)*grpdata.J);
    	for (int j = 0; j < grpdata.J; ++j)
    	{
    		fscanf(input_file,"%lf",&grpdata.cost[i][j]); 
    	}
    }
	fclose(input_file); 

/*	make Cmax & P[j]	*/
    double* r_check;
    r_check =(double*)malloc_e(sizeof(double)*grpdata.N);
    int r_sum;
    int ch_sum;
   	int a;
    srand((unsigned int) time(0));

/*	make Cmax	*/	
	grpdata.Cmax = 0.0;
    for (int i = 0; i < grpdata.I; ++i)
    {
	    r_sum = 0;
	    for (int j = 0; j < grpdata.N; ++j)
	    {
	    	r_check[j] = -1;
	    }
	    
		while(r_sum < grpdata.N)
	    {
	    	int ch_sum = 0;
	    	a = (int)((rand() / ((double) RAND_MAX + 1.0)) * grpdata.J);
	    	for (int k = 0; k < grpdata.N; ++k)
	    	{
	    		if (grpdata.cost[i][a] == r_check[k])
	    		{
			    	break;
	    		}
	    		else
	    		{
	    			++ch_sum;
	    		}
	    	}	

		   	if (ch_sum==grpdata.N)
		   	{
			    ++r_sum;
			    grpdata.Cmax += grpdata.cost[i][a];
			    r_check[r_sum-1] = grpdata.cost[i][a];
	    	}
	   	}
	}
    
    

/*	make P[J]*/
	double P_sum;
	grpdata.P = (double*)malloc_e(sizeof(double)*grpdata.J);
	for (int j = 0; j < grpdata.J; ++j)
	{
	    r_sum = 0;
	    P_sum = 0.0;
	    for (int i = 0; i < grpdata.N; ++i)
	    {
	    	r_check[i] = -1;
	    }

		while(r_sum < grpdata.N)
	    {
	    	ch_sum = 0;
	    	a = (int)((rand() / ((double) RAND_MAX + 1.0)) * grpdata.I);
	    	for (int k = 0; k < grpdata.N; ++k)
	    	{
	    		if (grpdata.res[a]*grpdata.cost[a][j] == r_check[k])
	    		{
			    	break;
	    		}

	    		else
	    		{
	    			++ch_sum;
	    		}
	    	}	

		   	if (ch_sum==grpdata.N)
		   	{
			    ++r_sum;
			    r_check[r_sum-1] = grpdata.res[a]*grpdata.cost[a][j];
			    P_sum += grpdata.res[a]*grpdata.cost[a][j];
	    	}

	   	}
	   	grpdata.P[j] = P_sum;
	}

/*	make sol	*/
	vdata.sol = (int**)malloc(sizeof(int*)*grpdata.I);
	for (int i = 0; i < grpdata.I; ++i)
	{
		vdata.sol[i] = (int*)malloc(sizeof(int)*grpdata.N);
	}

	printf("I = %d\n",grpdata.I);
	printf("res[0] = %lf\n",grpdata.res[0]);
	printf("cost[0][0] = %lf\n",grpdata.cost[0][0]);
	printf("res[I-1] = %lf\n",grpdata.res[grpdata.I-1]);
	printf("cost[I-1][J-1] = %lf\n",grpdata.cost[grpdata.I-1][grpdata.J-1]);
	printf("Cmax = %lf\n",grpdata.Cmax);
	printf("P[0] = %lf\n",grpdata.P[0]);
	printf("P[J-1] = %lf\n",grpdata.P[grpdata.J-1]);


/* subgrident method */
	vdata.start_time = cpu_time();
	Subgrad(&grpdata,&vdata);
	vdata.finish_time = cpu_time();
	
	printf("CPU time: %.3f seconds\n",vdata.finish_time - vdata.start_time); 
	return 0;
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

void Subgrad(Grpdata* grpdata_ptr,Vdata* vdata_ptr){
	int t = 0;
	double lmd = 0.0;
	double* mu;
	mu = (double*)malloc(sizeof(double)*grpdata_ptr->J);
	double s = 1;
	double** cf;
	cf = (double**)malloc(sizeof(double*)*grpdata_ptr->I);
	for (int i = 0; i < grpdata_ptr->I; ++i)
	{
		cf[i] = (double*)malloc(sizeof(double)*grpdata_ptr->J);
	}
	int* idx;
	idx = (int*)malloc(sizeof(int)*grpdata_ptr->J);
	double* sum_Gi;
	sum_Gi = (double*)malloc(sizeof(double)*grpdata_ptr->J);
	double* d;
	d = (double*)malloc(sizeof(double)*(grpdata_ptr->J + 1));
	double sum_C;
	double sum_G;
	double sum_muG;
	double d_nol;

	while(s > NUM_EPSILON)
	{
		for (int j = 0; j < grpdata_ptr->J; ++j)
		{
			sum_Gi[j] = 0.0;
		}
		sum_C = 0.0;
		sum_G = 0.0;
		sum_muG = 0.0;
		d_nol = 0.0;

		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			for (int j = 0; j < grpdata_ptr->J; ++j)
			{
				idx[j] = j;
				cf[i][j] = (1 + mu[j])*grpdata_ptr->res[i]*grpdata_ptr->cost[i][j] -lmd*grpdata_ptr->cost[i][j];
			}
		
			k_select(cf[i],idx,grpdata_ptr->N,grpdata_ptr->J);
			for (int n = 0;n < grpdata_ptr->N; ++n)
			{
				vdata_ptr->sol[i][n] = idx[n];
			}
		}

/*	calculation of ub & d	 */
		for (int i = 0; i < grpdata_ptr->I; ++i)
		{
			for (int n = 0;  n < grpdata_ptr->N; ++n)
			{
				sum_C += grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
				sum_G += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
				sum_Gi[vdata_ptr->sol[i][n]] += grpdata_ptr->res[i]*grpdata_ptr->cost[i][vdata_ptr->sol[i][n]];
			}
		}
		d[0] = grpdata_ptr->Cmax - sum_C;
		d_nol += d[0]*d[0];
		for (int j = 1; j < grpdata_ptr->J+1; ++j)
		{
			d[j] = sum_Gi[j] - grpdata_ptr->P[j];
			sum_muG += mu[j-1]*d[j];
			d_nol += d[j]*d[j];
		}
		vdata_ptr->ub = sum_G + lmd*d[0] + sum_muG;
		printf("ub:%lf\n",vdata_ptr->ub);
/*	recalculation of lmd,mu 	*/
		d_nol = (double)sqrt(d_nol);

		if (lmd-s*d[0] > 0.0)
		{
			lmd = lmd-s*d[0];
		}else{
			lmd = 0.0;
		}

		for (int j = 1; j < grpdata_ptr->J+1; ++j)
		{
			if(mu[j-1]-s*d[j]/d_nol > 0.0)
			{
				mu[j-1] = mu[j-1]-s*d[j]/d_nol;
			}else{
				mu[j-1]= 0.0;
			}
		}
		s = pow(0.9,t);
		t++;
	}
}

void k_select(double* cf,int* idx,int N,int J){
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
				swap_d(i,j,cf);
				swap_i(i,j,idx);
				i++; j--;
			}
		}
		if(i <= k) left = j + 1;
		if(k <= j) right = i - 1;
	}

	return;
}

void swap_d(int i,int j,double* A)
{
	int temp;
	temp = A[i];
	A[i] = A[j];
	A[j] = temp;
	return;
}

void swap_i(int i,int j,int* A)
{
	int temp;
	temp = A[i];
	A[i] = A[j];
	A[j] = temp;
	return;
}