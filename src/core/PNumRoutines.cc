/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "PNumRoutines.h"
#include <iostream>
using namespace std;

int ncom; //Global variables communicate with f1dim.
double *pcom,*xicom;
FunctFunctor *nrfunc;
FunctFunctor *Funct_Glob;
DerivFunctor *Deriv_Glob;
int N_Glob;

double my_f(const gsl_vector * x, void *params){
	double *v;
	v = (double*)calloc(N_Glob+1,sizeof(double));
	for (int i = 0; i < N_Glob; i++) v[i+1]=gsl_vector_get(x,i);
	double y0 = (*Funct_Glob)(v);
	free(v);
        return y0;
}
void my_df(const gsl_vector * x, void * p, gsl_vector * g){
	double *v;
	double *xi;
 	v = (double*)calloc(N_Glob+1,sizeof(double));
	xi = (double*)calloc(N_Glob+1,sizeof(double));
	for (int i = 0; i < N_Glob; i++) v[i+1]=gsl_vector_get(x,i);
	(*Deriv_Glob)(v,xi);
        for (int i = 0; i < N_Glob; i++)
              gsl_vector_set(g,i,xi[i+1]);
	free(v);
	free(xi);
}
void my_fdf (const gsl_vector * x, void * p, double * f, gsl_vector * g) {
	double *v;
	double *xi;
 	v = (double*)calloc(N_Glob+1,sizeof(double));
	xi = (double*)calloc(N_Glob+1,sizeof(double));
	for (int i = 0; i < N_Glob; i++) v[i+1]=gsl_vector_get(x,i);
	*f = (*Funct_Glob)(v);	
	(*Deriv_Glob)(v,xi);
        for (int i = 0; i < N_Glob; i++)
              gsl_vector_set(g,i,xi[i+1]);
	free(v);
	free(xi);
}


void PNumRoutines::nr_svd(double **a, int m, int n, double w[], double **v){
	gsl_matrix *a_g = gsl_matrix_alloc(n,m);
        for (int i = 0; i < n; i++)
         for (int j = 0; j < m; j++)
           gsl_matrix_set(a_g, i, j, a[j+1][i+1]);
        gsl_vector *S = gsl_vector_alloc(m);
        gsl_matrix *V = gsl_matrix_alloc(m,m);
        gsl_vector *work = gsl_vector_alloc(m);
        gsl_linalg_SV_decomp (a_g, V, S, work);
        for(int j = 0; j < m; j++) w[j+1] = gsl_vector_get(S,j);
	for(int j = m; j < n; j++) w[j+1] = 0.0;
        gsl_vector *tau = gsl_vector_alloc(m);//should be min(m,n) we assume m<n always
	gsl_linalg_QR_decomp (a_g,tau);
	gsl_matrix *Q = gsl_matrix_alloc(n,n);
	gsl_matrix *R = gsl_matrix_alloc(n,m);
	gsl_linalg_QR_unpack (a_g,tau,Q,R);
	 for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++)
                v[i+1][j+1] = gsl_matrix_get(Q,i,j);
	/*
	gsl_matrix *at = gsl_matrix_alloc(m,n);
	 for (int i = 0; i < m; i++)
         for (int j = 0; j < n; j++)
           gsl_matrix_set(at, i, j, a[i+1][j+1]);
	
	double temp;
	for (int i = 0; i < m;i++){
	for(int j=0;j<n;j++) {
	double temp=0.0;
	for(int k=0;k<n;k++) temp+=gsl_matrix_get(at,i,k)*gsl_matrix_get(Q,k,j);
		cout<<temp<<",";
	}
	cout<<"::"<<endl;
	}
	cout<<"DONE"<<endl;
	gsl_matrix_free(at);
	*/	
        gsl_matrix_free (V);
        gsl_vector_free(S);
        gsl_vector_free (work);
        gsl_matrix_free(a_g);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);
	gsl_vector_free(tau);
}


void PNumRoutines::nr_multimin( double p[], int n, double ftol, int *iter, double *fret, FunctFunctor *Funct, DerivFunctor *Deriv){
	N_Glob = n;
	int iter_g = 0;
	int status;
	Funct_Glob = Funct;
	Deriv_Glob = Deriv;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	// Position of the minimum (1,2).
	double par[2] = { 1.0, 2.0 };

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;
	
	my_func.f = &my_f;
	my_func.df = &my_df;
	my_func.fdf = &my_fdf;
	my_func.n = n;
	//my_func.params = &par;
	
	// Starting point, x = (5,7)
	x = gsl_vector_alloc (n);
	for (int i = 0; i < n; i++) gsl_vector_set(x, i, p[i+1]);

	//T = gsl_multimin_fdfminimizer_conjugate_pr;
	T = gsl_multimin_fdfminimizer_steepest_descent;
//	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, n);
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.1, 0.5);//with SDescent
//	gsl_multimin_fdfminimizer_set(s,&my_func,x,0.0001,0.0001);
	
	do
		{
		iter_g++;
		double fp=s->f;
		status = gsl_multimin_fdfminimizer_iterate (s);
	
		if (status)
		break;
	
		status = gsl_multimin_test_gradient (s->gradient, 1e-2);
		double fpnew=s->f;
		 if (2.0*fabs(fpnew-fp) <= ftol*(fabs(fpnew)+fabs(fp)+1E-10)){
			status = GSL_SUCCESS;
		}	
		
		if (status == GSL_SUCCESS)
		cout<<"Minimum value is:"<<s->f<<endl;
	
		}
	while (status == GSL_CONTINUE && iter_g < 200);
	cout<<gsl_strerror (status)<<endl;
	cout<<"iter_g:"<<iter_g<<endl;
        for (int i = 0; i < n; i++) p[i+1]=gsl_vector_get(s->x, i);
	*iter = iter_g;
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
}

void PNumRoutines::nr_inverse(double **a, int N, double **y){
	gsl_matrix *a_g = gsl_matrix_alloc(N,N);
	for (int i = 0; i < N; i++)
         for (int j = 0; j < N; j++)
           gsl_matrix_set(a_g, i, j, a[i+1][j+1]);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(N);
	gsl_matrix *a_g_inv = gsl_matrix_alloc(N,N);
	gsl_linalg_LU_decomp(a_g, p, &s);
	int status = gsl_linalg_LU_invert(a_g, p,a_g_inv);
	cout << status << endl;
	for (int i = 0; i < N; i++)
         for (int j = 0; j < N; j++)
           y[i+1][j+1] = gsl_matrix_get(a_g_inv,i,j);
	gsl_matrix_free(a_g);
	gsl_matrix_free(a_g_inv);
	gsl_permutation_free(p);

}

void PNumRoutines::nr_inverse(double **a, int N, double **y, int& status){
	gsl_matrix *a_g = gsl_matrix_alloc(N,N);
	for (int i = 0; i < N; i++)
         for (int j = 0; j < N; j++)
           gsl_matrix_set(a_g, i, j, a[i+1][j+1]);
	int s;
	gsl_permutation * p = gsl_permutation_alloc(N);
	gsl_matrix *a_g_inv = gsl_matrix_alloc(N,N);
	gsl_linalg_LU_decomp(a_g, p, &s);
	status = gsl_linalg_LU_invert(a_g, p,a_g_inv);
//	cout << status << endl;
	if( status == -1)
		return;
	for (int i = 0; i < N; i++)
         for (int j = 0; j < N; j++)
           y[i+1][j+1] = gsl_matrix_get(a_g_inv,i,j);
	gsl_matrix_free(a_g);
	gsl_matrix_free(a_g_inv);
	gsl_permutation_free(p);

}
