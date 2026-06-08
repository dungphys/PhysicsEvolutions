/*
	Multidimensional integration routines
*/

#include "nintlib.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <string>
#include <vector>
#include <sys/time.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <functional>

using namespace std;
using namespace std::placeholders;

/*****************************************************************************/
void dec2base(int b, int x, int res[], int ndigit)
/* Converse a decimal integer x to number in base b<10
   Results is stored in res[ ]	
*/
{
	if (x==0) {
		for (int i=0;i<ndigit;i++) res[i] = 0;
	}
	else {
		int r, xcopy=x;
		for (int i=ndigit-1;i>=0;i--) {
			r = xcopy%b;
			xcopy = xcopy/b;
			res[i] = r;
		}
	}
	//for (int i=0;i<ndigit;i++) cout << res[i];
	//cout << endl;
}

/*****************************************************************************/
void gauss_legendre_nd ( double func ( double x[], size_t dim_num, void *params ), size_t dim_num, void *params, double lo[], double hi[],
  int order, double xtab[], double weight[], int *eval_num, double *res)
/*
	This routine performs the multidimensional integration using the fixed-point Gaussian-Legendre quadrature
	Input:
		dim_num: number of dimensions
		func(): the function to be integrated
		lo[]: lower limits
		hi[]: upper limits
		order: order of the Gaussian-Legendre quadrature
		xtab[]: abscissae of the corresponding order
		weight[]: weights of the corresponding abscissae
	Output:
		res: the numerical estimation 
		eval_num: number of function evaluations
*/
{
	/*------ check dimension and quadrature order ------*/
	if (dim_num<1) {
		cerr << "G-L quadrature error!" << endl;
		cerr << "	Integration dimension (dim_num) should be at least 1:" << endl;
		exit(1);
	}
	/**/
	if (order<1) {
		cerr << "G-L quadrature error!" << endl;
		cerr << "	The quadrature order should be at least 1:" << endl;
		exit(1);
	}
	
	int dim, i;
	int index_num = (int) pow(order, (int) dim_num);
	
	/*array of index for loops*/
	int **index_array = new int*[index_num];
	for (i = 0; i < index_num; ++i) {
		index_array[i] = new int[(int) dim_num];
		dec2base(order, i, index_array[i], dim_num);
	}
	/* weight and function blocks*/
	double weight_block, func_block, normalization = 1;
	double *x_eval = new double[(int)dim_num];
	*eval_num = 0; *res = 0;
	for (i = 0; i < index_num; ++i) {
		weight_block = 1;
		for ( dim=0;dim<(int)dim_num;dim++) {
			x_eval[dim] = 0.5*(hi[dim]-lo[dim])*xtab[index_array[i][dim]] + 0.5*(hi[dim]+lo[dim]); 
			weight_block = weight_block*weight[index_array[i][dim]];
		}
		func_block = func(x_eval, dim_num, params);
		*eval_num = (*eval_num) + 1;
		*res = (*res) + (func_block*weight_block);
	}
	
	/*volume factor*/
	for (dim=0;dim<(int)dim_num;dim++) normalization = normalization*(0.5*(hi[dim]-lo[dim]));
	*res = normalization*(*res);
	
	/* free alllocating memory*/
	for(int i = 0; i < index_num; ++i) delete [] index_array[i];
	delete [] index_array; 
	delete [] x_eval;
}

void MC_nd (double func ( double x[], size_t dim_num, void *params ), size_t dim_num,  void *params, double lo[], double hi[], 
	std::string algorithm, int ncalls, double *res)
/*
	Integration routine using the MC technique with GSL
	Possible algorithms:
	- algorithm = "vegas": gsl_monte_vegas
	- algorithm = "plain": gsl_monte_plain
	- algorithm = "miser": gsl_monte_miser
	Note: Require GSL library to be installed 
*/
{
	double err;
	const gsl_rng_type *T_rng;
	gsl_rng *r_rng;
	gsl_rng_env_setup ();
	
	T_rng = gsl_rng_default;
  	r_rng = gsl_rng_alloc (T_rng);
  	
  	struct timeval tv; // Seed generation based on time
    	gettimeofday(&tv,0);
    	unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    	gsl_rng_set(r_rng, mySeed);
	
	gsl_monte_function G;
	G.f = func;
	G.dim = dim_num;
	G.params = params;
	
	if ( algorithm.compare("plain")==0 ) {
		gsl_monte_plain_state *s = gsl_monte_plain_alloc (dim_num);
    		gsl_monte_plain_integrate (&G, lo, hi, dim_num, ncalls, r_rng, s,
                               res, &err);
    		gsl_monte_plain_free (s);
	}
	else if ( algorithm.compare("miser") ==0 ) {
		gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim_num);
		gsl_monte_miser_integrate (&G, lo, hi, dim_num, ncalls, r_rng, s,
                               res, &err);
    		gsl_monte_miser_free (s);
	}
	else if ( algorithm.compare("vegas")==0 ) {
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_num);
    		gsl_monte_vegas_integrate (&G, lo, hi, dim_num, 200000, r_rng, s,
                               res, &err);
    		//cout << "vegas warm-up" << *res << " " << err << endl;
    		do {
        		gsl_monte_vegas_integrate (&G, lo, hi, dim_num, ncalls/5, r_rng, s,
                                   res, &err);
      	}
    		while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    		gsl_monte_vegas_free (s);
    		
	}
	gsl_rng_free (r_rng);
}

