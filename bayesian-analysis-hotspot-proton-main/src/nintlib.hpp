#ifndef nintlib_hpp
#define nintlib_hpp
#include <string>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

void dec2base(int b, int x, int res[], int ndigit);
void gauss_legendre_nd ( double func ( double x[], size_t dim_num, void *params ), size_t dim_num, void *params, double lo[], double hi[], 
  int order, double xtab[], double weight[], int *eval_num, double *res );
void MC_nd (double func ( double x[], size_t dim_num, void *params ), size_t dim_num, void *params, double lo[], double hi[], 
	std::string algorithm, int ncalls, double *res);

#endif
