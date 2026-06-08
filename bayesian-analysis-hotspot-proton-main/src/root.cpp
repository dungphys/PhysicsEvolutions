
#include "root.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;

double function(double x, void *params)
{
	struct root_params *p = (struct root_params *) params;

  	double Nh = p->Nh;
  	return ((x/(1-std::exp(-x))) - Nh);
}

double function_deriv (double x, void *params)
{
	return ((1-std::exp(-x)-x*std::exp(-x))/SQR(1-std::exp(-x)));
}

double root_solver(double Nh_mean)
{
	int status;
	int iter = 0, max_iter = 150;
	double r = Nh_mean;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	
	double x_lo = 0.01, x_hi = 2*Nh_mean+3;
  	gsl_function F;
  	struct root_params params = {Nh_mean};
  	
  	F.function = &function;
  	F.params = &params;

  	T = gsl_root_fsolver_brent;
  	s = gsl_root_fsolver_alloc (T);
  	gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  	
  	do
    	{
      	iter++;
      	status = gsl_root_fsolver_iterate (s);
      	r = gsl_root_fsolver_root (s);
      	x_lo = gsl_root_fsolver_x_lower (s);
      	x_hi = gsl_root_fsolver_x_upper (s);
      	status = gsl_root_test_interval (x_lo, x_hi, 0, 0.0001);

      	/*if (status == GSL_SUCCESS) { 
      		printf ("Converged:\n");
      		printf ("lambda = %.7f\n",r);
      	}*/
    	}
  	while (status == GSL_CONTINUE && iter < max_iter);

  	gsl_root_fsolver_free (s);
  	return r;
}
