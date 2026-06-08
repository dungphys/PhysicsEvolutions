
#ifndef root_hpp
#define root_hpp

#include <gsl/gsl_roots.h>

inline double SQR(double x) { return (x*x); }		/**/
inline double CUBE(double x) { return (x*x*x); }		/**/

struct root_params
{
    double Nh;
};

double function (double x, void *params);
double function_deriv (double x, void *params);
//void function_fdf (double x, void *params,
//                    double *y, double *dy);
double root_solver(double Nh_mean);



#endif /*root_hpp*/
