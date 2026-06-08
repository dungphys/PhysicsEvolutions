#ifndef INTERP_H
#define INTERP_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include "utils.h"

/* Spline interpolation */
extern gsl_interp_accel *acc;
extern gsl_spline *spline;
void initSpline(int nr);
void freeSpline();
double csplineInterpS(double r, std::vector<double> S_in, Lattice &lattice);

/* Linear interpolation */
double linearInterpS(double r, std::vector<double> S_in, Lattice &lattice);

#endif // INTERP_H