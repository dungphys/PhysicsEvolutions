/**/
#define _USE_MATH_DEFINES

#include "nintlib.hpp"
#include "root.hpp"
#include "aux.hpp"
#include "cross_sections.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <sys/time.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <vector> // For std::vector
#include <regex>

int n_coh = 8, n_icoh=8;
double t_coh[8] = {0.02, 0.08, 0.14, 0.21, 0.30, 0.41, 0.58, 0.90};
double t_icoh[8] = {0.10, 0.29, 0.52, 0.78, 1.12, 1.55, 2.21, 3.30};

using namespace std;

int main(int argc, char **argv)
{
	gsl_set_error_handler(&ErrHandler);
	
	model_params params;
	xsections dsigma;
	norm_factors norm;
	
	params.p1 = 0;
	params.p2 = 0;
	params.lambda_rH = 0;
	params.lambda_R = 0;
	params.sigma_Qs = 0;
	
	double Q2 = 0.1;
	double W = 78;
	double gev2_to_nb = 0.389379e6;
	
	//ofstream outfile;	
	/* norm factor */
	
	//double ct, cl, pol_factor, gev2_to_nb = 0.389379e6;
	//double t, C, DC;
	double aq;
	
	aq = norm.AQ('c');
	cout << aq << endl;
	
	
	
	//outfile.close();	
	return 0;
}
