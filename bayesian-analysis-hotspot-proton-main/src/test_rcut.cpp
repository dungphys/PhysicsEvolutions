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
#include <gsl/gsl_spline.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <vector> // For std::vector
#include <regex>

using namespace std;

int main(int argc, char **argv)
{
	gsl_set_error_handler(&ErrHandler);
	
	model_params params;
	xsections dsigma;
	norm_factors norm;
	
	
	double Q2 = 0.1;
	double W = 78;
	double gev2_to_nb = 0.389379e6;
	
	/* read rcut*/
	vector<double> vecX, vecY;
	double x, y;
	char iname[32];
	sprintf(iname, "rmax_Q2=%.2lf.txt", Q2);
	std::ifstream infile(iname);
	
	while (infile >> x >> y)
	{
    		vecX.push_back(x);
    		vecY.push_back(y);
	}
	int nt = vecX.size();
	double interp_t[nt], interp_rmax[nt];
    
    // Using copy() function to copy elements
    copy(vecX.begin(),vecX.end(),interp_t);
    copy(vecY.begin(),vecY.end(),interp_rmax);
	
	gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    	gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, nt);
	gsl_spline_init (spline, interp_t, interp_rmax, nt);
	
	int code_Nh[2] = {0,1};
	double sigma_Qs[2] = {0, 0.6};
	
	A = 0.213; B = -0.0157;
	MASS_QUARK = mC;
	CHARGE_QUARK = 2./3;
	
	params.p0 = 3.0;
	params.p1 = 0;
	params.p2 = 0;
	params.lambda_H = 0;
	params.lambda_p = 0;
	
	params.R0 = 2.12;
	params.rH0 = 0.54;
	params.gmu0 = 3.34;
	params.m2 = 0.044;
	
	
	double dt = 0.1, t0 = 0.02, t;
	double coh, ic;
	ofstream outfile0, outfile1;
	std::string name;
	
	
	for (int i=0;i<2;i++) {
		params.sigma_Qs = sigma_Qs[i];
		for (int j=0;j<2;j++) {	
			name = "dsigma_test_w_rcut_";
			name += std::to_string(i);
			name += "_";
			name += std::to_string(j);
			name += ".txt";
			outfile0.open(name);
			
			name = "dsigma_test_wo_rcut_";
			name += std::to_string(i);
			name += "_";
			name += std::to_string(j);
			name += ".txt";
			outfile1.open(name);
			
			cout << "Files opened, run for different t" << endl;
			
			for (int k=0;k<21;k++) {
				t = t0 + k*dt;
				cout << "t = " << t << endl;
				RCUT = 299; //GeV-1 
				coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh[j], Q2, W, -t, params);
				cout << "coherent nocut done! " << coh << endl;
				ic = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh[j], Q2, W, -t, params);
				cout << "incoherent nocut done! " << ic-coh << endl;
				outfile0 << t << "\t" << coh  << "\t" << ic - coh << endl;
				
				RCUT = gsl_spline_eval (spline, t, acc); //GeV-1
				cout << RCUT << endl;
				coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh[j], Q2, W, -t, params);
				cout << "coherent done! " << coh << endl;
				ic = (gev2_to_nb)*dsigma.dsigma_ci_rcut(code_Nh[j], Q2, W, -t, params);
				cout << "incoherent done! " << ic-coh << endl;
				outfile1 << t << "\t" << coh  << "\t" <<  ic - coh << endl;
			}
			outfile0.close();
			outfile1.close();
		}
	}	
	
	gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
	return 0;
}
