/**/
#define _USE_MATH_DEFINES

#include "root.hpp"
#include "aux.hpp"
#include "cross_sections.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_spline.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector> // For std::vector
#include <fstream>
#include <limits>

using namespace std;

void rmax_read(char *name, vector<double> &vecX, vector<double> &vecY)
{
	double x, y;
	std::ifstream rmax_infile;
	rmax_infile.open(name);
	if ( rmax_infile.is_open() ) { // always check whether the file is open
		cout << "rcut file opened successfully!" << endl;
	}
	else {
		cout << "Cannot open rcut file! Check again !" << endl;
		exit(1);
	}
	
	while (rmax_infile >> x >> y)
	{
    		vecX.push_back(x);
    		vecY.push_back(y);
	}
	rmax_infile.close();
}

int main(int argc, char* argv[])
{
	if (argc!=7) {
		cout << "Syntax should be: ./run infile outfile scenario wf Nh Qs" << endl;
		exit(1);
	}
	
	gsl_set_error_handler(&ErrHandler);
	
	model_params params;
	xsections dsigma;
	//norm_factors norm;
	
	/* Random number init */
	gsl_rng_env_setup();
    	struct timeval tv; // Seed generation based on time
    	gettimeofday(&tv,0);
    	unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    	T_rng = gsl_rng_default; // Generator setup
    	r_rng = gsl_rng_alloc (T_rng);
    	gsl_rng_set(r_rng, mySeed);
    	
    	/* code Nh */
    	int code_Nh;
    	if (strcmp(argv[5],"with_Nh_fluct")==0) code_Nh = 1;
    	else code_Nh = 0;
    	
    	/* J/psi wf */
    	if (strcmp(argv[4],"rel")==0) {A = 0.213; B = -0.0157;}
    	else {A = 0.211; B = 0;}
    	MASS_QUARK = mC;
	CHARGE_QUARK = eC;
    	double gev2_to_nb = 0.389379e6;
	
	/*init params*/
     params.sigma_Qs = 0;
     params.p0 = 3.; params.p1 = 0.; params.p2 = 0.;
     params.lambda_p = 0.; params.lambda_H = 0.;
     	
    	std::string line;
     std::ifstream infile;
     std::string col1;
     double col2;
	infile.open(argv[1]);
	if ( infile.is_open() ) { // always check whether the file is open
		cout << "infile opened successfully!" << endl;
	}
	else {
		cout << "Cannot open infile! Check again !" << endl;
		exit(1);
	}
     	
     while (getline(infile, line)) {
		std::stringstream ss(line);
		ss >> col1 >> col2;
		if ( col1.compare("gmu0")==0 )
			params.gmu0 = col2;
		else if ( col1.compare("m2")==0 )
			params.m2 = col2;
		else if ( col1.compare("sigma_Qs")==0 )
			params.sigma_Qs = col2;
		else if ( col1.compare("R0")==0 )
			params.R0 = col2;
		else if ( col1.compare("lambda_p")==0 )
			params.lambda_p = col2;
		else if ( col1.compare("lambda_H")==0 )
			params.lambda_H = col2;
		else if ( col1.compare("rH0")==0 )
			params.rH0 = col2;
		else if ( col1.compare("p0")==0 )
			params.p0 = col2;
		else if ( col1.compare("p1")==0 )
			params.p1 = col2;
		else if ( col1.compare("p2")==0 )
			params.p2 = col2;	
     }
     cout << "Params initialization done!" << endl;	
     
     /* open datafile to read */
	ifstream datafile_coh, datafile_incoh;
	datafile_coh.open("/home/anh-dung/projects/current_projects/hot_spot_NR/bayes/bin/expr_data/HERA_coherent_data.txt");
	datafile_incoh.open("/home/anh-dung/projects/current_projects/hot_spot_NR/bayes/bin/expr_data/HERA_incoherent_data.txt");
	if ( datafile_coh.is_open() && datafile_coh.is_open() ) { // always check whether the file is open
		cout << "Data files opened successfully!" << endl;
	}
	else {
		cout << "Cannot open data file! Check again !" << endl;
		exit(1);
	}
	
	
	getline(datafile_coh, line);	/* skip the first line containing labels*/
	getline(datafile_incoh, line);	/* skip the first line containing labels*/
	
	/* First scenario: No evolution */
	if ( (strcmp(argv[3],"no_evolution_W78")==0) || (strcmp(argv[3],"no_evolution_W55")==0) ) {
		double Q2 = 0.1, W = 78;
		std::string s = argv[3];
		int p1 = s.find("W");
		int p2 = s.length();
		W = stod(s.substr(p1 + 1, p2 - p1 - 1));
		cout << "W = " << W << endl;
		/**/
		/* read rcut*/
		vector<double> vec_t, vec_rmax;
		char iname[120];
		sprintf(iname, "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/rmax_Q2=%.2lf.txt", Q2);
		rmax_read(iname, vec_t, vec_rmax);
		int nt = vec_t.size();
		double interp_t[nt], interp_rmax[nt];
    
    		copy(vec_t.begin(),vec_t.end(),interp_t);
    		copy(vec_rmax.begin(),vec_rmax.end(),interp_rmax);
	
		gsl_interp_accel *acc
      		= gsl_interp_accel_alloc ();
    		gsl_spline *spline
      		= gsl_spline_alloc (gsl_interp_cspline, nt);
		gsl_spline_init (spline, interp_t, interp_rmax, nt);
		/* read t arrays for running from dataset at fixed W=78 and Q2=0.1 */
		vector<double> t_coh, t_incoh;
		double col_Q2, col_W, col_t;
		while (getline(datafile_coh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			if ((col_Q2==Q2)&&(col_W==W))
				t_coh.push_back(col_t);
     	}
     	/**/
     	while (getline(datafile_incoh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			if ((col_Q2==Q2)&&(col_W==W))
				t_incoh.push_back(col_t);
     	}
		
		ofstream outfile;
		outfile.open(argv[2]);
		outfile << "Process\tQ2\tW\tt\tdsigma/dt" << endl;
		double coh, ci;
		/* run coherent */
		for (int i=0; i<(int)t_coh.size(); i++) {
			cout << "Coherent, t = " << t_coh[i] << endl;
			if (strcmp(argv[4],"rel")==0)
				RCUT = gsl_spline_eval (spline, t_coh[i], acc);
			else RCUT = 399;
			coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params);
			outfile << "Coherent" << "\t" << Q2 << "\t" << W << "\t" << t_coh[i] << "\t" << coh << endl;
		}
		/* run incoherent */
		for (int i=0; i<(int)t_incoh.size(); i++) {
			cout << "Incoherent, t = " << t_incoh[i] << endl;
			if (strcmp(argv[4],"rel")==0) {
				RCUT = gsl_spline_eval (spline, t_incoh[i], acc);
				coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params);
				ci = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_incoh[i], params);
			}
			else {
				RCUT = 399;
				coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params);
				ci = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_incoh[i], params);
			}
			outfile << "Incoherent" << "\t" << Q2 << "\t" << W << "\t" << t_incoh[i] << "\t" << ci - coh << endl;
		}		
		outfile.close();
		gsl_spline_free (spline);
    		gsl_interp_accel_free (acc);
     }	
	/* Second scenario: With evolution*/
	else {
		vector<double> Q2_coh, Q2_incoh;
		vector<double> W_coh, W_incoh;
		vector<double> t_coh, t_incoh;
		double col_Q2, col_W, col_t;
		while (getline(datafile_coh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			Q2_coh.push_back(col_Q2);
			W_coh.push_back(col_W);
			t_coh.push_back(col_t);
     	}
     	/**/
     	while (getline(datafile_incoh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			Q2_incoh.push_back(col_Q2);
			W_incoh.push_back(col_W);
			t_incoh.push_back(col_t);
     	}
     	char iname[120];
		ofstream outfile;
		outfile.open(argv[2]);
		outfile << "Process\tQ2\tW\tt\tdsigma/dt" << endl;
		double coh, ci;
		double *interp_t, *interp_rmax;
		vector<double> vec_t, vec_rmax;
		gsl_interp_accel *acc
      			= gsl_interp_accel_alloc ();
      	gsl_spline *spline;
		/* run coherent */
		for (int i=0; i<(int)t_coh.size(); i++) {
			/* rcut */
			sprintf(iname, "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/rmax_Q2=%.2lf.txt", Q2_coh[i]);
			rmax_read(iname, vec_t, vec_rmax);
			interp_t = (double*) malloc(vec_t.size()*sizeof(double));
			interp_rmax = (double*) malloc(vec_t.size()*sizeof(double));
			copy(vec_t.begin(),vec_t.end(),interp_t);
    			copy(vec_rmax.begin(),vec_rmax.end(),interp_rmax);
    			
      		spline = gsl_spline_alloc (gsl_interp_cspline, vec_t.size());
			gsl_spline_init (spline, interp_t, interp_rmax, vec_t.size());			
			RCUT = gsl_spline_eval (spline, t_coh[i], acc);
			/* run coh dsigma*/
			coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_coh[i], W_coh[i], -t_coh[i], params);
			outfile << "Coherent" << "\t" << Q2_coh[i] << "\t" << W_coh[i] << "\t" << t_coh[i] << "\t" << coh << endl;
		}
		/* run incoherent */
		for (int i=0; i<(int)t_incoh.size(); i++) {
			/* rcut */
			sprintf(iname, "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/rmax_Q2=%.2lf.txt", Q2_incoh[i]);
			rmax_read(iname, vec_t, vec_rmax);
			interp_t = (double*) malloc(vec_t.size()*sizeof(double));
			interp_rmax = (double*) malloc(vec_t.size()*sizeof(double));
			copy(vec_t.begin(),vec_t.end(),interp_t);
    			copy(vec_rmax.begin(),vec_rmax.end(),interp_rmax);
      		spline = gsl_spline_alloc (gsl_interp_cspline, vec_t.size());
			gsl_spline_init (spline, interp_t, interp_rmax, vec_t.size());			
			RCUT = gsl_spline_eval (spline, t_coh[i], acc);
			/* run coh dsigma*/
			coh = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_incoh[i], W_incoh[i], -t_incoh[i], params);
			ci = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2_incoh[i], W_incoh[i], -t_incoh[i], params);
			outfile << "Incoherent" << "\t" << Q2_incoh[i] << "\t" << W_incoh[i] << "\t" << t_incoh[i] << "\t" << ci - coh << endl;
		}		
		
	}
	
	datafile_coh.close();
	datafile_incoh.close();
	gsl_rng_free(r_rng);
    	
	return 0;
}
