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


int main(int argc, char **argv)
{
	if (argc!=6) {
		cout << "Syntax should be: ./run scenario wf Nh Qs" << endl;
		exit(1);
	}
	
	gsl_set_error_handler(&ErrHandler);
	
	model_params params_MAP, params_MED, params_add;
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
    	if (strcmp(argv[3],"with_Nh_fluct")==0) code_Nh = 1;
    	else code_Nh = 0;
    	
    	/* J/psi wf */
    	if (strcmp(argv[2],"rel")==0) {A = 0.213; B = -0.0157;}
    	else {A = 0.211; B = 0;}
    	MASS_QUARK = mC;
	CHARGE_QUARK = eC;
    	double gev2_to_nb = 0.389379e6;
	
	/*init params*/
     params_MAP.sigma_Qs = 0;
     params_MAP.p0 = 3.; params_MED.p1 = 0.; params_MED.p2 = 0.;
     params_MAP.lambda_p = 0.; params_MAP.lambda_H = 0.;
     
     params_MED.sigma_Qs = 0;
     params_MED.p0 = 3.; params_MAP.p1 = 0.; params_MED.p2 = 0.;
     params_MED.lambda_p = 0.; params_MED.lambda_H = 0.;
     
     std::string line;
     std::ifstream infile;
     std::string col1;
     double col2, col3, col4, col5, col6, col7;
     char iname[150];
     sprintf(iname,"/home/anh-dung/projects/current_projects/hot_spot_NR/bayes/params_bayes_output_%s_%s_%s_%s_%s.txt",argv[1],argv[3],argv[4],argv[2],argv[5]);
     
	infile.open(iname);
	if ( infile.is_open() ) { // always check whether the file is open
		cout << "infile opened successfully!" << endl;
	}
	else {
		cout << "Cannot open infile! Check again !" << endl;
		exit(1);
	}
	
     getline(infile, line); //skip 1st line	
     while (getline(infile, line)) {
		std::stringstream ss(line);
		ss >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7;
		cout << col1 << endl;
		if ( col1.compare("gmu0")==0 ) {
			params_MAP.gmu0 = col2;
			params_MED.gmu0 = col5;
			cout << "gmu0 extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("m2")==0 ) {
			params_MAP.m2 = col2;
			params_MED.m2 = col5;
			cout << "m2 extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("sigma_Qs")==0 ){
			params_MAP.sigma_Qs = col2;
			params_MED.sigma_Qs = col5;
			cout << "sigma_s extracted" << endl;
			cout << col2 << endl;
		}
		else if ( ( col1.compare("R0")==0 ) ) {
			params_MAP.R0 = col2;
			params_MED.R0 = col5;
			cout << "R0 extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("lambda_p")==0 )	{
			params_MAP.lambda_p = col2;
			params_MED.lambda_p = col5;
			cout << "l_p extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("lambda_H")==0 ) {
			params_MAP.lambda_H = col2;
			params_MED.lambda_H = col5;
			cout << "l_H extracted" << endl;
			cout << col2 << endl;
		}
		else if ( ( col1.compare("rH0")==0) ) {
			params_MAP.rH0 = col2;
			params_MED.rH0 = col5;
			cout << "rH0 extracted" << endl;
			cout << col2 << endl;
		}
		else if (( col1.compare("p0")==0 ) ) {
			params_MAP.p0 = col2;
			params_MED.p0 = col5;
			cout << "p0 extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("p1")==0 ) {
			params_MAP.p1 = col2;
			params_MED.p1 = col5;
			cout << "p1 extracted" << endl;
			cout << col2 << endl;
		}
		else if ( col1.compare("p2")==0 ) {
			params_MAP.p2 = col2;
			params_MED.p2 = col5;
			cout << "p2 extracted" << endl;
			cout << col2 << endl;
		}
     }
     params_add.sigma_Qs = params_MAP.sigma_Qs;
     params_add.gmu0 = params_MAP.gmu0;
     params_add.m2 = params_MAP.m2;
     params_add.rH0 = params_MAP.rH0;
     params_add.R0 = params_MAP.R0;
     params_add.lambda_p = params_MAP.lambda_p;
     params_add.lambda_H = params_MAP.lambda_H;
     params_add.p0 = params_MAP.p0;
     params_add.p1 = params_MAP.p1;
     params_add.p2 = params_MAP.p2;  
     if (strcmp(argv[2],"rel")==0)
     	{params_add.gmu0 = 4.306; params_add.R0 = 2.048;
     	 params_add.m2 = 0.069; params_add.rH0 = 0.578;
     	 params_add.p0 = 5.595; params_add.sigma_Qs = 0.571;}
     if (strcmp(argv[2],"nonrel")==0)
     	{params_add.gmu0 = 2.946; params_add.R0 = 2.382;
     	 params_add.m2 = 0.063; params_add.rH0 = 0.002;
     	 params_add.p0 = 3.303; params_add.sigma_Qs = 0;}
     if (strcmp(argv[2],"nonrel")==0)
     	{params_MAP.gmu0 = 2.946; params_MAP.R0 = 2.382;
     	 params_MAP.m2 = 0.063; params_MAP.rH0 = 0.2;
     	 params_MAP.p0 = 3.303; params_MAP.sigma_Qs = 0;}
     cout << "Params initialization done!" << endl;
	 /* open datafile to read */
	 
	 double col_Q2, col_W, col_t;
	 vector<double> t_coh, t_incoh;
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
	
	if ( (strcmp(argv[1],"no_evolution_W78")==0) || (strcmp(argv[1],"no_evolution_W55")==0) ) {
		double Q2 = 0.1, W = 78;
		std::string s = argv[1];
		int p1 = s.find("W");
		int p2 = s.length();
		W = std::stod(s.substr(p1 + 1, p2 - p1 - 1));
		cout << "W = " << W << endl;
		/**/
		while (getline(datafile_coh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			if ((col_Q2 == Q2) && (col_W==W)) 
				t_coh.push_back(col_t);
     	}
     	/**/
     	while (getline(datafile_incoh, line)) {
			std::stringstream ss(line);
			ss >> col_Q2 >> col_W >> col_t;
			if ((col_Q2 == Q2) && (col_W==W)) 
				t_incoh.push_back(col_t);
     	}
		
		/* read rcut*/
		vector<double> vec_t, vec_rmax;
		char iname[150];
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
		
		double t, dt = 0.04, tmin = 0.001;
		
     	ofstream outfile_coh, outfile_incoh;
     	ofstream outfile_chi2_coh, outfile_chi2_incoh;
     	char oname_coh[150], oname_incoh[150];
     	char for_chi2_coh[150], for_chi2_incoh[150];  
     	sprintf(oname_coh,"dsigma_bayes_fit_%s_%s_%s_%s_coh_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5]);
     	sprintf(oname_incoh,"dsigma_bayes_fit_%s_%s_%s_%s_incoh_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5]);
		outfile_coh.open(oname_coh);
		outfile_incoh.open(oname_incoh);
		/**/
		sprintf(for_chi2_coh,"dsigma_chi2_bayes_%s_%s_%s_%s_coh_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5]);
     	sprintf(for_chi2_incoh,"dsigma_chi2_bayes_%s_%s_%s_%s_incoh_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5]);
		outfile_chi2_coh.open(for_chi2_coh);
		outfile_chi2_incoh.open(for_chi2_incoh);
		outfile_chi2_coh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)\tdsigma/dt(add)" << endl;
		outfile_chi2_incoh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)\tdsigma/dt(add)" << endl;
		outfile_coh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)\tdsigma/dt(add)" << endl;
		outfile_incoh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)\tdsigma/dt(add)" << endl;
		double coh_MAP, ci_MAP, coh_MED, ci_MED, ci_add, coh_add;
		for (int i=0; i<200; i++) {
			//cout << "Incoherent, t = " << t_incoh[i] << endl;
			t = tmin + i*dt;
			cout << "t = " << t << endl;
			if (strcmp(argv[2],"rel")==0) {
				RCUT = gsl_spline_eval (spline, t, acc);
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t, params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t, params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t, params_MED);
				outfile_coh << Q2 << "\t" << W << "\t" << t << "\t" << coh_MAP << "\t" << coh_MED << "\t" << coh_add << endl;
				outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED << "\t" << ci_add - coh_add <<  endl;
			}
			else {
				RCUT = 399;
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t, params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t, params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t, params_MED);
				outfile_coh << Q2 << "\t" << W << "\t" << t << "\t" << coh_MAP << "\t" << coh_MED << "\t" << coh_add << endl;
				outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED << "\t" << ci_add - coh_add << endl;
			}
		}
		/**/
		cout << "For chi2:" << endl;
		for (int i=0; i<(int)t_coh.size(); i++) {
			//cout << "Incoherent, t = " << t_incoh[i] << endl;
			cout << "t = " << t_coh[i] << endl;
			if (strcmp(argv[2],"rel")==0) {
				RCUT = gsl_spline_eval (spline, t_coh[i], acc);
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_coh[i], params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_coh[i], params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_coh[i], params_MED);
				outfile_chi2_coh << Q2 << "\t" << W << "\t" << t_coh[i] << "\t" << coh_MAP << "\t" << coh_MED << "\t" << coh_add << endl;
				//outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED <<  endl;
			}
			else {
				RCUT = 399;
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_coh[i], params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_coh[i], params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_coh[i], params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_coh[i], params_MED);
				outfile_chi2_coh << Q2 << "\t" << W << "\t" << t_coh[i] << "\t" << coh_MAP << "\t" << coh_MED << "\t" << coh_add << endl;
				//outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED <<  endl;
			}
		}
		
		for (int i=0; i<(int)t_incoh.size(); i++) {
			//cout << "Incoherent, t = " << t_incoh[i] << endl;
			cout << "t = " << t_incoh[i] << endl;
			if (strcmp(argv[2],"rel")==0) {
				RCUT = gsl_spline_eval (spline, t_incoh[i], acc);
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_incoh[i], params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_incoh[i], params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t_incoh[i], params_MED);
				outfile_chi2_incoh << Q2 << "\t" << W << "\t" << t_incoh[i] << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED << "\t" << ci_add - coh_add << endl;
				//outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED <<  endl;
			}
			else {
				RCUT = 399;
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_MAP);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_incoh[i], params_MAP);
				coh_add = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_add);
				ci_add = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_incoh[i], params_add);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t_incoh[i], params_MED);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t_incoh[i], params_MED);
				outfile_chi2_incoh << Q2 << "\t" << W << "\t" << t_incoh[i] << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED << "\t" << ci_add - coh_add << endl;
				//outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED <<  endl;
			}
		}
				
		outfile_chi2_coh.close();
		outfile_chi2_incoh.close();
		outfile_coh.close();
		outfile_incoh.close();
		gsl_spline_free (spline);
    		gsl_interp_accel_free (acc);
    	}
    	/* Second scenario: With evolution*/
	else {
		vector<double> Q2_coh, Q2_incoh;
		vector<double> W_coh, W_incoh;
		
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
     	
		ofstream outfile_coh, outfile_incoh;
     	char oname_coh[150], oname_incoh[150]; 
     	sprintf(oname_coh,"dsigma_bayes_fit_%s_%s_%s_%s_coh.txt",argv[1],argv[2],argv[3],argv[4]);
     	sprintf(oname_incoh,"dsigma_bayes_fit_%s_%s_%s_%s_incoh.txt",argv[1],argv[2],argv[3],argv[4]);
		outfile_coh.open(oname_coh);
		outfile_incoh.open(oname_incoh);
		outfile_coh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)" << endl;
		outfile_incoh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)" << endl;
		double coh_MAP, ci_MAP, coh_MED, ci_MED;
		
		double *interp_t, *interp_rmax;
		vector<double> vec_t, vec_rmax;
		gsl_interp_accel *acc
      			= gsl_interp_accel_alloc ();
      	gsl_spline *spline;
      	
      	double t, dt = 0.05, tmin = 0.001;
		/* run coherent */
		for (int i=0; i<(int)Q2_coh.size();i++) {
			/* rcut */
			sprintf(iname, "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/rmax_Q2=%.2lf.txt", Q2_coh[i]);
			rmax_read(iname, vec_t, vec_rmax);
			interp_t = (double*) malloc(vec_t.size()*sizeof(double));
			interp_rmax = (double*) malloc(vec_t.size()*sizeof(double));
			copy(vec_t.begin(),vec_t.end(),interp_t);
    			copy(vec_rmax.begin(),vec_rmax.end(),interp_rmax);
    			
      		spline = gsl_spline_alloc (gsl_interp_cspline, vec_t.size());
			gsl_spline_init (spline, interp_t, interp_rmax, vec_t.size());
			for (int j=0; j<62; j++) {
				t = tmin + j*dt;				
				RCUT = gsl_spline_eval (spline, t, acc);
				/* run coh dsigma*/
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_coh[i], W_coh[i], -t, params_MAP);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_coh[i], W_coh[i], -t, params_MED);
				outfile_coh << Q2_coh[i] << "\t" << W_coh[i] << "\t" << t << "\t" << coh_MAP << "\t" << coh_MED << endl;
				//outfile_incoh << Q2_coh[i] << "\t" << W_coh[i] << "\t" << t << "\t" << coh_MAP << "\t" << coh_MED << endl;
			}
			gsl_spline_free (spline);
		}
		/* run incoherent */
		for (int i=0; i<(int)Q2_incoh.size();i++) {
			/* rcut */
			sprintf(iname, "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/rmax_Q2=%.2lf.txt", Q2_incoh[i]);
			rmax_read(iname, vec_t, vec_rmax);
			interp_t = (double*) malloc(vec_t.size()*sizeof(double));
			interp_rmax = (double*) malloc(vec_t.size()*sizeof(double));
			copy(vec_t.begin(),vec_t.end(),interp_t);
    			copy(vec_rmax.begin(),vec_rmax.end(),interp_rmax);
    			
      		spline = gsl_spline_alloc (gsl_interp_cspline, vec_t.size());
			gsl_spline_init (spline, interp_t, interp_rmax, vec_t.size());
			for (int j=0; j<62; j++) {
				t = tmin + j*dt;			
				RCUT = gsl_spline_eval (spline, t, acc);
				/* run coh dsigma*/
				coh_MAP = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_incoh[i], W_incoh[i], -t, params_MAP);
				coh_MED = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2_incoh[i], W_incoh[i], -t, params_MED);
				ci_MAP = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2_incoh[i], W_incoh[i], -t, params_MAP);
				ci_MED = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2_incoh[i], W_incoh[i], -t, params_MED);
				outfile_coh << Q2_incoh[i] << "\t" << W_incoh[i] << "\t" << t << "\t" << ci_MAP - coh_MAP << "\t" << ci_MED - coh_MED << endl;
				//outfile_incoh << Q2_coh[i] << "\t" << W_coh[i] << "\t" << t << "\t" << coh_MAP << "\t" << coh_MED << endl;
			}
			gsl_spline_free (spline);
		}
		gsl_interp_accel_free (acc);
		outfile_coh.close();		
		outfile_incoh.close();
	}
	datafile_coh.close();
	datafile_incoh.close();
	gsl_rng_free(r_rng);	
	return 0;
}
