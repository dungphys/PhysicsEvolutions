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
	if (argc!=4) {
		cout << "Syntax should be: ./run_test scenario Nh Qs" << endl;
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
    	if (strcmp(argv[2],"with_Nh_fluct")==0) code_Nh = 1;
    	else code_Nh = 0;
    	
    	/* J/psi wf */
    	//if (strcmp(argv[2],"rel")==0) {A = 0.213; B = -0.0157;}
    	//else {A = 0.211; B = 0;}
    	MASS_QUARK = mC;
	CHARGE_QUARK = eC;
    	double gev2_to_nb = 0.389379e6;
	
	/*init params*/
     params.sigma_Qs = 0.0;
     params.p0 = 3.0; params.p1 = 0.; params.p2 = 0.;
     params.lambda_p = 0.; params.lambda_H = 0.;
     params.gmu0 = pow(43.222,0.25);
     params.m2 = 0.0484;
     params.R0 = 1.817;
     //params.rH0 = 0.7;
     params.rH0 = 0.837;
      
     cout << "Params initialization done!" << endl;
	 /* open datafile to read */
	 
	 std::string line;
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
		
		double t, dt = 0.05, tmin = 0.001;
		
     	ofstream outfile_coh, outfile_incoh;
     	//ofstream outfile_chi2_coh, outfile_chi2_incoh;
     	char oname_coh[150], oname_incoh[150];
     	//char for_chi2_coh[150], for_chi2_incoh[150];  
     	sprintf(oname_coh,"test_dsigma_R_vs_NR_%s_%s_%s_coh_2a.txt",argv[1],argv[2],argv[3]);
     	sprintf(oname_incoh,"test_dsigma_R_vs_NR_%s_%s_%s_incoh_2a.txt",argv[1],argv[2],argv[3]);
		outfile_coh.open(oname_coh);
		outfile_incoh.open(oname_incoh);
		/**/
		//sprintf(for_chi2_coh,"test_dsigma_chi2_bayes_%s_%s_%s_%s_coh.txt",argv[1],argv[2],argv[3],argv[4]);
     	//sprintf(for_chi2_incoh,"test_dsigma_chi2_bayes_%s_%s_%s_%s_incoh.txt",argv[1],argv[2],argv[3],argv[4]);
		//outfile_chi2_coh.open(for_chi2_coh);
		//outfile_chi2_incoh.open(for_chi2_incoh);
		//outfile_chi2_coh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)" << endl;
		//outfile_chi2_incoh << "Q2\tW\tt\tdsigma/dt(MAP)\tdsigma/dt(MED)" << endl;
		outfile_coh << "Q2\tW\tt\tdsigma/dt(R)\tdsigma/dt(NR)" << endl;
		outfile_incoh << "Q2\tW\tt\tdsigma/dt(R)\tdsigma/dt(NR)" << endl;
		double coh_NR, ci_NR, coh_R, ci_R;
		for (int i=0; i<150; i++) {
			//cout << "Incoherent, t = " << t_incoh[i] << endl;
			t = tmin + i*dt;
			cout << "t = " << t << endl;
			
			// rel
			RCUT = gsl_spline_eval (spline, t, acc); A = 0.213; B = -0.0157;
				coh_R = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params);
				ci_R = (gev2_to_nb)*dsigma.dsigma_ci_rcut_red(code_Nh, Q2, W, -t, params);
			
			// nonrel
			RCUT = 399; A = 0.211; B = 0.0;
				coh_NR = (gev2_to_nb)*dsigma.dsigma_coherent(code_Nh, Q2, W, -t, params);
				ci_NR = (gev2_to_nb)*dsigma.dsigma_ci(code_Nh, Q2, W, -t, params);
			outfile_coh << Q2 << "\t" << W << "\t" << t << "\t" << coh_R << "\t" << coh_NR << endl;
			outfile_incoh << Q2 << "\t" << W << "\t" << t << "\t" << ci_R - coh_R << "\t" << ci_NR - coh_NR <<  endl;
			
		}
		
		outfile_coh.close();
		outfile_incoh.close();
		gsl_spline_free (spline);
    		gsl_interp_accel_free (acc);
    	}
    
	datafile_coh.close();
	datafile_incoh.close();
	gsl_rng_free(r_rng);	
	return 0;
}
