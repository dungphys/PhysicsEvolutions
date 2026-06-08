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
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector> // For std::vector
#include <fstream>
#include <limits>

using namespace std;

int main(int argc, char* argv[])
{
	gsl_set_error_handler(&ErrHandler);

	aux_functions aux;
	double A = 0.123, B = -0.0157;
	double m = 0.3; 
	double eps, Delta;
	//char quark = 'c'; 
	
	double Z0, Z1, Z2;
	double Q2, Q2max = 30, Q2min=0, dQ2 = 1;
	double t, tmax = 5, tmin=0.0, dt = 0.2;
	int nQ = (int) std::round((Q2max-Q2min)/dQ2) + 1;
	int nt = (int) std::round((tmax-tmin)/dt) + 1;
	double trans, longi, ratio, tmp1, tmp2, tmp3, tmp4;
	
	/* Q2 dependence */
	ofstream outfile_Q2;
	ofstream outfile_t;
	ifstream infile;
	char file_name[100];
	
	std::vector<double> t_data, R_data, error_up, error_down;
	string line;
	
	infile.open ("/home/ad_le/current_projects/hot_spot/data/ratio_data_t.txt");
	if ( infile.is_open() ) { // always check whether the file is open
		cout << "Datafile opened successfully! \n now read ..." << endl;
	}
	else {
		cout << "Error opening file" << endl;
		return 1;
	}
	
	for (int i=0;i<1;i++) getline(infile, line); // skip first nskip=1 lines;
	
	while (getline(infile, line)) {
		std::stringstream ss(line);
		/**/
		getline(ss, line, '\t');
		t_data.push_back(stod(line));
		/**/
		getline(ss, line, '\t');
		R_data.push_back(stod(line));
		/**/
		getline(ss, line, '\t'); tmp1 = stod(line);
		getline(ss, line, '\t'); tmp2 = stod(line);
		getline(ss, line, '\t'); tmp3 = stod(line);
		getline(ss, line, '\t'); tmp4 = stod(line);
		error_up.push_back(sqrt(tmp1*tmp1 + tmp3*tmp3));
		error_down.push_back(sqrt(tmp2*tmp2 + tmp4*tmp4));
	} 
      infile.close();
      
      double *t_hera = t_data.data();
      double *R_hera = R_data.data();
      double *std_up = error_up.data();
      double *std_down = error_down.data();
      int datasize = t_data.size();
      for (int i=0;i<datasize;i++) {
      	cout << t_hera[i] << "\t" << R_hera[i] << "\t" << std_up[i] << "\t" << std_down[i] << endl;
      }
	//Q2 = (double *) malloc(nQ*sizeof(double));
	//t = (double *) malloc(nt*sizeof(double));
	for (int i=0;i<nt;i++) {
		t = -(tmin + i*dt);
		Delta = std::sqrt(-t);
		sprintf(file_name,"ratio_t=%.2lf_m=%.2lf.txt",t,m);
		outfile_Q2.open(file_name);
		cout << "t = " << t << endl;
		for (int j=0;j<nQ;j++) {
			Q2 = j*dQ2;
			eps = std::sqrt(0.25*Q2 + SQR(mC));
			Z0 = aux.Zeta(Delta, SQR(m), eps); 
			Z1 = aux.Zeta_1(Delta, SQR(m), eps);
			Z2 = aux.Zeta_2(Delta, SQR(m), eps);
			trans = A*Z0 + (B/(SQR(mC)))*( 3.5*Z0 + (SQR(mC) - 0.125*t)*Z1 - ((Q2 + 2*SQR(mC))/(2*eps))*Z2 );
			longi = A*Z0 + (B/(SQR(mC)))*( 4.5*Z0 + (SQR(mC) - 0.125*t)*Z1 - ((Q2)/(4*eps))*Z2 );
			ratio = SQR(longi/trans)*(Q2/(4*SQR(mC)));
			outfile_Q2 << Q2 << "\t"<< ratio << "\t" << (Q2/(4*SQR(mC))) << endl; 
		}
		outfile_Q2.close();	
	}
	/* t dependence */
	Q2 = 6.8;
	eps = std::sqrt(0.25*Q2 + SQR(mC));
	sprintf(file_name,"ratio_Q2=%.2f_m=%.2f.txt",Q2,m);
	outfile_t.open(file_name);
	//double chi2_m;
	//double m_min = 0.02, m_max=0.6, dm = 0.01;
	//int nm = (int) std::round((m_max-m_min)/dm) + 1;
	for (int j=0;j<datasize;j++) {
		Delta = std::sqrt(t_hera[j]);
		Z0 = aux.Zeta(Delta, SQR(m), eps); 
		Z1 = aux.Zeta_1(Delta, SQR(m), eps);
		Z2 = aux.Zeta_2(Delta, SQR(m), eps);
		trans = A*Z0 + (B/(SQR(mC)))*( 3.5*Z0 + (SQR(mC) + 0.125*t_hera[j])*Z1 - ((Q2 + 2*SQR(mC))/(2*eps))*Z2 );
		longi = A*Z0 + (B/(SQR(mC)))*( 4.5*Z0 + (SQR(mC) + 0.125*t_hera[j])*Z1 - ((Q2)/(4*eps))*Z2 );
		ratio = SQR(longi/trans)*(Q2/(4*SQR(mC)));
		outfile_t << t_hera[j] << "\t"<< ratio << "\t" << (Q2/(4*SQR(mC))) << endl;
	}
	
	outfile_t.close();
	free(t_hera); free(R_hera);
	free(std_down); free(std_up);
	/* m dependence */
	/*int nm = 16;
	nt = 3;
	Q2 = 6.8;
	eps = std::sqrt(0.25*Q2 + SQR(mC));
	
	for (int i=0;i<nt;i++) {
		t = -0.02*std::pow(10,i);
		Delta = std::sqrt(-t);
		sprintf(file_name,"ratio_Q2=%.2lf_t=%.2lf.txt",Q2,-t);
		file.open(file_name);
		for (int j=0;j<nm;j++) {
			m = 0.08 + j*0.02;
			Z0 = aux.Zeta(Delta, SQR(m), eps); 
			Z1 = aux.Zeta_1(Delta, SQR(m), eps);
			Z2 = aux.Zeta_2(Delta, SQR(m), eps);
			trans = A*Z0 + (B/(SQR(mC)))*( 3.5*Z0 + (SQR(mC) - 0.25*t)*Z1 - ((Q2 + 2*SQR(mC))/(2*eps))*Z2 );
			longi = A*Z0 + (B/(SQR(mC)))*( 4.5*Z0 + (SQR(mC) - 0.25*t)*Z1 - ((Q2)/(4*eps))*Z2 );
			ratio = SQR(longi/trans)*(Q2/(4*SQR(mC)));
			file << m << "\t"<< ratio << "\t" << (Q2/(4*SQR(mC))) << endl;
		}
		file.close();
	}
	*/
	return 0;
}
