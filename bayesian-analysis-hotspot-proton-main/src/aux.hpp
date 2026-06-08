/*
 	This defines auxiliary functions and parameters to use
*/

#ifndef aux_hpp
#define aux_hpp

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <vector> // For std::vector

extern const gsl_rng_type * T_rng;
extern gsl_rng * r_rng;
//extern double Q2_INPUT;
//extern std::vector<double> EXPR_DATA[7];

const int NCUT = 36;
const int NC = 3; 	/* number of color */ 
const int NQ = 3; 	/* number of quark */
const double Echarge = 0.30282212;	/* elementary charge in term of natural unit */  
const int Nh_fixed = 3;
const double mC = 1.275; /* mass of c quark (GeV, from NRQCD analysis) */
const double mB = 4.183; /* mass of b quark (GeV, from PDG 2024) */
const double eC = 2./3; /* charge of c quark */
const double eB = -1./3; /* charge of b quark */
const double alpha_em = 1/137.035; 			/* fine structure constant */ 
const double Jpsi_ee_decay_width = 5.55e-6;		/* GeV */
const double upsilon1S_ee_decay_width = 1.291078e-6;	/* GeV */
const double mJpsi = 3.096;	// (GeV)
const double mUpsilon = 9.46;	// (GeV)
const double mProton = 0.981;	// (GeV) 
//const double A_rel = 0.213;
//const double B_rel = -0.0157;

extern double A, B;
extern double MASS_QUARK, MASS_VM, CHARGE_QUARK;
extern double RCUT;


/* fundamental Casimir */
inline double Cf(void) {return (((NC*NC)-1.)/(2.*NC)); }	

/* for error handle */
#ifndef LINEINFO
#define LINEINFO __FILE__ << ":" << __LINE__
#endif

void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno);
                        
/* Parameters of the model to be tuned: xP-dependent scheme */
struct model_params
{
	int NPARAMS; 		/* number of active params (#0) */
	/*
		- NPARAMS = 5 if no evolution and without Qs fluctuation
		- NPARAMS = 6 if no evolution and with Qs fluctuation
		- NPARAMS = 7 if evolution in size/Nh and without Qs fluct
		- NPARAMS = 8 if evolution in size/Nh and with Qs fluct
	*/
	/* for hot spot size */
	double rH0;		
	double lambda_H;	/* zero if no hot spot size evolution */
	/* for proton size */
	double R0;		
	double lambda_p; 	/* zero if no proton size evolution */
	/* for mean number of hot spot*/
	double p0;		/*identical to Nhmean in case of no Nh fluct*/
	double p1;		/* zero in case of no Nh fluct */
	double p2;		/* zero in case of no Nh fluct */
	/* overall constant g*mu0 encodes the field strenght of a hot spot */
	double gmu0;
	/* width of Qs fluctuation */
	double sigma_Qs;	/* zero in case of no Qs fluct */
	/* IR cutoff */
	double m2;		
};

class norm_factors
{
	public:
		double AQ(char quark); 
		double CT(double Q2); // for transverse part 
		double CL(double Q2); // for longitudinal part
};
/* Parameters for Mi integrations */
struct int_par_Mi_1 {
	double r;
	double Delta;
	double epsilon;
	double Q2;
	double m2;
	int pol;
};

struct int_par_Mi_2 {
	double Delta;
	double epsilon;
	double Q2;
	double m2;
	int pol;
};

/* integrand for MC integration of Ki */
struct int_par_Ki {
	double R;
	double Delta;
	double m2;
	double epsilon;
	double Q2;
	int pol; //pol.state (L=0, T=1)
};

struct int_par_Ki_gen {
	double R1;
	double R2;
	double Delta;
	double m2;
	double epsilon;
	double Q2;
	double Nh_factor;
	double sigma_Qs_factor;
};

double factorial(int n);
double truncated_Poisson(int n, double lambda);
double wf(double r, double Delta ,double Q2, double eps, int part, int pol);
double Mi_integrand_eta(double eta, void *params);
double Mi_integrand(double r, void *params);
double X_Delta(double Delta, double eps, double Q2, int pol);
double X_k(double k, double Delta, double eps, double Q2, int pol);
double X_Delta_cut(double Delta, double eps, double Q2, int pol);
double X_k_cut(double k, double Delta, double eps, double Q2, int pol);
double Ki_integrand(double x[], size_t dim_num, void *params);
double Ki_integrands_tot(double x[], size_t dim_num, void *params);
double Ki_integrand_rcut(double x[], size_t dim_num, void *params);
double Ki_integrands_tot_rcut(double x[], size_t dim_num, void *params);

class aux_functions
{		
	public:
		double Mi(double Delta, double Q2, double eps, double m2, int pol);
		double sum_Nh(double Rp, double Delta, double lambda);
		double F_lambda(double lambda);
		double Ki_rcut(double R, double Delta, double m2, double eps, double Q2, int pol);
		double Ki_rcut_gen(double R1, double R2, double Delta, double m2, double eps, double Q2, double Nh_factor, double sigma_Qs_factor);
		double Ki(double R, double Delta, double m2, double eps, double Q2, int pol);
		double Ki_gen(double R1, double R2, double Delta, double m2, double eps, double Q2, double Nh_factor, double sigma_Qs_factor);
};

class evolution
{
	public:
		/* evolution of the mean number of hot spot */
		double meanNh(double x, double p0, double p1, double p2);
		/* evolution of the target size */
		double radius_proton(double x, double R0, double lp);
		/* evolution of the */
		double radius_hotspot(double x, double rh0, double lh);
};

#endif /*aux_hpp*/
