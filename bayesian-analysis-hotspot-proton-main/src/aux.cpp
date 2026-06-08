/*
 	This defines auxiliary functions and parameters to use
*/
#define _USE_MATH_DEFINES

#include "nintlib.hpp"
#include "root.hpp"
#include "aux.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector> // For std::vector
#include <functional>

using namespace std;
using namespace std::placeholders;

const gsl_rng_type * T_rng;
gsl_rng * r_rng;

double A, B;
double MASS_QUARK, MASS_VM, CHARGE_QUARK;
double RCUT;

int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
   
    // Errors related to convergence of integrals are handled when
    // gsl_integration functions are called
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    if (gsl_errno == 14 or gsl_errno == 18 or gsl_errno == 11)
        return;

    // 15: underflows
    if (gsl_errno == 15 ) return;
    // 16: overflows
    if (gsl_errno == 16 ) return;


    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}
/*----------------------------------------------------------------------------------------------------------------*/
/*--------------------------       Define auxiliary functions      -----------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------*/


/* truncated Poisson distribution for number fluctuation*/
double factorial(int n)
{
	if (n==1) return 1;
	else return (n*factorial(n-1)); 	
}

double truncated_Poisson(int n, double lambda)
{
	return (std::pow(lambda,n)/((std::exp(lambda)-1)*gsl_sf_gamma(n+1)));
}

/* Defined M_{TL} functions */
double wf(double r, double Delta ,double Q2, double eps, int part, int pol)
{
	gsl_sf_result eval;
	gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off ();
	
	gsl_set_error_handler(old_error_handler);
	int code1 = gsl_sf_bessel_K1_e(r*eps,&eval);
	double besselK1 = (code1==GSL_SUCCESS) ? eval.val : 0.0;
	int code0 = gsl_sf_bessel_K0_e(r*eps,&eval);
	double besselK0 = (code0==GSL_SUCCESS) ? eval.val : 0.0;
	double f2 = (B/SQR(MASS_QUARK))*(1./8)*SQR(r*Delta)*besselK0;		
	double f1 = A*besselK0 + B/SQR(MASS_QUARK)*((4.5-pol)*besselK0  + SQR(MASS_QUARK*r)*besselK0 + (Q2/(4*eps) + pol)*r*besselK1) + f2;
	return ((part-1)*f2 - (part-2)*f1);
}

double Mi_integrand_eta(double eta, void *params)
{
	int_par_Mi_1* par = (int_par_Mi_1*) params;	
	double r = par->r, Q2 = par->Q2, Delta = par->Delta, eps = par->epsilon, m2 = par->m2;
	int pol = par->pol;
	gsl_sf_result eval;
	gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off ();
	
	gsl_set_error_handler(old_error_handler);
	double var = std::sqrt(-SQR(Delta)*SQR(eta) + 0.25*SQR(Delta) + m2);
	int codeK1 = gsl_sf_bessel_K1_e(r*var,&eval);
	double besselK1 = (codeK1==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ0 = gsl_sf_bessel_J0_e(eta*Delta*r,&eval);
	double besselJ0 = (codeJ0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ2 = gsl_sf_bessel_Jn_e(2,eta*Delta*r,&eval);
	double besselJ2 = (codeJ2==GSL_SUCCESS) ? eval.val : 0.0;
	double u = r*besselK1/var;
	//cout << wf(r,Delta,Q2,eps,1,pol) << endl;
	return (u*(besselJ0*wf(r,Delta,Q2,eps,1,pol) - besselJ2*wf(r,Delta,Q2,eps,2,pol)));
}

double Mi_integrand(double r, void *params)
{
	int_par_Mi_2* par = (int_par_Mi_2*) params;	
	double Q2 = par->Q2, Delta = par->Delta, eps = par->epsilon, m2 = par->m2;
	int pol = par->pol;
	int_par_Mi_1 p;
	p.r = r; p.Q2 = Q2; p.Delta = Delta; p.m2 = m2; 
	p.epsilon = eps; p.pol = pol;
	
	double Psi0 = 2*gsl_atanh(Delta/sqrt(SQR(Delta) + 4*m2))/(Delta*sqrt(SQR(Delta) + 4*m2));
	gsl_sf_result eval;
	gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off ();
	
	gsl_set_error_handler(old_error_handler);
	int codeJ0 = gsl_sf_bessel_J0_e(0.5*Delta*r,&eval);
	double besselJ0 = (codeJ0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ2 = gsl_sf_bessel_Jn_e(2,0.5*Delta*r,&eval);
	double besselJ2 = (codeJ2==GSL_SUCCESS) ? eval.val : 0.0;
	
	double term0 = Psi0*(besselJ0*wf(r,Delta,Q2,eps,1,pol) - besselJ2*wf(r,Delta,Q2,eps,2,pol));
	//cout << Psi0 << endl;
	//cout << term0 << endl;
	
	gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(1000);
	double res, err;
	gsl_function f;
	f.params = &p; f.function = Mi_integrand_eta;
		
	int status = gsl_integration_qag(&f, 0, 0.5, 0, 1e-3, 1000, GSL_INTEG_GAUSS61, wspace, &res, &err);
	gsl_integration_workspace_free(wspace);
		
	if(status)
         cerr << "The Mi_integrand integration failed at " << LINEINFO <<", result " << res
                        << ", relerr " << std::abs(res-err)/res << endl;
     //cout << r*(term0 - res) << endl;
     return (r*(term0 - res));	  
}

double aux_functions::Mi(double Delta, double Q2, double eps, double m2, int pol)
{
	int_par_Mi_2 params;
	params.Q2 = Q2; params.Delta = Delta; params.m2 = m2; 
	params.epsilon = eps; params.pol = pol;
	
	gsl_integration_workspace *wspace = gsl_integration_workspace_alloc(1000);
	double res, err;
	gsl_function f;
	f.params = &params; f.function = Mi_integrand;
	
	int status = gsl_integration_qag(&f, 0, RCUT, 0, 1e-3, 1000, GSL_INTEG_GAUSS61, wspace, &res, &err);
	gsl_integration_workspace_free(wspace);
	
	if(status)
         cerr << "The Mi integration failed at " << LINEINFO <<", result " << res
                        << ", relerr " << std::abs(res-err)/res << endl;
	
	return (2*M_PI*res);
}

/* Defined sums over Nh*/
double aux_functions::sum_Nh(double Rp, double Delta, double lambda)
{
	double sum = 0;
	for (int i=1;i<=NCUT;i++) {
		sum = sum + (truncated_Poisson(i, lambda)*exp(SQR(Rp*Delta)/(2*i)));
	}
	return sum;
}

double aux_functions::F_lambda(double lambda)
{
	return (gsl_sf_expint_Ei(lambda) - log(lambda) - M_EULER)/(exp(lambda) - 1);
}

/* X functions */
double X_Delta(double Delta, double eps, double Q2, int pol)
{
	double term1, term2, term3, term4;
	double var = SQR(eps) + 0.25*SQR(Delta);
	term1 = (A + (4.5-pol)*(B/(SQR(MASS_QUARK))))/var;
	term2 = 4*B*(1+SQR(Delta)/(8*SQR(MASS_QUARK)))*(SQR(eps) - 0.25*SQR(Delta))/CUBE(var);
	term3 = (2*B/(SQR(MASS_QUARK)))*(0.25*Q2 + pol*SQR(eps))/SQR(var);
	term4 = (B/SQR(2*MASS_QUARK))*(SQR(Delta*Delta)/CUBE(var));  
	//cout << term1 << "\t" << term2 << "\t" << term3 << "\t" << term4 << endl;
	return (term1 + term2 - term3 - term4);
}

double X_k(double k, double Delta, double eps, double Q2, int pol)
{
	double term1, term2, term3, term4;
	double var = SQR(eps) + SQR(k);
	term1 = (A + (4.5-pol)*(B/(SQR(MASS_QUARK))))/var;
	term2 = 4*B*(1+SQR(Delta)/(8*SQR(MASS_QUARK)))*(SQR(eps) - SQR(k))/CUBE(var);
	term3 = (2*B/(SQR(MASS_QUARK)))*(0.25*Q2 + pol*SQR(eps))/SQR(var);
	term4 = (B/SQR(MASS_QUARK))*(SQR(Delta*k)/CUBE(var));  
	//cout << eps << "\t" << var << endl;
	//cout << B/ << "\t" << term2 << "\t" << term3 << "\t" << term4 << endl;
	return (term1 + term2 - term3 - term4);
}

double X_Delta_cut(double Delta, double eps, double Q2, int pol)
{	
	gsl_sf_result eval;
	gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off ();
	
	gsl_set_error_handler(old_error_handler);
	int codeJ0 = gsl_sf_bessel_J0_e(0.5*Delta*RCUT,&eval);
	double J0 = (codeJ0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ1 = gsl_sf_bessel_J1_e(0.5*Delta*RCUT,&eval);
	double J1 = (codeJ1==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ2 = gsl_sf_bessel_Jn_e(2,0.5*Delta*RCUT,&eval);
	double J2 = (codeJ2==GSL_SUCCESS) ? eval.val : 0.0;
	int codeK0 = gsl_sf_bessel_K0_e(eps*RCUT,&eval);
	double K0 = (codeK0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeK1 = gsl_sf_bessel_K1_e(eps*RCUT,&eval);
	double K1 = (codeK1==GSL_SUCCESS) ? eval.val : 0.0;
	
	double P = SQR(eps) + SQR(Delta/2), M = SQR(eps) - SQR(Delta/2);
	double mQ2 = SQR(MASS_QUARK);
	double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0;
	
	t1 = ((A + ((4.5-pol)*B)/(mQ2))/P)*(0.5*Delta*RCUT*J1*K0 - eps*RCUT*J0*K1 + 1);
	t2 = B*(1. + SQR(Delta)/(8*mQ2))*( (4*M/CUBE(P)) + (2*SQR(RCUT)/SQR(P))*(-M*J0*K0 + eps*Delta*K1*J1) );
	t3 = -B*(1. + SQR(Delta)/(8*mQ2))*( (4*M*RCUT/CUBE(P)) + CUBE(RCUT)/P )*(eps*J0*K1 - 0.5*Delta*J1*K0);
	t4 = -(B/mQ2)*((Q2/(4*eps)) + pol*eps)*( (2*eps/SQR(P)) + (2*eps*RCUT/SQR(P))*(-eps*J0*K1 + 0.5*Delta*J1*K0) + (SQR(RCUT)/P)*(0.5*Delta*J1*K1 - eps*J0*K0) );
	t5 = -(B*SQR(Delta)/(8*mQ2))*( (2*SQR(Delta)/CUBE(P)) - (CUBE(RCUT)/eps)*J2*K1 );
	t6 = (B*SQR(Delta)/(8*mQ2))*( Delta*SQR(RCUT)/(eps*SQR(P)) )*(Delta*eps*J0*K0 + M*J1*K1);
	t7 = (B*SQR(Delta)/(8*mQ2))*( 2*SQR(Delta)*RCUT/CUBE(P) )*(eps*J0*K1 - 0.5*Delta*J1*K0);
	t8 = (B*SQR(Delta)/(8*mQ2))*( Delta*CUBE(RCUT)/(2*eps*P) )*(0.5*Delta*J0*K1 + eps*J1*K0);
	return (t1+t2+t3+t4+t5+t6+t7+t8);
}

double X_k_cut(double k, double Delta, double eps, double Q2, int pol)
{
	gsl_sf_result eval;
	gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off ();
	
	gsl_set_error_handler(old_error_handler);
	int codeJ0 = gsl_sf_bessel_J0_e(k*RCUT,&eval);
	double J0 = (codeJ0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ1 = gsl_sf_bessel_J1_e(k*RCUT,&eval);
	double J1 = (codeJ1==GSL_SUCCESS) ? eval.val : 0.0;
	int codeJ2 = gsl_sf_bessel_Jn_e(2,k*RCUT,&eval);
	double J2 = (codeJ2==GSL_SUCCESS) ? eval.val : 0.0;
	int codeK0 = gsl_sf_bessel_K0_e(eps*RCUT,&eval);
	double K0 = (codeK0==GSL_SUCCESS) ? eval.val : 0.0;
	int codeK1 = gsl_sf_bessel_K1_e(eps*RCUT,&eval);
	double K1 = (codeK1==GSL_SUCCESS) ? eval.val : 0.0;
	
	double P = SQR(eps) + SQR(k), M = SQR(eps) - SQR(k);
	double mQ2 = SQR(MASS_QUARK);
	double t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0;
	
	t1 = ((A + ((4.5-pol)*B)/(mQ2))/P)*( k*RCUT*J1*K0 - eps*RCUT*J0*K1 + 1);
	t2 = B*(1. + SQR(Delta)/(8*mQ2))*( (4*M/CUBE(P)) + (2*SQR(RCUT)/SQR(P))*(-M*J0*K0 + 2*eps*k*K1*J1) );
	t3 = -B*(1. + SQR(Delta)/(8*mQ2))*( (4*M*RCUT/CUBE(P)) + CUBE(RCUT)/P )*(eps*J0*K1 - k*J1*K0);
	t4 = -(B/mQ2)*((Q2/(4*eps)) + pol*eps)*( (2*eps/SQR(P)) + (2*eps*RCUT/SQR(P))*(-eps*J0*K1 + k*J1*K0) + (SQR(RCUT)/P)*(k*J1*K1 - eps*J0*K0) );
	t5 = -(B*SQR(Delta)/(8*mQ2))*( (8*SQR(k)/CUBE(P)) - (CUBE(RCUT)/eps)*J2*K1 );
	t6 = (B*SQR(Delta)/(8*mQ2))*( 2*k*SQR(RCUT)/(eps*SQR(P)) )*(2*k*eps*J0*K0 + M*J1*K1);
	t7 = (B*SQR(Delta)/(8*mQ2))*( 8*SQR(k)*RCUT/CUBE(P) )*(eps*J0*K1 - k*J1*K0);
	t8 = (B*SQR(Delta)/(8*mQ2))*( k*CUBE(RCUT)/(eps*P) )*(k*J0*K1 + eps*J1*K0);
	
	return (t1+t2+t3+t4+t5+t6+t7+t8);
}

/* Define K functions */

double Ki_integrand(double x[], size_t dim_num, void *params)
/*
	(x[0],x[1],x[2],x[3]) = (k1,k2,theta_1,theta_2)
	theta_1 = angle(k1, Delta), theta_2 = angle(k2, Delta)
*/
{
	int_par_Ki* p = (int_par_Ki*) params;
	double k1 = x[0], k2 = x[1], theta_1 = x[2], theta_2 = x[3];
	double R = p->R, Delta = p->Delta, m2 = p->m2, eps = p->epsilon, Q2 = p->Q2;
	int pol = p->pol; 
	double k_delta_1, k_delta_2, Dk1, Dk2; 
	double k1sq = SQR(k1), k2sq = SQR(k2), Deltasq = SQR(Delta);
	double expfactor = std::exp(-R*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	k_delta_1 = SQR(k1sq + 0.25*Deltasq + m2) - k1sq*Deltasq*std::cos(theta_1)*std::cos(theta_1);
	k_delta_2 = SQR(k2sq + 0.25*Deltasq + m2) - k2sq*Deltasq*std::cos(theta_2)*std::cos(theta_2);
	Dk1 = X_Delta(Delta, eps, Q2, pol) - X_k(k1, Delta, eps, Q2, pol);
	Dk2 = X_Delta(Delta, eps, Q2, pol) - X_k(k2, Delta, eps, Q2, pol);
	
	return ((k1*k2*expfactor*Dk1*Dk2)/(k_delta_1*k_delta_2));
}

double Ki_integrands_tot(double x[], size_t dim_num, void *params)
/*
	(x[0],x[1],x[2],x[3]) = (k1,k2,theta_1,theta_2)
	theta_1 = angle(k1, Delta), theta_2 = angle(k2, Delta)
*/
{
	norm_factors norm;
	int_par_Ki_gen* p = (int_par_Ki_gen*) params;
	double R1 = p->R1, R2 = p->R2, Delta = p->Delta, m2 = p->m2;
	double eps = p->epsilon, Q2 = p->Q2;
	double Nh_factor = p->Nh_factor, sigma_Qs_factor = p->sigma_Qs_factor;
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double k1 = x[0], k2 = x[1], theta_1 = x[2], theta_2 = x[3];
	double k_delta_1, k_delta_2, Dk1_T, Dk2_T, Dk1_L, Dk2_L; 
	double k1sq = SQR(k1), k2sq = SQR(k2), Deltasq = SQR(Delta);
	double expfactor_1 = std::exp(-R1*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	double expfactor_2 = std::exp(-R2*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	k_delta_1 = SQR(k1sq + 0.25*Deltasq + m2) - k1sq*Deltasq*std::cos(theta_1)*std::cos(theta_1);
	k_delta_2 = SQR(k2sq + 0.25*Deltasq + m2) - k2sq*Deltasq*std::cos(theta_2)*std::cos(theta_2);
	Dk1_T = X_Delta(Delta, eps, Q2, 1) - X_k(k1, Delta, eps, Q2, 1);
	Dk1_L = X_Delta(Delta, eps, Q2, 0) - X_k(k1, Delta, eps, Q2, 0);
	Dk2_T = X_Delta(Delta, eps, Q2, 1) - X_k(k2, Delta, eps, Q2, 1);
	Dk2_L = X_Delta(Delta, eps, Q2, 0) - X_k(k2, Delta, eps, Q2, 0);
	
	double C1 = sigma_Qs_factor*Nh_factor*expfactor_1*( SQR(ct)*Dk1_T*Dk2_T + SQR(cl)*Dk1_L*Dk2_L );
	double C2 = (1-Nh_factor)*expfactor_2*( SQR(ct)*Dk1_T*Dk2_T + SQR(cl)*Dk1_L*Dk2_L );
	
	return (k1*k2*(C1 + C2)/(k_delta_1*k_delta_2));
}


double Ki_integrand_rcut(double x[], size_t dim_num, void *params)
{
	int_par_Ki* p = (int_par_Ki*) params;
	//int code = p->pol_code;
	double k1 = x[0], k2 = x[1], theta_1 = x[2], theta_2 = x[3];
	double R = p->R, Delta = p->Delta, m2 = p->m2;
	double eps = p->epsilon, Q2 = p->Q2;
	int pol = p->pol;
	double k_delta_1, k_delta_2, Dk1, Dk2; 
	double k1sq = SQR(k1), k2sq = SQR(k2), Deltasq = SQR(Delta);
	double expfactor = std::exp(-R*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	k_delta_1 = SQR(k1sq + 0.25*Deltasq + m2) - k1sq*Deltasq*std::cos(theta_1)*std::cos(theta_1);
	k_delta_2 = SQR(k2sq + 0.25*Deltasq + m2) - k2sq*Deltasq*std::cos(theta_2)*std::cos(theta_2);
	Dk1 = X_Delta_cut(Delta, eps, Q2, pol) - X_k_cut(k1, Delta, eps, Q2, pol);
	Dk2 = X_Delta_cut(Delta, eps, Q2, pol) - X_k_cut(k2, Delta, eps, Q2, pol);

	return ((k1*k2*expfactor*Dk1*Dk2)/(k_delta_1*k_delta_2));
}

double Ki_integrands_tot_rcut(double x[], size_t dim_num, void *params)
{
	norm_factors norm;
	int_par_Ki_gen* p = (int_par_Ki_gen*) params;
	double R1 = p->R1, R2 = p->R2, Delta = p->Delta, m2 = p->m2;
	double eps = p->epsilon, Q2 = p->Q2;
	double Nh_factor = p->Nh_factor, sigma_Qs_factor = p->sigma_Qs_factor;
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double k1 = x[0], k2 = x[1], theta_1 = x[2], theta_2 = x[3];
	double k_delta_1, k_delta_2, Dk1_T, Dk2_T, Dk1_L, Dk2_L; 
	double k1sq = SQR(k1), k2sq = SQR(k2), Deltasq = SQR(Delta);
	double expfactor_1 = std::exp(-R1*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	double expfactor_2 = std::exp(-R2*(k1sq + k2sq + 2*k1*k2*std::cos(theta_1-theta_2))); 
	k_delta_1 = SQR(k1sq + 0.25*Deltasq + m2) - k1sq*Deltasq*std::cos(theta_1)*std::cos(theta_1);
	k_delta_2 = SQR(k2sq + 0.25*Deltasq + m2) - k2sq*Deltasq*std::cos(theta_2)*std::cos(theta_2);
	Dk1_T = X_Delta_cut(Delta, eps, Q2, 1) - X_k_cut(k1, Delta, eps, Q2, 1);
	Dk1_L = X_Delta_cut(Delta, eps, Q2, 0) - X_k_cut(k1, Delta, eps, Q2, 0);
	Dk2_T = X_Delta_cut(Delta, eps, Q2, 1) - X_k_cut(k2, Delta, eps, Q2, 1);
	Dk2_L = X_Delta_cut(Delta, eps, Q2, 0) - X_k_cut(k2, Delta, eps, Q2, 0);
	
	double C1 = sigma_Qs_factor*Nh_factor*expfactor_1*( SQR(ct)*Dk1_T*Dk2_T + SQR(cl)*Dk1_L*Dk2_L );
	double C2 = (1-Nh_factor)*expfactor_2*( SQR(ct)*Dk1_T*Dk2_T + SQR(cl)*Dk1_L*Dk2_L );
	
	return (k1*k2*(C1 + C2)/(k_delta_1*k_delta_2));
}


double aux_functions::Ki_rcut(double R, double Delta, double m2, double eps, double Q2, int pol)
{
	size_t dim_num = 4;
	double lo[dim_num] = {0., 0., 0., 0.}, hi[dim_num] = {499., 499., 2*M_PI, 2*M_PI};
	int_par_Ki params;
	params.R = R; params.Delta = Delta; params.m2 = m2;
	params.epsilon = eps; params.Q2 = Q2;
	params.pol = pol;
	int ncalls = 1000000;
	double res;
	
	MC_nd (&Ki_integrand_rcut, dim_num, &params, lo, hi, "vegas" , ncalls, &res);
	
	return res;
}

double aux_functions::Ki_rcut_gen(double R1, double R2, double Delta, double m2, double eps, double Q2, double Nh_factor, double sigma_Qs_factor)
{
	size_t dim_num = 4;
	double lo[dim_num] = {0., 0., 0., 0.}, hi[dim_num] = {599., 599., 2*M_PI, 2*M_PI};
	int_par_Ki_gen params;
	params.R1 = R1; params.R2 = R2; params.Delta = Delta; params.m2 = m2;
	params.epsilon = eps; params.Q2 = Q2;
	params.sigma_Qs_factor = sigma_Qs_factor;
	params.Nh_factor = Nh_factor;
	int ncalls = 1500000;
	double res;
	
	MC_nd (&Ki_integrands_tot_rcut, dim_num, &params, lo, hi, "vegas" , ncalls, &res);
	
	return res;
}

double aux_functions::Ki(double R, double Delta, double m2, double eps, double Q2, int pol)
{
	size_t dim_num = 4;
	double lo[dim_num] = {0., 0., 0., 0.}, hi[dim_num] = {599., 599., 2*M_PI, 2*M_PI};
	int_par_Ki params;
	params.R = R; params.Delta = Delta; params.m2 = m2;
	params.epsilon = eps; params.Q2 = Q2;
	params.pol = pol;
	int ncalls = 1200000;
	double res;
	
	MC_nd (&Ki_integrand, dim_num, &params, lo, hi, "vegas" , ncalls, &res);
	
	return res;
}

double aux_functions::Ki_gen(double R1, double R2, double Delta, double m2, double eps, double Q2, double Nh_factor, double sigma_Qs_factor)
{
	size_t dim_num = 4;
	double lo[dim_num] = {0., 0., 0., 0.}, hi[dim_num] = {599., 599., 2*M_PI, 2*M_PI};
	int_par_Ki_gen params;
	params.R1 = R1; params.R2 = R2; params.Delta = Delta; params.m2 = m2;
	params.epsilon = eps; params.Q2 = Q2;
	params.sigma_Qs_factor = sigma_Qs_factor;
	params.Nh_factor = Nh_factor;
	int ncalls = 1500000;
	double res;
	
	MC_nd (&Ki_integrands_tot, dim_num, &params, lo, hi, "vegas" , ncalls, &res);
	
	return res;
}

/*-----------------------------------------------------------------------------------------------------------------*/
/*--------------------------       Define prefactors CL and CT      -----------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------------------*/

/* Prefactors*/
double  norm_factors::AQ(char quark)
{
	if (quark == 'c') {
		return std::sqrt(Jpsi_ee_decay_width*mJpsi*mJpsi/(16*M_PI*eC*eC*alpha_em*alpha_em));
	}
	else if (quark == 'b') { return std::sqrt(upsilon1S_ee_decay_width*mUpsilon*mUpsilon/(16*M_PI*eB*eB*alpha_em)); }
	else {cout << "Only J/psi (c) or Upsilon (b)\n"; return 0;}
	
}

double  norm_factors::CT(double Q2)
{
	return (-(std::sqrt(2*MASS_QUARK*NC))*CHARGE_QUARK*Echarge/(4*M_PI)); 
}

double  norm_factors::CL(double Q2)
{
	return (-(std::sqrt(2*NC*Q2/MASS_QUARK))*CHARGE_QUARK*Echarge/(8*M_PI)); 
}

/*----------------------------------------------------------------------------------------------------------------*/
/*--------------------------              Evolutions               -----------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------*/

double evolution::meanNh(double x, double p0, double p1, double p2)
/* Mean number of hotspots */
{
	return (p0*std::pow(1./x,p1)*(1+p2*std::sqrt(x)));
}

double evolution::radius_proton(double W, double R0, double lp)
/* Proton radius */
{
	return R0*std::sqrt( 1.0 + lp*log(W) );
}

double evolution::radius_hotspot(double x, double rh0, double lh)
/* hotspot radius */
{
	return rh0*std::pow(1./x,lh);
}

/**/
