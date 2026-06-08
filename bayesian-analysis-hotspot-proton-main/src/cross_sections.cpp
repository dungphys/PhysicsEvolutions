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

using namespace std;

double xsections::dsigma_coherent(int code_Nh, double Q2, double W, double t, model_params mparams)
{
	aux_functions aux;
	evolution evol;
	norm_factors norm;

	double Delta = std::sqrt(-t);
	double xP = (SQR(mJpsi) + Q2 - t)/(Q2 + SQR(W) - SQR(mProton));
	double eps = std::sqrt(0.25*Q2 + SQR(MASS_QUARK));
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double MT = aux.Mi(Delta, Q2, eps, mparams.m2, 1);
	double ML = aux.Mi(Delta, Q2, eps, mparams.m2, 0); 
	double prefactor = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0*Cf()/M_PI);

	/* evolution */
	double Nhmean = evol.meanNh(xP, mparams.p0, mparams.p1, mparams.p2);
	double Rp = evol.radius_proton(W, mparams.R0, mparams.lambda_p);
	double rH = evol.radius_hotspot(xP, mparams.rH0, mparams.lambda_H);
	
	/* without Nh fluctuation */
	if (code_Nh == 0) {
		Nhmean = std::round(Nhmean);
		double exp_factor = std::exp(SQR(rH)*t + SQR(Rp)*t)*std::exp(-SQR(Rp)*t/(Nhmean));
		return (exp_factor*prefactor*(SQR(ct*MT) + SQR(cl*ML)));
	}
	/* with Nh fluctuation */
	else if (code_Nh == 1) {
		double lambda = root_solver(Nhmean);
		//cout << "lbd = " << lambda;
		double exp_factor = std::exp(SQR(rH)*t + SQR(Rp)*t);
		double sum_factor = SQR(aux.sum_Nh(Rp, Delta, lambda));
		return (prefactor*exp_factor*sum_factor*(SQR(ct*MT) + SQR(cl*ML)));
	}
	else {
		cout << "code only 0 or 1\n";
		return 1;
	}		
}

double xsections::dsigma_ci(int code_Nh, double Q2, double W, double t, model_params mparams)
{
	aux_functions aux;
	evolution evol;
	norm_factors norm;
	
	double Delta = std::sqrt(-t);
	double xP = (SQR(mJpsi) + Q2 - t)/(Q2 + SQR(W) - SQR(mProton));
	double eps = std::sqrt(0.25*Q2 + SQR(MASS_QUARK));
	
	/* evolution */
	double Nhmean = evol.meanNh(xP, mparams.p0, mparams.p1, mparams.p2);
	double Rp = evol.radius_proton(W, mparams.R0, mparams.lambda_p);
	double rH = evol.radius_hotspot(xP, mparams.rH0, mparams.lambda_H);

	double prefactor_DC = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0*Cf()/M_PI);
	double prefactor_C = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0/M_PI)*(Cf()/NC);
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double MT = aux.Mi(Delta, Q2, eps, mparams.m2, 1);
	double ML = aux.Mi(Delta, Q2, eps, mparams.m2, 0);
	//double KT1 = aux.Ki(SQR(rH), Delta, mparams.m2, eps, Q2, 1);
	//double KT2 = aux.Ki(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 1);
	//double KL1 = aux.Ki(SQR(rH), Delta, mparams.m2, eps, Q2, 0);
	//double KL2 = aux.Ki(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 0);
	double DC1, DC2, C;//C1, C2;
	double esigma = std::exp(SQR(mparams.sigma_Qs));
	
	//cout << KT1 << "\t" << KL1 << endl;
	//cout << KT2 << "\t" << KL2 << endl; 
	 	
	if (code_Nh == 0) {
		Nhmean = std::round(Nhmean);
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*(1./Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*((Nhmean - 1.)/Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		//cout << DC1 << "\t"<<  DC2 << endl;
		//C1 = prefactor_C*(1./Nhmean)*(SQR(ct)*KT1 + SQR(cl)*KL1);
		//C2 = prefactor_C*((Nhmean - 1.)/Nhmean)*((SQR(ct)*KT2 + SQR(cl)*KL2));
		C = prefactor_C*aux.Ki_gen(SQR(rH), SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 1./Nhmean, esigma);
	}
	else {
		double lambda = root_solver(Nhmean);
		//cout << lambda << aux.F_lambda(lambda) << endl;
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*aux.F_lambda(lambda)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*(1-aux.F_lambda(lambda))*(SQR(ct*MT) + SQR(cl*ML));
		//C1 = prefactor_C*(aux.F_lambda(lambda))*(SQR(ct)*KT1 + SQR(cl)*KL1);
		//C2 = prefactor_C*(1-aux.F_lambda(lambda))*(SQR(ct)*KT2 + SQR(cl)*KL2);
		C = prefactor_C*aux.Ki_gen(SQR(rH), SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, aux.F_lambda(lambda), esigma);
	}
	
	return (esigma*DC1 + DC2 + C);
}

double xsections::dsigma_ci_rcut(int code_Nh, double Q2, double W, double t, model_params mparams)
{
	aux_functions aux;
	evolution evol;
	norm_factors norm;
	
	double Delta = std::sqrt(-t);
	double xP = (SQR(mJpsi) + Q2 - t)/(Q2 + SQR(W) - SQR(mProton));
	double eps = std::sqrt(0.25*Q2 + SQR(MASS_QUARK));
	
	/* evolution */
	double Nhmean = evol.meanNh(xP, mparams.p0, mparams.p1, mparams.p2);
	double Rp = evol.radius_proton(W, mparams.R0, mparams.lambda_p);
	double rH = evol.radius_hotspot(xP, mparams.rH0, mparams.lambda_H);

	double prefactor_DC = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0*Cf()/M_PI);
	double prefactor_C = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0/M_PI)*(Cf()/NC);
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double MT = aux.Mi(Delta, Q2, eps, mparams.m2, 1);
	double ML = aux.Mi(Delta, Q2, eps, mparams.m2, 0);
	double KT1 = aux.Ki_rcut(SQR(rH), Delta, mparams.m2, eps, Q2, 1);
	double KT2 = aux.Ki_rcut(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 1);
	double KL1 = aux.Ki_rcut(SQR(rH), Delta, mparams.m2, eps, Q2, 0);
	double KL2 = aux.Ki_rcut(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 0);
	double DC1, DC2, C1, C2;
	double esigma = std::exp(SQR(mparams.sigma_Qs));
	 	
	if (code_Nh == 0) {
		Nhmean = std::round(Nhmean);
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*(1./Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*((Nhmean - 1.)/Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		C1 = prefactor_C*(1./Nhmean)*(SQR(ct)*KT1 + SQR(cl)*KL1);
		C2 = prefactor_C*((Nhmean - 1.)/Nhmean)*(SQR(ct)*KT2 + SQR(cl)*KL2);
	}
	else {
		double lambda = root_solver(Nhmean);
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*aux.F_lambda(lambda)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*(1-aux.F_lambda(lambda))*(SQR(ct*MT) + SQR(cl*ML));
		C1 = prefactor_C*(aux.F_lambda(lambda))*(SQR(ct)*KT1 + SQR(cl)*KL1);
		C2 = prefactor_C*(1-aux.F_lambda(lambda))*(SQR(ct)*KT2 + SQR(cl)*KL2);
	}
	
	return (esigma*DC1 + DC2 + esigma*C1 + C2);
}

double xsections::dsigma_ci_rcut_red(int code_Nh, double Q2, double W, double t, model_params mparams)
{
	aux_functions aux;
	evolution evol;
	norm_factors norm;
	
	double Delta = std::sqrt(-t);
	double xP = (SQR(mJpsi) + Q2 - t)/(Q2 + SQR(W) - SQR(mProton));
	double eps = std::sqrt(0.25*Q2 + SQR(MASS_QUARK));
	
	/* evolution */
	double Nhmean = evol.meanNh(xP, mparams.p0, mparams.p1, mparams.p2);
	double Rp = evol.radius_proton(W, mparams.R0, mparams.lambda_p);
	double rH = evol.radius_hotspot(xP, mparams.rH0, mparams.lambda_H);

	double prefactor_DC = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0*Cf()/M_PI);
	double prefactor_C = (9./(16*M_PI))*SQR(mparams.gmu0*mparams.gmu0/M_PI)*(Cf()/NC);
	double ct = norm.CT(Q2), cl = norm.CL(Q2);
	double MT = aux.Mi(Delta, Q2, eps, mparams.m2, 1);
	double ML = aux.Mi(Delta, Q2, eps, mparams.m2, 0);
	//double KT1 = aux.Ki_rcut(SQR(rH), Delta, mparams.m2, eps, Q2, 1);
	//double KT2 = aux.Ki_rcut(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 1);
	//double KL1 = aux.Ki_rcut(SQR(rH), Delta, mparams.m2, eps, Q2, 0);
	//double KL2 = aux.Ki_rcut(SQR(rH) + SQR(Rp), Delta, mparams.m2, eps, Q2, 0);
	double DC1, DC2, C;
	double esigma = std::exp(SQR(mparams.sigma_Qs));
	 	
	if (code_Nh == 0) {
		Nhmean = std::round(Nhmean);
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*(1./Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*((Nhmean - 1.)/Nhmean)*(SQR(ct*MT) + SQR(cl*ML));
		C = prefactor_C*aux.Ki_rcut_gen(SQR(rH), SQR(Rp)+SQR(rH), Delta, mparams.m2, eps, Q2, 1./Nhmean, esigma);
	}
	else {
		double lambda = root_solver(Nhmean);
		DC1 = prefactor_DC*std::exp(SQR(rH)*t)*aux.F_lambda(lambda)*(SQR(ct*MT) + SQR(cl*ML));
		DC2 = prefactor_DC*std::exp((SQR(Rp)+SQR(rH))*t)*(1-aux.F_lambda(lambda))*(SQR(ct*MT) + SQR(cl*ML));
		C = prefactor_C*aux.Ki_rcut_gen(SQR(rH), SQR(Rp)+SQR(rH), Delta, mparams.m2, eps, Q2, aux.F_lambda(lambda), esigma);
	}
	
	return (esigma*DC1 + DC2 + C);
}









