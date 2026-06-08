/*
	testing
*/

#include "nintlib.hpp"
#include "root.hpp"
#include "aux.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <string>
#include <vector>
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
#include <functional>
#include <fstream>

using namespace std;
using namespace std::placeholders;

int main(void)
{	
	aux_functions aux;
	double m2 = 0.044, Q2=0.1;
	double eps = sqrt(0.25*Q2 + SQR(mC));
	double R2 = 2.5;
	double t, Zi, Ki, Delta;
	ofstream f;
	f.open("test_Zi_t.txt");
	for (int i=0; i<101;i++) {
		t = i*0.001;
		Delta = std::sqrt(t);
		Zi = aux.Zi(Delta, m2, eps);
		Ki = aux.Ki(R2, Delta, m2, eps);
		f << t << "\t" << Zi << "\t" << Ki << endl;
	}
	f.close();
	return 0;
}
