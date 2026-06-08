#ifndef INIT_CONDITION_H
#define INIT_CONDITION_H

#include <vector>
#include <cmath>
#include "utils.h"

typedef enum
{
    GBW,  // (modified) Golec-Biernat-Wusthoff
    STEP, // Step function
    MV,   // (modified) McLerran-Venugopalan
} InitCond;

/**
 * GBW profle: S(r, Y=0) = exp( -(r²Q²/4)^gamma )
 */
double gbw_gamma(double r, std::vector<double> &params);

/**
 * Step function: S(r, Y=0) = 1 if r < r0, else 0
 */
double step(double r, std::vector<double> &params);

/**
 * MV profile: S(r, Y=0) = exp( -(r²Q²/4)^gamma * ln(e + 1/(r*lQCD)) )
 */
double mv_gamma(double r, std::vector<double> &params);

void getInitProfile(Lattice &lattice, InitCond key, std::vector<double> &params);

#endif // INIT_CONDITION_H