#ifndef SOLVER_H
#define SOLVER_H

#include "utils.h"
#include "interp.h"

typedef enum
{
    LO,        // loBK
    RC_PARENT, // rcBK parent-dipole kernel
    RC_MIN,    // rcBK smallest-dipole kernel
    RC_BAL,    // rcBK Balitsky kernel
    RC_KW      // rcBK Kovchegov-Weigert
} KernelType;

/* QCD running coupling */
double alpha_rc(double r);
/* kernel evaluation */
double evaluateKernel(KernelType type, double r, double r1, double r2);
/* individual kernel implementations */
double kernel_lo(double r, double r1, double r2);
double kernel_rcParent(double r, double r1, double r2);
double kernel_rcMin(double r, double r1, double r2);
double kernel_rcBal(double r, double r1, double r2);
double kernel_rcKW(double r, double r1, double r2);

/* core */
std::vector<double> computeRHS(
    KernelType type,
    std::vector<double> S_in,  
    Lattice &lattice);
double integrateKernel(
    KernelType type, 
    int i,
    std::vector<double> S_in, 
    Lattice &lattice);
void stepRK4(KernelType type, Lattice &lattice);

#endif // SOLVER_H