#ifndef SOLVER_H
#define SOLVER_H

#include "utils.h"

typedef enum
{
    CANONICAL, // standard FKPP
    ZELDOVICH, // ignition-type nonlinearity
    POWER_LAW, // power-law nonlinearity
    WEAK_ALLE, // weak Allee effect
    QUADRATIC  // quadratic nonlinearity
} ReactionType;
/**/
double reactionFunction(
    double u,
    std::vector<double> &params,
    ReactionType type);
/**/
void advanceFKPP_ExplicitEuler(
    Lattice &lattice,
    std::vector<double> &u_new,
    ReactionType reacType,
    std::vector<double> &params);
/**/
void solveTridiagonal(
    std::vector<double> &a,
    std::vector<double> &b,
    std::vector<double> &c,
    std::vector<double> &rhs,
    std::vector<double> &u,
    int nx);
/**/
void advanceFKPP_ImplicitEuler(
    Lattice &lattice,
    std::vector<double> &u_new,
    ReactionType reacType,
    std::vector<double> &params,
    int max_iters);

#endif