#include "init_condition.h"

/* Golec-Biernat-Wusthoff dipole S-matrix element */
double gbw_gamma(double r, std::vector<double> &params)
{
    double Qs0 = params[0], gamma = params[1];
    double tmp = std::pow(0.5 * r * Qs0, 2 * gamma);
    return std::exp(-tmp);
}

/* Step-like dipole S-matrix element */
double step(double r, std::vector<double> &params)
{
    return (r < params[0]) ? 1.0 : 0.0;
}

/* McLerran-Venugopalan dipole S-matrix element */
double mv_gamma(double r, std::vector<double> &params)
{
    double Qs0 = params[0], gamma = params[1];
    double tmp = std::pow(0.5 * r * Qs0, 2 * gamma);
    double ln_tmp = std::log(M_E + 1 / (r * LQCD));
    return std::exp(-tmp * ln_tmp);
}

void getInitProfile(Lattice &lattice, InitCond key, std::vector<double> &params)
{
    int NR = lattice.getNR();
    double val;
    for (int i = 0; i < NR; i++)
    {
        double r = lattice.getR(i);
        switch (key)
        {
        case GBW:
            val = gbw_gamma(r, params);
            break;
        case STEP:
            val = step(r, params);
            break;
        case MV:
            val = mv_gamma(r, params);
            break;
        default:
            std::cerr << "Unknown Initial Condition: " << key << std::endl;
            exit(EXIT_FAILURE);
        }
        lattice.setSolution(i, val);
    }
}
