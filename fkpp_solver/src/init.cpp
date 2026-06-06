#include "init.h"

/* Step function*/
void heaviside(Lattice &lattice)
{
    int NX = lattice.getNX();
    for (int i = 0; i < NX; i++)
    {
        if (lattice.getX(i) > 0.)
            lattice.setSolution(i, 0);
        else
            lattice.setSolution(i, 1);
    }
    std::cout << "IC: Heaviside step function generated successfully!\n"
              << std::endl;
}

/* Exponetial function */
int exponential(Lattice &lattice, double slope)
{
    int NX = lattice.getNX();

    if (slope <= 0)
    {
        std::cerr << "Error: init. slope should be positive" << std::endl;
        return -1;
    }

    for (int i = 0; i < NX; i++)
    {
        double temp = -expm1(-exp(-slope * lattice.getX(i)));
        lattice.setSolution(i, temp);
    }
    std::cout << "IC: exponential step function generated successfully!\n"
              << std::endl;
    return 0;
}

/* Logistic function*/
int logistic(Lattice &lattice, double slope)
{
    int NX = lattice.getNX();
    if (slope <= 0)
    {
        std::cerr << "Error: init. slope should be positive" << std::endl;
        return -1;
    }
    for (int i = 0; i < NX; i++)
    {
        double temp = 1 / (1 + exp(-slope * lattice.getX(i)));
        lattice.setSolution(i, temp);
    }
    std::cout << "IC: logistic step function generated successfully!\n"
              << std::endl;
    return 0;
}
