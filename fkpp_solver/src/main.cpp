#include "utils.h"
#include "init.h"
#include "solver.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

std::string formatDouble(double val, int decimals = 3)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimals) << val;
    return oss.str();
}

int main()
{
    // Example usage of utility functions
    /** Initialize the lattice */
    Lattice lattice(-400, 100000, 0.2, 50.0, 0.01, 0.5);

    /* Initialize the solution */
    std::vector<double> u_new(lattice.getNX());
    ReactionType reacType = CANONICAL;
    std::vector<double> params = {1.0}; // r = 1.

    heaviside(lattice);
    std::string filename = "../bin/solution_" + formatDouble(lattice.getT(0)) + ".txt";
    lattice.writeSolutionToFiles(filename);
    lattice.writeGridsToFiles("../bin/grid_x.txt", "../bin/grid_t.txt");

    std::cout << "Starting iteration..." << std::endl;
    // advance the solution using the explicit Euler method with an exponential initial condition
    for (int n = 0; n < lattice.getNT(); n++)
    {
        advanceFKPP_ImplicitEuler(lattice, u_new, reacType, params, 10);

        for (int i = 0; i < lattice.getNX(); ++i)
            lattice.setSolution(i, u_new[i]);
        // save the solution every 10 steps
        if (((n + 1) % 10) == 0)
        {
            filename = "../bin/solution_" + formatDouble(lattice.getT(n + 1)) + ".txt";
            lattice.writeSolutionToFiles(filename);
        }
        std::cout << "Time step " << n + 1 << "/" << lattice.getNT() << " completed." << std::endl;
    }
    return 0;
}