#ifndef UTILS_H
#define UTILS_H

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

/* ===== Lattice Structure ===== */
class Lattice
{
public:
    Lattice(double xmin, double xmax, double dx,
            double tmax, double dt, double D = 1.0);
    /** Write spatial and temporal grids to files */
    void writeGridsToFiles(const std::string &fileX_name,
                           const std::string &fileT_name);
    /** Write solution to file */
    void writeSolutionToFiles(const std::string &file_name);
    /** Get spatial coordinate at index i*/
    double getX(int i) const { return X[i]; }
    /** Get temporal coordinate at index n*/
    double getT(int n) const { return T[n]; }
    /** Get solution at index i*/
    double getSolution(int i) const { return solution[i]; }
    double getDiffusionCoef() const { return diffusion_coef; }
    double getDX() const { return dx; }
    double getDT() const { return dt; }
    int getNX() const { return NX; }
    int getNT() const { return NT; }
    /** Set solution value at index i */
    void setSolution(int i, double value) { solution[i] = value; }

private:
    std::vector<double> X;
    std::vector<double> T;
    std::vector<double> solution;
    int NX, NT;

    double xmin, xmax, dx;
    double tmax, dt;
    double diffusion_coef;

    // Additional private members can be added here
};

#endif // UTILS_H