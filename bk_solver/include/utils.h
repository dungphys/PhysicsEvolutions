#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <algorithm>

/* ===== QCD Physical Constants ===== */
extern const double ALPHA_S;      // fixed QCD coupling
extern const double ALPHA_FROZEN; // frozen coupling at large r
extern const double LQCD;         // QCD scale parameter
extern const double RHO_CUTOFF;   // Ultraviolet dipole cutoff
extern const int NC;              // Number of colors
extern const int NF;              // Number of active flavors

/* ===== parameter ===== */
extern double C2;


/* ===== Lattice Class ===== */
class Lattice
{
public:
    Lattice(double rmin, double rmax, double dr,
            double dtheta, double ymax, double dY);
    ~Lattice() = default;
    /* === Grid Accessors === */
    double getR(int i) const;
    double getTheta(int j) const;
    double getY(int n) const;
    double getSolution(int i) const;

    int getIndexR(double r) const;
    int getIndexTheta(double theta) const;
    int getIndexY(double y) const;

    /* === Grid Parameters === */
    double getDR() const { return dr; }
    double getDTheta() const { return dtheta; }
    double getDY() const { return dY; }

    /* === Grid Sizes === */
    int getNR() const { return NR; }
    int getNTheta() const { return NTheta; }
    int getNY() const { return NY; }

    /* === Utilities ===*/
    void setSolution(int i, double value);
    void writeSolutionToFile(const std::string &filename) const;
    void printGridInfo() const;

private:
    /* ===== Grid arrays ===== */
    std::vector<double> R;
    std::vector<double> Theta;
    std::vector<double> Y;
    std::vector<double> solution;

    /* ===== Grid parameters ===== */
    int NR, NTheta, NY;
    double rmin, rmax;
    double dr;
    double dtheta;
    double ymax, dY;

    void validateGridParameters() const;
};

#endif // UTILS_H