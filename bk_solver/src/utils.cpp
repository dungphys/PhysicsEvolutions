/**/
#include "utils.h"

/* ===== QCD Physical Constants ===== */
const double ALPHA_S = 0.18;     // fixed QCD coupling
const double ALPHA_FROZEN = 0.7; // frozen coupling at large r
const double LQCD = 0.241;       // QCD scale parameter
const double RHO_CUTOFF = 1e-42; // Ultraviolet dipole cutoff
const int NC = 3;                // Number of colors
const int NF = 3;                // Number of active flavors

/* ===== parameter ===== */
double C2;

Lattice::Lattice(double rmin, double rmax, double dr,
                 double dtheta, double ymax, double dY)
    : rmin(rmin), rmax(rmax), dr(dr), dtheta(dtheta), ymax(ymax), dY(dY)
{
    validateGridParameters();

    NR = (int)std::round(std::log(rmax / rmin) / dr) + 1;
    NTheta = (int)std::round(2.0 * M_PI / dtheta);
    NY = (int)std::round(ymax / dY) + 1;

    R.resize(NR);
    Theta.resize(NTheta);
    Y.resize(NY);
    solution.resize(NR, 0.0);

    for (int i = 0; i < NR; i++)
        R[i] = rmin * std::exp(i * dr);
    for (int i = 0; i < NTheta; i++)
        Theta[i] = i * dtheta;
    for (int n = 0; n < NY; n++)
        Y[n] = n * dY;

    // std::cout << "Lattice: NR=" << NR << ", NTheta=" << NTheta
    //           << ", NY=" << NY << std::endl;
}

void Lattice::validateGridParameters() const
{
    if (rmax <= rmin)
        throw std::invalid_argument("rmax > rmin required");
    if (dr <= 0 || dtheta <= 0 || dY <= 0)
        throw std::invalid_argument("Grid steps must be positive");
}

/* ========================================================= */
/* ===                   Grid Accessors                   === */
/* ========================================================= */

double Lattice::getR(int i) const
{
    if (i < 0 || i >= NR)
        throw std::out_of_range("R index out of bounds");
    return R[i];
}

double Lattice::getTheta(int i) const
{
    if (i < 0 || i >= NTheta)
        throw std::out_of_range("Theta index out of bounds");
    return Theta[i];
}

double Lattice::getY(int n) const
{
    if (n < 0 || n >= NY)
        throw std::out_of_range("Y index out of bounds");
    return Y[n];
}

double Lattice::getSolution(int i) const
{
    if (i < 0 || i >= (int)solution.size())
    {
        throw std::out_of_range("Solution index out of bounds");
    }
    return solution[i];
}


int Lattice::getIndexR(double r) const
{
    int idx = (int)std::round(std::log(r / rmin) / dr);
    return std::max(0, std::min(NR - 1, idx));
}

int Lattice::getIndexTheta(double theta) const
{
    double normalized = std::fmod(theta, 2.0 * M_PI);
    if (normalized < 0)
        normalized += 2.0 * M_PI;
    int idx = (int)std::round(normalized / dtheta);
    return idx % NTheta;
}

int Lattice::getIndexY(double y) const
{
    int idx = (int)std::round(y / dY);
    return std::max(0, std::min(NY - 1, idx));
}


/* ========================================================= */
/* ===                     Utilities                     === */
/* ========================================================= */

void Lattice::setSolution(int i, double value)
{
    if (i < 0 || i >= (int)solution.size())
        throw std::out_of_range("Solution index out of bounds");
    solution[i] = value;
}

void Lattice::writeSolutionToFile(const std::string &filename) const
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
        throw std::runtime_error("Cannot open: " + filename);

    for (int i = 0; i < NR; i++)
    {
        outfile << std::scientific << std::setprecision(6)
                << R[i] << "\t"
                << 1 - solution[i] << std::endl;
        /// outfile << "\n";
    }
    outfile.close();
    std::cout << "Solution written to " << filename << std::endl;
}

void Lattice::printGridInfo() const
{
    std::cout << "\n===== Lattice Info =====" << std::endl;
    std::cout << "R: [" << rmin << ", "
              << rmax << "], NR=" << NR << ", dR=" << dr << std::endl;
    std::cout << "NTheta=" << NTheta << ", dTheta=" << dtheta << std::endl;
    std::cout << "Y: [0, " << ymax << "], NY=" << NY << ", dY=" << dY << std::endl;
    std::cout << "========================\n"
              << std::endl;
}