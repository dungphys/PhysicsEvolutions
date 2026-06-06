#include "utils.h"

Lattice::Lattice(double xmin, double xmax, double dx,
                 double tmax, double dt, double D)
    : xmin(xmin), xmax(xmax), dx(dx),
      tmax(tmax), dt(dt), diffusion_coef(D)
{
    // Calculate number of points
    NX = static_cast<int>(std::round((xmax - xmin) / dx)) + 1;
    NT = static_cast<int>(std::round(tmax / dt)) + 1;

    // Fill spatial grid
    X.resize(NX);
    for (int i = 0; i < NX; ++i)
        X[i] = xmin + i * dx;

    // Fill temporal grid
    T.resize(NT);
    for (int n = 0; n < NT; ++n)
        T[n] = n * dt;

    // Initialize solution to zero
    solution.resize(NX, 0.0);
}

void Lattice::writeGridsToFiles(const std::string &fileX_name,
                                const std::string &fileT_name)
{
    // Write spatial grid to file
    std::ofstream fileX(fileX_name);
    if (!fileX.is_open())
        throw std::runtime_error("Could not open file: " + fileX_name);
    fileX << std::scientific << std::setprecision(4);
    for (int i = 0; i < NX; ++i)
        fileX << X[i] << "\n";
    std::cout << "Spatial grid written to " << fileX_name << std::endl;
    fileX.close();

    // Write temporal grid to file
    std::ofstream fileT(fileT_name);
    if (!fileT.is_open())
        throw std::runtime_error("Could not open file: " + fileT_name);
    fileT << std::scientific << std::setprecision(4);
    for (int n = 0; n < NT; ++n)
        fileT << T[n] << "\n";
    std::cout << "Temporal grid written to " << fileT_name << std::endl;
    fileT.close();
}

void Lattice::writeSolutionToFiles(const std::string &file_name)
{
    std::ofstream file(file_name);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + file_name);
    file << std::scientific << std::setprecision(4);
    for (int i = 0; i < NX; ++i)
        file << solution[i] << "\n";
    std::cout << "Solution written to " << file_name << std::endl;
    file.close();
}