#include "interp.h"

gsl_spline* spline = nullptr;
gsl_interp_accel* acc = nullptr;
/* ====================== SPLINE CONFIG ====================== */
void initSpline(int nr)
{
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, nr);
}
void freeSpline()
{
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

double csplineInterpS(double r, 
    std::vector<double> S_in, 
    Lattice &lattice)
/* cubic spline interpolation */
{
    int nr = lattice.getNR();
    std::vector<double> r_grid(nr);
    std::vector<double> S_values = S_in;
    for (int i = 0; i < nr; i++)
    {
        r_grid[i] = lattice.getR(i);
    }
    // Bounds checking
    if (r < r_grid[0])
        return 1; // S(r~0)=1:  color transparency
    if (r > r_grid[nr - 1])
        return 0; // S(r>>Qs)=0:  saturation

    // Initialize spline with r_grid and S values
    int status = gsl_spline_init(spline, r_grid.data(), S_values.data(), nr);

    if (status != GSL_SUCCESS)
    {
        std::cerr << "Error initializing spline: "
                  << gsl_strerror(status) << std::endl;
    }

    // Evaluate spline at r
    double result = gsl_spline_eval(spline, r, acc);

    return result;
}

double linearInterpS(double r,
    std::vector<double> S_in,  
    Lattice &lattice)
{
    int nr = lattice.getNR();
    std::vector<double> r_grid(nr);
    std::vector<double> S_values = S_in;
    for (int i = 0; i < nr; i++)
    {
        //std::cout << i << " " << nr << std::endl;
        r_grid[i] = lattice.getR(i);
    }
    // Bounds checking
    if (r < r_grid[0])
        return 1; // S(r~0)=1:  color transparency
    if (r > r_grid[nr - 1])
        return 0; // S(r>>Qs)=0:  saturation

    // Find grid index closest to r
    int idx = lattice.getIndexR(r);
    if (lattice.getR(idx) > r)
        idx = idx - 1;

    double r_i = r_grid[idx];
    double r_ip1 = r_grid[idx + 1];

    double S_i = S_in[idx];
    double S_ip1 = S_in[idx + 1];

    // Linear interpolation formula
    double alpha = (r - r_i) / (r_ip1 - r_i);
    return S_i + alpha * (S_ip1 - S_i);
}