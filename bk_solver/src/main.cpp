/**/
#include "utils.h"
#include "init_condition.h"
#include "interp.h"
#include "solver.h"

std::string formatDouble(double val, int decimals = 3)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimals) << val;
    return oss.str();
}

int main()
{
    try
    {
        double rmin = 1e-30, rmax = 1e3, dr = 0.18;
        double dtheta = 2.0 * M_PI / 80;
        double ymax = 5.0, dY = 0.02;
        std::string fname;

        Lattice lattice(rmin, rmax, dr, dtheta, ymax, dY);
        lattice.printGridInfo();
        initSpline(lattice.getNR());
        // ===== INITIAL CONDITION =====
        std::vector<double> params(2);
        params[0] = 1.0;
        params[1] = 1.0;
        C2 = 0.5;
        InitCond key = MV;
        getInitProfile(lattice, key, params);
        fname = "../bin/bk_solution_" + formatDouble(lattice.getY(0)) + ".txt";
        lattice.writeSolutionToFile(fname);

        // ===== TIME EVOLUTION =====
        int nsteps = (int)(ymax / dY);

        KernelType type = LO;
        for (int n = 0; n < nsteps; n++)
        {
            stepRK4(type, lattice);
            double Y = (n + 1) * dY;
            std::cout << "Y = " << std::fixed << std::setprecision(2)
                      << Y << " / " << ymax << std::endl;
            fname = "../bin/bk_solution_" + formatDouble(lattice.getY(n+1)) + ".txt";
            lattice.writeSolutionToFile(fname);
        }
        freeSpline();
        std::cout << "\nDone! Solution: solution_running_coupling.txt" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}