#include "solver.h"

/* ===== KERNEL EVALUATION ===== */

double alpha_rc(double r)
/* QCD Running coupling at scale r² */
{
    double denom = (11 * NC - 2 * NF) * std::log(4 * C2 / (r * r * LQCD * LQCD));
    if (denom == 0)
        return ALPHA_FROZEN;
    return std::min({(12 * M_PI / denom), ALPHA_FROZEN});
}

std::string getKernelName(KernelType kernel_type)
/* Get particular kernel name based on kernel type */
{
    switch (kernel_type)
    {
    case LO:
        return "Leading-order kernel";
    case RC_PARENT:
        return "Running-coupling kernel with parent-dipole prescription";
    case RC_MIN:
        return "Running-coupling kernel with smallest-dipole prescription";
    case RC_BAL:
        return "Running-coupling kernel with Balitsky prescription";
    case RC_KW:
        return "Running-coupling kernel with Kovchegov-Weigert prescription";
    default:
        throw std::runtime_error("Unknown kernel type");
    }
}

double kernel_lo(double r, double r1, double r2)
/* Leading-order kernel */
{
    double r_sq = r * r, r1_sq = r1 * r1, r2_sq = r2 * r2;
    double constant = NC * ALPHA_S / (2 * M_PI * M_PI);
    double rPart = r_sq / (r1_sq * r2_sq);
    return constant * rPart;
}

double kernel_rcParent(double r, double r1, double r2)
/* Parent-dipole rc kernel */
{
    double r_sq = r * r, r1_sq = r1 * r1, r2_sq = r2 * r2;
    double rcConstant = NC * alpha_rc(r) / (2 * M_PI * M_PI);
    double rPart = r_sq / (r1_sq * r2_sq);
    return rcConstant * rPart;
}

double kernel_rcMin(double r, double r1, double r2)
/* Smallest-dipole rc kernel */
{
    double r_sq = r * r, r1_sq = r1 * r1, r2_sq = r2 * r2;
    double smallest = std::min({r, r1, r2});
    double rcConstant = NC * alpha_rc(smallest) / (2 * M_PI * M_PI);
    double rPart = r_sq / (r1_sq * r2_sq);
    return rcConstant * rPart;
}

double kernel_rcBal(double r, double r1, double r2)
/* Balitsky rc kernel */
{
    double r_sq = r * r, r1_sq = r1 * r1, r2_sq = r2 * r2;
    double rcConstant = NC * alpha_rc(r) / (2 * M_PI * M_PI);
    double alpha_ratio = alpha_rc(r1) / alpha_rc(r2);
    double rPart_1 = r_sq / (r1_sq * r2_sq);
    double rPart_2 = (1 / r1_sq) * (alpha_ratio - 1);
    double rPart_3 = (1 / r2_sq) * ((1 / alpha_ratio) - 1);
    return rcConstant * (rPart_1 + rPart_2 + rPart_3);
}

double kernel_rcKW(double r, double r1, double r2)
/* KW rc kernel */
{
    double r_sq = r * r, r1_sq = r1 * r1, r2_sq = r2 * r2;
    if ((r1_sq + r2_sq - r_sq) == 0)
        return 1;
    double A = (r1_sq + r2_sq) - 4 * r1_sq * r2_sq / (r1_sq + r2_sq - r_sq);
    double constant = NC / (2 * M_PI * M_PI);
    double R, R_sq;
    if (r1 == r2)
    {
        R_sq = r1 * r2 * std::exp(-A / (2 * r1_sq));
        R = std::sqrt(R_sq);
    }
    else
    {
        double power = A / (r1_sq - r2_sq);
        R_sq = r1 * r2 * std::pow(r2 / r1, power);
        R = std::sqrt(R_sq);
    }

    double rPart_1 = alpha_rc(r1) / r1_sq + alpha_rc(r2) / r2_sq;
    double alpha_ratio = alpha_rc(r1) * alpha_rc(r2) / alpha_rc(R);
    double rPart_2 = alpha_ratio * (r1_sq + r2_sq - r_sq) / (r1_sq * r2_sq);
    return constant * (rPart_1 + rPart_2);
}

double evaluateKernel(
    KernelType type, 
    double r, double r1, double r2)
/* Evaluate kernel based on selected type */
{
    switch (type)
    {
    case LO:
        return kernel_lo(r, r1, r2);
    case RC_PARENT:
        return kernel_rcParent(r, r1, r2);
    case RC_MIN:
        return kernel_rcMin(r, r1, r2);
    case RC_BAL:
        return kernel_rcBal(r, r1, r2);
    case RC_KW:
        return kernel_rcKW(r, r1, r2);
    default:
        throw std::runtime_error("Unknown kernel type");
    }
}

/* ===== CORE SOLVER ===== */
std::vector<double> computeRHS(
    KernelType type, 
    std::vector<double> S_in,
    Lattice &lattice
)
/* Compute right-hand-sider integral */
{
    int NR = lattice.getNR();
    std::vector<double> dS(NR);

    /* For each spatial point r_i */
    for (int i = 0; i < NR; i++)
    {   
        //std::cout << i << std::endl;
        dS[i] = integrateKernel(type, i, S_in, lattice);
    }

    return dS;
}

double integrateKernel(
    KernelType type, 
    int i, 
    std::vector<double> S_in,
    Lattice &lattice)
/* integral kernel at index i */
{
    double r = lattice.getR(i);
    double dr = lattice.getDR();
    int NR = lattice.getNR();
    int NTheta = lattice.getNTheta();
    double integral = 0.0, angular_integral;
    double S_r1, S_r2, S_r = S_in[i];

    // Loop over r1 on the grid
    for (int i1 = 0; i1 < NR; i1++)
    {
        double r1 = lattice.getR(i1);
        if (r1 <= RHO_CUTOFF)
            continue;
        // Angular integration over theta1
        S_r1 = S_in[i1];
        angular_integral = 0.0;
        for (int m = 0; m < NTheta; m++)
        {
            double theta1 = (2.0 * M_PI * m) / NTheta;
            // Law of cosines
            double r2_sq = r * r + r1 * r1 - 2.0 * r * r1 * std::cos(theta1);
            double r2 = std::sqrt(r2_sq);

            if (r2 <= RHO_CUTOFF)
                continue;

            // Kernel value
            double K_val = evaluateKernel(type, r, r1, r2);
            //S_r2 = csplineInterpS(r2, S_in, lattice);7
            S_r2 = linearInterpS(r2, S_in, lattice);
            // S only depends on r, not theta
            // Integrand: K(r, r1, r2, C2) * [S(r1)*S(r2) - S(r)]
            double integrand = K_val * (S_r1 * S_r2 - S_r);

            angular_integral += integrand;
        }

        angular_integral *= (2.0 * M_PI) / NTheta;

        // Add r1 dr1 weight for radial integration
        integral += (angular_integral * r1 * r1 * (std::exp(dr) - 1));
    }

    return integral;
}

void stepRK4(KernelType type, Lattice &lattice)
{
    int NR = lattice.getNR();
    int dY = lattice.getDY();
    std::vector<double> k1, k2, k3, k4;
    std::vector<double> S_temp(NR);
    std::vector<double> S_in(NR);

    for (int i = 0; i < NR; i++) 
        S_in[i] = lattice.getSolution(i);  
    std::cout << "Zeroth step ok!" << std::endl;

    k1 = computeRHS(type, S_in, lattice);
    std::cout << "First step ok! " << k1[400] << std::endl;

    for (int i = 0; i < NR; i++)
        S_temp[i] = S_in[i] + 0.5 * dY * k1[i];
    k2 = computeRHS(type, S_temp, lattice);
    std::cout << "Second step ok! " << k2[400] << std::endl;

    for (int i = 0; i < NR; i++)
        S_temp[i] = S_in[i] + 0.5 * dY * k2[i];
    k3 = computeRHS(type, S_temp, lattice);
    std::cout << "Third step ok!" << k3[400] << std::endl;

    for (int i = 0; i < NR; i++)
        S_temp[i] = S_in[i] + dY * k3[i];
    k4 = computeRHS(type, S_temp, lattice);
    std::cout << "Fourth step ok!" << k4[400] << std::endl;

    for (int i = 0; i < NR; i++)
    {
        S_in[i] += (dY / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        // S[i] = std::max(0.0, std::min(1.0, S[i]));
    }
    for (int i=0;i<NR;i++) 
    {
        lattice.setSolution(i, S_in[i]);
    }
}