#include "solver.h"

/*-----------------------------------------------------*/
/*                   reaction kernels                  */
/*-----------------------------------------------------*/

double reactionFunction(
    double u,
    std::vector<double> &params,
    ReactionType type)
{
    switch (type)
    {
    case CANONICAL:
    {
        double r = params[0];
        return r * u * (1 - u);
    }
    case ZELDOVICH:
    {
        double r = params[0];
        double thres = params[1];
        if (u <= thres)
            return 0;
        else
            return r * u * (1 - u);
    }
    case POWER_LAW:
    {
        double r = params[0];
        int k = (int)params[1];
        return r * pow(u, k) * (1 - u);
    }
    case WEAK_ALLE:
    {
        double r = params[0];
        double delta = params[1];
        return r * pow(u, 1 + delta) * (1 - u);
    }
    case QUADRATIC:
    {
        double r = params[0];
        int n = (int)params[1];
        return r * u * (1 - pow(u, n));
    }
    default:
    {
        std::cerr << "Unknown reaction type: " << type << std::endl;
        exit(EXIT_FAILURE);
    }
    }
}

/*--------------------------------------------------------------*/
/*                        explicit solver                      */
/*--------------------------------------------------------------*/

void advanceFKPP_ExplicitEuler(
    Lattice &lattice,
    std::vector<double> &u_new,
    ReactionType reacType,
    std::vector<double> &params)
{
    int NX = lattice.getNX();
    double dX = lattice.getDX();
    double dT = lattice.getDT();
    double diffusion_coef = lattice.getDiffusionCoef();
    double alpha = diffusion_coef * dT / (dX * dX);
    double u_prev, u_curr, u_next;
    double diffusion, reaction;

    // Boundary conditions at index 0
    u_prev = 1;
    u_curr = lattice.getSolution(0);
    u_next = lattice.getSolution(1);

    diffusion = alpha * (u_prev - 2.0 * u_curr + u_next);
    reaction = dT * reactionFunction(u_curr, params, reacType);
    u_new[0] = u_curr + diffusion + reaction;

    // Boundary conditions at index NX-1
    u_prev = lattice.getSolution(NX - 2);
    u_curr = lattice.getSolution(NX - 1);
    u_next = 0;

    diffusion = alpha * (u_prev - 2.0 * u_curr + u_next);
    reaction = dT * reactionFunction(u_curr, params, reacType);
    u_new[NX - 1] = u_curr + diffusion + reaction;

    // Interior points
    for (int i = 1; i < NX - 1; i++)
    {
        u_prev = lattice.getSolution(i - 1);
        u_curr = lattice.getSolution(i);
        u_next = lattice.getSolution(i + 1);

        diffusion = alpha * (u_prev - 2.0 * u_curr + u_next);
        reaction = dT * reactionFunction(u_curr, params, reacType);

        u_new[i] = u_curr + diffusion + reaction;
    }
}

/*--------------------------------------------------------------*/
/*                        implicit solver                      */
/*--------------------------------------------------------------*/
/* *
 * Solve tridiagonal system using Thomas algorithm (LU decomposition)
 * System: a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = rhs[i]
 *
 * @param a     subdiagonal (lower)
 * @param b     main diagonal
 * @param c     superdiagonal (upper)
 * @param rhs   right-hand side
 * @param u     solution (output)
 * @param nx    system size
 */
void solveTridiagonal(
    std::vector<double> &a,
    std::vector<double> &b,
    std::vector<double> &c,
    std::vector<double> &rhs,
    std::vector<double> &u,
    int nx)
{
    // Forward elimination
    std::vector<double> c_star(nx), d_star(nx);

    c_star[0] = c[0] / b[0];
    d_star[0] = rhs[0] / b[0];

    for (int i = 1; i < nx - 1; i++)
    {
        double denom = b[i] - a[i] * c_star[i - 1];
        c_star[i] = c[i] / denom;
        d_star[i] = (rhs[i] - a[i] * d_star[i - 1]) / denom;
    }

    // Last row (no superdiagonal)
    d_star[nx - 1] = (rhs[nx - 1] - a[nx - 1] * d_star[nx - 2]) /
                     (b[nx - 1] - a[nx - 1] * c_star[nx - 2]);

    // Back substitution
    u[nx - 1] = d_star[nx - 1];
    for (int i = nx - 2; i >= 0; i--)
    {
        u[i] = d_star[i] - c_star[i] * u[i + 1];
    }
}

void advanceFKPP_ImplicitEuler(
    Lattice &lattice,
    std::vector<double> &u_new,
    ReactionType reacType,
    std::vector<double> &params,
    int max_iters)
{
    int nx = lattice.getNX();
    double D = lattice.getDiffusionCoef();
    double dt = lattice.getDT();
    double dx = lattice.getDX();
    double alpha = D * dt / (dx * dx);

    std::vector<double> a(nx), b(nx), c(nx), rhs(nx);

    // Interior points
    for (int i = 1; i < nx - 1; i++)
    {
        a[i] = -alpha;
        b[i] = 1.0 + 2.0 * alpha;
        c[i] = -alpha;
    }

    // Boundary conditions
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    /**/
    a[nx - 1] = 0.0;
    b[nx - 1] = 1.0;
    c[nx - 1] = 0.0;

    // Initial guess for nonlinear iteration
    for (int i = 0; i < nx; i++)
        u_new[i] = lattice.getSolution(i);
    // Nonlinear Picard iteration
    for (int iter = 0; iter < max_iters; iter++)
    {
        // Update the right-hand side with the reaction term
        for (int i = 0; i < nx; i++)
        {
            double reaction = dt * reactionFunction(u_new[i], params, reacType);
            rhs[i] = lattice.getSolution(i) + reaction;
        }
        // Solve system
        std::vector<double> u_old = u_new;
        solveTridiagonal(a, b, c, rhs, u_new, nx);
        // Check convergence
        double error = 0.0;
        for (int i = 0; i < nx; i++)
            error += std::abs(u_new[i] - u_old[i]) / u_old[i];

        if (error < 1e-5)
            break;
    }

    // Update solution
    for (int i = 0; i < nx; i++)
        lattice.setSolution(i, u_new[i]);
}