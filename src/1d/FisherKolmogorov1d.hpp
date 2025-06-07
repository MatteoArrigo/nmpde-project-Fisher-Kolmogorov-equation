/**
 *
 * This class implements the numerical solution of the Fisher-Kolmogorov equation using the finite element method (FEM) in one spatial dimension. 
 * In this context, the equation is used to simulate the concentration of misfolded proteins in the human brain.
 * 
 * Main features:
 * - 1D mesh generation with uniform subdivisions.
 * - Finite element space defined with polynomial degree r.
 * - Time-stepping using a fixed time step deltat.
 * - Weak formulation of the Fisher-Kolmogorov equation.
 * - Assembly of the system matrices and vectors.
 * - Solution of the linearized system using a conjugate gradient solver.
 * - Linearization via Newton's method to handle nonlinearity.
 * - Output of results to files for post-processing.
 * 
 * Class structure:
 * - Internal classes for the definition of physical coefficients and initial conditions.
 * - Public methods for initialization, solution, and output.
 * - Protected methods for system assembly, linear solution, and Newton's method management.
 * - Data members for mesh management, finite element space, degree of freedom handler, system matrices, and vectors.
 * 
 * Main constructor parameters:
 * @param N Number of mesh subdivisions.
 * @param r Polynomial degree of the finite elements.
 * @param T Final simulation time.
 * @param deltat Time step.
 * 
 * Example usage:
 * @code
 * FisherKolmogorov1D solver(N, r, T, deltat);
 * solver.setup();
 * solver.solve();
 * @endcode
 *
 */
#ifndef FISHIER_KOLMOGORV_1D_HPP
#define FISHIER_KOLMOGORV_1D_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>


using namespace dealii;

inline const std::string output_directory = "output";

// Class representing the non-linear diffusion problem.
class FisherKolmogorov1D {
   public:
    // Physical dimension (1D, 2D, 3D)
    static constexpr unsigned int dim = 1;

    // Function for the mu_0 coefficient.
    class FunctionD : public Function<dim> {
       public:
        virtual double value(
            const Point<dim> & /*p*/,
            const unsigned int /*component*/ = 0) const override {
            return 0.0001 * 4;
        }
    };

    // Function for the mu_1 coefficient.
    class FunctionAlpha : public Function<dim> {
       public:
        virtual double value(
            const Point<dim> & /*p*/,
            const unsigned int /*component*/ = 0) const override {
            return 1.0;
        }
    };

    // Function for initial conditions.
    class FunctionU0 : public Function<dim> {
       public:
        virtual double value(
            const Point<dim> &p,
            const unsigned int /*component*/ = 0) const override {
            if (p[0] == 0.0)
                return 0.1;
            else
                return 0.0;
        }
    };

    /**
     * Initializes the Fisher-Kolmogorov 1D problem with the specified mesh size,
     * polynomial degree, final time, and time step.
     *
     * @param N_ Number of segments in the 1D mesh.
     * @param r_ Polynomial degree for the finite element space.
     * @param T_ Final simulation time.
     * @param deltat_ Time step size.
     */
    FisherKolmogorov1D(const unsigned &N_, const unsigned int &r_, const double &T_,
                  const double &deltat_)
        : T(T_),
          N(N_),
          r(r_),
          deltat(deltat_)
          {}

    // Initialization.
    void setup();

    // Solve the problem.
    void solve();

   protected:
    // Assemble the tangent problem.
    void assemble_system();

    // Solve the linear system associated to the tangent problem.
    void solve_linear_system();

    // Solve the problem for one time step using Newton's method.
    void solve_newton();

    // Output.
    void output(const unsigned int &time_step) const;

    ///////////////////////////////////////////////////////////////////
    //                    Problem Definition                         //
    ///////////////////////////////////////////////////////////////////

    // mu_0 coefficient.
    FunctionD d;

    // mu_1 coefficient.
    FunctionAlpha alpha;

    // Initial conditions.
    FunctionU0 u_0;

    // Current time.
    double time;

    // Final time.
    const double T;

    ///////////////////////////////////////////////////////////////////
    //                      Discretization                           //
    ///////////////////////////////////////////////////////////////////

    // Number of triangulation segments in the mesh
    const unsigned N;

    // Polynomial degree.
    const unsigned int r;

    // Time step.
    const double deltat;

    // Mesh.
    Triangulation<dim> mesh;

    // Finite element space.
    std::unique_ptr<FiniteElement<dim>> fe;

    // Quadrature formula.
    std::unique_ptr<Quadrature<dim>> quadrature;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // Sparsity pattern.
    SparsityPattern sparsity;

    // Jacobian matrix.
    SparseMatrix<double> jacobian_matrix;

    // Residual vector.
    Vector<double> residual_vector;

    // Increment of the solution between Newton iterations.
    Vector<double> delta;

    // System solution (including ghost elements).
    Vector<double> solution;

    // System solution at previous time step.
    Vector<double> solution_old;
};

#endif /* FISHIER_KOLMOGORV_1D_HPP */