/**
 *
 * This header defines the FisherKolmogorov3D class, which implements a parallel finite element solver for the Fisher-Kolmogorov equation in three dimensions.
 * The solver supports various types of diffusion tensors (isotropic, radial, circumferential) and is designed for use with distributed meshes
 * and MPI parallelism.
 *
 * The parameters for the problem and the simulation are configurable through SimulationParams structure, 
 * which can be initialized from a CSV file containing key-value pairs.
 *
 * Key features:
 * - Support for non-linear, possibly anisotropic diffusion tensors.
 * - Homogeneous Neumann boundary conditions.
 * - Implicit Euler time integration (theta-method with theta=1, not configurable).
 * - Newton's method for solving the non-linear system at each time step.
 * - Parallel assembly and solution using Trilinos and deal.II distributed DoF handling.
 * - Output of results for visualization.
 *
 * Main components:
 * - SimulationParams: Structure holding all simulation parameters, with a constructor for initialization from a map.
 * - readCSVData: Utility function to read simulation parameters from a CSV file.
 * - FisherKolmogorov3D: Main class encapsulating the problem setup, assembly, solution, and output.
 *   - FunctionD: TensorFunction representing the (possibly anisotropic) diffusion tensor D(x), configurable by diffusion type.
 *   - FunctionAlpha: Scalar function for the reaction coefficient alpha, possibly spatially varying.
 *   - FunctionU0: Initial condition function, with support for seeding in a region around a specified point.
 *   - setup(): Initializes mesh, DoF handler, finite element, and vectors.
 *   - solve(): Advances the solution in time using Newton's method at each step.
 *   - assemble_system(): Assembles the Jacobian and residual for the current Newton iteration.
 *   - solve_linear_system(): Solves the linearized system for the Newton update.
 *   - solve_newton(): Performs Newton iterations for a single time step.
 *   - output(): Writes solution data for visualization (made a posteriori with Paraview).
 *
 * Usage:
 * 1. Construct a SimulationParams object (e.g., via readCSVData).
 * 2. Instantiate FisherKolmogorov3D with the parameters.
 * 3. Call setup() to initialize.
 * 4. Call solve() to run the simulation.
 *
 * Notes:
 * - The code is designed for extensibility; diffusion types and initial seeding can be customized.
 * - The implementation assumes homogeneous Neumann boundary conditions.
 * - The mesh and polynomial degree are specified in the simulation parameters.
 * - The code is parallelized using MPI and supports distributed memory execution.
 * 
 */
#ifndef FISHER_KOLMOGORV_3D_HPP
#define FISHER_KOLMOGORV_3D_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>

#include "DiffusionType.hpp"
#include "SimulationParams.hpp"

using namespace dealii;

// Class representing the non-linear diffusion problem.
class FisherKolmogorov3D {
   public:
    // Physical dimension (3D).)
    static constexpr unsigned int dim = 3;

    // Function representing the diffusion tensor D(x).
    class FunctionD : public TensorFunction<2, dim, double> {
       public:
        FunctionD(const SimulationParams &params)
            : TensorFunction<2, dim>(),
              d_axn(params.d_axn),
              d_ext(params.d_ext),
              diffusion_type(params.diffusion_type),
              seed_point(params.seed_point) {}

        virtual Tensor<2, dim, double> value(
            const Point<dim> &p) const override {
            Tensor<2, dim> matrix;
            Tensor<1, dim> n = compute_normal_vector(p);

            for (unsigned int i = 0; i < dim; i++) {
                matrix[i][i] = d_ext;
            }

            for (unsigned int i = 0; i < dim; ++i) {
                for (unsigned int j = 0; j < dim; ++j) {
                    matrix[i][j] += d_axn * n[i] * n[j];
                }
            }

            return matrix;
        }

       private:
       // Compute the normal vector based on the diffusion type.
        Tensor<1, dim> compute_normal_vector(const Point<dim> &p) const {
            Tensor<1, dim> n;
            Tensor<1, dim> t;  // Circumferential

            switch (diffusion_type) {
                case ISOTROPIC: {
                    for (unsigned int i = 0; i < dim; ++i) n[i] = 0.0;
                    break;
                }

                case RADIAL: {
                    n = p - seed_point;
                    double norm = n.norm();
                    if (norm > 1e-12) n /= norm;
                    break;
                }

                case CIRCUMFERENTIAL: {
                    t = p - seed_point;
                    double norm_t = t.norm();
                    if (norm_t > 1e-12) t /= norm_t;

                    n[0] = -t[1];
                    n[1] = t[0];
                    if (dim == 3) n[2] = 0;
                    double norm = n.norm();
                    if (norm > 1e-12) n /= norm;
                    break;
                }
            }
            return n;
        }

        const double d_axn;
        const double d_ext;
        const DiffusionType diffusion_type;
        const Point<dim> seed_point;
    };

    // Function for the reaction coefficient alpha.
    class FunctionAlpha : public Function<dim> {
       public:
        FunctionAlpha(const SimulationParams &params) : alpha(params.alpha) {}

        virtual double value(
            const Point<dim> & /*p*/,
            const unsigned int /*component*/ = 0) const override {
            return alpha;
        }

       private:
        const double alpha;
    };

    // Function for initial conditions.
    class FunctionU0 : public Function<dim> {
       public:
        FunctionU0(const SimulationParams &params) : 
            seed_point(params.seed_point),
            seed_radius(params.seed_radius),
            initial_seed(params.initial_seed) {}

        virtual double value(
            const Point<dim> &p,
            const unsigned int /*component*/ = 0) const override {

                // Spherical seed condition
            if ((p - seed_point).norm_square() <= seed_radius * seed_radius) {
                return initial_seed;
            } else {
                return 0.0;
            }
        }

       private:
        const Point<dim> seed_point;
        const double seed_radius;
        const double initial_seed;
    };

    // Constructor. The input parameters are read from a SimulationParams object.
    // The mesh is read from a file, and the output directory is set.
    // The polynomial degree and time step are also set from the parameters.
    FisherKolmogorov3D(SimulationParams params)
        : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
          mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
          pcout(std::cout, mpi_rank == 0),
          simulation_params(params),
          D(params),
          alpha(params),
          u_0(params),
          T(params.T),
          mesh_file_name(params.mesh_file_name),
          output_directory(params.output_directory),
          r(params.degree),
          deltat(params.deltat),
          mesh(MPI_COMM_WORLD) {}

    // Initialization.
    void setup();

    // Solve the problem.
    void solve();

    const auto &get_mesh() const { return mesh; }
    const auto &get_dof_handler() const { return dof_handler; }

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
    //                        MPI parallel                           //
    ///////////////////////////////////////////////////////////////////

    // Number of MPI processes.
    const unsigned int mpi_size;

    // This MPI process.
    const unsigned int mpi_rank;

    // Parallel output stream.
    ConditionalOStream pcout;

    ///////////////////////////////////////////////////////////////////
    //                    Problem Definition                         //
    ///////////////////////////////////////////////////////////////////

    const SimulationParams simulation_params;

    // mu_0 coefficient.
    FunctionD D;

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

    // Path to the mesh file.
    const std::string mesh_file_name;

    // Path to output directory.
    const std::string output_directory;

    // Polynomial degree.
    const unsigned int r;

    // Time step.
    const double deltat;

    // Mesh.
    parallel::fullydistributed::Triangulation<dim> mesh;

    // Finite element space.
    std::unique_ptr<FiniteElement<dim>> fe;

    // Quadrature formula.
    std::unique_ptr<Quadrature<dim>> quadrature;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;

    // DoFs relevant to the current process (including ghost DoFs).
    IndexSet locally_relevant_dofs;

    // Jacobian matrix.
    TrilinosWrappers::SparseMatrix jacobian_matrix;

    // Residual vector.
    TrilinosWrappers::MPI::Vector residual_vector;

    // Increment of the solution between Newton iterations.
    TrilinosWrappers::MPI::Vector delta_owned;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::Vector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::Vector solution;

    // System solution at previous time step.
    TrilinosWrappers::MPI::Vector solution_old;
};

#endif /* FISHIER_KOLMOGORV_3D_HPP */