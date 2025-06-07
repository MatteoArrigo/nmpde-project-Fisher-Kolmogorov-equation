#include <fstream>
#include "FisherKolmogorov3d.hpp"

/*
 * The main function initializes the MPI environment, reads simulation parameters from a CSV file,
 * prints them, creates a FisherKolmogorov3D problem instance with those parameters, sets up the problem,
 * and then solves it.
 */

// Main function.
int main(int argc, char *argv[]) {
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
    
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    double total_start_time = MPI_Wtime();

    // Determine parameter file path
    std::string param_file_path = "../data/params.csv";
    if (argc > 1) {
        param_file_path = argv[1];
    }
    SimulationParams simulation_params = readCSVData(param_file_path);
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    simulation_params.printParameters();
    
    double setup_start_time = MPI_Wtime();
    FisherKolmogorov3D problem(simulation_params);
    problem.setup();
    double setup_end_time = MPI_Wtime();
    double setup_time = setup_end_time - setup_start_time;

    double solve_start_time = MPI_Wtime();
    problem.solve();
    double solve_end_time = MPI_Wtime();
    double solve_time = solve_end_time - solve_start_time;

    double total_end_time = MPI_Wtime();
    double total_time = total_end_time - total_start_time;

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
        std::ofstream outfile;
        outfile.open("../data/timing_results.csv", std::ios::app);
        outfile << mpi_size << "," 
                << problem.get_mesh().n_global_active_cells() << ","  // Number of elements
                << problem.get_dof_handler().n_dofs() << ","          // Number of DoFs
                << setup_time << ","
                << solve_time << ","
                << total_time << "\n";
        outfile.close();

        std::cout << "Execution Time with " << mpi_size << " tasks:" << std::endl;
        std::cout << "  Number of elements: " << problem.get_mesh().n_global_active_cells() << std::endl;
        std::cout << "  Number of DoFs: " << problem.get_dof_handler().n_dofs() << std::endl;
        std::cout << "  Setup Time: " << setup_time << " seconds" << std::endl;
        std::cout << "  Solve Time: " << solve_time << " seconds" << std::endl;
        std::cout << "  Total Time: " << total_time << " seconds" << std::endl;
    }
}
