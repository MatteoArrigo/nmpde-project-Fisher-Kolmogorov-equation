#ifndef SIMULATION_PARAMS_H
#define SIMULATION_PARAMS_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include "DiffusionType.hpp"
#include <deal.II/base/point.h>


/**
 * 
 * This file provides functionality to:
 *   - Read simulation parameters from a CSV file and store them as key-value pairs in a map.
 *   - Construct a SimulationParams object from the parsed configuration map.
 *   - Print the loaded simulation parameters in a formatted, right-aligned table.
 *
 * Main Components:
 * ----------------
 * 1. readCSVData:
 *    - Reads a CSV file where each line contains a key-value pair separated by a comma (with no spaces).
 *    - Strips leading and trailing whitespace from both keys and values.
 *    - Stores the pairs in a std::map<std::string, std::string>.
 *    - Returns a SimulationParams object constructed from the parsed map.
 *
 * 2. SimulationParams Constructor:
 *    - Initializes all simulation parameters from the provided configuration map.
 *    - Converts string values to appropriate types (e.g., double, unsigned long, enum).
 *    - Expects the following keys to be present in the map:
 *        - "mesh_file_name": Path to the mesh file (std::string).
 *        - "output_dir": Output directory for simulation results (std::string).
 *        - "degree": Polynomial degree for the simulation (unsigned long).
 *        - "T": Final simulation time (double).
 *        - "deltat": Time step size (double).
 *        - "d_axn": Axonal diffusion coefficient (double).
 *        - "d_ext": Extracellular diffusion coefficient (double).
 *        - "diffusion_type": Type of diffusion (integer, cast to DiffusionType enum).
 *              Possible values: 0 (ISOTROPIC), 1 (CIRCUMFERENTIAL), 2 (RADIAL).
 *  *           See DiffusionType.hpp for more details.       
 *        - "seed_point_x", "seed_point_y", "seed_point_z": Coordinates of the initial seed point (double).
 *        - "alpha": Reaction rate parameter (double).
 *        - "initial_seed": Initial seed value (double).
 *
 * 3. printParameters: prints all simulation parameters in a human-readable, right-aligned table format.
 *
 * Usage Example:
 * --------------
 * @code
 * SimulationParams params = readCSVData("config.csv");
 * params.printParameters();
 * @endcode
 *
 * @note
 * - The CSV file must not contain a header row.
 * - Each line should be of the form: key,value
 * - Whitespace around keys and values is ignored.
 * - All required keys must be present in the CSV file.
 *
 */


struct SimulationParams {
    const std::string mesh_file_name;
    const std::string output_directory;
    const unsigned int degree;
    const double T;
    const double deltat;
    const double d_axn;
    const double d_ext;
    const double alpha;
    const DiffusionType diffusion_type;
    const dealii::Point<3> seed_point;
    const double seed_radius;
    const double initial_seed;

    // Constructor to initialize from a key-value map
    SimulationParams(const std::map<std::string, std::string> &configMap);

    // Function to print parameters in a formatted way
    void printParameters() const;
};

// Function to read a CSV file and store key-value pairs in a map
SimulationParams readCSVData(const std::string &filename);

#endif /* SIMULATION_PARAMS_H */