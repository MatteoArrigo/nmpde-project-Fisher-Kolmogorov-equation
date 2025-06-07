#include "SimulationParams.hpp"

// Function to read a CSV file and store key-value pairs in a map
SimulationParams readCSVData(const std::string &filename) {
    std::map<std::string, std::string> config;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string key, value;

        if (std::getline(ss, key, ',') && std::getline(ss, value)) {
            // Strip whitespace from key and value before adding to the map
            auto strip = [](std::string &s) {
                size_t start = s.find_first_not_of(" \t\r\n");
                size_t end = s.find_last_not_of(" \t\r\n");
                if (start == std::string::npos) {
                    s.clear();
                } else {
                    s = s.substr(start, end - start + 1);
                }
            };
            strip(key);
            strip(value);
            config[key] = value;
        }
    }

    return SimulationParams(config);
}

// Constructor to initialize SimulationParams from a configuration map
SimulationParams::SimulationParams(const std::map<std::string, std::string> &configMap)
    : mesh_file_name(configMap.at("mesh_file_name")),
      output_directory(configMap.at("output_dir")),
      degree(std::stoul(configMap.at("degree"))),
      T(std::stod(configMap.at("T"))),
      deltat(std::stod(configMap.at("deltat"))),
      d_axn(std::stod(configMap.at("d_axn"))),
      d_ext(std::stod(configMap.at("d_ext"))),
      alpha(std::stod(configMap.at("alpha"))),
      diffusion_type(static_cast<DiffusionType>(std::stoi(configMap.at("diffusion_type")))),
      seed_point{
          std::stod(configMap.at("seed_point_x")),
          std::stod(configMap.at("seed_point_y")),
          std::stod(configMap.at("seed_point_z")),
      },
      seed_radius(std::stod(configMap.at("seed_radius"))),
      initial_seed(std::stod(configMap.at("initial_seed")))
    {}
    
// Print parameters as a right-aligned table
void SimulationParams::printParameters() const {
    const int label_width = 20;
    std::cout << "Parameters\n"
                << std::right << std::setw(label_width) << "mesh_file_name: " << mesh_file_name << "\n"
                << std::right << std::setw(label_width) << "output_directory: " << output_directory << "\n"
                << std::right << std::setw(label_width) << "degree: " << degree << "\n"
                << std::right << std::setw(label_width) << "T: " << T << "\n"
                << std::right << std::setw(label_width) << "deltat: " << deltat << "\n"
                << std::right << std::setw(label_width) << "d_axn: " << d_axn << "\n"
                << std::right << std::setw(label_width) << "d_ext: " << d_ext << "\n"
                << std::right << std::setw(label_width) << "alpha: " << alpha << "\n"
                << std::right << std::setw(label_width) << "diffusion_type: " << DiffusionName(diffusion_type) << "\n"
                << std::right << std::setw(label_width) << "seed_point: "
                << "(" << seed_point[0] << ", " << seed_point[1] << ", " << seed_point[2] << ")\n"
                << std::right << std::setw(label_width) << "seed_radius: " << seed_radius << "\n"
                << std::right << std::setw(label_width) << "initial_seed: " << initial_seed << "\n\n";
}