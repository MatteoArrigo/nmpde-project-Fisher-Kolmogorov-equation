cmake_minimum_required(VERSION 3.12.0)
project(FisherKolmogorovBrain LANGUAGES CXX C)

include(cmake-common.cmake)

add_executable(simulation1d src/1d/simulation1d.cpp src/1d/FisherKolmogorov1d.cpp)
deal_ii_setup_target(simulation1d)

add_executable(simulation3d src/3d/simulation3d.cpp src/3d/FisherKolmogorov3d.cpp src/3d/SimulationParams.cpp)
deal_ii_setup_target(simulation3d)