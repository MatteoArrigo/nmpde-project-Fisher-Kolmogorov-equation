#include <fstream>

#include "FisherKolmogorov1d.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{

  const unsigned int degree = 1;

  const double T      = 20.0;
  const double deltat = 0.1;

  FisherKolmogorov1D problem(200, degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}