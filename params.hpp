#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"
#include <fstream>
template<typename ntype>
class simpars
{
public:
  int nx, ny, nz; /* nx*ny*nz particelle */
  double T, P; // temperature and pressure
  int Np; // numero di particelle
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps, mc_step; // savemeasure=steps at which save measures, totsteps = simulations steps, outstps steps print something on current simulation status
  double rho, rc; // density
  int simtype; // simulation type (see below)
  int seed; // -1 means random
  pvector<double, 3> L; // box
  double sigma, epsilon, mass; // Lennard-Jones parameters
  double dt, deltra, vmax; // parameter of MC moves


  std::string path; // path to save data 

  simpars()
    {
      std::ifstream ri;
      path = "simulations/simulation1";
      totsteps = 409600*2;
      mc_step = 10;
    }
};
#endif
