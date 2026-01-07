#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"
#include <fstream>
template<typename ntype>
class simpars
{
public:
  ntype T_min, T_max; // temperature and pressure
  int N_T; // number of temperatures in PT
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps, mc_step; // savemeasure=steps at which save measures, totsteps = simulations steps, outstps steps print something on current simulation status
  int seed; // -1 means random
  int L; // box size
  double sigma, epsilon, mass; // Lennard-Jones parameters
  double dt, deltra, vmax; // parameter of MC moves


  std::string path; // path to save data 

  simpars()
    {
      std::ifstream ri;
      path = "simulations/simulation1";
      totsteps = 409600*2;
      mc_step = 10;
      T_min = 0.265;
      T_max = 0.45;
      N_T = 35;
      L = 24;
    }
};
#endif
