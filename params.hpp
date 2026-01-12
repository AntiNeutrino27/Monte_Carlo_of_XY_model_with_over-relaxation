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
  long int maxadjstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps, mc_step; // savemeasure=steps at which save measures, totsteps = simulations steps, outstps steps print something on current simulation status
  int seed; // -1 means random
  int L; // box size
  ntype theta_max; // initial maximum angle for trial moves in metropolis


  std::string path; // path to save data 

  simpars()
    {
      std::ifstream ri;
      path = "simulations/simulation1";
      totsteps = 409600*2;
      save_mgl_snapshot = 100;
      mc_step = 10;
      T_min = 0.265;
      T_max = 0.45;
      N_T = 35;
      L = 24;
      seed = 1;
      maxadjstps = totsteps;
      adjstps = 2000;
      savemeasure = 100;
      theta_max = 3.14159265359 * 0.05;
    }
};
#endif
