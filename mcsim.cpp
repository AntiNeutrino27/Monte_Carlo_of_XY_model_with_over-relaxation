#include "./sim.hpp"
#include "./randnumgen.hpp"
#include "./params.hpp"
#include <iostream>
#include <filesystem>
#include "./enviorment.hpp"

using simp = simpars<double>;

int main(int argc, char **argv)
{

  simp pars;
  if (pars.seed < 0)
    rng.rseed();
  else
    rng.seed(pars.seed);
  
  mcsim<double, XYmodel<double>> mc;
  mc.prepare_initial_conf(); 
  mcsim<double, XYmodel<double>> mc_clone = mc.clone_system();
  mc_clone.prepare_initial_conf(); 
  mc.run();
  mc_clone.run();


  return 0;
}
