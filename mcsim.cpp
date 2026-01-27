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
  
  std::string path = pars.path;
  std::string full_path;

  for(int i = 0; i < 10; i++) {
    full_path = path + std::to_string(i);
    mcsim<double, XYmodel<double>> mc(full_path);
    mc.prepare_initial_conf(); 
    mcsim<double, XYmodel<double>> mc_clone = mc.clone_system();
    mc_clone.prepare_initial_conf(); 
    mc.run();
    mc_clone.run();
  }
  

  return 0;
}
