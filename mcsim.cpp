#include "./sim.hpp"
int main(int argc, char **argv)
{
  mcsim<double, particleLJ<double>> mc;
  mc.init_rng(0); // inizializzo il generatore di numeri casuali
  mc.prepare_initial_conf(); // creo la configurazione iniziale
  mc.run(); // lancio la simulazione MC
  return 0;
}
