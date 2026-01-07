#include "./sim.hpp"

#include <iostream>
#include <filesystem>
#include "./enviorment.hpp"

int main(int argc, char **argv)
{
  mcsim<double, XYmodel<double>> mc;
  mc.init_rng(0); // inizializzo il generatore di numeri casuali
  mc.prepare_initial_conf(); // creo la configurazione iniziale
  mc.run(); // lancio la simulazione MC

  return 0;
}
