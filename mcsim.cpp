#include "./sim.hpp"
#include "./randnumgen.hpp"
#include "./params.hpp"
#include <iostream>
#include <filesystem>
#include "./enviorment.hpp"
#include <thread>
#include <vector>

using simp = simpars<double>;

int main(int argc, char **argv)
{
    simp pars;

    unsigned seed_base;
    if (pars.seed < 0) {
        std::random_device rd;
        seed_base = rd();  // true entropy
    } else {
        seed_base = pars.seed;
    }

    std::string path = pars.path;
    int Nsim = 10;

    std::vector<std::thread> threads;

    for (int i = 0; i < Nsim; ++i) {
        threads.emplace_back([i, &pars, &path, seed_base]() {
            unsigned thread_seed = seed_base + i;
            randnumgen<double,std::mt19937_64> rng_local;  
            rng_local.seed(thread_seed);
            std::string full_path = path + std::to_string(i);

            mcsim<double, XYmodel<double>> mc(full_path, rng_local);
            mc.prepare_initial_conf();

            auto mc_clone = mc.clone_system();
            mc_clone.prepare_initial_conf();

            mc.run();
            mc_clone.run();
        });
    }

    for (auto &t : threads)
        t.join();

    return 0;
}
