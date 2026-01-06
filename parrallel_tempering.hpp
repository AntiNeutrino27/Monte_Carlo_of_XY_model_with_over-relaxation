#ifndef _PT_
#define _PT_

#include "./pvector.hpp"
#include "./particle.hpp"
#include <vector>
#include <algorithm>


template<typename ntype, typename model_type>
class ModelPT
{
public:
    ntype T_min, T_max;
    int N_T;
    std::vector<ntype> temperatures;
    std::vector<model_type> replicas;
    std::vector<ntype> thetas;
    ModelPT(ntype T_min, ntype T_max, int N_T_, int Lbox): T_min(T_min), T_max(T_max), N_T(N_T)
    {
        temperatures.resize(N_T);
        replicas.resize(N_T);
        thetas.resize(N_T);

        // Initialize temperatures with geometric spacing
        for(int i=0; i<N_T; i++){
            temperatures[i] = T_min * pow(T_max/T_min, ntype(i)/(N_T-1.0));
        }

        // Initialize replicas
        model1 = model_type(Lbox);
        replicas[0] = model1;
        for(int i=1; i<N_T; i++){
            replicas[i] = model1.make_replica();
        }
    }

    ModelPT<ntype> make_copy() const
    {  
        int Lbox = replicas[0].get_L();
        ModelPT<ntype> copy(T_min, T_max, N_T, Lbox);

        // copying bonds
        for(int i=0; i<N_T; i++){
            copy.replicas[i] = replicas[i].make_replica();
        }
        return copy;
    }

    void over_relaxation_sweep()
    {
        for(int i=0; i<N_T; i++){
            replicas[i].or_sweep();
        }
    }

    std::vector<long int> metropolis_sweep()
    {
        std::vector<int> rej_counts;
        rej_counts.resize(N_T);
        for(int i=0; i<N_T; i++){
            rej_counts[i] = replicas[i].metropolis_sweep(temperatures[i], thetas[i]);
        }
        return rej_counts;
    }

    void pt_sweep()
    {
        for(int i=0; i<N_T-1; i++){
            ntype E1 = replicas[i].total_energy();
            ntype E2 = replicas[i+1].total_energy();
            ntype dBeta = 1.0/temperatures[i+1] - 1.0/temperatures[i];
            ntype delta = (E2 - E1) * dBeta;
            ntype xi = rng.ranf();
            if(delta < 0.0 || xi < exp(-delta)){
                std::swap(replicas[i], replicas[i+1]);
            }
        }
    }
};

#endif