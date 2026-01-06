#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include "./enviorment.hpp"
#include "./parrallel_tempering.hpp"
#include <iomanip> // for setprecision()

template<typename ntype, typename model_type>
class sim
{
  using simp = simpars<ntype>;
protected:
  simp pars;  // parameters of the simulation
  ModelPT<ntype, model_type> PT_model;

  void save_mgl_snapshot(long int t) 
  {
    for (int i=0; i<PT_model.N_T; i++)
    {
      model_type model = PT_model.replicas[i];
      std::fstream f;
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/cnf-" + std::to_string(t) + ".mgl";
      f.open(s, std::ios::out|std::ios::trunc);
      for (int j=0; j < model.get_N(); j++)
      {
        f << model.parts[j].s(0) << " " << model.parts[j].s(1) << "\n";
      }
      f.close();
    }
  }
 
public:
  void prepare_initial_conf(void) 
    {

    }
 
  void init_rng(int n) 
    {
      if (n < 0)
        rng.rseed();
      else
        rng.seed(n);
    }

  void run(void) 
    {
      // intentionally void
    }; 
};

template<typename ntype, typename particle_type>
class mcsim: public sim<ntype, particle_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<ntype, particle_type>;
  using bc::PT_model, bc::pars, bc::save_mgl_snapshot;
  
  // counters used to calculate acceptance rates
  long int tot_tra;
  std::vector<long int> rej_counts;

  void calc_acceptance_and_adjust(void)
    {
      for(int i=0; i<PT_model.N_T; i++)
      {
        ntype acc_rate = 1.0 - ((ntype)rej_counts[i])/((ntype)tot_tra);
        if(acc_rate > 0.5)
        {
          PT_model.thetas[i] *= 1.1;
        }
        else
        {
          PT_model.thetas[i] /= 1.1;
        }
      }
      restore_rej_count();
      tot_tra = 0;
    }

  void init_measures(void)
  {
      // // open files in writing mode to truncate them to 0
      // std::fstream f;
      // f.open("energy.dat", std::ios::out|std::ios::trunc);
      // f.close();
      // if (pars.simtype==1) // simtype == 1 means an NPT simluations
      //   {
      //     std::fstream f;
      //     f.open("density.dat", std::ios::out|std::ios::trunc);
      //     f.close();
      //   }

      // // init your own measures here
  }

 void save_measures(long int t)
  {
      // std::fstream f;
      // f.open("energy.dat", std::ios::out|std::ios::app);
      // // save potential energy per particle
      // f << t << " " << totenergy()/pars.Np << "\n";
      // f.close();
      // if (pars.simtype==1) // 1 means NPT simulation, save density in this case
      //   {
      //     f.open("density.dat", std::ios::out|std::ios::app);
      //     // save potential energy per particle
      //     f << t << " " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
      //     f.close();
      //   }
      // // you can add your own measures here
  }

  void restore_rej_count()
  {
      for(int i=0; i<PT_model.N_T; i++)
      {
        rej_counts[i] = 0;
      }
  }

  void add_rej_counts(std::vector<long int> loc_rej_counts)
  {
      for(int i=0; i<PT_model.N_T; i++)
      {
        rej_counts[i] += loc_rej_counts[i];
      }
  }

 public:
  void prepare_initial_conf(void) 
  {
    rej_counts.resize(PT_model.N_T);
    restore_rej_count();
  }

  void run(void) 
  {
      // loop over MC steps
      int i, t, ip;
      ntype iv;
      tot_tra = 0;
      restore_rej_count();
      init_measures();
      for (t = 0; t < pars.totsteps; t++)
        {
          PT_model.over_relaxation_sweep();
          if (t % pars.mc_p == 0)
            {
              std::vector<long int> loc_rej_counts(PT_model.N_T);
              loc_rej_counts = PT_model.metropolis_sweep();
              PT_model.pt_sweep();
              tot_tra += PT_model.replicas[ip].get_N();
              add_rej_counts(loc_rej_counts);
            }
  
          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }

          if (t > 0 && pars.save_mgl_snapshot > 0 && 
              t % pars.save_mgl_snapshot == 0)
            {
              save_mgl_snapshot(t);
            }

          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && 
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }
  } 
};
#endif
