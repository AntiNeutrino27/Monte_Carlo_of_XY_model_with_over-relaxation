#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include "./parrallel_tempering.hpp"
#include <iomanip> // for setprecision()
#include <iostream>

namespace fs = std::filesystem;

template<typename ntype, typename model_type>
class sim
{
  using simp = simpars<ntype>;
protected:
  simp pars;  // parameters of the simulation
  ModelPT<ntype, model_type> PT_model; // parallel tempering model

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

  bool create_dir(std::string dir)
  {
    std::error_code ec;

    if (fs::exists(dir, ec))
      return fs::is_directory(dir, ec);

    return fs::create_directories(dir, ec);
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

template<typename ntype, typename model_type>
class mcsim: public sim<ntype, model_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<ntype, model_type>;
  using bc::PT_model, bc::pars, bc::save_mgl_snapshot, bc::create_dir;
  
  // counters used to calculate acceptance rates
  long int tot_met, tot_PT;
  std::vector<long int> rej_counts_met, rej_counts_PT;

  void calc_acceptance_and_adjust(void)
  {
      for(int i=0; i<PT_model.N_T; i++)
      {
        ntype acc_rate = 1.0 - ((ntype)rej_counts_met[i])/((ntype)tot_met);
        if(acc_rate > 0.5 && PT_model.thetas[i]*1.1 < M_PI/2.)
        {
          PT_model.thetas[i] *= 1.1;
        }
        else
        {
          PT_model.thetas[i] /= 1.1;
        }
      }
      restore_rej_count(1);
      tot_met = 0;
  }

  void init_measures(void)
  {
    for(int i=0; i<PT_model.N_T; i++)
    {
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/energy.dat";
      std::fstream f;
      f.open(s,  std::ios::out|std::ios::trunc);
      f.close();
    }

    for(int i=0; i<PT_model.N_T; i++)
    {
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/theta.dat";
      std::fstream f;
      f.open(s,  std::ios::out|std::ios::trunc);
      f.close();
    }

    if (!pars.use_pt)
      return;

    for(int i=0; i<PT_model.N_T; i++)
    {
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/PT_acc_rate.dat";
      std::fstream f;
      f.open(s,  std::ios::out|std::ios::trunc);
      f.close();
    }

  }

 void save_measures(long int t)
  {
    for (int i=0; i<PT_model.N_T; i++)
    {
      model_type model = PT_model.replicas[i];
      ntype totenergy = model.total_energy();
      std::fstream f;
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/energy.dat";
      f.open(s, std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << totenergy << "\n";
      f.close();
    }

    for (int i=0; i<PT_model.N_T; i++)
    {
      ntype theta = PT_model.thetas[i];
      std::fstream f;
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/theta.dat";
      f.open(s, std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << theta << "\n";
      f.close();
    }

    if (!pars.use_pt)
      return;

    for (int i=0; i<PT_model.N_T; i++)
    {
      ntype acc_rate = 1.0 - ((ntype)rej_counts_PT[i])/((ntype)tot_PT);
      std::fstream f;
      std::string s;
      s = pars.path + "/temp_" + std::to_string(i) + "/PT_acc_rate.dat";
      f.open(s, std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << acc_rate << "\n";
      f.close();
    }
    restore_rej_count(2);
    tot_PT = 0;

  }

  void restore_rej_count(int s) // if s=1 reset metropolis swap counts, if s=2 reset PT counts
  {
    if (s == 1)
    {
      for(int i=0; i<pars.N_T; i++)
      {
        rej_counts_met[i] = 0;
      }
    } 
    else if (s == 2)
    {
      for(int i=0; i<pars.N_T; i++)
      {
        rej_counts_PT[i] = 0;
      }
    }
      
  }

  void add_rej_counts(std::vector<long int> loc_rej_counts, int s)
  {
    if (s == 1)
    {
      for(int i=0; i<PT_model.N_T; i++)
      {
        rej_counts_met[i] += loc_rej_counts[i];
      }
    } else if (s == 2)
    {
      for(int i=0; i<PT_model.N_T; i++)
      {
        rej_counts_PT[i] += loc_rej_counts[i];
      }
    }
  }

  bool create_save_dir()
  {
    return create_dir(pars.path);
  }

  bool create_temp_dirs()
  {
    std::string base_dir = pars.path;
    int N_T = pars.N_T;
    for (int i=0; i<N_T; i++)
    {
      std::string temp_dir = base_dir + "/temp_" + std::to_string(i);
      if (!create_dir(temp_dir))
        return false;
    }
    return true;
  }

  bool write_readme()
  {
    std::string filename = pars.path + "/README.txt";
    std::fstream f;
    f.open(filename, std::ios::out|std::ios::trunc);
    if (!f.is_open())
      return false; 
    f << "Parameters of the simulation:\n";
    f << "-----------------------------------\n";
    f << "date and time: " << __DATE__ << " " << __TIME__ << "\n";
    f << "Box size: " << pars.L << "\n";
    f << "Number of particles per replica: " << PT_model.replicas[0].get_N() << "\n";
    f << "Use overrelaxation: " << (pars.use_or ? "yes" : "no") << "\n";
    f << "Initial maximum angle for trial moves in metropolis: " << pars.theta_max << "\n";
    f << "Use parallel tempering: " << (pars.use_pt ? "yes" : "no") << "\n";
    f << "Number of temperatures: " << pars.N_T << "\n";
    f << "Temperature range: " << pars.T_min << " to " << pars.T_max << "\n";
    f << "Number of MC steps: " << pars.totsteps << "\n";
    f << "MC steps period: " << pars.mc_step << "\n";
    f << "Save MGL snapshot every: " << pars.save_mgl_snapshot << "\n";
    f << "Save measures every: " << pars.savemeasure << "\n";
    f << "Adjust steps every: " << pars.adjstps << "\n";
    f << "Maximum adjustment steps: " << pars.maxadjstps << "\n";
    f << "Random seed: " << pars.seed << "\n";
    f << "-----------------------------------\n";

    f.close();
    return true;
  }

 public:
  void prepare_initial_conf(void) 
  {
    if(!pars.use_pt)
    {
      pars.T_max = pars.T_min;
      pars.N_T = 1;
    }
    PT_model.init(pars.T_min, pars.T_max, pars.N_T, pars.L, pars.theta_max);
    create_save_dir();
    write_readme();
    create_temp_dirs();

    rej_counts_met.resize(pars.N_T);
    rej_counts_PT.resize(pars.N_T);
    restore_rej_count(1);
    restore_rej_count(2);
  }

  void run(void) 
  {
      // loop over MC steps
      int t;
      int step = 500;
      tot_met = 0;
      tot_PT = 0;
      std::cout << "Starting MC simulation with Parallel Tempering...\n";
      restore_rej_count(1);
      restore_rej_count(2);
      std::cout << "Preparing initial measures...\n";
      init_measures();
      std::cout << "Running simulation...\n";
      std::cout << pars.totsteps << " MC steps to be performed.\n";
      for (t = 0; t < pars.totsteps; t++)
      {
        if(pars.use_or && t % 10 == 0)
          PT_model.over_relaxation_sweep();

        if (t % pars.mc_step == 0)
        {
          std::vector<long int> loc_rej_counts_met(PT_model.N_T);
          std::vector<long int> loc_rej_counts_PT(PT_model.N_T);
          loc_rej_counts_met = PT_model.metropolis_sweep();
          if(pars.use_pt)
            loc_rej_counts_PT = PT_model.pt_sweep();
          tot_met += PT_model.replicas[0].get_N();
          tot_PT += 1;
          add_rej_counts(loc_rej_counts_met, 1);
          add_rej_counts(loc_rej_counts_PT, 2);
        }
  
        if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
        {
          save_measures(t);
        }

        if (t > 0 && pars.save_mgl_snapshot > 0 && t % pars.save_mgl_snapshot == 0)
        {
          save_mgl_snapshot(t);
        }

        if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 && t < pars.maxadjstps)
        {
          calc_acceptance_and_adjust();
        }
      }
  } 
};
#endif
