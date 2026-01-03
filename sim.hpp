#ifndef _SIMCLASS_
#define _SIMCLASS_
/*
 * X = done 
 * * = todo
 * E = exercise
 * OE = optional exercise
 *
 *  ===== base class (common stuff for MC and MD) =====  [ DD/MM/2025 ] 
 *
 *   [X]  calcenergyi()
 *   [X]  calctotenergy()
 *   [X]  pbc() [ to apply pbc ]
 *   [X]  save_mgl_snapshot() save configurations in mgl format
 *   [X]  prepare_initial_conf() 
 *   [X]  init_rng()
 *   [X]  run() method (main loop) 
 *
 *  ===== MC class =====
 *
 *  NTV ensemble [ 25/11/2025 ]
 *
 *   [X]  alpha()
 *   [X]  acc()
 *   [X]  move_NTV()
 *
 *   [X]  calc_acceptance and adjust()
 *   [X]  init_measures() truncate measure files
 *   [X]  save_measures() append measures
 *
 *   [OE]  calculate tail corrections to energy and pressure
 *
 *  optional exercises:
 *
 *   // in "sim" base class
 *   [OE]  save_snapshot() save configurations 
 *   [OE]  save_restart()
 *   
 *   NPT ensemble [ 09/12/2025 ]
 *
 *   // in "mcsim" derived class
 *
 *   [E] init_measures() and save_measures(): add density saving
 *
 *   [E] modify run() method to attempt box moves too
 *   
 *   [ 09/12/2025 ] EXERCISES FOR NPT ensemble:
 *
 *   [E] restore_all_pars (restore positions of all particles)
 *
 *   [E] store_all_pars (store position of all particles)
 *
 *   [E] alpha_box() i.e. trial box move
 *   
 *   [E] acc_box() i.e. acceptance of trial box move
 *
 *** MD class [ DD/MM/2025 ]
 *
 *   [ ]  init_measures() truncate measure files
 *
 *   [ ]  calcK() calculate the total kinetic energy
 *   
 *   [ ]  save_measures() save measures
 *   
 *   [ ]  gauss() (gauss stochastic variable using Box-Muller)
 *   
 *   [ ]  calc_forces()  <---
 *   
 *   [ ]  run() (implement velocity verlet as a symmetric factorization of propagator) <---
 *   
 *   [ ]  prepare_initial_conf() (initialize velocities) <---
 *   
 *
 */
#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <iomanip> // for setprecision()

template<typename ntype, typename particle_type>
class sim
{
  using simp = simpars<ntype>;
protected:
  simp pars;  // parametri per composizione
  std::vector<particle_type> parts;
  // if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt=0) 
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < pars.Np; j++) // pars.Np è il numero totale di particelle
        {
          if (opt==1 && i >= j)
            continue;
          if (i==j)
            continue;
          // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
          enei += parts[i].vij(parts[j], pars.L); 
          // pars.L è un vettore con i lati del box 
        }
      return enei;
    }
  
  ntype totenergy(void)
    {
      ntype ene=0.0;

      for (auto i=0; i < static_cast<int>(parts.size()); i++)
        {
          ene+=calcenergyi(i, 1);
        }
      return ene;
    }

  void pbc(int i) 
    {
      auto Dr = parts[i].r;
      Dr = mulcw(pars.L,rint(divcw(Dr,pars.L))); // L*rint(Dr/L)
      parts[i].r = parts[i].r - Dr;
    }

  void save_mgl_snapshot(long int t) 
    {
       std::fstream f;
       std::string s;

       s = "cnf-" + std::to_string(t) + ".mgl";
       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << parts[i].r(0) << " " << parts[i].r(1) << " " <<
             parts[i].r(2) << " @ " << pars.sigma*0.5 << "\n";
         }
       f.close();
    }
 
public:
  void prepare_initial_conf(void) 
    {
      // SC 
      int ix, iy, iz;
      int cc=0;
      pars.Np = pars.nx*pars.ny*pars.nz;
      parts.resize(pars.Np);
      ntype vcell = pow(pars.sigma,3.0);
      ntype rhomax = 1.0/vcell;
      ntype sf;
      sf = cbrt(rhomax/pars.rho);
      pars.L ={ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
      std::cout << "Box size is " << pars.L << "\n";
      ntype clen = sf*pars.sigma;
      pars.L *= clen;
      std::cout << "Initial density is " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
      for (ix = 0; ix < pars.nx; ix++)
        for (iy = 0; iy < pars.ny; iy++)
          for (iz = 0; iz < pars.nz; iz++)
            {
              parts[cc].r = {ix*clen, iy*clen, iz*clen};
              parts[cc].r -= pars.L*0.5;
              parts[cc].set_sigma(pars.sigma);
              parts[cc].set_epsilon(pars.epsilon);
              parts[cc].set_rcut(pars.rc);
              cc++;
            }
      // ...or BCC or FCC lattice
 
    }

  void set_sim_type(int type)
    {
      pars.simtype = type;
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
  using bc::calcenergyi, bc::parts, bc::pars, 
        bc::totenergy, bc::save_mgl_snapshot, bc::pbc;
  
  // counters used to calculate acceptance rates
  long int tot_tra, tot_vol, tra_rej, vol_rej;

  void alpha(int i)
    {
      // displacement pars.deltra
      // move particle i by using as max displacemente paramater pars.deltra
      // and apply PBC 
      pvector<ntype,3> delr;
      delr ={pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5)};
      parts[i].tra_move(delr);
      pbc(i);
    }

  void acc(int i, ntype eno)
    {
      // ACCEPT OR REJECT TRIAL MOVE ACCORDING TO METROPOLIS CRITERIUM
      // We assume kB=1 (Boltzmann constants) so that beta=1/T. 
      // Calculate new energy and accept move with probability
      // min {1, exp(-\beta deltaU}, i.e.
      // if deltaU < 0 accept the move
      // otherwise generate a random number xi in [0,1] 
      // and if xi < exp(-beta*deltaU) accept the move otherwise reject it 
      ntype enn = calcenergyi(i);
      ntype delu= enn-eno;
      ntype xi = rng.ranf();
      if (delu > 0.0 && xi >= exp(-delu/pars.T))
        {
          // reject move
          tra_rej++;
          parts[i].restore();
        }
    }
 
  void move_NTV(int i)
    {
      // 1) calculate energy of particle i
      ntype eno;
      eno = calcenergyi(i);
      // 2) store position of particle i
      parts[i].store();
      // 3) trial move 
      alpha(i);
      // 4) acceptance 
      acc(i, eno);
    }

  void move_box(void)// <------------------------------------------------
    {
      ntype DG, fact; // DG is \Delta G as discussed during lectures (see pdf of lectures)
      store_all_pars(); // store all particle positions before attempting a box move
      alpha_box(DG, fact);
      acc_box(DG, fact);
    }

  void alpha_box(ntype& DG, ntype& fact)
    {
      // trial box move
      //
      // 1) choose randomly a new volume
      //
      // 2) scale all particle positions and box size by factor "fact"
      // 
      // 3) calculate new total interaction energy
      //
      // 4) calculate \DeltaG (see pdf of lectures)

      // wrote your code here
      ntype tot_energy_old = totenergy();
      ntype deltaV = pars.vmax * (rng.ranf()-0.5);
      ntype V = pars.L(0)*pars.L(1)*pars.L(2);
      fact = pow((V + deltaV)/V, 1./3);
      for(int i = 0; i<pars.Np; i++){
        parts[i].r *= fact;
      }
      ntype tot_eneryg_new = totenergy();
      ntype deltaU = tot_eneryg_new - tot_eneryg_new;
      DG = pars.P * deltaV + deltaU - pars.Np * pars.T * 1/3 * log(fact);

    }
              
  void acc_box(ntype DG, ntype fact)
    {
      // accept or reject box move
      //
      // 1) calculate e^{-\Delta G} (see pdf of lectures)
      // 2) generate a random numner \xi and check wheter to accept the box trial move      
      // 3) if move is rejected:
      //    i)   restore all particle positions thourgh method restore_all_pars()
      //    ii)  restore box size
      //    iii) update counter vol_rej of rejected box moves for calculating acceptance ratio 

      // EXERCISE: write your code here
      ntype xi = rng.ranf();
      if (DG > 0.0 && xi >= exp(-DG/pars.T))
        {
          restore_all_pars();
          vol_rej++;
        }else{
          pars.L *= fact;
        }

    }
 

 void restore_all_pars()
   {
     // restore all particle positions 
     // EXERCISE: write your code here
     for(int i =0; i<pars.Np; i++){
        parts[i].restore();
     }
   }


 void store_all_pars()
   {
     // store all particle position 
     // EXERCISE: write your code here
     for(int i =0; i<pars.Np; i++){
        parts[i].store();
     }
   }

  void calc_acceptance_and_adjust(void)
    {
      // CALC ACCEPTANCE RATES AND 
      // ADJUST pars.deltra 
      ntype r_tra, r_vol;
      if (tot_tra > 0)
        {
          // estimate acceptance rate
          r_tra = ((double)(tot_tra - tra_rej))/tot_tra;
          std::cout << "rate tra: " << r_tra << " deltra=" << pars.deltra << "\n";
          if (r_tra > 0.5)
            {
              pars.deltra *= 1.1;
            }
          else
            {
              pars.deltra /= 1.1;
            }
          tot_tra=tra_rej=0;
        }

      // adjust maximum volume "displacement" in alpha_box
      // so that acceptance rate of box move is around 0.5
      if (tot_vol > 0) // <------------------------------------------------
        {
          r_vol = ((ntype)(tot_vol - vol_rej)) / tot_vol;
          std::cout << "rate vol: " << r_vol << " vmax=" << pars.vmax << "\n";
          if (r_vol > 0.5) 
            {
              pars.vmax *= 1.1;
           }
          else
            {
              pars.vmax /= 1.1;
            }
          tot_vol=vol_rej=0.0;
        }
    }

  void init_measures(void)
    {
      // open files in writing mode to truncate them to 0
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::trunc);
      f.close();
      if (pars.simtype==1) // simtype == 1 means an NPT simluations
        {
          std::fstream f;
          f.open("density.dat", std::ios::out|std::ios::trunc);
          f.close();
        }

      // init your own measures here
    }

 void save_measures(long int t)
    {
      std::fstream f;
      f.open("energy.dat", std::ios::out|std::ios::app);
      // save potential energy per particle
      f << t << " " << totenergy()/pars.Np << "\n";
      f.close();
      if (pars.simtype==1) // 1 means NPT simulation, save density in this case
        {
          f.open("density.dat", std::ios::out|std::ios::app);
          // save potential energy per particle
          f << t << " " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
          f.close();
        }
      // you can add your own measures here
    }

 public:
  void run(void) 
    {
      // loop over MC steps
      int i, t, ip;
      ntype iv;
      tot_tra = tra_rej = tra_rej = vol_rej = 0;
      init_measures();
      for (t = 0; t < pars.totsteps; t++)
        {
          // ATTEMPT TO MOVE ALL PARTICLES:
          // Np move attempts = 1 MC 
          if(pars.simtype==1){
            iv = rng.ranf();
            if(iv<1./pars.Np){
              move_box();
              tot_vol++;
              continue;
            }
          }
          for (i=0; i < pars.Np; i++)
            {
              // EXERCISE: choose between box move and particle move, write your own code here
              ip = rng.ranf()*pars.Np;
              move_NTV(ip);
              tot_tra++;
            }
        
          if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
            {
              save_measures(t);
            }

          if (t > 0 && t % pars.outstps == 0) 
            {
              std::cout << "Step #" << t << "\n";
              // NOTE: to compare with Johnson paper data the internal energy must be calculated (i.e. by including
              // the kinetic contribution)
              std::cout << "total energy per particle is " << totenergy()/pars.Np << "\n";
              if(pars.simtype==1){
                std::cout << "density is " << pars.Np/(pars.L(0)*pars.L(1)*pars.L(2)) << "\n";
              }
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
