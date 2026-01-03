#ifndef _PARTICLE_
#define _PARTICLE_
/*
 * X = done 
 * * = todo
 * E = exercise
 * OE = optional exercise 
 *
 * particle class
 *
 * MC [ 25/11/2025 ]
 *
 * [X] store/restore methods
 * 
 * [X] tra_move()
 * 
 * MD [ DD/MM/2025 ]
 * 
 * [ ] expiLq  <----
 *
 * [ ] expiLp  <----
 *
 * particleLJ class 
 *
 * [X] vij()
 * 
 * [X] set_vcut()
 *
 * [X] fij()  <----
 *
 * particleSS class
 *
 * [X] vij()
 * 
 * [OE] set_vcut()
 * 
 * [OE] fij)
 *
 */
#include "./pvector.hpp"
template<typename ntype>
class particle
{
  // ADD later forces and velocities
  pvector<ntype,3> rold; // to store particle's position 
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, m;
  pvector<ntype,3> r; 

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void);
 
  void tra_move(pvector<ntype,3> delr)
    {
      r+=delr; 
    }

  void store()
    {
     rold = r;
    }

  void restore()
    {
      r = rold;
    }
  void set_sigma(ntype sig)
    {
      sigma = sig;
    }

  void set_epsilon(ntype eps)
    {
      epsilon = eps;
    }

  void set_rcut(ntype rcut)
    {
      rc = rcut;
    }

  particle()
    {
      sigma=1.0;
      epsilon=1.0;
      rc=2.5;
      vcut = 0.0;
    }
};

template <typename ntype>
class particleLJ: public particle<ntype>
{
  using pt=particle<ntype>;
protected:
  using pt::vcut;
public:
  using pt::r, pt::epsilon, pt::rc, pt::sigma;
  ntype vij(particleLJ P, pvector<ntype,3> L)
    {
      ntype ene;
      pvector<ntype,3> Dr;
     
      Dr = r - P.r;
      // MINIMUM IMAGE CONVENTION 
      Dr = Dr - mulcw(L,rint(divcw(Dr,L))); // Dr - L*rint(Dr/L)
      ntype rsq, rn = Dr.norm();
      rsq = rn*rn; 
      if (rsq < rc*rc) // interaction potential cut-off
        ene = 4.0*epsilon*(pow(sigma/rn,12.0)-pow(sigma/rn,6));
      else
        ene = 0.0;
      return ene;

    }
  pvector<ntype,3> fij(particleLJ P, pvector<ntype,3> L, ntype &vij, ntype& vijs, ntype &wij)
    {
      // L is a vector with box dimensions
      // vij will be the interaction potential between i and j
      // vijs will be the shifted interaction potential
      // wij is the virial (i.e. rij*fij) for calculating the pressure
      ntype fij, srij2, srij6, srij12;
      ntype rijsq, epsilon24 = epsilon*24.0;
      ntype epsilon4 = epsilon*4.0;
      pvector<ntype,3> Dr, fijv;
      Dr = r - P.r;
      
      // minimum image convention
      Dr = Dr - mulcw(L,rint(divcw(Dr,divcw(L)))); // Dr - L*rint(Dr/L)
      ntype rn = Dr.norm();
      if (rn >= rc)
        {
          wij=vij=vijs=0;
          fijv = {0,0,0};
          return fijv;
        }
      rijsq = rn*rn; 
      srij2  = sigma*sigma / rijsq;
      srij6  = srij2 * srij2 * srij2;
      srij12 = srij6 * srij6;
      vij = srij12 - srij6;

      // virial (for calculating pressure)
      wij = vij + srij12;
      wij *= epsilon24;

      // modulus of the force divided by modulus of rij
      // fij = 24*epsilon*(2*(sigma/r)^12 - (sigma/r)^6)
      fij =  wij / rijsq;

      /* force between two atoms */
      fijv  = fij * Dr;          
      vij = epsilon4 * vij;
      vijs = vij - vcut;
      return fijv;
    }
  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void)
    {
      ntype rijsq, srij2, srij6, srij12;
      ntype epsilon4=epsilon*4.0;
      rijsq = rc*rc;
      srij2 = sigma*sigma / rijsq;
      srij6  = srij2 * srij2 * srij2;
      srij12 = srij6*srij6;
      vcut = epsilon4*(srij12 - srij6);
    }
};

template <typename ntype>
class particleSS: public particle<ntype>
{
  using pt=particle<ntype>;
protected:
  using pt::vcut;
public:
  using pt::r, pt::epsilon, pt::rc, pt::sigma;
 
  ntype vij(particleSS P, pvector<ntype,3> L)
    {
      ntype ene;
      pvector<ntype,3> Dr;
      Dr = r - P.r;
      // MINIMUM IMAGE CONVENTION 
      Dr = Dr - mulcw(L,rint(divcw(Dr,L))); // Dr - L*rint(Dr/L)
      ntype rsq, rn = Dr.norm();
      rsq = rn*rn; 
      if (rsq < rc*rc) // interaction potential cut-off
        ene = 4.0*epsilon*(pow(sigma/rn,12.0));
      else
        ene = 0.0;
      return ene;
    }

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void)
    {
      // write your own code here
      vcut  =  4.0*epsilon*(pow(sigma/rc,12.0));
    }
};
#endif 
