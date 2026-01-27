#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
#include "./randnumgen.hpp"
template<typename ntype>
class particle
{
  pvector<ntype,2> sold; 
  randnumgen<ntype, std::mt19937_64> *rng = nullptr;
public:
  pvector<ntype,2> s;
 
  void tra_move_metro(ntype dtheta)
    {
      s(0)= s(0)*cos(dtheta) - s(1)*sin(dtheta);
      s(1)= s(0)*sin(dtheta) + s(1)*cos(dtheta);
      s /= s.norm();
    }

  void ov_rel_step(pvector<ntype,2> H)
    {
      s= (2*(s*H)/(H*H))*H - s;
      s /= s.norm();
    }

  void store()
    {
     sold = s;
    }

  void restore()
    {
      s = sold;
    }

  void show_old()
    {
      std::cout << "sold=(" << sold.get(0) << "," << sold.get(1) << ")\n";
    }

  particle(randnumgen<ntype, std::mt19937_64> *rng_instance) : rng(rng_instance)
    {
      ntype theta = 2.*rng->ranf()*M_PI;
      s = {cos(theta), sin(theta)};
      store();
    }

  particle () = default;
};
#endif 
