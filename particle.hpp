#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
template<typename ntype>
class particle
{
  pvector<ntype,2> sold; 
public:
  pvector<ntype,2> s;
 
  void tra_move_metro(ntype dtheta)
    {
      s(0)= s(0)*cos(dtheta) - s(1)*sin(dtheta);
      s(1)= s(0)*sin(dtheta) + s(1)*cos(dtheta);
    }

  void ov_rel_step(pvector<ntype,2> H)
    {
      s-= (2*(s*H)/(H*H))*H;
    }

  void store()
    {
     sold = s;
    }

  void restore()
    {
      s = sold;
    }

  particle()
    {
      theta = rng.ranf()*M_PI;
      s = {cos(theta), sin(theta)};
    }
};
#endif 
