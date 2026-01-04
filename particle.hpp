#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
template<typename ntype>
class particle
{
  pvector<ntype,2> sold; 
public:
  pvector<ntype,2> s;
 
  void tra_move_metro(pvector<ntype,2> ds)
    {
      s+=ds;
      s/=norm(s); 
    }

  void tra_move_rel(pvector<ntype,2> H)
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
      s = {sin(theta), cos(theta)};
    }
};
#endif 
