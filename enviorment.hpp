#ifndef _ENVIORMENT_
#define _ENVIORMENT_

#include "./pvector.hpp"
#include "./particle.hpp"
#include<vector>
#include<array>
template<typename ntype>
class XYmodel
{
public:
 int L;
 std::vector<particle<ntype>> parts; // particles
 std::vector<pvector<ntype,3>> J; // couplings for each particle to + direction
    XYmodel(int Lbox): L(Lbox)
    {
        int N = L*L*L;
        parts.resize(N);
        J.resize(N);
        for(int i=0; i<N; i++){
            parts[i]= particle<ntype>();
            J[i] = {rng.randn(), rng.randn(), rng.randn()};
        }
    }

    int d3_d1(int ix, int iy, int iz)
    {
        ix%=L;
        iy%=L;
        iz%=L;
        return ix + L*iy + L*L*iz;
    }

    std::array<int,3> d1_3d(int index)
    {
        index%= (L*L*L);
        int iz = index/(L*L);
        int iy = (index - iz*L*L)/L;
        int ix = index - iz*L*L - iy*L;
        return {ix, iy, iz};
    }

    pvector<ntype, 2> calc_H(int index){
        auto [ix, iy, iz] = d1_3d(index);
        pvector<ntype,2> si = parts[index].s;
        pvector<ntype,2> H = {0.0, 0.0};

        pvector<ntype,2> sxn = parts[d3_d1(ix-1, iy, iz)].s;
        ntype Jxn = J[d3_d1(ix-1, iy, iz)].get(0);

        pvector<ntype,2> sxp = parts[d3_d1(ix+1, iy, iz)].s;
        ntype Jxp = J[index].get(0);

        pvector<ntype,2> syn = parts[d3_d1(ix, iy-1, iz)].s;
        ntype Jyn = J[d3_d1(ix, iy-1, iz)].get(1);

        pvector<ntype,2> syp = parts[d3_d1(ix, iy+1, iz)].s;
        ntype Jyp = J[index].get(1);

        pvector<ntype,2> szn = parts[d3_d1(ix, iy, iz-1)].s;
        ntype Jzn = J[d3_d1(ix, iy, iz-1)].get(2);

        pvector<ntype,2> szp = parts[d3_d1(ix, iy, iz+1)].s;
        ntype Jzp = J[index].get(2);

        pvector<ntype,2> Hx = Jxn*sxn + Jxp*sxp;
        pvector<ntype,2> Hy = Jyn*syn + Jyp*syp;
        pvector<ntype,2> Hz = Jzn*szn + Jzp*szp;

        H = Hx + Hy + Hz;

        return H;
    }

    ntype energy_particle(int index)
    {
        pvector<ntype,2> H = calc_H(index);
        ntype e = - (parts[index].s)*H;
        return e;
    }

    ntype total_energy()
    {
        ntype E = 0.0;
        int N = L*L*L;
        for(int i=0; i<N; i++){
            E += energy_particle(i);
        }
        return E/2.0; // each bond counted twice
    }

};


#endif
