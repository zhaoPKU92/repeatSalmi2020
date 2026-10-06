// Exact thermal redistribution kernel with energy-transfer fitted quadrature.
// See docs/energy_bin_integration.md for the approximation and validation scope.
#include "salmi/redistribution.hpp"
#include "salmi/compton.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
namespace salmi {
Redistribution Redistribution::assemble_bins(const PhotonGrid& g,double th,int naz,int ng,int np,int energy_order,int rate_order) {
    if(g.size()==0 || naz<4)throw std::invalid_argument("Invalid redistribution grid/quadrature");
    ThermalKernel kernel(th,ng,np);Redistribution r;r.grid=g;r.theta=th;
    int A=int(g.mu.size()),E=int(g.x.size()),N=g.size();
    r.probability.assign(std::size_t(N)*N,0);r.rate.assign(N,0);
    std::vector<double> cosphi(naz);
    for(int k=0;k<naz;++k)cosphi[k]=std::cos(2*cgs::pi*(k+.5)/naz);
    if(energy_order<4 || rate_order<32)throw std::invalid_argument("Insufficient bin quadrature order");
    auto qe=gauss_legendre(energy_order,0,1);
    // Energy pairs and angular pairs are symmetric. Detailed balance is imposed
    // pairwise using the exact continuum identity, without adjusting total rates.
    bool mirror=true;
    for(int a=0;a<A;++a)if(g.mu[a]!=-g.mu[A-1-a] || g.dmu[a]!=g.dmu[A-1-a])mirror=false;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int e=0;e<E;++e)for(int f=e;f<E;++f)for(int a=0;a<A;++a)for(int b=0;b<=a;++b) {
        if(mirror && a+b>A-1)continue;
        double x=g.x[e],y=g.x[f],av=0;
        double aa=g.mu[a]*g.mu[b],bb=std::sqrt((1-g.mu[a]*g.mu[a])*(1-g.mu[b]*g.mu[b]));
        for(double cp:cosphi) {
            double eta=aa+bb*cp;
            if(e!=f) {
                double lo=.5*(std::sqrt(g.dx[e]*g.dx[e]+4*x*x)-g.dx[e]);
                double hi=e+1==f?y:lo+g.dx[e],val=0;
                double rest=y/(1+y*(1-eta)),closest=std::max(lo,std::min(hi,rest));
                if((minimum_electron_gamma(closest,y,eta)-1)/th>700)continue;
                for(std::size_t q=0;q<qe.x.size();++q) {
                    double xx=lo+(hi-lo)*qe.x[q];
                    double kv=kernel(xx,y,eta);if(kv==0)continue;
                    double fit=(y-xx)/(y-x)*std::expm1((xx-y)/th)/std::expm1((x-y)/th);
                    val+=qe.w[q]*(hi-lo)*(xx/y)*kv*fit;
                }
                av+=val/g.dx[e]*(y/x);
            }
            else {
                // Resolve the elastic forward peak inside the diagonal bin.
                // Split at x=y; evaluating this cusp only at its maximum biases
                // the total cross section even when energy/angle symmetry holds.
                double lo=.5*(std::sqrt(g.dx[e]*g.dx[e]+4*x*x)-g.dx[e]);
                double hi=lo+g.dx[e],v=0;
                for(int side=0;side<2;++side) {
                    double l=side?x:lo,h=side?hi:x;
                    for(std::size_t j=0;j<qe.x.size();++j) {
                        double xx=l+(h-l)*qe.x[j];
                        v+=qe.w[j]*(h-l)*(xx/y)*kernel(xx,y,eta);
                    }
                }
                av+=v/g.dx[e];
            }
        }
        av*=2*cgs::pi/naz;
        double reverse=av*std::exp((x-y)/th);
        auto put=[&](int aa0,int bb0) {
            int i=g.state(e,aa0),j=g.state(f,bb0);
            r.probability[std::size_t(i)*N+j]=(x/y)*av*g.dx[e]*g.dmu[aa0];
            r.probability[std::size_t(j)*N+i]=(y/x)*reverse*g.dx[f]*g.dmu[bb0];
        };
        put(a,b);if(a!=b)put(b,a);
        if(mirror && a+b<A-1) {
            put(A-1-a,A-1-b);
            if(a!=b)put(A-1-b,A-1-a);
        }
    }

    // Constrain the zeroth moment at each incoming energy and angular pair.
    // Inelastic rates are limited only if they exceed the independently
    // integrated full angular cross section; the remainder is elastic.
    auto qr=gauss_legendre(rate_order,0,1);
    int limited=0;double worst=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:limited) reduction(max:worst)
#endif
    for(int a=0;a<A;++a)for(int b=0;b<=a;++b) {
        std::vector<double> total(E),factor(E,1),sum(E);
        double aa=g.mu[a]*g.mu[b],bb=std::sqrt((1-g.mu[a]*g.mu[a])*(1-g.mu[b]*g.mu[b]));
        double lower=.5*(std::sqrt(g.dx.front()*g.dx.front()+4*g.x.front()*g.x.front())-g.dx.front());
        double upper=.5*(std::sqrt(g.dx.back()*g.dx.back()+4*g.x.back()*g.x.back())+g.dx.back());
        for(int f=0;f<E;++f) {
            double y=g.x[f];
            for(double cp:cosphi) {
                double eta=aa+bb*cp;
                double width=std::max(.08,12*std::sqrt(2*th*(1-eta))+6*y*(1-eta));
                for(int side=0;side<2;++side) {
                    double lo=side?std::log(y):std::max(std::log(lower),std::log(y)-width);
                    double hi=side?std::min(std::log(upper),std::log(y)+width):std::log(y);
                    for(std::size_t q=0;q<qr.x.size();++q) {
                        double xx=std::exp(lo+(hi-lo)*qr.x[q]);
                        total[f]+=qr.w[q]*(hi-lo)*xx*(xx/y)*kernel(xx,y,eta)*(2*cgs::pi/naz);
                    }
                }
            }
            for(int e=0;e<E;++e)if(e!=f)
                sum[f]+=r.probability[std::size_t(g.state(e,a))*N+g.state(f,b)]/g.dmu[a];
            if(sum[f]>total[f]) {
                factor[f]=total[f]/sum[f];++limited;worst=std::max(worst,1-factor[f]);
            }
        }
        for(int e=0;e<E;++e)for(int f=0;f<E;++f)if(e!=f) {
            double z=std::min(factor[e],factor[f]);
            r.probability[std::size_t(g.state(e,a))*N+g.state(f,b)]*=z;
            if(a!=b)r.probability[std::size_t(g.state(e,b))*N+g.state(f,a)]*=z;
        }
        for(int f=0;f<E;++f) {
            double out=0;
            for(int e=0;e<E;++e)if(e!=f)
                out+=r.probability[std::size_t(g.state(e,a))*N+g.state(f,b)]/g.dmu[a];
            double elastic=std::max(0.,total[f]-out);
            r.probability[std::size_t(g.state(f,a))*N+g.state(f,b)]=elastic*g.dmu[a];
            if(a!=b)r.probability[std::size_t(g.state(f,b))*N+g.state(f,a)]=elastic*g.dmu[b];
        }
    }
    std::cerr<<"Moment quadrature theta="<<th<<" limited angular energy rows="<<limited<<" maximum factor reduction="<<worst<<std::endl;
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)r.rate[j]+=r.probability[std::size_t(i)*N+j];
    return r;
}
}
