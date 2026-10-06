#include "salmi/redistribution.hpp"
#include "salmi/compton.hpp"
#include "salmi/constants.hpp"
#include "salmi/radiation.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <string>
namespace salmi {
int redistribution_integration_version() {
    const char* value=std::getenv("SALMI_COMPTON_BIN_INTEGRATION");
    if(!value || std::string(value)=="midpoint")return 0;
    if(std::string(value)=="moment_v1")return 1;
    throw std::invalid_argument("Unknown SALMI_COMPTON_BIN_INTEGRATION; use midpoint or moment_v1");
}
PhotonGrid PhotonGrid::logarithmic(double xmin,double xmax,int ne,int nmu_half) {
    if(!(xmin>0 && xmax>xmin) || ne<4 || nmu_half<1)
        throw std::invalid_argument("Invalid photon grid");
    PhotonGrid g;auto edges=log_grid(xmin,xmax,ne+1);
    auto q=gauss_legendre(nmu_half,0,1);
    for(int i=0;i<nmu_half;++i){g.mu.push_back(-q.x[nmu_half-1-i]);g.dmu.push_back(q.w[nmu_half-1-i]);}
    for(int i=0;i<nmu_half;++i){g.mu.push_back(q.x[i]);g.dmu.push_back(q.w[i]);}
    for(int e=0;e<ne;++e) {
        g.x.push_back(std::sqrt(edges[e]*edges[e+1]));g.dx.push_back(edges[e+1]-edges[e]);
        for(double w:g.dmu)g.weight.push_back(g.x.back()*g.x.back()*g.dx.back()*w);
    }return g;
}
Redistribution Redistribution::assemble(const PhotonGrid& g,double th,int naz,int ng,int np) {
    if(redistribution_integration_version()==1)return assemble_bins(g,th,naz,ng,np);
    if(g.size()==0 || naz<4)throw std::invalid_argument("Invalid redistribution grid/quadrature");
    ThermalKernel kernel(th,ng,np);Redistribution r;r.grid=g;r.theta=th;
    int A=int(g.mu.size()),E=int(g.x.size()),N=g.size();
    r.probability.assign(std::size_t(N)*N,0);r.rate.assign(N,0);
    std::vector<double> cosphi(naz);
    for(int k=0;k<naz;++k)cosphi[k]=std::cos(2*cgs::pi*(k+.5)/naz);
    auto qe=gauss_legendre(8,0,1);
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
            if(e!=f)av+=kernel(x,y,eta);
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
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)r.rate[j]+=r.probability[std::size_t(i)*N+j];
    return r;
}
std::vector<double> Redistribution::bose_einstein(double chemical) const {
    if(!(chemical>=0))throw std::invalid_argument("Negative dimensionless chemical potential");
    std::vector<double> n(grid.size());
    for(std::size_t e=0;e<grid.x.size();++e)for(std::size_t a=0;a<grid.mu.size();++a) {
        double u=grid.x[e]/theta+chemical;
        n[grid.state(int(e),int(a))]=u<700?1/std::expm1(u):0;
    }return n;
}
CollisionTerms Redistribution::collision(const std::vector<double>& n,bool induced) const {
    int N=grid.size();if(int(n.size())!=N)throw std::invalid_argument("Occupation size mismatch");
    for(double v:n)if(!(v>=0) || !std::isfinite(v))throw std::invalid_argument("Invalid occupation");
    CollisionTerms r;r.gain.assign(N,0);r.loss.assign(N,0);r.net.assign(N,0);
    std::vector<double> photons(N);for(int j=0;j<N;++j)photons[j]=grid.weight[j]*n[j];
    if(!induced)r.loss=rate;
    for(int i=0;i<N;++i) {
        double gain=0,enhancement=1+n[i];
        for(int j=0;j<N;++j) {
            double p=probability[std::size_t(i)*N+j];
            gain+=p*photons[j];
            if(induced)r.loss[j]+=p*enhancement;
        }
        r.gain[i]=gain/grid.weight[i];
    }
    for(int i=0;i<N;++i)r.net[i]=(induced?1+n[i]:1)*r.gain[i]-n[i]*r.loss[i];
    return r;
}
ComptonSlab solve_compton_slab(const Redistribution& k,const std::vector<double>& tau,
    const std::vector<double>& bottom,const std::vector<double>& top,bool induced,double tol,int maxit) {
    int N=k.grid.size(),D=int(tau.size()),A=int(k.grid.mu.size());
    if(D<2 || int(bottom.size())!=N || int(top.size())!=N || !(tol>0) || maxit<1)
        throw std::invalid_argument("Invalid Compton slab input");
    for(int d=0;d<D;++d)if(!std::isfinite(tau[d]) || (d && !(tau[d]>tau[d-1])))
        throw std::invalid_argument("Nonincreasing slab depth");
    for(double v:bottom)if(!(v>=0) || !std::isfinite(v))throw std::invalid_argument("Invalid bottom boundary");
    for(double v:top)if(!(v>=0) || !std::isfinite(v))throw std::invalid_argument("Invalid top boundary");
    ComptonSlab r;r.tau=tau;r.occupation.assign(std::size_t(D)*N,0);
    for(int d=0;d<D;++d)for(int i=0;i<N;++i)
        r.occupation[std::size_t(d)*N+i]=.5*(bottom[i]+top[i]);
    std::vector<double> extinction(std::size_t(D)*N),source(extinction.size()),next(extinction.size());
    for(int it=0;it<maxit;++it) {
        for(int d=0;d<D;++d) {
            auto begin=r.occupation.begin()+std::size_t(d)*N;
            std::vector<double> local(begin,begin+N);auto c=k.collision(local,induced);
            for(int i=0;i<N;++i) {
                auto p=std::size_t(d)*N+i;extinction[p]=c.loss[i];
                source[p]=c.loss[i]>0?(induced?1+local[i]:1)*c.gain[i]/c.loss[i]:0;
            }
        }
        for(int i=0;i<N;++i) {
            double mu=k.grid.mu[i%A];int start=mu>0?D-1:0,end=mu>0?-1:D,step=mu>0?-1:1;
            next[std::size_t(start)*N+i]=mu>0?bottom[i]:top[i];
            for(int d=start+step;d!=end;d+=step) {
                std::size_t p=std::size_t(d)*N+i,u=std::size_t(d-step)*N+i;
                double dt=std::abs(tau[d]-tau[d-step])*.5*(extinction[p]+extinction[u])/std::abs(mu);
                auto w=short_characteristic_weights(dt);
                next[p]=w.transmission*next[u]+w.upstream*source[u]+w.local*source[p];
            }
        }
        double err=0,maxval=*std::max_element(next.begin(),next.end());
        for(std::size_t p=0;p<next.size();++p)
            err=std::max(err,std::abs(next[p]-r.occupation[p])/std::max({next[p],maxval*1e-14,1e-300}));
        r.history.push_back(err);r.iterations=it+1;r.residual=err;
        r.occupation.swap(next);
        if(err<tol) {
            r.flux.assign(D,0);r.photon_flux.assign(D,0);
            for(int d=0;d<D;++d)for(int i=0;i<N;++i) {
                double val=k.grid.weight[i]*k.grid.mu[i%A]*r.occupation[std::size_t(d)*N+i];
                r.photon_flux[d]+=val;r.flux[d]+=val*k.grid.x[i/A];
            }return r;
        }
    }throw std::runtime_error("Compton slab failed to converge");
}
}
