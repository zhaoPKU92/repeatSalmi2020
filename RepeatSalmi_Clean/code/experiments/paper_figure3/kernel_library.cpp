#include "paper.hpp"
#include "salmi/constants.hpp"
#include "salmi/stopping.hpp"
#include "salmi/beam.hpp"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace salmi {
namespace {
std::size_t pair_template_budget() {
    static const std::size_t bytes=[](){
        const char* raw=std::getenv("SALMI_PAIR_CACHE_MB");
        if(!raw)return std::size_t(512)*1024*1024;
        std::size_t used=0;std::string value(raw);unsigned long long mb=std::stoull(value,&used);
        if(used!=value.size() || mb<1 || mb>8192)throw std::invalid_argument("SALMI_PAIR_CACHE_MB must be an integer1..8192");
        return std::size_t(mb)*1024*1024;
    }();
    return bytes;
}
}
double SparseScattering::symmetric_value(int i,std::size_t z) const {
    if(!shared_template || !symmetric_kernel.empty())return symmetric_kernel[z];
    const auto& p=*shared_template;int j=p.column[z],A=int(grid.mu.size()),E=int(grid.x.size());
    double down=(1-temperature_fraction)*p.low[z]+temperature_fraction*p.high[z];
    if(down==0)return 0;
    // down is itself one of the four candidates. In the common case it
    // already proves that the original cutoff is passed, without divisions.
    if(!(down>1e-28)) {
        double boltz=boltzmann_pair[(i/A)*E+j/A];
        double up=down*(grid.weight[j]/grid.weight[i])*boltz;
        double gain=down*(grid.weight[j]/grid.weight[i]);
        double largest=std::max({gain,up,down,up*grid.weight[i]/grid.weight[j]});
        if(!(largest>1e-28))return 0;
    }
    return float(down*(grid.x[j/A]/grid.x[i/A])/(grid.dx[i/A]*grid.dmu[i%A]));
}
double SparseScattering::loss(int i,std::size_t z) const {
    if(symmetric_kernel.empty() && !shared_template)return loss_coefficient[z];
    int A=int(grid.mu.size()),E=int(grid.x.size()),j=shared_template?shared_template->column[z]:column[z],e=i/A,f=j/A;
    double base=symmetric_value(i,z)*(grid.x[f]/grid.x[e])*grid.dx[f]*grid.dmu[j%A];
    return e<f?base*boltzmann_pair[e*E+f]:base;
}
CollisionTerms SparseScattering::terms(const double* n) const {
    int N=grid.size();CollisionTerms c;c.gain.assign(N,0);c.loss.assign(N,0);c.net.assign(N,0);
    if(paired) {
        int A=int(grid.mu.size()),E=int(grid.x.size());
        const auto& columns=shared_template?shared_template->column:column;
        std::vector<double> inverse_energy(E),phase_weight(N);
        std::vector<int> energy(N),reflected(N);
        for(int e=0;e<E;++e)inverse_energy[e]=1/grid.x[e];
        for(int i=0;i<N;++i) {
            int e=i/A,a=i%A;energy[i]=e;reflected[i]=e*A+A-1-a;
            phase_weight[i]=grid.x[e]*grid.dx[e]*grid.dmu[a];
        }
        for(int i=0;i<N;++i)for(std::size_t z=start[i];z<start[i+1];++z) {
            int j=columns[z],e=energy[i],f=energy[j];double av=symmetric_kernel.empty()?symmetric_value(i,z):double(symmetric_kernel[z]),boltz=boltzmann_pair[e*E+f];
            if(av==0)continue;
            double ij=av*phase_weight[j]*inverse_energy[e];
            double ji=av*phase_weight[i]*inverse_energy[f];
            auto add=[&](int u,int v) {
                c.gain[u]+=ij*n[v];c.loss[u]+=ij*boltz*(1+n[v]);
                if(u!=v){c.gain[v]+=ji*boltz*n[u];c.loss[v]+=ji*(1+n[u]);}
            };
            add(i,j);int im=reflected[i],jm=reflected[j];
            if(std::min(im,jm)!=i || std::max(im,jm)!=j)add(im,jm);
        }
        for(int i=0;i<N;++i)c.net[i]=(1+n[i])*c.gain[i]-n[i]*c.loss[i];
        return c;
    }
    for(int i=0;i<N;++i) {
        double g=0,l=0;
        if(symmetric_kernel.empty()) {
            for(std::size_t z=start[i];z<start[i+1];++z){double v=n[column[z]];g+=gain_coefficient[z]*v;l+=loss_coefficient[z]*(1+v);}
        }else {
            int A=int(grid.mu.size()),E=int(grid.x.size()),e=i/A;
            for(std::size_t z=start[i];z<start[i+1];++z) {
                int j=column[z],f=j/A;double v=n[j];
                double base=double(symmetric_kernel[z])*(grid.x[f]/grid.x[e])*grid.dx[f]*grid.dmu[j%A],boltz=boltzmann_pair[e*E+f];
                g+=base*(e>f?boltz:1)*v;l+=base*(e<f?boltz:1)*(1+v);
            }
        }
        c.gain[i]=g;c.loss[i]=l;c.net[i]=(1+n[i])*g-n[i]*l;
    }return c;
}
std::vector<double> SparseScattering::linearized_gain(const double* departure,const double* B,const double* mask) const {
    int N=grid.size(),A=int(grid.mu.size()),E=int(grid.x.size());std::vector<double> result(N,0);
    // Thin absolute-field layers have no LTE departure targets. In particular,
    // hot skins have dense redistribution matrices; do not traverse them only
    // to discard every contribution. This preserves every nonzero operation.
    if(mask && std::none_of(mask,mask+E,[](double value){return value!=0;}))return result;
    if(!paired) {
        for(int i=0;i<N;++i)if(!mask || mask[i/A]!=0)for(std::size_t z=start[i];z<start[i+1];++z) {
            int j=column[z];result[i]+=gain(i,z)*(1+B[i/A])/(1+B[j/A])*departure[j];
        }
        return result;
    }
    const auto& columns=shared_template?shared_template->column:column;
    for(int i=0;i<N;++i)for(std::size_t z=start[i];z<start[i+1];++z) {
        int j=columns[z],e=i/A,f=j/A;
        if(mask && mask[e]==0 && mask[f]==0)continue;
        double av=symmetric_kernel.empty()?symmetric_value(i,z):double(symmetric_kernel[z]),boltz=boltzmann_pair[e*E+f];
        if(av==0)continue;
        double ij=av*(grid.x[f]/grid.x[e])*grid.dx[f]*grid.dmu[j%A]*(1+B[e])/(1+B[f]);
        double ji=av*(grid.x[e]/grid.x[f])*grid.dx[e]*grid.dmu[i%A]*boltz*(1+B[f])/(1+B[e]);
        auto add=[&](int u,int v) {
            if(!mask || mask[e]!=0)result[u]+=ij*departure[v];
            if(u!=v && (!mask || mask[f]!=0))result[v]+=ji*departure[u];
        };
        add(i,j);int im=e*A+A-1-i%A,jm=f*A+A-1-j%A;
        if(std::min(im,jm)!=i || std::max(im,jm)!=j)add(im,jm);
    }
    return result;
}
// Derivative of the same temperature-interpolated Compton energy exchange,
// at fixed photon occupation. Float coefficient rounding is not differentiated;
// this is a smooth Jacobian preconditioner, not a change to the collision terms.
double compton_cooling_temperature_response(const SparseScattering& k,const double* n,double step) {
 if(!k.paired || !k.shared_template || !(step>0))throw std::invalid_argument("Thermal response requires shared paired kernel storage");
 const auto& g=k.grid;const auto& p=*k.shared_template;int A=g.mu.size(),E=g.x.size(),N=g.size();
 double theta=cgs::kb*k.temperature/cgs::mec2,total=0,compensation=0;
 auto accumulate=[&](double v){double y=v-compensation,t=total+y;compensation=(t-total)-y;total=t;};
 for(int i=0;i<N;++i)for(std::size_t z=k.start[i];z<k.start[i+1];++z) {
  int j=p.column[z],e=i/A,f=j/A;if(e==f)continue;
  double av=k.symmetric_value(i,z);if(av==0)continue;
  double dav=((p.high[z]-p.low[z])/step)*(g.x[f]/g.x[e])/(g.dx[e]*g.dmu[i%A]);
  double phase=g.x[f]*g.dx[f]*g.dmu[j%A];double rate=av*phase/g.x[e],drate=dav*phase/g.x[e];
  double b=k.boltzmann_pair[e*E+f],db=b*std::abs(g.x[f]-g.x[e])/theta;
  auto add=[&](int u,int v){double forward=n[v]*(1+n[u]),reverse=n[u]*(1+n[v]);double derivative=drate*(forward-b*reverse)-rate*db*reverse;accumulate(g.weight[u]*(g.x[e]-g.x[f])*derivative);};
  add(i,j);int im=e*A+A-1-i%A,jm=f*A+A-1-j%A;if(std::min(im,jm)!=i || std::max(im,jm)!=j)add(im,jm);
 }
 return 2*cgs::pi*occupation_intensity_x(1,1)*cgs::kappa_e*total;
}
KernelLibrary::KernelLibrary(const PhotonGrid& g,const std::string& p,double step):grid_(g),directory_(p),log_step_(step) {
    if(!(step>0 && step<=.5))throw std::invalid_argument("Kernel temperature spacing must be in (0,.5]");
    std::cout<<"Reusable pair-template cache budget MB="<<pair_template_budget()/(1024*1024)<<std::endl;
}
const Redistribution& KernelLibrary::table(int index) {
    auto found=tables_.find(index);if(found!=tables_.end())return *found->second;
    double T=1e6*std::exp(index*log_step_),theta=cgs::kb*T/cgs::mec2;
    std::ostringstream path;path<<directory_<<"/thermal_v3_method"<<redistribution_integration_version()<<"_e"<<grid_.x.size()<<"_a"<<grid_.mu.size()<<"_k"<<index<<".bin";
    int N=grid_.size();double signature[]={double(N),grid_.x.front(),grid_.x.back(),theta,log_step_,16,24,double(redistribution_integration_version())};
    auto r=std::make_shared<Redistribution>();r->grid=grid_;r->theta=theta;
    bool valid=false;
    {
        std::ifstream in(path.str(),std::ios::binary);double check[8]={};
        if(in.read(reinterpret_cast<char*>(check),sizeof(check))) {
            valid=true;for(int i=0;i<8;++i)if(check[i]!=signature[i])valid=false;
            if(valid) {
                r->probability.resize(std::size_t(N)*N);r->rate.resize(N);
                in.read(reinterpret_cast<char*>(r->probability.data()),r->probability.size()*sizeof(double));
                in.read(reinterpret_cast<char*>(r->rate.data()),r->rate.size()*sizeof(double));valid=bool(in);
            }
        }
    }
    if(!valid) {
        std::cout<<"Kernel T="<<T<<" K"<<std::endl;
        *r=Redistribution::assemble(grid_,theta,16,24,16);
        std::ofstream out(path.str(),std::ios::binary);
        if(!out)throw std::runtime_error("Create kernel cache directory before running");
        out.write(reinterpret_cast<const char*>(signature),sizeof(signature));
        out.write(reinterpret_cast<const char*>(r->probability.data()),r->probability.size()*sizeof(double));
        out.write(reinterpret_cast<const char*>(r->rate.data()),r->rate.size()*sizeof(double));
    }
    tables_[index]=r;return *r;
}
SparseScattering KernelLibrary::at(double T,int compact) {
    if(compact<0 || compact>4)throw std::invalid_argument("Kernel storage must be 0, 1, 2, 3 or 4");
    if(!(T>1e5 && T<8e9))throw std::invalid_argument("Atmosphere temperature outside kernel range");
    double z=std::log(T/1e6)/log_step_;int lower=int(std::floor(z));
    double f=z-lower,theta=cgs::kb*T/cgs::mec2;
    SparseScattering s;s.grid=grid_;s.temperature=T;int N=grid_.size(),A=int(grid_.mu.size());
    s.paired=compact>=2;
    if(s.paired)for(int a0=0;a0<A;++a0)if(grid_.mu[a0]!=-grid_.mu[A-1-a0] || grid_.dmu[a0]!=grid_.dmu[A-1-a0])throw std::invalid_argument("Paired storage requires a reflection-symmetric angular grid");
    if(compact) {
        int E=int(grid_.x.size());s.boltzmann_pair.resize(E*E);
        for(int e=0;e<E;++e)for(int f=0;f<E;++f)s.boltzmann_pair[e*E+f]=std::exp(-std::abs(grid_.x[e]-grid_.x[f])/theta);
    }
    auto append=[&](int i,int j,double down) {
            int low=std::min(i,j),high=std::max(i,j);
            if(down==0)return;
            double boltz=compact?s.boltzmann_pair[(low/A)*grid_.x.size()+high/A]:
                std::exp((grid_.x[low/A]-grid_.x[high/A])/theta);
            double up=down*(grid_.weight[high]/grid_.weight[low])*
                boltz;
            double pij=low==i?down:up,pji=low==i?up:down;
            double gain=pij*(grid_.weight[j]/grid_.weight[i]);
            double largest=s.paired?std::max({gain,pji,pij,pji*grid_.weight[i]/grid_.weight[j]}):std::max(gain,pji);
            if(largest>1e-28) {
                s.column.push_back(j);
                if(compact)s.symmetric_kernel.push_back(float(down*(grid_.x[high/A]/grid_.x[low/A])/(grid_.dx[low/A]*grid_.dmu[low%A])));
                else {s.gain_coefficient.push_back(gain);s.loss_coefficient.push_back(pji);}
            }
    };
    if(s.paired) {
        // The canonical nonzero pair layout depends on the two bracketing
        // tables, not on the atmosphere depth. Cache it once instead of
        // scanning an N*N dense array at every depth and outer iteration.
        auto found=pairs_.find(lower);
        if(found==pairs_.end()) {
            const auto& a=table(lower);const auto& b=table(lower+1);
            auto p=std::make_shared<ScatteringPairTemplate>();
            for(int i=0;i<N;++i) {
                p->start.push_back(p->column.size());int im=(i/A)*A+A-1-i%A;
                for(int j=i;j<N;++j) {
                    std::size_t idx=std::size_t(i)*N+j;
                    double lo=a.probability[idx],hi=b.probability[idx];if(lo==0 && hi==0)continue;
                    // Interpolated downward rates cannot exceed max(lo,hi),
                    // and detailed balance factors are <= 1. This bound is
                    // ten times below the existing actual-rate cutoff, so
                    // discarded tails cannot enter any matrix in this bin.
                    // The wide margin also covers floating-point roundoff.
                    double upper=std::max(lo,hi)*std::max(1.,grid_.weight[j]/grid_.weight[i]);
                    if(upper<1e-29)continue;
                    int jm=(j/A)*A+A-1-j%A;
                    if(std::make_pair(std::min(im,jm),std::max(im,jm))<std::make_pair(i,j))continue;
                    p->column.push_back(j);p->low.push_back(lo);p->high.push_back(hi);
                }
            }
            p->start.push_back(p->column.size());p->column.shrink_to_fit();p->low.shrink_to_fit();p->high.shrink_to_fit();
            pairs_[lower]=p;found=pairs_.find(lower);
        }
        auto p=found->second;p->used=++pair_clock_;
        if(compact>=3) {
            s.shared_template=p;s.temperature_fraction=f;s.start=p->start;
            if(compact==4) {
                // Retain the shared column layout, but evaluate the identical
                // rounded coefficient once per depth rather than per inner
                // transfer iteration. Zero entries remain explicit. Use a
                // temporary vector so symmetric_value still computes entries
                // until the entire cache is valid.
                std::vector<float> values(p->column.size());
                for(int i=0;i<N;++i)for(std::size_t k=p->start[i];k<p->start[i+1];++k)
                    values[k]=float(s.symmetric_value(i,k));
                s.symmetric_kernel=std::move(values);
            }
        }else {
            s.column.reserve(p->column.size());s.symmetric_kernel.reserve(p->column.size());
            for(int i=0;i<N;++i) {
                s.start.push_back(s.column.size());
                for(std::size_t k=p->start[i];k<p->start[i+1];++k)append(i,p->column[k],(1-f)*p->low[k]+f*p->high[k]);
            }
        }
        std::size_t bytes=0;for(const auto& entry:pairs_)bytes+=entry.second->bytes();
        while(pairs_.size()>1 && bytes>pair_template_budget()) {
            auto victim=pairs_.end();
            for(auto it=pairs_.begin();it!=pairs_.end();++it)if(it->first!=lower && (victim==pairs_.end() || it->second->used<victim->second->used))victim=it;
            bytes-=victim->second->bytes();pairs_.erase(victim);
        }
    }else {
        const auto& a=table(lower);const auto& b=table(lower+1);
        for(int i=0;i<N;++i) {
            s.start.push_back(s.column.size());
            for(int j=0;j<N;++j) {
                std::size_t idx=std::size_t(std::min(i,j))*N+std::max(i,j);
                append(i,j,(1-f)*a.probability[idx]+f*b.probability[idx]);
            }
        }
    }
    if(compact<3)s.start.push_back(s.column.size());
    s.column.shrink_to_fit();s.symmetric_kernel.shrink_to_fit();s.gain_coefficient.shrink_to_fit();s.loss_coefficient.shrink_to_fit();
    // Bound dense table residency during the final many-angle calculation.
    while(tables_.size()>2 && tables_.size()*std::size_t(N)*N*sizeof(double)>512u*1024u*1024u) {
        auto victim=tables_.begin();while(victim!=tables_.end() && (victim->first==lower || victim->first==lower+1))++victim;
        if(victim==tables_.end())break;tables_.erase(victim);
    }
    return s;
}
}
