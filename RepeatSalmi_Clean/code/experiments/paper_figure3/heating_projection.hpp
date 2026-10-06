#pragma once
#include "salmi/stopping.hpp"
#include "salmi/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace salmi {
// Positive nodal-hat projection of the continuous track loss measure and its
// terminal thermalization. Partition of unity preserves energy/momentum and
// moves the terminal contribution continuously between neighboring nodes.
// This is a numerical projection choice, not a claimed author recipe.
inline void project_track_hats(const Track& tr,const std::vector<double>& m,
 double norm,std::vector<double>& Q,std::vector<double>& ram){
 const int D=int(m.size());if(D<2 || Q.size()!=m.size() || ram.size()!=m.size())throw std::runtime_error("Invalid hat deposition grid");
 std::vector<double> width(D);
 for(int d=0;d<D;++d)width[d]=d==0?.5*(m[1]-m[0]):d==D-1?.5*(m[d]-m[d-1]):.5*(m[d+1]-m[d-1]);
 auto add=[&](int d,double energy,double momentum){Q[d]+=norm*energy/width[d];ram[d]+=norm/cgs::c*momentum/width[d];};
 const double gx[]={-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526};
 const double gw[]={.3478548451374538,.6521451548625461,.6521451548625461,.3478548451374538};
 for(std::size_t j=1;j<tr.tau.size();++j){
  double x0=tr.tau[j-1]/cgs::kappa_e,x1=tr.tau[j]/cgs::kappa_e;
  double left=std::max(x0,m.front()),right=std::min(x1,m.back());if(right<=left)continue;
  int lo=std::min(D-2,std::max(0,int(std::upper_bound(m.begin(),m.end(),left)-m.begin())-1));
  while(left<right && lo<D-1){
   double end=std::min(right,m[lo+1]);
   double ga=tr.gamma[j-1]+(tr.gamma[j]-tr.gamma[j-1])*(left-x0)/(x1-x0);
   double gb=tr.gamma[j-1]+(tr.gamma[j]-tr.gamma[j-1])*(end-x0)/(x1-x0);
   double de=std::max(0.,ga-gb),fa=(left-m[lo])/(m[lo+1]-m[lo]),fb=(end-m[lo])/(m[lo+1]-m[lo]);
   fa=std::max(0.,std::min(1.,fa));fb=std::max(0.,std::min(1.,fb));
   if(de>0){
    double high_energy=de*.5*(fa+fb);
    double dp=de*(ga+gb)/(std::sqrt((ga-1)*(ga+1))+std::sqrt((gb-1)*(gb+1)));
    double weighted=0,total=0;
    for(int k=0;k<4;++k){double u=.5*(1+gx[k]),g=ga-u*de,p=std::sqrt((g-1)*(g+1)),weight=gw[k]*g/p;total+=weight;weighted+=weight*(fa+(fb-fa)*u);}
    double high_momentum=dp*weighted/total;
    add(lo,de-high_energy,dp-high_momentum);add(lo+1,high_energy,high_momentum);
   }
   left=end;if(end==m[lo+1])++lo;
  }
 }
 double stop=tr.tau.back()/cgs::kappa_e;
 if(stop>=m.front() && stop<=m.back()){
  int lo=std::min(D-2,std::max(0,int(std::upper_bound(m.begin(),m.end(),stop)-m.begin())-1));
  double f=(stop-m[lo])/(m[lo+1]-m[lo]),e=tr.cutoff-1,p=std::sqrt(e*(e+2));
  add(lo,e*(1-f),p*(1-f));add(lo+1,e*f,p*f);
 }
}
}
