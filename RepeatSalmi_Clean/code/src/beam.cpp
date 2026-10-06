#include "salmi/beam.hpp"
#include <cmath>
#include <stdexcept>
namespace salmi {
PowerLaw::PowerLaw(double lo,double hi,double d):minimum(lo),maximum(hi),delta(d) {
    if(!(lo>1 && hi>lo) || !std::isfinite(d)) throw std::invalid_argument("Invalid power law");
    normalization=1/power_integral(lo,hi,-d);
}
double PowerLaw::density(double g) const {
    if(g<minimum || g>maximum) return 0;return normalization*std::pow(g,-delta);
}
double PowerLaw::mean_gamma() const {return normalization*power_integral(minimum,maximum,1-delta);}
Quadrature PowerLaw::quadrature(int panels,int order) const {
    if(panels<1) throw std::invalid_argument("Invalid beam panels");
    Quadrature out;double l=std::log(minimum),step=std::log(maximum/minimum)/panels;
    for(int j=0;j<panels;++j) {
        auto q=gauss_legendre(order,l+j*step,l+(j+1)*step);
        for(int i=0;i<order;++i) {
            double g=std::exp(q.x[i]);out.x.push_back(g);out.w.push_back(q.w[i]*g*density(g));
        }
    }return out;
}
}
