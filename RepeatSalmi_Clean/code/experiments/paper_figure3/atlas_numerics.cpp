#include "atlas_numerics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace salmi {
// The numerical prescriptions of deriv/integ in the public ATLAS utilities.
// Integration is expressed with interval-local curvature to avoid subtracting
// large polynomial terms on a mass grid spanning twelve decades.
std::vector<double> atlas_derivative(const std::vector<double>& x,const std::vector<double>& f) {
 const int n=x.size();if(n<2 || f.size()!=x.size())throw std::invalid_argument("Invalid derivative grid");
 std::vector<double> result(n);
 result[0]=(f[1]-f[0])/(x[1]-x[0]);result[n-1]=(f[n-1]-f[n-2])/(x[n-1]-x[n-2]);
 for(int j=1;j<n-1;++j) {
   double scale=std::max({std::abs(f[j-1]),std::abs(f[j]),std::abs(f[j+1])})/std::abs(x[j]);if(scale==0)scale=1;
   double a=(f[j]-f[j-1])/(x[j]-x[j-1])/scale,b=(f[j+1]-f[j])/(x[j+1]-x[j])/scale;
   a=a/(std::sqrt(1+a*a)+1);b=b/(std::sqrt(1+b*b)+1);
   result[j]=scale*(a+b)/(1-a*b);
 }
 return result;
}
std::vector<double> atlas_integral(const std::vector<double>& x,const std::vector<double>& f,double start) {
 const int n=x.size();if(n<2 || f.size()!=x.size())throw std::invalid_argument("Invalid integral grid");
 std::vector<double> curvature(n,0),result(n,start);
 for(int j=3;j<n-1;++j) {
   double left=(f[j]-f[j-1])/(x[j]-x[j-1]),right=(f[j+1]-f[j])/(x[j+1]-x[j]);
   curvature[j]=(right-left)/(x[j+1]-x[j-1]);
 }
 for(int i=0;i<n-1;++i) {
   double c=0;
   if(i>=3 && i<n-2 && curvature[i]!=0) {
     double wt=std::abs(curvature[i+1])/(std::abs(curvature[i])+std::abs(curvature[i+1]));
     c=curvature[i+1]+wt*(curvature[i]-curvature[i+1]);
   }
   double h=x[i+1]-x[i];result[i+1]=result[i]+h*(.5*(f[i]+f[i+1])-c*h*h/6);
 }
 return result;
}
}
