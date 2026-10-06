#pragma once
#include <vector>
namespace salmi {
class AndersonMixing {
    struct Sample {std::vector<double> x,f;};
    std::vector<Sample> history_;
    int memory_;
    double last_norm_=0;
public:
    explicit AndersonMixing(int memory=4):memory_(memory){}
    std::vector<double> step(const std::vector<double>& x,const std::vector<double>& f,
        const std::vector<double>& scale,double damping,int positive_count);
};
}
