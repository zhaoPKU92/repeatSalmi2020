#include "salmi/anderson.hpp"
#include <algorithm>
#include <cmath>
namespace salmi {
std::vector<double> AndersonMixing::step(const std::vector<double>& x,const std::vector<double>& f,
    const std::vector<double>& scale,double beta,int positive_count) {
    std::vector<double> normalized_f(scale.size());
    double current_norm=0;for(std::size_t i=0;i<scale.size();++i){normalized_f[i]=f[i]/scale[i];current_norm+=normalized_f[i]*normalized_f[i];}
    if(last_norm_>0 && current_norm>16*last_norm_)history_.clear();
    last_norm_=current_norm;
    int n=int(x.size()),m=std::min(memory_,int(history_.size()));
    std::vector<double> answer(n);for(int i=0;i<n;++i)answer[i]=x[i]+beta*f[i];
    if(m) {
        int start=int(history_.size())-m;
        // Only normalized correction differences are needed for the Gram
        // matrix. Reconstruct full differences during the final update rather
        // than allocating two full history-sized temporary arrays.
        std::vector<std::vector<double>> normalized_df(m,std::vector<double>(scale.size()));
        auto normalize_row=[&](int j) {
            const auto& before=history_[start+j];
            const auto& after_f=j==m-1?f:history_[start+j+1].f;
            for(std::size_t i=0;i<scale.size();++i)normalized_df[j][i]=(after_f[i]-before.f[i])/scale[i];
        }
        ;
        if(scale.size()>20000){
            #pragma omp parallel for schedule(static)
            for(int j=0;j<m;++j)normalize_row(j);
        }else for(int j=0;j<m;++j)normalize_row(j);
        std::vector<std::vector<double>> matrix(m,std::vector<double>(m+1));double trace=0;
        auto gram_row=[&](int j) {
            for(int k=j;k<m;++k) {
                double sum=0;for(std::size_t i=0;i<scale.size();++i)sum+=normalized_df[j][i]*normalized_df[k][i];
                matrix[j][k]=sum;matrix[k][j]=sum;
            }
            for(std::size_t i=0;i<scale.size();++i)matrix[j][m]+=normalized_df[j][i]*normalized_f[i];
        }
        ;
        if(scale.size()>20000){
            #pragma omp parallel for schedule(static)
            for(int j=0;j<m;++j)gram_row(j);
        }else for(int j=0;j<m;++j)gram_row(j);
        for(int j=0;j<m;++j)trace+=matrix[j][j];
        bool valid=trace>1e-280;
        for(int j=0;j<m;++j)matrix[j][j]+=1e-12*trace;
        for(int j=0;j<m && valid;++j) {
            int pivot=j;for(int k=j+1;k<m;++k)if(std::abs(matrix[k][j])>std::abs(matrix[pivot][j]))pivot=k;
            std::swap(matrix[j],matrix[pivot]);double p=matrix[j][j];
            if(!(std::abs(p)>1e-300)){valid=false;break;}
            for(int k=j;k<=m;++k)matrix[j][k]/=p;
            for(int row=0;row<m;++row)if(row!=j) {
                double c=matrix[row][j];for(int k=j;k<=m;++k)matrix[row][k]-=c*matrix[j][k];
            }
        }
        if(valid) {
            auto candidate=answer;
            for(int j=0;j<m;++j) {
                double gamma=matrix[j][m];if(!std::isfinite(gamma)||std::abs(gamma)>1e5){valid=false;break;}
                const auto& before=history_[start+j];
                const auto& after_x=j==m-1?x:history_[start+j+1].x;
                const auto& after_f=j==m-1?f:history_[start+j+1].f;
                for(int i=0;i<n;++i)candidate[i]-=gamma*((after_x[i]-before.x[i])+beta*(after_f[i]-before.f[i]));
            }
            for(int i=0;i<positive_count && valid;++i) {
                if(!std::isfinite(candidate[i])){valid=false;break;}
                if(candidate[i]<0) {
                    // Project only the affected components. A global step
                    // bound lets one negligible Wien-tail value suppress all
                    // acceleration. Keep the paired LTE departure consistent.
                    if(n==2*positive_count)candidate[positive_count+i]-=candidate[i];
                    candidate[i]=0;
                }
            }
            if(valid)answer.swap(candidate);
        }
    }
    history_.push_back({x,f});if(int(history_.size())>memory_)history_.erase(history_.begin());
    return answer;
}
}
