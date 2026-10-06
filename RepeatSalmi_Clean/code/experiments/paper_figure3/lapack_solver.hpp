#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
namespace salmi {
// Optional LP64 OpenBLAS/LAPACK runtime. The default remains the existing
// standalone Jacobi solver. DGELSD solves the same minimum-norm least-squares
// problem, with the SAME relative singular-value cutoff (1e-6).
inline std::vector<double> lapack_temperature_solve(
    const std::vector<std::vector<double>>& rows,const std::vector<double>& rhs) {
#ifdef _WIN32
    struct Runtime {
        using Solve=void(*)(int*,int*,int*,double*,int*,double*,int*,double*,double*,int*,double*,int*,int*,int*);
        HMODULE library;Solve solve;
        Runtime(){
            const char* path=std::getenv("SALMI_LAPACK_DLL");
            if(!path || !*path)throw std::runtime_error("Missing SALMI_LAPACK_DLL");
            library=LoadLibraryExA(path,nullptr,0x00000100|0x00001000);
            if(!library)throw std::runtime_error("Cannot load temperature LAPACK DLL");
            solve=reinterpret_cast<Solve>(GetProcAddress(library,"dgelsd_"));
            auto config=reinterpret_cast<const char*(*)()>(GetProcAddress(library,"openblas_get_config"));
            auto threads=reinterpret_cast<void(*)(int)>(GetProcAddress(library,"openblas_set_num_threads"));
            if(!solve || !config || !threads)throw std::runtime_error("Missing LP64 OpenBLAS entry points");
            if(std::string(config()).find("USE64BITINT")!=std::string::npos)
                throw std::runtime_error("ILP64 DLL is incompatible with this LP64 interface");
            threads(1);
        }
    };
    static Runtime runtime;
    int m=int(rows.size()),n=m?int(rows.front().size()):0,nrhs=1;
    if(!m || !n || int(rhs.size())!=m)throw std::runtime_error("Invalid temperature least-squares shape");
    int lda=m,ldb=std::max(m,n),rank=0,info=0,lwork=-1;
    std::vector<double> a(std::size_t(m)*n),b(ldb,0),s(std::min(m,n));
    for(int i=0;i<m;++i){
        if(int(rows[i].size())!=n)throw std::runtime_error("Ragged temperature matrix");
        b[i]=rhs[i];
        for(int j=0;j<n;++j)a[std::size_t(j)*m+i]=rows[i][j];
    }
    // More than the documented LIWORK bound for these atmosphere dimensions;
    // include its logarithmic growth so larger grids remain supported.
    int minmn=std::min(m,n),levels=std::max(0,int(std::log2(std::max(1,minmn)))+1);
    std::vector<int> iwork(std::size_t(minmn)*(3*levels+11)+100);
    double rcond=1e-6,query=0;
    runtime.solve(&m,&n,&nrhs,a.data(),&lda,b.data(),&ldb,s.data(),&rcond,&rank,&query,&lwork,iwork.data(),&info);
    if(info || !std::isfinite(query) || query<1 || query>std::numeric_limits<int>::max())
        throw std::runtime_error("Temperature DGELSD workspace query failed");
    lwork=int(std::ceil(query));std::vector<double> work(lwork);
    runtime.solve(&m,&n,&nrhs,a.data(),&lda,b.data(),&ldb,s.data(),&rcond,&rank,work.data(),&lwork,iwork.data(),&info);
    if(info)throw std::runtime_error("Temperature DGELSD did not converge");
    b.resize(n);
    for(double x:b)if(!std::isfinite(x))throw std::runtime_error("Nonfinite temperature DGELSD solution");
    return b;
#else
    throw std::runtime_error("The optional LAPACK DLL interface requires Windows");
#endif
}
}
