//
// Created by weld on 05/12/23.
//

#ifndef SIMD_DDM_FIT_HPP
#define SIMD_DDM_FIT_HPP

#include "delays.hpp"
#include "ddm.hpp"
#include "utils.hpp"

#include <functional>
#include <iostream>


template<typename T, typename Callable>
class Fit {
public:
    Fit(DDM<T> &ddm, utils::Options &options, Callable fn)
    :  delays(ddm.delays), ddm(ddm), options(options), fn(fn) {
        this->ROI = findROI();
        this->parameters = fftwf_alloc_real((ROI*ROI)*3);
    };

    void process();

protected:
    int ROI;

    const Delays<T>& delays;
    const DDM<T>& ddm;
    utils::Options options;

    Callable fn;

    void fit();

    int findROI();
    float *parameters;
private:
    // to be implemented
    virtual void smooth() = 0;
    virtual void save() = 0;
};


template<typename T, typename Callable>
class QuadraticSmoothingFit : Fit<T, Callable> {
// gpoy's routine
public:
    QuadraticSmoothingFit(DDM<T> &ddm, utils::Options &options, std::function<T (T, T, T, T)> fn)
    : Fit<T, Callable>(ddm, options, fn) {};

    using Fit<T,Callable>::process;
private:

    void smooth();
    void save();
};


#endif //SIMD_DDM_FIT_HPP
