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
    virtual void fit() = 0;

    Fit(DDM<T> &ddm, utils::Options &options, Callable fn)
    : ddm(ddm), options(options), fn(fn), delays(ddm.delays) {};

protected:
    const Delays<T>& delays;
    const DDM<T>& ddm;
    utils::Options options;

    Callable fn;
};


template<typename T, typename Callable>
class Fit2D : Fit<T, Callable> {
// gpoy's routine
public:
    void fit();

    Fit2D(DDM<T> &ddm, utils::Options &options, std::function<T (T,T,T,T)> fn)
    : Fit<T, Callable>(ddm, options, fn) {
        ROI = findROI();
    };

private:
    // region of interest
    T ROI;
    T findROI();
};


#endif //SIMD_DDM_FIT_HPP
