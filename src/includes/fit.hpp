//
// Created by weld on 05/12/23.
//

#ifndef SIMD_DDM_FIT_HPP
#define SIMD_DDM_FIT_HPP

#include "delays.hpp"
#include "ddm.hpp"
#include "utils.hpp"

#include <iostream>

template<typename T>
class Fit {
public:
    Fit(DDM<T> &ddm, utils::Options &options) :
        ddm(ddm),
        options(options) {

        this->ROI = findROI();
        this->parameters.resize(ROI*ROI*3);
    };

    void process();

protected:
    int ROI;

    const DDM<T>& ddm;
    utils::Options options;

    void fit();

    int findROI();
    std::vector<float> parameters;

private:
    // to be implemented
    virtual void smooth() = 0;
    virtual void save() = 0;
};


template<typename T>
class QuadraticSmoothingFit : Fit<T> {
// gpoy's routine
public:
    QuadraticSmoothingFit(DDM<T> &ddm, utils::Options &options)
    : Fit<T>(ddm, options) {};

    using Fit<T>::process;

private:
    void smooth();
    void save();
};


#endif //SIMD_DDM_FIT_HPP
