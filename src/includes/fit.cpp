#include "fit.h"

#include <iostream>

#include "curve_fit.hpp"
#include "utils.hpp"

// ddm already have delays, not necessary to take it in the constructor
// make const what needs to be const

template<typename T, typename Callable>
T Fit2D<T,Callable>::findROI() {

    auto ddm_buffer = this->ddm.exposeDdmBuffer();
    std::vector<double> times(this->ddm.delays.getTime().begin(),
                              this->ddm.delays.getTime().end());

    int Ntau = this->options.Ntau;
    float delay_max = this->options.delayMax;
    float frequency_threshold = this->options.frequencyThreshold;
    int ddm_size = this->ddm.ddmSize();

    double relaxation_frequency = this->delays.getSamplingTime()*frequency_threshold;

    int idx_center = this->ddm.ddm_width/2;

    std::vector<double> frequencies;
    std::vector<double> kx_along_tau;
    for(int ikx = 0; ikx < idx_center; ikx++) {
        kx_along_tau.clear();
        int tmp_idx = idx_center*this->ddm.ddm_width + idx_center + ikx;
        for(int t = 0; t < Ntau; t++) {
            kx_along_tau.push_back(ddm_buffer[t * ddm_size + tmp_idx]);
        }

        double A = ddm_buffer[(Ntau-1)*ddm_size + tmp_idx];
        double B = 0.0;
        double f = 1./this->delays.getTime()[Ntau-1];

        auto res = curve_fit(this->fn, {A, B, f}, times, kx_along_tau);
        frequencies.push_back(res[2]);
    }

    return utils::closest_index(frequencies.begin(), frequencies.end(), relaxation_frequency);
}

template<typename T, typename Callable>
void Fit2D<T,Callable>::fit() {
    std::cout << std::endl << "FIT() : " << std::endl;
    std::cout << this->ROI << std::endl;
};


template class Fit<float,std::function<float(float,float,float,float)> >;
template class Fit<double,std::function<float(float,float,float,float)> >;

template class Fit2D<float, std::function<float(float,float,float,float)> >;
template class Fit2D<double,std::function<float(float,float,float,float)> >;