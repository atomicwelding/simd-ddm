#include "fit.hpp"

#include <iostream>

#include "curve_fit.hpp"
#include "utils.hpp"

// ddm already have delays, not necessary to take it in the constructor
// make const what needs to be const

template<typename T, typename Callable>
T Fit2D<T,Callable>::findROI() {
    /**
     * Not working properly.
     */
    const auto& ddmBuffer = this->ddm.exposeDdmBuffer();
    const auto ddmSize = this->ddm.ddmSize();

    // i'll need to modify curve_fit to support both float and double
    const std::vector<double> times(this->ddm.delays.getTime().begin(),
                                    this->ddm.delays.getTime().end());

    int idxCenter = this->ddm.ddm_width/2 + 1;
    std::cout << std::endl << idxCenter << std::endl;

    std::vector<double> frequencies;
    std::vector<double> kxsAlongDelays;
    for(int ikx = 0; ikx < idxCenter; ikx++) {
        kxsAlongDelays.clear();

        auto temp = idxCenter*this->ddm.ddm_width + idxCenter + ikx;
        for(int it = 0; it < times.size(); it++) {
            kxsAlongDelays.push_back(ddmBuffer[it*ddmSize + temp]);
        }

        auto A = kxsAlongDelays.back();
        auto f = 1./times.back();

        auto res = curve_fit(this->fn, {A,0.0,f}, times, kxsAlongDelays);

        frequencies.push_back(res[2]);
    }

    float relaxationFrequency = this->delays.getSamplingTime() * this->options.frequencyThreshold;
    return utils::closest_index(frequencies.begin(), frequencies.end(), relaxationFrequency);
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