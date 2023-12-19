#include "fit.hpp"

#include <iostream>
#include <algorithm>
#include <tinytiffwriter.h>



#include "curve_fit.hpp"
#include "utils.hpp"
#include "timer.hpp"
// ddm already have delays, not necessary to take it in the constructor
// make const what needs to be const

template<typename T, typename Callable>
void Fit<T,Callable>::process() {
    Timer timer;
    std::cout << "* Fit... ROI : " << ROI << "\n" << std::flush;

    std::cout << "    -> Fitting ..." << std::flush;
    timer.start();
    this->fit();
    std::cout << "             " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "    -> Smoothing ..." << std::flush;
    timer.start();
    this->smooth();
    std::cout << "           " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "    -> Saving ..." << std::flush;
    timer.start();
    this->save();
    std::cout << "              " << timer.elapsedSec() << "s" << std::endl;
};

template<typename T, typename Callable>
int Fit<T,Callable>::findROI() {
    const auto& ddmBuffer = this->ddm.exposeDdmBuffer();
    const auto ddmSize = this->ddm.ddmSize();

    // i'll need to modify curve_fit to support both float and double
    const std::vector<double> times(this->ddm.delays.getTime().begin(),
                                    this->ddm.delays.getTime().end());

    int idxCenter = this->ddm.ddm_width/2 + 1;

    std::vector<double> frequencies;
    std::vector<double> kxsAlongDelays;
    for(int ikx = 0; ikx < this->ddm.ddm_width/2; ikx++) {
        kxsAlongDelays.clear();

        auto temp = idxCenter*this->ddm.ddm_width + idxCenter + ikx;
        for(int it = 0; it < times.size(); it++) {
            kxsAlongDelays.push_back(ddmBuffer[it*ddmSize + temp]);
        }

        auto A = kxsAlongDelays.back();


        // Estimation of the relaxation frequency
        auto imin = *std::min_element(kxsAlongDelays.begin(), kxsAlongDelays.end());
        auto imax = *std::max_element(kxsAlongDelays.begin(), kxsAlongDelays.end());
        auto ithresh = imin+(imax-imin)*(1-std::exp(-1));

        int t;
        for(t = 0; t < times.size(); t++) {
            if(kxsAlongDelays[t]>ithresh)
                break;
        }
        auto f = 1./times[t];

        auto res = curve_fit(this->fn, {A,0.0,f}, times, kxsAlongDelays);

        frequencies.push_back(res[2]);
    }

    float relaxationFrequency = this->delays.getSamplingTime() * this->options.frequencyThreshold;
    int HalfROI = utils::closest_index(frequencies.begin(), frequencies.end(), relaxationFrequency);
    return 2*HalfROI;
}

template<typename T, typename Callable>
void Fit<T,Callable>::fit() {
    // wishful thinking
    // roi = roi();
    // fit();
    // smooth();  ---> average vs azimuthal

    // shortcut ptrs
    auto ddmbuf = this->ddm.exposeDdmBuffer();

    float* As = parameters;
    float* Bs  = &parameters[(ROI*ROI)];
    float* raw_fs = &parameters[(ROI*ROI)* 2];

    auto indexes = this->delays.getIndex();
    const std::vector<double> times(this->ddm.delays.getTime().begin(),
                                    this->ddm.delays.getTime().end());

    int idxCenter = this->ddm.ddm_width/2 + 1;
    int roim = idxCenter - ROI/2;
    int roip = idxCenter + ROI/2;

    // values
    std::vector<double> ksAlongDelays;
    double imin, imax, ithresh;

    for(int iky = 0; iky < ROI; iky++) {
        for(int ikx = 0; ikx < ROI; ikx++) {

            ksAlongDelays.clear();
            for (int it = 0; it < times.size(); it++) {
                int temp = it * this->ddm.ddmSize()
                           + (iky + roim) * this->ddm.ddm_width
                           + (ikx + roim);
                ksAlongDelays.push_back(ddmbuf[temp]);
            }

            auto A = ksAlongDelays.back();
            auto B = 0.0;

            // freq estimation
            imin = *std::min_element(ksAlongDelays.begin(), ksAlongDelays.end());
            imax = *std::max_element(ksAlongDelays.begin(), ksAlongDelays.end());
            ithresh = imin + (imax - imin ) * (1 - std::exp(-1));

            int t;
            for(t = 0; t < times.size(); t++) {
                if(ksAlongDelays[t]>ithresh)
                    break;
            }

            auto f = 1./times[t];
            auto params_fitted = curve_fit(this->fn, {A,B,f}, times, ksAlongDelays);

            As[iky * ROI + ikx] = params_fitted[0];
            Bs[iky * ROI + ikx] = params_fitted[1];
            raw_fs[iky * ROI + ikx] = params_fitted[2];
        }
    }
};

template<typename T, typename Callable>
void QuadraticSmoothingFit<T,Callable>::smooth() {
    std::cout << std::endl << "DEBUG SMOOTHING" << std::endl;


    std::cout << "END DEBUGGING SMOOTHING" << std::endl;
}

template<typename T, typename Callable>
void QuadraticSmoothingFit<T,Callable>::save() {
    TinyTIFFWriterFile* fit_tiff = TinyTIFFWriter_open("fit.tif", 32,
                                                       TinyTIFFWriter_Float,1,
                                                       this->ROI, this->ROI,
                                                       TinyTIFFWriter_Greyscale);
    if(!fit_tiff) {
        std::cout << "Can't write fitting parameters into tif!" << std::endl;
    }

    for(int param = 0; param < 3; param++)
        TinyTIFFWriter_writeImage(fit_tiff,
                                  &this->parameters[param*this->ROI*this->ROI]);
    TinyTIFFWriter_close(fit_tiff);
}

template class Fit<float,std::function<float(float,float,float,float)> >;
template class Fit<double,std::function<float(float,float,float,float)> >;

template class QuadraticSmoothingFit<float, std::function<float(float,float,float,float)> >;
template class QuadraticSmoothingFit<double,std::function<float(float,float,float,float)> >;
