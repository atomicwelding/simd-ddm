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
    const auto ddmSize = this->ddm.ddm_size;

    // i'll need to modify curve_fit to support both float and double
    const std::vector<double> times(this->ddm.delays.getTime().begin(),
                                    this->ddm.delays.getTime().end());

    int ix_center = this->ddm.ddm_width/2 + 1;
    int iy_center = this->ddm.ddm_height/2 + 1;

    std::vector<double> frequencies;
    std::vector<double> kxsAlongDelays;
    for(int ikx = 0; ikx < this->ddm.ddm_width/2; ikx++) {
        kxsAlongDelays.clear();

        auto temp = iy_center*this->ddm.ddm_width + ix_center + ikx;
        for(int it = 0; it < times.size(); it++)
            kxsAlongDelays.push_back(ddmBuffer[it*ddmSize + temp]);

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

    int HalfROI = utils::closest_index(
			frequencies.begin(), frequencies.end(), this->options.frequencyThreshold);
    return 2*HalfROI;
}

template<typename T, typename Callable>
void Fit<T,Callable>::fit() {
    // shortcut ptrs
    auto ddmbuf = this->ddm.exposeDdmBuffer();

    float* As = parameters;
    float* Bs  = &parameters[(ROI*ROI)];
    float* raw_fs = &parameters[(ROI*ROI)* 2];

    auto indexes = this->delays.getIndex();
    const std::vector<double> times(this->ddm.delays.getTime().begin(),
                                    this->ddm.delays.getTime().end());

    int idxCenter = this->ddm.ddm_width/2 + 1;
    int ROIStart = idxCenter - ROI/2;
    // values
    std::vector<double> ksAlongDelays;
    double imin, imax, ithresh;

    for(int iky = 0; iky < ROI; iky++) {
        for(int ikx = 0; ikx < ROI; ikx++) {

            ksAlongDelays.clear();
            for (int it = 0; it < times.size(); it++) {
                int temp = it * this->ddm.ddm_size
                           + (iky + ROIStart) * this->ddm.ddm_width
                           + (ikx + ROIStart);
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

    printf("Smooth() : not implemented yet\n");

    std::cout << "END DEBUGGING SMOOTHING";
}

template<typename T, typename Callable>
void QuadraticSmoothingFit<T,Callable>::save() {
   //  TinyTIFFWriterFile* tif = TinyTIFFWriter_open(
			// (this->options.pathOutput+"_ddm_fit.tif").c_str(), 32,
			// TinyTIFFWriter_Float,1, this->ROI, this->ROI, TinyTIFFWriter_Greyscale);
	TIFF* tif = TIFFOpen((this->options.pathOutput+"_ddm_fit.tif").c_str(), "w");
    if(!tif)
        std::cout << "Can't write fit parameters into tif!" << std::endl;

    for(int param = 0; param < 3; param++)
        // TinyTIFFWriter_writeImage(tif, &this->parameters[param*this->ROI*this->ROI]);
		utils::libTIFFWriter_writeImage(
				tif, &this->parameters[param*this->ROI*this->ROI], this->ROI, this->ROI);
    // TinyTIFFWriter_close(tif);
	TIFFClose(tif);
}


// static std::vector<std::vector<int>> generateR(int di) {
//
// }

template class Fit<float,std::function<float(float,float,float,float)> >;
template class Fit<double,std::function<float(float,float,float,float)> >;

template class QuadraticSmoothingFit<float, std::function<float(float,float,float,float)> >;
template class QuadraticSmoothingFit<double,std::function<float(float,float,float,float)> >;
