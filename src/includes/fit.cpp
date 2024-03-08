#include "fit.hpp"

#include <iostream>
#include <algorithm>
#include <tinytiffwriter.h>

#include "utils.hpp"
#include "timer.hpp"
#include "fit_models.hpp"

// ddm already have delays, not necessary to take it in the constructor
// make const what needs to be const

template<typename T>
void Fit<T>::process() {
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

template<typename T>
int Fit<T>::findROI() {
    const auto& ddmBuffer = this->ddm.exposeDdmBuffer();
    const auto ddm_size = this->ddm.ddm_size;
    const auto &tau_vals = this->ddm.delays.getTime();

    FitSolver<SingleExpFitModel> fit_solver(tau_vals.size());
    fit_solver.fit_data.tau_vals.assign(tau_vals.begin(), tau_vals.end());

    int ix_center = this->ddm.ddm_width/2 + 1;
    int iy_center = this->ddm.ddm_height/2 + 1;
    int offset;

    std::vector<double> frequencies(this->ddm.ddm_width/2);
    for(int ikx = 0; ikx < this->ddm.ddm_width/2; ikx++) {
        offset = iy_center*this->ddm.ddm_width + ix_center + ikx;
        for(int it = 0; it < tau_vals.size(); it++)
            fit_solver.fit_data.ISF_vals[it] = ddmBuffer[it*ddm_size + offset];
        fit_solver.optimize();
        frequencies[ikx] = fit_solver.get_fit_param(2);
    }

    int HalfROI = utils::closest_index(
			frequencies.begin(), frequencies.end(), this->options.frequencyThreshold);
    return 2*HalfROI+1;
}

template<typename T>
void Fit<T>::fit() {
    const auto &ddmBuffer = this->ddm.exposeDdmBuffer();
    const auto &tau_vals = this->ddm.delays.getTime();
    const auto ddm_size = this->ddm.ddm_size;
    const auto ddm_width = this->ddm.ddm_width;
    const auto ddm_height = this->ddm.ddm_height;

    int ix_start = ddm_width/2 - ROI/2 + 1;
    int iy_start = ddm_height/2 - ROI/2 + 1;
    int offset;

    FitSolver<SingleExpFitModel> fit_solver(tau_vals.size());
    fit_solver.fit_data.tau_vals.assign(tau_vals.begin(), tau_vals.end());
    #pragma omp parallel for schedule(nonmonotonic:dynamic) private(offset) \
            firstprivate(fit_solver)
    for(int diy = 0; diy < ROI; diy++) {
        for(int dix = 0; dix < ROI; dix++) {
            offset = (iy_start+diy)*ddm_width + (ix_start+dix);
            for(int it = 0; it < tau_vals.size(); it++)
                fit_solver.fit_data.ISF_vals[it] = ddmBuffer[it*ddm_size + offset];
            fit_solver.optimize();

            parameters[diy*ROI + dix] = fit_solver.get_fit_param(0);
            parameters[ROI*ROI + diy*ROI + dix] = fit_solver.get_fit_param(1);
            parameters[ROI*ROI*2 + diy*ROI + dix] = fit_solver.get_fit_param(2);
        }
    }
}

template<typename T>
void QuadraticSmoothingFit<T>::smooth() {
    std::cout << std::endl << "DEBUG SMOOTHING" << std::endl;

    printf("Smooth() : not implemented yet\n");

    std::cout << "END DEBUGGING SMOOTHING";
}

template<typename T>
void QuadraticSmoothingFit<T>::save() {
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

template class Fit<float>;
template class Fit<double>;

template class QuadraticSmoothingFit<float>;
template class QuadraticSmoothingFit<double>;