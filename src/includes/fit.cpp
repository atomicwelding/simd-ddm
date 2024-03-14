#include "fit.hpp"
#include "timer.hpp"
#include "fit_models.hpp"

#include <iostream>
#include <algorithm>

// TODO: make const what needs to be const

void Fit::process() {
    Timer timer;
    std::cout << "* Fit... ROI : " << ROI << "\n" << std::flush;

    std::cout << "    -> Fitting ..." << std::flush;
    timer.start();
    fit();
    std::cout << "             " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "    -> Smoothing ..." << std::flush;
    timer.start();
    smooth();
    std::cout << "           " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "    -> Saving ..." << std::flush;
    timer.start();
    save();
    std::cout << "              " << timer.elapsedSec() << "s" << std::endl;
};

int Fit::findROI() {
    const auto& ddm_buffer = ddm.get_ddm_buffer();
    const auto& lag_times = ddm.get_lag_times();
    const auto ddm_size = ddm.ddm_size;
    const auto n_lags = ddm.n_lags;

    int ix_center = ddm.ddm_width/2 + 1;
    int iy_center = ddm.ddm_height/2 + 1;
    int offset, lag_idx, dix;

    FitSolver<SingleExpFitModel> fit_solver(n_lags);
    for(lag_idx=0; lag_idx<n_lags; lag_idx++)
        fit_solver.fit_data.tau_vals[lag_idx] = lag_times[lag_idx];

    std::vector<double> frequencies(ddm.ddm_width/2);
    for(dix = 0; dix < ddm.ddm_width/2; dix++) {
        offset = iy_center*ddm.ddm_width + ix_center + dix;
        for(lag_idx=0; lag_idx<n_lags; lag_idx++)
            fit_solver.fit_data.ISF_vals[lag_idx] = ddm_buffer[lag_idx*ddm_size + offset];
        fit_solver.optimize();
        frequencies[dix] = fit_solver.get_fit_param(2);
    }

    int HalfROI = utils::closest_index(
			frequencies.begin(), frequencies.end(), options.max_decay_freq);
    return 2*HalfROI+1;
}


void Fit::fit() {
    const auto& ddm_buffer = ddm.get_ddm_buffer();
    const auto& lag_times = ddm.get_lag_times();    
    const auto n_lags = ddm.n_lags;
    const auto ddm_size = ddm.ddm_size;
    const auto ddm_width = ddm.ddm_width;
    const auto ddm_height = ddm.ddm_height;

    int ix_start = ddm_width/2 - ROI/2 + 1;
    int iy_start = ddm_height/2 - ROI/2 + 1;
    int offset, lag_idx, dix, diy;

    FitSolver<SingleExpFitModel> fit_solver(n_lags);
    for(lag_idx=0; lag_idx<n_lags; lag_idx++)
        fit_solver.fit_data.tau_vals[lag_idx] = lag_times[lag_idx];

    #pragma omp parallel for schedule(nonmonotonic:dynamic) private(offset) \
            firstprivate(fit_solver)
    for(diy = 0; diy < ROI; diy++) {
        for(dix = 0; dix < ROI; dix++) {
            offset = (iy_start+diy)*ddm_width + (ix_start+dix);
            for(lag_idx=0; lag_idx<n_lags; lag_idx++)
                fit_solver.fit_data.ISF_vals[lag_idx] = ddm_buffer[lag_idx*ddm_size + offset];
            fit_solver.optimize();

            parameters[diy*ROI + dix] = fit_solver.get_fit_param(0);
            parameters[ROI*ROI + diy*ROI + dix] = fit_solver.get_fit_param(1);
            parameters[ROI*ROI*2 + diy*ROI + dix] = fit_solver.get_fit_param(2);
        }
    }
}


void QuadraticSmoothingFit::smooth() {
    std::cout << std::endl << "DEBUG SMOOTHING" << std::endl;

    printf("Smooth() : not implemented yet\n");

    std::cout << "END DEBUGGING SMOOTHING";
}


void QuadraticSmoothingFit::save() {
   //  TinyTIFFWriterFile* tif = TinyTIFFWriter_open(
			// (options.pathOutput+"_ddm_fit.tif").c_str(), 32,
			// TinyTIFFWriter_Float,1, ROI, ROI, TinyTIFFWriter_Greyscale);
	TIFF* tif = TIFFOpen((options.output_path+"_ddm_fit.tif").c_str(), "w");
    if(!tif)
        std::cout << "Can't write fit parameters into tif!" << std::endl;

    for(int param = 0; param < 3; param++)
        // TinyTIFFWriter_writeImage(tif, &parameters[param*ROI*ROI]);
		utils::libTIFFWriter_writeImage(
				tif, &parameters[param*ROI*ROI], ROI, ROI);
    // TinyTIFFWriter_close(tif);
	TIFFClose(tif);
}