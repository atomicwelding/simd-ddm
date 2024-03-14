#include "ddm.hpp"
#include "timer.hpp"

#include <gsl/gsl_sf_lambert.h>
#include <omp.h>

#include <ranges>
#include <iostream>

DDM::DDM(
		Stack& stack, utils::Options& options) :
    raw_ddm_width(stack.aoi_width/2+1),
    raw_ddm_height(stack.aoi_height),
    ddm_width(2*(stack.aoi_width/2)+1),
    ddm_height(2*(stack.aoi_height/2)+1),
    Nt(options.N_frames),
    n_lags(options.N_lags),
    max_lag_shift(std::floor(options.max_lag_time / stack.mean_sampling_time())),
	raw_ddm_size(raw_ddm_width*raw_ddm_height),
	ddm_size(ddm_width*ddm_height),
    options(options),
	stack(stack) {

    // initialize buffers
    stack_fft = utils::allocate_complex_float_array(raw_ddm_size*Nt);
    raw_ddm_buffer = utils::allocate_float_array(raw_ddm_size*n_lags);
    ddm_buffer = utils::allocate_float_array(ddm_size*n_lags);

    compute_lags();
    compute_FFT();
}


DDM::~DDM() {

    fftwf_free(stack_fft);
    fftwf_free(raw_ddm_buffer);
    fftwf_free(ddm_buffer);
}


void DDM::compute_lags() {

    std::string mode = options.log_lags? "logarithmic" : "linear";
    if(mode == "linear")
        lag_shifts = utils::linspace<int>(1, max_lag_shift, n_lags);
    else if(mode == "logarithmic") {
        double max_shift = max_lag_shift; // Largest lag shift (in double for floating point divisions)

        // The first Nlin lag shifts should be linearly distributed to avoid index repetitions.
        // We estimate Nlin from an analytical formula found with Mathematica and adjust it if
        // necessary.
        int n_lin = std::round(-n_lags/gsl_sf_lambert_Wm1(-n_lags/(max_shift*std::exp(1)))-0.5);
        while(std::round(n_lin*std::pow(max_shift/n_lin, 1./(n_lags-n_lin)))==n_lin)
            n_lin++;

        lag_shifts = utils::linspace<int>(1, n_lin-1, n_lin-1);
        auto log_lag_shifts = utils::logspace<int>(n_lin, max_lag_shift, n_lags-n_lin+1);
        lag_shifts.insert(lag_shifts.end(), log_lag_shifts.begin(), log_lag_shifts.end());
    }

    lag_times.resize(n_lags);
    for(int lag_idx=0; lag_idx<n_lags; lag_idx++)
        lag_times[lag_idx] = lag_shifts[lag_idx]*stack.mean_sampling_time();
    
    if(n_lags > options.N_frames )
        throw std::runtime_error(
				"Error : the file you want to process is shorter than the desired max lag index");
}

void DDM::compute_FFT() {
    
    Timer timer;
    timer.start();

    std::cout << "* Calculating spatial DFT...      " << std::flush;

    fftwf_plan plan = fftwf_plan_dft_r2c_2d(
			stack.aoi_height, stack.aoi_width, stack.images, stack_fft, FFTW_ESTIMATE);

    #pragma omp parallel for
    for(int it=0; it<Nt; it++)
        fftwf_execute_dft_r2c(
            plan, &stack.images[stack.image_size*it], &stack_fft[raw_ddm_size*it]);

    fftwf_destroy_plan(plan);

    std::cout << timer.elapsedSec() << "s" << std::endl;
}

void DDM::ddm_shift() {

    Timer timer;

    std::cout << "* Mirroring DDM images...         " << std::flush;
    timer.start();

    int Nx_half = raw_ddm_width-1;
    int Ny_half = raw_ddm_height/2;

    auto raw_idx = [&](int lag_idx, int iy, int ix) -> size_t {
        return raw_ddm_size*lag_idx + raw_ddm_width*iy + ix;
    };
    auto idx = [&](int lag_idx, int iy, int ix) -> size_t {
        return ddm_size*lag_idx + ddm_width*iy + ix;
    };

    // See header doc comment for block defs
    int lag_idx, ix, iy;
    for(lag_idx=0; lag_idx<n_lags; lag_idx++) {

        // Copy 1 -> 3'
        for(iy=0; iy<Ny_half+1; iy++)
            for(ix=0; ix<Nx_half+1; ix++)
                ddm_buffer[idx(lag_idx, iy+Ny_half, ix+Nx_half)] =
                    raw_ddm_buffer[raw_idx(lag_idx, iy, ix)];

        // Copy 4 -> 2'
        for(iy=raw_ddm_height-Ny_half; iy<raw_ddm_height; iy++)
            for(ix=0; ix<Nx_half+1; ix++)
                ddm_buffer[idx(lag_idx, iy+Ny_half-raw_ddm_height, ix+Nx_half)] =
                    raw_ddm_buffer[raw_idx(lag_idx, iy, ix)];

        // Mirror 3'->1' and 2'->4'
        for(iy=0; iy<ddm_height; iy++)
            for(ix=Nx_half+1; ix<ddm_width; ix++)
                ddm_buffer[idx(lag_idx, ddm_height-1-iy, ddm_width-1-ix)] =
                    ddm_buffer[idx(lag_idx, iy, ix)];
    }

    std::cout << "" << timer.elapsedSec() << "s" << std::endl;
}


void DDM::save() {
   
    Timer timer;
    std::cout << "* Writing DDM files...            " << std::flush;
    timer.start();

	TIFF* tif_file = TIFFOpen((options.output_path+"_ddm.tif").c_str(), "w");
    if(!tif_file)
        throw std::runtime_error("Can't write files");

    for(int frame = 0; frame < n_lags; frame++)
		utils::libTIFFWriter_writeImage(
				tif_file, &ddm_buffer[frame*ddm_size], ddm_width, ddm_height);
    TIFFClose(tif_file);

    std::ofstream lag_file;
    lag_file.open(options.output_path  + "_lags.dat");
    lag_file << "lag_shift" << "," << "lag_time" << "\n";
    for(int i = 0; i < lag_times.size(); i++)
        lag_file << lag_shifts[i] << "," << lag_times[i] << "\n";
    lag_file.close();

    timer.stop();
    std::cout << timer.elapsedSec() << "s" << std::endl;	
}