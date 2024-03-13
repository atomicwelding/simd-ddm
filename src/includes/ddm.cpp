#include "ddm.hpp"
#include "utils.hpp"
#include "timer.hpp"

#include <gsl/gsl_sf_lambert.h>
#include <tinytiffwriter.h>
#include <omp.h>
#include <fftw3.h>

#include <ranges>
#include <iostream>


DDM::DDM(
		Stack& stack, utils::Options& options) :
    raw_ddm_width(stack.aoi_width/2+1),
    raw_ddm_height(stack.aoi_height),
	raw_ddm_size(raw_ddm_width*raw_ddm_height),
    ddm_width(2*(stack.aoi_width/2)+1),
    ddm_height(2*(stack.aoi_height/2)+1),
	ddm_size(ddm_width*ddm_height),
    Nt(options.loadNframes),
    n_lags(options.Ntau),
    max_lag_shift(std::floor(options.delayMax / stack.mean_sampling_time())),
    options(options),
	stack(stack) {

    // initialize buffers
#ifdef __AVX512F__
    stack_fft = (fftwf_complex*) _mm_malloc(2*Nt*raw_ddm_size*sizeof(float), 64);
    raw_ddm_buffer = (float*) _mm_malloc(n_lags*raw_ddm_size*sizeof(float), 64);
    ddm_buffer = (float*) _mm_malloc(n_lags*ddm_size*sizeof(float), 64);
#else
    stack_fft  = fftwf_alloc_complex(Nt*raw_ddm_size);
    raw_ddm_buffer = fftwf_alloc_real(n_lags*raw_ddm_size);
    ddm_buffer = fftwf_alloc_real(n_lags*ddm_size);
#endif

    compute_lags();
    compute_FFT();
}


DDM::~DDM() {

    fftwf_free(stack_fft);
    fftwf_free(raw_ddm_buffer);
    fftwf_free(ddm_buffer);
}


void DDM::compute_lags() {

    std::string mode = options.doLogScale? "logarithmic" : "linear";
    if(mode == "linear")
        lag_shifts = utils::linspace<int>(1, max_lag_shift, n_lags);
    else if(mode == "logarithmic") {
        double max_shift = max_lag_shift; // Largest lag shift (in double for floating point divisions)

        // The first Nlin lag shifts should be linearly distributed to avoid index repetitions.
        // We estimate Nlin from an analytical formula found with Mathematica and adjust it if necessary.
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
    
    if(n_lags > options.loadNframes )
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
    for(int it=0; it<Nt; it++) {
        fftwf_execute_dft_r2c(
            plan, &stack.images[it*stack.image_size], &stack_fft[it*raw_ddm_size]);
    }   

    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);

    std::cout << timer.elapsedSec() << "s" << std::endl;
}

void DDM::ddm_shift() {

    // See header doc comment for block defs
    Timer timer;

    std::cout << "* Mirroring DDM images...         " << std::flush;
    timer.start();

    int Nx_Half = raw_ddm_width-1;
    int Ny_Half = raw_ddm_height/2;

    auto raw_idx = [&](int i, int j, int k) -> int {
        return i*raw_ddm_size + j*raw_ddm_width + k;
    };
    auto idx = [&](int i, int jbis, int kbis) -> int {
        return i*ddm_size + jbis*ddm_width + kbis;
    };

    for(int i = 0; i < n_lags; i++) {

        // Copy 1 -> 3'
        for(int j=0; j<Ny_Half+1; j++) {
            for(int k=0; k<Nx_Half+1; k++) {
                int jbis = j + Ny_Half;
                int kbis = k + Nx_Half;
                ddm_buffer[idx(i,jbis,kbis)] = raw_ddm_buffer[raw_idx(i,j,k)];
            }
        }

        // Copy 4 -> 2'
        for(int j=raw_ddm_height-Ny_Half; j<raw_ddm_height; j++) {
            for(int k=0; k<Nx_Half+1; k++) {
                int jbis = j + Ny_Half - raw_ddm_height;
                int kbis = k + Nx_Half;
                ddm_buffer[idx(i,jbis,kbis)] = raw_ddm_buffer[raw_idx(i,j,k)];
            }
        }

        // Mirror 3'->1' and 2'->4'
        for(int jbis = 0; jbis<ddm_height; jbis++) {
            for(int kbis = Nx_Half+1; kbis<ddm_width; kbis++) {
                int jbisbis = ddm_height - 1 - jbis;
                int kbisbis = ddm_width - 1 - kbis;
                ddm_buffer[idx(i, jbisbis, kbisbis)] = ddm_buffer[idx(i, jbis, kbis)];
            }
        }
    }

    std::cout << "" << timer.elapsedSec() << "s" << std::endl;
}


void DDM::save() {
   
    Timer timer;
    std::cout << "* Writing DDM files...            " << std::flush;
    timer.start();

	TIFF* tif_file = TIFFOpen((options.pathOutput+"_ddm.tif").c_str(), "w");
    if(!tif_file)
        throw std::runtime_error("Can't write files");

    for(int frame = 0; frame < n_lags; frame++)
		utils::libTIFFWriter_writeImage(
				tif_file, &ddm_buffer[frame*ddm_size], ddm_width, ddm_height);
    TIFFClose(tif_file);

    std::ofstream lag_file;
    lag_file.open(options.pathOutput  + "_lags.dat");
    lag_file << "lag_shift" << "," << "lag_time" << "\n";
    for(int i = 0; i < lag_times.size(); i++)
        lag_file << lag_shifts[i] << "," << lag_times[i] << "\n";
    lag_file.close();

    timer.stop();
    std::cout << timer.elapsedSec() << "s" << std::endl;	
}