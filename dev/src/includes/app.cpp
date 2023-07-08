#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

#include <tinytiffwriter.h>
#include <omp.h>
#include <fftw3.h>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"

App::App(Options& options) : options(&options) {}
App::~App()= default;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

void App::run() {
    if(utils::stoe(this->options->encoding) != Mono12Packed) {
        std::cout << "Encoding not supported yet" << std::endl;
        return;
    }

	Timer timer;

    std::cout << "* Loading images..." << std::flush;
	timer.start();
    auto* stack = new Stack<float>(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->do_normalize);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    /*
     * TODO:
     *      add inplace fourier transform
     *      refactor code (naming convention), move it to stack
     */

    std::cout << "* Creating FFTW plan with " << omp_get_max_threads() << " threads..." << std::flush;
    int r = fftwf_init_threads();
    if(r == 0) {
        // throw une vraie erreur
        std::cerr << std::endl << "ERROR CAN'T SPAWN THREAD" << std::endl;
        return;
    }

	timer.start();
    fftwf_plan_with_nthreads(omp_get_max_threads());
	int fft_size = stack->aoi_height*(stack->aoi_width/2+1);
    fftwf_complex* stack_fft  = fftwf_alloc_complex(fft_size * this->options->loadNframes);

    int rank = 2;
    int n_in[] = {stack->aoi_height, stack->aoi_width};
    int n_out[] = {stack->aoi_height, stack->aoi_width/2+1};
    int idist = n_in[0] * n_in[1];
    int odist = n_out[0] * n_out[1];
    fftwf_plan plan = fftwf_plan_many_dft_r2c(rank, n_in, this->options->loadNframes,
                                              stack->N_images_buffer, n_in,
                                              1, idist,
                                              stack_fft, n_out,
                                              1, odist,
                                              FFTW_ESTIMATE);

	std::cout << "  " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Performing DFT..." << std::flush;
	timer.start();
    fftwf_execute(plan);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Computing DDM differences..." << std::flush;
	timer.start();

	int tau_max = this->options->tauMax;
	int N_frames =  this->options->loadNframes;
	std::vector<float> ddm(tau_max * fft_size);

	float mean_weight = 1. / ( N_frames - tau_max );

    // parallelize on every tau
    #pragma omp parallel for
    for(int tau = 0; tau < tau_max; tau++) {
		// stack indices
		int t, pix;

		// shortcut ptrs
		fftwf_complex *i1, *i2;
		float* ddm_cur = &ddm[tau * fft_size];

		// Update of ddm averages
		for(pix = 0; pix < fft_size; pix++)
			ddm_cur[pix] = 0;
		for(t = tau_max; t < N_frames; t++) { // for all times
			i1 =  &stack_fft[t *fft_size];
			i2 =  &stack_fft[(t-tau-1) * fft_size];

			for(pix = 0; pix < fft_size; pix++) // for all pixels
				ddm_cur[pix] +=
					std::pow(i1[pix][REAL]-i2[pix][REAL], 2.) + 
					std::pow(i1[pix][IMAG]-i2[pix][IMAG], 2.);
        }
		for(pix = 0; pix < fft_size; pix++)
			ddm_cur[pix] *= mean_weight;
    }
	std::cout << "          " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Writing files ..." << std::flush;
    timer.start();
    TinyTIFFWriterFile* tif = TinyTIFFWriter_open(this->options->pathOutput.c_str(), 32, TinyTIFFWriter_Float,
                                                  1, n_out[1], n_out[0],
                                                  TinyTIFFWriter_Greyscale);
    if(!tif) {
        std::cout << "[ERROR] CANNOT WRITE " << std::endl;
        return;
    }

    for(int frame = 0; frame < this->options->tauMax; frame++) {
        float* data = &ddm[frame * fft_size];
        TinyTIFFWriter_writeImage(tif, data);

    }
    timer.stop();
    std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;


    std::cout << "* Cleaning ..." << std::endl;
    delete stack;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    fftw_free(stack_fft);
    TinyTIFFWriter_close(tif);
};
