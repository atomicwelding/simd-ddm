#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

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
        std::cerr << "ERROR CAN'T SPAWN THREAD" << std::endl;
        return;
    }

	timer.start();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    fftwf_complex* stack_fft  = fftwf_alloc_complex(stack->image_size * this->options->loadNframes);
    fftwf_plan plan = fftwf_plan_dft_r2c_2d(stack->aoi_width,
                                            stack->aoi_height*this->options->loadNframes,
                                            stack->N_images_buffer,
                                            stack_fft,
                                            FFTW_ESTIMATE);
	std::cout << " " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Performing DFT..." << std::flush;
	timer.start();
    fftwf_execute(plan);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Computing DDM differences..." << std::flush;
	timer.start();

	std::vector<float> ddm(this->options->tau * stack->image_size);

	int tau_max = this->options->tau;
	float mean_weight = 1. / ( this->options->loadNframes - tau_max );

    // chaque tau est traité en //
    #pragma omp parallel for
    for(int tau = 0; tau < tau_max; tau++) {
		// stack indices
		int t, pix;

		// shortcut ptrs
		fftwf_complex *i1, *i2;
		float* ddm_cur = &ddm[tau * stack->image_size];

		// Update of ddm averages
		for(t = tau_max; t < this->options->loadNframes; t++) { // for all times
			i1 =  &stack_fft[t * stack->image_size];
			i2 =  &stack_fft[(t-tau-1) * stack->image_size];

			for(pix = 0; pix < stack->image_size; pix++)
				ddm_cur[pix] +=
					std::pow(i1[pix][REAL]-i2[pix][REAL], 2.) + 
					std::pow(i1[pix][IMAG]-i2[pix][IMAG], 2.);
        }
		for(pix = 0; pix < stack->image_size; pix++)
			ddm_cur[pix] *= mean_weight;
    }
	std::cout << "          " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Cleaning ..." << std::endl;
    delete stack;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    fftw_free(stack_fft);
};
