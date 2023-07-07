#include <omp.h>
#include <vector>
#include <tinytiffwriter.h>


#include <fftw3.h>
#include <iostream>

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
    fftwf_complex* out  = fftwf_alloc_complex(stack->image_size * this->options->loadNframes);

    int rank = 2;
    int n[] = {stack->aoi_height, stack->aoi_width};
    int idist = n[0] * n[1];
    int odist = idist;
    fftwf_plan plan = fftwf_plan_many_dft_r2c(rank, n, this->options->loadNframes,
                                              stack->N_images_buffer, n,
                                              1, idist,
                                              out, n,
                                              1, odist,
                                              FFTW_ESTIMATE);

	std::cout << "  " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Performing DFT..." << std::flush;
	timer.start();
    fftwf_execute(plan);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Computing DDM differences..." << std::flush;
	timer.start();

	std::vector<float> ddm(this->options->tauMax * stack->image_size);

	int tau_max = this->options->tauMax;

    // parallelize on every tau
    #pragma omp parallel for
    for(int tau = 0; tau < tau_max; tau++) {
        // flat indices
        int idx_k, idx_kp, idx_out;
		// stack indices
		int t, pix;
		// tmp vars for ddm calculation 
		float sqr_diff;

		// Update of ddm averages with online ergodic mean
		for(t = tau_max; t < this->options->loadNframes; t++) { // for all times
   	     	for(pix = 0; pix < stack->image_size; pix++) { // for all pixels

				idx_out = pix + tau * stack->image_size;
				idx_k =  pix + t * stack->image_size;
				idx_kp = pix + (t-tau-1) * stack->image_size;

				sqr_diff = 
					utils::sqr(out[idx_k][REAL]-out[idx_kp][REAL]) + 
					utils::sqr(out[idx_k][IMAG]-out[idx_kp][IMAG]);
                ddm[idx_out] = ( (t-tau_max) * ddm[idx_out] + sqr_diff ) / (t-tau_max+1.);
            }
        }
    }
	std::cout << "          " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Writing files ..." << std::flush;
    timer.start();
    TinyTIFFWriterFile* tif = TinyTIFFWriter_open(this->options->pathOutput.c_str(), 32, TinyTIFFWriter_Float,
                                                  1, stack->aoi_width, stack->aoi_height,
                                                  TinyTIFFWriter_Greyscale);
    if(!tif) {
        std::cout << "[ERROR] CANNOT WRITE " << std::endl;
        return;
    }

    for(int frame = 0; frame < this->options->tauMax; frame++) {
        float* data = &ddm[frame * stack->image_size];
        TinyTIFFWriter_writeImage(tif, data);

    }
    timer.stop();
    std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;


    std::cout << "* Cleaning ..." << std::endl;
    delete stack;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    fftw_free(out);
    TinyTIFFWriter_close(tif);
};
