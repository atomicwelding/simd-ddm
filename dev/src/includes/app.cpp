#include <vector>
#include <iostream>

#include <tinytiffwriter.h>
#include <omp.h>
#include <fftw3.h>
#include <cmath>
#include <functional>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm.hpp"
#include "fitting.hpp"

App::App(utils::Options& options) : options(&options) {}
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

    int rank = 2;
    int n_in[] = {stack->aoi_height, stack->aoi_width};
    int n_out[] = {stack->aoi_height, stack->aoi_width/2+1};
    int fft_size = n_out[0] * n_out[1];
    fftwf_complex* stack_fft  = fftwf_alloc_complex(fft_size * this->options->loadNframes);
    fftwf_plan plan = fftwf_plan_many_dft_r2c(rank, n_in, this->options->loadNframes,
                                              stack->images, n_in,
                                              1, stack->image_size,
                                              stack_fft, n_out,
                                              1, fft_size,
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
    float* ddm = fftwf_alloc_real(tau_max * fft_size);
#ifdef __AVX2__
	DDM::ddm_loop_avx(ddm, stack_fft, fft_size, *(this->options));
#else
    DDM::ddm_loop_autovec(ddm, stack_fft, fft_size, *(this->options*));
#endif

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

    if(this->options->do_fit) {
        timer.start();
        std::cout << "* Fitting ..." << std::flush;

        auto exp_to_fit  = [](double tau, double A, double B, double f) -> double {
                    return A*(1-std::exp(-tau*f))+B;
        };

        fit::fit_routine(exp_to_fit, stack, ddm, tau_max, fft_size);

        timer.stop();
        std::cout << "                           " << timer.elapsedSec() << "s" << std::endl;
    }

    std::cout << "* Cleaning ..." << std::endl;

    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    fftwf_free(stack_fft);
    delete stack;
    TinyTIFFWriter_close(tif);
};


