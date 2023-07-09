#include <complex>
#include <vector>
#include <cmath>
#include <iostream>

#include <immintrin.h>
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
    float* ddm = fftwf_alloc_real(tau_max * fft_size);
	// ddm_loop_autovec(ddm, stack_fft, fft_size);
	ddm_loop_avx(ddm, stack_fft, fft_size);

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

void App::ddm_loop_autovec(float* ddm, const fftwf_complex* stack_fft, const int fft_size) {

	int tau_max = this->options->tauMax;
	int N_frames =  this->options->loadNframes;

	for(int i=0; i<fft_size*tau_max; i++)
		ddm[i] = 0;

    // parallelize on every tau, and let the compiler vectorize the code as much as
	// possible
    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int tau = 0; tau < tau_max; tau++) {
		// stack indices
		int t, pix;

		// shortcut ptrs
		const fftwf_complex *i1, *i2;
		float* ddm_cur = &ddm[tau * fft_size];

		// Update of ddm averages
		for(t = tau_max; t < N_frames; t++) { // for all times
			i1 =  &stack_fft[t *fft_size];
			i2 =  &stack_fft[(t-tau-1) * fft_size];

			for(pix = 0; pix < fft_size; pix++) // for all pixels
				ddm_cur[pix] +=
					std::pow(i1[pix][REAL]-i2[pix][REAL], 2.) + 
					std::pow(i1[pix][IMAG]-i2[pix][IMAG], 2.);
        }
    }

	float mean_weight = 1. / ( fft_size * (N_frames-tau_max) );
	for(int i=0; i<fft_size*tau_max; i++)
		ddm[i] *= mean_weight;
}

void App::ddm_loop_avx(float* ddm, const fftwf_complex* stack_fft, const int fft_size) {

#ifdef __AVX2__
	int tau_max = this->options->tauMax;
	int N_frames =  this->options->loadNframes;

	for(int i=0; i<fft_size*tau_max; i++)
		ddm[i] = 0;

	#pragma omp parallel for schedule(nonmonotonic:dynamic)
	for(int tau=0; tau<tau_max; tau++) {
		int t, pix;
		const float *i1, *i2;
		float* ddm_cur = &ddm[tau*fft_size];

		// AVX and SSE registers for calculations
		__m256 avx_a, avx_b, avx_c;
		__m128 sse_a, sse_b, sse_c;
		const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

		for(t=tau_max; t<N_frames; t++) {
			i1 = static_cast<const float*>(&stack_fft[t*fft_size][0]);
			i2 = static_cast<const float*>(&stack_fft[(t-tau-1)*fft_size][0]);
			for(pix=0; pix<fft_size; pix+=4) {
				avx_a = _mm256_load_ps(&i1[2*pix]);
				avx_b = _mm256_load_ps(&i2[2*pix]);
				avx_c = _mm256_sub_ps(avx_a, avx_b);
				avx_a = _mm256_permutevar8x32_ps(avx_c, perm);

				sse_a = _mm_load_ps(&ddm_cur[pix]);
				sse_b = _mm256_extractf128_ps(avx_a, 0);
				sse_c = _mm_fmadd_ps(sse_b, sse_b, sse_a);
				sse_a = _mm256_extractf128_ps(avx_a, 1);
				sse_b = _mm_fmadd_ps(sse_a, sse_a, sse_c);
				_mm_store_ps(&ddm_cur[pix], sse_b);
			}
		}
	}

	float mean_weight = 1. / ( fft_size * (N_frames-tau_max) );
	for(int i=0; i<fft_size*tau_max; i++)
		ddm[i] *= mean_weight;

#else
	throw std::exception("AVX intrinsics not supported on your CPU platform");
#endif
}
