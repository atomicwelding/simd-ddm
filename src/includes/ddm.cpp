#include "ddm.hpp"
#include "utils.hpp"
#include "timer.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __AVX512F__
#include <immintrin.h>
#endif

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
    n_lags(options.Ntau),
    max_lag_shift(std::floor(options.delayMax / stack.mean_sampling_time())),
    options(options),
	stack(stack) {    
    
    compute_lags();
    compute_FFT();
    compute_DDM();
    ddm_shift();
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
    std::cout << "* Creating FFTW plan with " << omp_get_max_threads() << " threads..." << std::flush;

    int r = fftwf_init_threads();
    if(r == 0)
        throw std::runtime_error("Can't spawn threads");

    timer.start();
    fftwf_plan_with_nthreads(omp_get_max_threads());

    int rank = 2;

    int n_in[] = {stack.aoi_height, stack.aoi_width};
    int n_out[] = {raw_ddm_height, raw_ddm_width};

#ifdef __AVX512F__
    this->stack_fft = (fftwf_complex*) _mm_malloc(
        2*raw_ddm_size*options.loadNframes*sizeof(float), 64);
#else
    this->stack_fft  = fftwf_alloc_complex(raw_ddm_size * options.loadNframes);
#endif
    fftwf_plan plan = fftwf_plan_many_dft_r2c(
			rank, n_in, options.loadNframes,
			stack.images, n_in, 1, stack.image_size,
			stack_fft, n_out, 1, raw_ddm_size, FFTW_ESTIMATE);
    std::cout << "  " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Performing DFT..." << std::flush;
    timer.start();
    fftwf_execute(plan);
    std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
}


void DDM::compute_DDM() {

    Timer timer;
    timer.start();    

#ifdef __AVX512F__
    std::cout << "* Computing DDM differences [AVX512]... " << std::flush;
    this->raw_ddm_buffer = (float*) _mm_malloc(n_lags*raw_ddm_size*sizeof(float), 64);
    ddm_loop_avx512();
#elif __AVX2__
    std::cout << "* Computing DDM differences [AVX2]...   " << std::flush;
    this->raw_ddm_buffer = fftwf_alloc_real(n_lags*raw_ddm_size);
    ddm_loop_avx2();
#else
    std::cout << "* Computing DDM differences [autovec]..." << std::flush;
    this->raw_ddm_buffer = fftwf_alloc_real(n_lags*raw_ddm_size);
    ddm_loop_autovec();
#endif

    std::cout << timer.elapsedSec() << "s" << std::endl;
}


void DDM::ddm_loop_avx512() {
    
#ifdef __AVX512F__
    // init
    for(int i = 0; i < raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    int chunk_size = max_lag_shift/2;
    float* unrolled_stack_fft = static_cast<float*>(&stack_fft[0][0]);

    for(int offset=max_lag_shift; offset<options.loadNframes-chunk_size; offset+=chunk_size) {
        #pragma omp parallel for schedule(nonmonotonic:dynamic) firstprivate(unrolled_stack_fft)
        for(int lag_idx = 0; lag_idx < n_lags; lag_idx++) {
            int it, buf_idx;
            const float *i1, *i2;

            float* ddm_cur = &raw_ddm_buffer[lag_idx * raw_ddm_size];

            // AVX512 and AVX2 registers for calculations
            __m512 m512_a, m512_b, m512_c;
            __m256 m256_a, m256_b, m256_c;
            const auto perm = _mm512_set_epi32(15,13,11,9,7,5,3,1,14,12,10,8,6,4,2,0);

            int max_it = std::min(offset+chunk_size, options.loadNframes);
            int lag_shift = lag_shifts[lag_idx];
            for(it = offset; it < max_it; it++) {
                i1 = &unrolled_stack_fft[it*raw_ddm_size*2];
                i2 = &unrolled_stack_fft[(it-lag_shift)*raw_ddm_size*2];
                for(buf_idx=0; buf_idx<2*raw_ddm_size; buf_idx+=16) {
                    m512_a = _mm512_load_ps(&i1[buf_idx]);
                    m512_b = _mm512_load_ps(&i2[buf_idx]);
                    m512_c = _mm512_sub_ps(m512_a, m512_b);
                    m512_a = _mm512_permutexvar_ps(perm, m512_c);

                    m256_a = _mm256_load_ps(&ddm_cur[buf_idx/2]);
                    m256_b = _mm512_extractf32x8_ps(m512_a, 0);
                    m256_c = _mm256_fmadd_ps(m256_b, m256_b, m256_a);
                    m256_a = _mm512_extractf32x8_ps(m512_a, 1);
                    m256_b = _mm256_fmadd_ps(m256_a, m256_a, m256_c);
                    _mm256_store_ps(&ddm_cur[buf_idx/2], m256_b);
                }
            }
        }
    }

    float mean_weight = 1. / ( 2*raw_ddm_size * (options.loadNframes-n_lags) );
    for(int i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX512");
#endif
}


void DDM::ddm_loop_avx2() {
    
#ifdef __AVX2__
    // init
    for(int i = 0; i < raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int lag_idx = 0; lag_idx < n_lags; lag_idx++) {
        int it, buf_idx;
        const float *i1, *i2;

        float* ddm_cur = &raw_ddm_buffer[lag_idx * raw_ddm_size];

        // AVX and SSE registers for calculations
        __m256 avx_a, avx_b, avx_c;
        __m128 sse_a, sse_b, sse_c;
        const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

        for(it = max_lag_shift; it < options.loadNframes; it++) {
            i1 = static_cast<const float*>(&stack_fft[it*raw_ddm_size][0]);
            i2 = static_cast<const float*>(&stack_fft[(it-lag_shifts[lag_idx])*raw_ddm_size][0]);
            for(buf_idx=0; buf_idx<2*raw_ddm_size; buf_idx+=8) {
                avx_a = _mm256_load_ps(&i1[buf_idx]);
                avx_b = _mm256_load_ps(&i2[buf_idx]);
                avx_c = _mm256_sub_ps(avx_a, avx_b);
                avx_a = _mm256_permutevar8x32_ps(avx_c, perm);

                sse_a = _mm_load_ps(&ddm_cur[buf_idx/2]);
                sse_b = _mm256_extractf128_ps(avx_a, 0);
                sse_c = _mm_fmadd_ps(sse_b, sse_b, sse_a);
                sse_a = _mm256_extractf128_ps(avx_a, 1);
                sse_b = _mm_fmadd_ps(sse_a, sse_a, sse_c);
                _mm_store_ps(&ddm_cur[buf_idx/2], sse_b);
            }
        }
    }

    float mean_weight = 1. / ( 2*raw_ddm_size * (options.loadNframes-n_lags) );
    for(int i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}


void DDM::ddm_loop_autovec() {
    
    // init
    for(int i = 0; i < raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int lag_idx = 0; lag_idx < n_lags; lag_idx++) {
        int it, pix;

        const fftwf_complex *i1, *i2;;
        float* ddm_cur = &raw_ddm_buffer[lag_idx * raw_ddm_size];

        //update of raw_ddm averages
        for(it = max_lag_shift; it < options.loadNframes; it++) {
            i1 = &stack_fft[it*raw_ddm_size];
            i2 = &stack_fft[(it-lag_shifts[lag_idx])*raw_ddm_size];

            for(pix = 0; pix < raw_ddm_size; pix++) // for all pixels
                ddm_cur[pix] +=
                        std::pow(i1[pix][REAL] - i2[pix][REAL], 2.) +
                        std::pow(i1[pix][IMAG] - i2[pix][IMAG], 2.);

        }
    }

    float mean_weight =
		1. / ( 2*raw_ddm_size * (options.loadNframes-n_lags) );
    for(int i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;
}


void DDM::ddm_shift() {

    // See header doc comment for block defs
    Timer timer;

    std::cout << "* Mirroring DDM images ..." << std::flush;
    timer.start();

    this->ddm_buffer = fftwf_alloc_real(n_lags*ddm_size);

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

    std::cout << "              " << timer.elapsedSec() << "s" << std::endl;
}


void DDM::save() {
   
    Timer timer;
    std::cout << "* Writing DDM files ..." << std::flush;
    timer.start();

	TIFF* tif_file = TIFFOpen((this->options.pathOutput+"_ddm.tif").c_str(), "w");
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
    std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;	
}