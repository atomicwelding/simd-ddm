#include "ddm.hpp"

#include <cmath>
#include <immintrin.h>
#include <vector>

#include <iostream>

#include "utils.hpp"

void DDM::ddm_loop_log_autovec(float* ddm,
                          const fftwf_complex* stack_fft,
                          const int fft_size,
                          std::vector<int> delays,
                          utils::Options &opt) {



    // init
    for(int i = 0; i < fft_size*opt.Ntau; i++)
        ddm[i] = 0;

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < opt.Ntau; idx_delay++) {
        int t, pix;

        const fftwf_complex *i1, *i2;;
        float* ddm_cur = &ddm[idx_delay * fft_size];

        //update of ddm averages
        for(t = delays[idx_delay]; t < opt.loadNframes; t++) {
            i1 = &stack_fft[t*fft_size];
            i2 = &stack_fft[(t-delays[idx_delay])*fft_size];

            for(pix = 0; pix < fft_size; pix++) // for all pixels
                ddm_cur[pix] +=
                        std::pow(i1[pix][REAL]-i2[pix][REAL], 2.) +
                        std::pow(i1[pix][IMAG]-i2[pix][IMAG], 2.);

        }
    }

    float mean_weight = 1. / ( 2*fft_size * (opt.loadNframes-opt.Ntau) );
    for(int i=0; i<fft_size*opt.Ntau; i++)
        ddm[i] *= mean_weight;
}

void DDM::ddm_loop_log_avx(float* ddm,
                       const fftwf_complex* stack_fft,
                       const int fft_size,
                       std::vector<int> delays,
                       utils::Options &opt) {
#ifdef __AVX2__
    // init
    for(int i = 0; i < fft_size*opt.Ntau; i++)
        ddm[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < opt.Ntau; idx_delay++) {
        int t, pix;
        const float *i1, *i2;

        float* ddm_cur = &ddm[idx_delay * fft_size];

        // AVX and SSE registers for calculations
        __m256 avx_a, avx_b, avx_c;
        __m128 sse_a, sse_b, sse_c;
        const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

        for(t = delays[idx_delay]; t < opt.loadNframes; t++) {
            i1 = static_cast<const float*>(&stack_fft[t*fft_size][0]);
            i2 = static_cast<const float*>(&stack_fft[(t-delays[idx_delay])*fft_size][0]);
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

    float mean_weight = 1. / ( 2*fft_size * (opt.loadNframes-opt.Ntau) );
    for(int i=0; i<fft_size*opt.Ntau; i++)
        ddm[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}

void DDM::ddm_loop_autovec(float* ddm,
                      const fftwf_complex* stack_fft,
                      const int fft_size,
                      utils::Options &opt) {
    int tau_max = opt.Ntau;
    int N_frames =  opt.loadNframes;

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

    float mean_weight = 1. / ( 2*fft_size * (N_frames-tau_max) );
    for(int i=0; i<fft_size*tau_max; i++)
        ddm[i] *= mean_weight;
}


void DDM::ddm_loop_avx(float* ddm,
                       const fftwf_complex* stack_fft,
                       const int fft_size,
                       utils::Options &opt) {

#ifdef __AVX2__
    int tau_max = opt.Ntau;
    int N_frames =  opt.loadNframes;

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

    float mean_weight = 1. / ( 2 * fft_size * (N_frames-tau_max) );
    for(int i=0; i<fft_size*tau_max; i++)
        ddm[i] *= mean_weight;
#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}

void DDM::ddmshift(float* raw_ddm,
              float* ddm,
              int raw_width, int raw_height,
              utils::Options &opt) {
    /**
     * Shifts zero frequencies to the center of the images,
     * from raw_ddm to ddm.
     * Further comments will describe operations refering the
     * following numerotation for subblocks of the image:
     *       #########
     *       # 1 # 2 #
     *       # 4 # 3 #
     *       #########
     */


    // need to modify to handle odd y dimension

    int width = (raw_width*2)-1;
    int height = (raw_height);

    int half_width = width/2;
    int half_height = height/2;

    // for each images
    for(int i = 0; i < opt.Ntau; i++) {
        int frame = i*width*height;
        int raw_frame = i*raw_width*raw_height;

        // move block 1's lobe to block 3
        for(int ydx = 0; ydx < half_height; ydx++) {
            for (int xdx = 0; xdx < raw_width; xdx++) {
                int cpy = frame + (ydx + half_height) * width + (xdx + half_width);
                int src = raw_frame + ydx * raw_width + xdx;
                ddm[cpy] = raw_ddm[src];
            }
        }

        // move block 4's lobe to block 2
        for(int ydx = 0; ydx < half_height; ydx++) {
            for(int xdx = 0; xdx < raw_width; xdx++) {
                int cpy = frame + ydx * width + (xdx + half_width);
                int src = raw_frame + (ydx + half_height) * raw_width + xdx;
                ddm[cpy] = raw_ddm[src];
            }
        }

        // mirror 3 -> 1
        for(int ydx = 0; ydx < half_height; ydx++) {
            for(int xdx = 1; xdx < half_width; xdx++) {
                int cpy = frame + (half_height-ydx)*width + (half_width-xdx);
                int src = frame + (ydx+half_height)*width + (xdx+half_width);
                ddm[cpy] = ddm[src];
            }
        }

        // mirror 2 -> 4
        for(int ydx = 0; ydx < half_height; ydx++){
            for(int xdx = 1; xdx < half_width; xdx++) {
                int cpy = frame + (half_height+ydx)*width + (half_width-xdx);
                int src = frame + (half_height-ydx)*width + (xdx+half_width);
                ddm[cpy] = ddm[src];
            }
        }
    }
}


