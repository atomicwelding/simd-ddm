#include "ddm.hpp"

#include <cmath>
#include <immintrin.h>
#include <vector>

#include <iostream>

#include "utils.hpp"

void DDM::ddm_loop_log_autovec(float* raw_ddm,
                          const fftwf_complex* stack_fft,
                          const int fft_size,
                          std::vector<int> delays,
                          utils::Options &opt) {



    // init
    for(int i = 0; i < fft_size*opt.Ntau; i++)
        raw_ddm[i] = 0;

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < opt.Ntau; idx_delay++) {
        int t, pix;

        const fftwf_complex *i1, *i2;;
        float* ddm_cur = &raw_ddm[idx_delay * fft_size];

        //update of raw_ddm averages
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
        raw_ddm[i] *= mean_weight;
}

void DDM::ddm_loop_log_avx(float* raw_ddm,
                       const fftwf_complex* stack_fft,
                       const int fft_size,
                       std::vector<int> delays,
                       utils::Options &opt) {
#ifdef __AVX2__
    // init
    for(int i = 0; i < fft_size*opt.Ntau; i++)
        raw_ddm[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < opt.Ntau; idx_delay++) {
        int t, pix;
        const float *i1, *i2;

        float* ddm_cur = &raw_ddm[idx_delay * fft_size];

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
        raw_ddm[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}

void DDM::ddm_loop_autovec(float* raw_ddm,
                      const fftwf_complex* stack_fft,
                      const int fft_size,
                      utils::Options &opt) {
    int tau_max = opt.Ntau;
    int N_frames =  opt.loadNframes;

    for(int i=0; i<fft_size*tau_max; i++)
        raw_ddm[i] = 0;

    // parallelize on every tau, and let the compiler vectorize the code as much as
    // possible
#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int tau = 0; tau < tau_max; tau++) {
        // stack indices
        int t, pix;

        // shortcut ptrs
        const fftwf_complex *i1, *i2;
        float* ddm_cur = &raw_ddm[tau * fft_size];

        // Update of raw_ddm averages
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
        raw_ddm[i] *= mean_weight;
}


void DDM::ddm_loop_avx(float* raw_ddm,
                       const fftwf_complex* stack_fft,
                       const int fft_size,
                       utils::Options &opt) {

#ifdef __AVX2__
    int tau_max = opt.Ntau;
    int N_frames =  opt.loadNframes;

    for(int i=0; i<fft_size*tau_max; i++)
        raw_ddm[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int tau=0; tau<tau_max; tau++) {
        int t, pix;
        const float *i1, *i2;
        float* ddm_cur = &raw_ddm[tau * fft_size];

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
        raw_ddm[i] *= mean_weight;
#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}

void DDM::ddmshift(const float* raw_ddm,
                   const int raw_width, const int raw_height,
                   float* ddm,
                   utils::Options &opt) {
    /**
     * After computing DDM using appropriate methods (e.g. DDM::ddm_loop_avx),
     * the "raw ddm images" looks like the following (due to symmetries of the fourier transform, we have "fft lobes"
     * in corners of the image) :
                        ┌───────────────────────┐
                        │     xx                │
                        │      x                │
                        │     xx                │
                        │    xx                 │
                        │xxxx                   │
                        │                       │
                        │    raw ddm image      │
                        │                       │
                        │                       │
                        │xxx                    │
                        │   xxx                 │
                        │     xx                │
                        │      x                │
                        │      x                │
                        └───────────────────────┘
      The width is half width of the real input signal + 1, while height is still the same.
      We then have raw_width = input_width//2 + 1, raw_height = input_height.

      We want to have, as a final result, a square image with the zero frequency centered. To achieve this,
      you need to provide a buffer of odd width and height that we call ddm.

      Let us represent our buffers :

                      ┌───┐    ┌────────┐
                      │ 1 │ -> │ 1'   2'│
                      │ 4 │ -> │ 4'   3'│
                      └───┘    └────────┘
                      beg buf   final buf

      Region 1,4 contains the aforementioned "fft lobes". We apply those operations :
          * Copy 1 -> 3'
          * Copy 4 -> 2'
          * Mirror 3'->1' and 2'->4'
     */

    int Nx_Half = raw_width-1;
    int Nx_Full = 2*Nx_Half + 1;

    int Ny_Half = raw_height/2;
    int Ny_Full = 2*Ny_Half+1;

    auto raw_idx = [&](int i, int j, int k) -> int {
        return i*raw_height*raw_width + j*raw_width + k;
    };

    auto idx = [&](int i, int jbis, int kbis) -> int {
        return i*Nx_Full*Ny_Full + jbis*Nx_Full + kbis;
    };

    for(int i = 0; i < opt.Ntau; i++) {

        // Copy 1 -> 3'
        for(int j = 0; j < (Ny_Half+1); j++) {
            for(int k = 0; k < Nx_Half+1; k++) {
                int jbis = j + Ny_Half;
                int kbis = k + Nx_Half;
                ddm[idx(i,jbis,kbis)] = raw_ddm[raw_idx(i,j,k)];
            }
        }

        // Copy 4 -> 2'
        for(int j = raw_height-Ny_Half; j < raw_height; j++) {
            for(int k = 0; k < Nx_Half+1; k++) {
                int jbis = j-raw_height+Ny_Half;
                int kbis = k+Nx_Half;
                ddm[idx(i,jbis,kbis)] = raw_ddm[raw_idx(i,j,k)];
            }
        }

        // Mirror 3'->1' and 2'->4'
        for(int jbis = 0; jbis < Ny_Full; jbis++) {
            for(int kbis = Nx_Half+1; kbis < Nx_Full; kbis++) {
                int jbisbis = Ny_Full - 1 - jbis;
                int kbisbis = Nx_Full - 1 - kbis;
                ddm[idx(i, jbisbis, kbisbis)] = ddm[idx(i, jbis, kbis)];
            }
        }

    }

}


