#include "ddm.hpp"
#include <stdio.h>
#include <immintrin.h>

void DDM::ddm_loop_avx_delays() {
#ifdef __AVX2__
    // init
    for(int i = 0; i < fft_size*options.Ntau; i++)
        raw_ddm_buffer[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < options.Ntau; idx_delay++) {
        int t, pix;
        const float *i1, *i2;

        float* ddm_cur = &raw_ddm_buffer[idx_delay * fft_size];

        // AVX and SSE registers for calculations
        __m256 avx_a, avx_b, avx_c;
        __m128 sse_a, sse_b, sse_c;
        const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

        for(t = delays.index.back(); t < options.loadNframes; t++) {
            i1 = static_cast<const float*>(&stack_fft[t*fft_size][0]);
            i2 = static_cast<const float*>(&stack_fft[(t-delays.index[idx_delay]-1)*fft_size][0]);
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

    float mean_weight = 1. / ( 2*fft_size * (options.loadNframes-options.Ntau) );
    for(int i=0; i<fft_size*options.Ntau; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}

void DDM::ddm_loop_autovec_delays() {
    // init
    for(int i = 0; i < raw_ddm_size() * options.Ntau; i++)
        raw_ddm_buffer[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < options.Ntau; idx_delay++) {
        int t, pix;

        const fftwf_complex *i1, *i2;;
        float* ddm_cur = &raw_ddm_buffer[idx_delay * fft_size];

        //update of raw_ddm averages
        for(t = delays.index.back(); t < options.loadNframes; t++) {
            i1 = &stack_fft[t*fft_size];
            i2 = &stack_fft[(t-delays.index[idx_delay]-1)*fft_size];

            for(pix = 0; pix < fft_size; pix++) // for all pixels
                ddm_cur[pix] +=
                        std::pow(i1[pix][REAL] - i2[pix][REAL], 2.) +
                        std::pow(i1[pix][IMAG] - i2[pix][IMAG], 2.);

        }
    }

    float mean_weight = 1. / ( 2*fft_size * (options.loadNframes-options.Ntau) );
    for(int i=0; i<fft_size*options.Ntau; i++)
        raw_ddm_buffer[i] *= mean_weight;
}

/*void DDM::ddm_loop_avx() {
#ifdef __AVX2__
    // init
    for(int i = 0; i < fft_size*options.Ntau; i++)
        raw_ddm_buffer[i] = 0;

#pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int idx_delay = 0; idx_delay < options.Ntau; idx_delay++) {
        int t, pix;
        const float *i1, *i2;

        float* ddm_cur = &raw_ddm_buffer[idx_delay * fft_size];

        // AVX and SSE registers for calculations
        __m256 avx_a, avx_b, avx_c;
        __m128 sse_a, sse_b, sse_c;
        const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

        for(t = delays.index[idx_delay]; t < options.loadNframes; t++) {
            i1 = static_cast<const float*>(&stack_fft[t*fft_size][0]);
            i2 = static_cast<const float*>(&stack_fft[(t-delays.index[idx_delay])*fft_size][0]);
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

    float mean_weight = 1. / ( 2*fft_size * (options.loadNframes-options.Ntau) );
    for(int i=0; i<fft_size*options.Ntau; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}
*/

/*void DDM::ddm_loop_autovec(float* raw_ddm,
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
*/