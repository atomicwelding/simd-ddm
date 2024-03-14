#include "ddm_diff.hpp"
#include "timer.hpp"

DDMDiff::DDMDiff(Stack &stack, utils::Options& options) :
    DDM(stack, options) {
        
    compute_DDM();
    ddm_shift();
}

void DDMDiff::compute_DDM() {

    Timer timer;
    timer.start();

#ifdef __AVX512F__
    std::cout << "* Computing DDM [diff,AVX512]...  " << std::flush;
    ddm_loop_avx512();
#elif __AVX2__
    std::cout << "* Computing DDM [diff,AVX2]...    " << std::flush;
    ddm_loop_avx2();
#else
    std::cout << "* Computing DDM [diff,autovec]... " << std::flush;
    ddm_loop_autovec();
#endif

    std::cout << timer.elapsedSec() << "s" << std::endl;
}

void DDMDiff::ddm_loop_avx512() {
    
#ifdef __AVX512F__
    // init
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    int chunk_size = max_lag_shift;
    float* unrolled_stack_fft = static_cast<float*>(&stack_fft[0][0]);

    for(int offset=max_lag_shift; offset<Nt; offset+=chunk_size) {
        #pragma omp parallel for schedule(nonmonotonic:dynamic)
        for(int lag_idx=0; lag_idx<n_lags; lag_idx++) {
            int it, buf_idx;
            const float *i1, *i2;

            float* ddm_cur = &raw_ddm_buffer[raw_ddm_size*lag_idx];

            // AVX512 and AVX2 registers for calculations
            __m512 m512_a, m512_b, m512_c;
            __m256 m256_a, m256_b, m256_c;
            const auto perm = _mm512_set_epi32(15,13,11,9,7,5,3,1,14,12,10,8,6,4,2,0);

            int max_it = std::min(offset+chunk_size, Nt);
            int lag_shift = lag_shifts[lag_idx];
            for(it=offset; it<max_it; it++) {
                i1 = &unrolled_stack_fft[raw_ddm_size*2*it];
                i2 = &unrolled_stack_fft[raw_ddm_size*2*(it-lag_shift)];
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

    float mean_weight = 1. / ( 2*raw_ddm_size * (Nt-max_lag_shift) );
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX512");
#endif
}


void DDMDiff::ddm_loop_avx2() {
    
#ifdef __AVX2__
    // init
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    float* unrolled_stack_fft = static_cast<float*>(&stack_fft[0][0]);

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int lag_idx=0; lag_idx<n_lags; lag_idx++) {
        int it, buf_idx;
        const float *i1, *i2;

        float* ddm_cur = &raw_ddm_buffer[raw_ddm_size*lag_idx];

        // AVX and SSE registers for calculations
        __m256 avx_a, avx_b, avx_c;
        __m128 sse_a, sse_b, sse_c;
        const auto perm = _mm256_set_epi32(7,5,3,1,6,4,2,0);

        int lag_shift = lag_shifts[lag_idx];
        for(it=max_lag_shift; it<Nt; it++) {
            i1 = &unrolled_stack_fft[raw_ddm_size*2*it];
            i2 = &unrolled_stack_fft[raw_ddm_size*2*(it-lag_shift)];
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

    float mean_weight = 1. / ( 2*raw_ddm_size * (options.N_frames-max_lag_shift) );
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;

#else
    throw std::runtime_error("Your cpu doesn't support AVX2");
#endif
}


void DDMDiff::ddm_loop_autovec() {
    
    // init
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] = 0;

    #pragma omp parallel for schedule(nonmonotonic:dynamic)
    for(int lag_idx=0; lag_idx<n_lags; lag_idx++) {
        int it, iq;

        const fftwf_complex *i1, *i2;;
        float* ddm_cur = &raw_ddm_buffer[raw_ddm_size*lag_idx];

        int lag_shift = lag_shifts[lag_idx];
        for(it=max_lag_shift; it<options.N_frames; it++) {
            i1 = &stack_fft[raw_ddm_size*it];
            i2 = &stack_fft[raw_ddm_size*(it-lag_shift)];

            for(iq=0; iq<raw_ddm_size; iq++) // for all pixels
                ddm_cur[iq] +=
                        std::pow(i1[iq][REAL] - i2[iq][REAL], 2.) +
                        std::pow(i1[iq][IMAG] - i2[iq][IMAG], 2.);
        }
    }

    float mean_weight =
		1. / ( 2*raw_ddm_size * (Nt-max_lag_shift) );
    for(size_t i=0; i<raw_ddm_size*n_lags; i++)
        raw_ddm_buffer[i] *= mean_weight;
}