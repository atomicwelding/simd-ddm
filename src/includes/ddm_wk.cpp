#include "ddm_wk.hpp"
#include "timer.hpp"

DDMWK::DDMWK(Stack &stack, utils::Options& options) :
    DDM(stack, options) {
        
    compute_DDM();
    ddm_shift();
}

void DDMWK::compute_DDM() {

    Timer timer;
    timer.start();

    std::cout << "* Computing DDM [wk,autovec]...   " << std::flush;

    ddm_loop_autovec();

    std::cout << timer.elapsedSec() << "s" << std::endl;
}

void DDMWK::ddm_loop_autovec() {

    const long N2 = std::pow(2.,std::ceil(std::log2(Nt+max_lag_shift)));

#ifdef __AVX512F__
    fftwf_complex* in = (fftwf_complex*) _mm_malloc(2*N2*sizeof(float), 64);
    fftwf_complex* out = (fftwf_complex*) _mm_malloc(2*N2*sizeof(float), 64);
#else
    fftwf_complex* in = fftwf_alloc_complex(2*N2);
    fftwf_complex* out = fftwf_alloc_complex(2*N2);
#endif

	fftwf_plan forward_plan = fftwf_plan_dft_1d(N2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_plan backward_plan = fftwf_plan_dft_1d(N2, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    #pragma omp parallel for schedule(nonmonotonic:dynamic) 
    for(int iq=0; iq<raw_ddm_size; iq++) {

#ifdef __AVX512F__
        fftwf_complex* data_in = (fftwf_complex*) _mm_malloc(2*N2*sizeof(float), 64);
        fftwf_complex* data_out = (fftwf_complex*) _mm_malloc(2*N2*sizeof(float), 64);
#else
        fftwf_complex* data_in = fftwf_alloc_complex(2*N2);
        fftwf_complex* data_out = fftwf_alloc_complex(2*N2);
#endif

        int it;
        for(it=0; it<Nt; it++) {
            data_in[it][0] = stack_fft[it*raw_ddm_size+iq][0];
            data_in[it][1] = stack_fft[it*raw_ddm_size+iq][1];
        }
        for(it=Nt; it<N2; it++) {
            data_in[it][0] = 0;
            data_in[it][1] = 0;
        }

        std::vector<float> asums(max_lag_shift+1,0);
        for(it=0; it<Nt-max_lag_shift; it++) {
            asums[max_lag_shift] +=
                    data_in[it][0]*data_in[it][0] +
                    data_in[it][1]*data_in[it][1] +
                    data_in[Nt-1-it][0]*data_in[Nt-1-it][0] +
                    data_in[Nt-1-it][1]*data_in[Nt-1-it][1];
        }
        for(it=max_lag_shift-1; it>0; it--) {
            asums[it] = asums[it+1] +
                    data_in[it][0]*data_in[it][0] +
                    data_in[it][1]*data_in[it][1] +
                    data_in[Nt-1-it][0]*data_in[Nt-1-it][0] +
                    data_in[Nt-1-it][1]*data_in[Nt-1-it][1];
        }
        for(int lag_idx=0; lag_idx<n_lags; lag_idx++)
            raw_ddm_buffer[lag_idx*raw_ddm_size+iq] =
                asums[lag_shifts[lag_idx]]/(Nt-lag_shifts[lag_idx]);

        fftwf_execute_dft(forward_plan, data_in, data_out);
        for(it=0; it<N2; it++) {
            data_out[it][0] =
                (data_out[it][0]*data_out[it][0]+data_out[it][1]*data_out[it][1])/N2;
            data_out[it][1] = 0;
        }
        fftwf_execute_dft(backward_plan, data_out, data_in);

        float mean_weight = 1. / ( 2*raw_ddm_size );
        for(long lag_idx=0; lag_idx<n_lags; lag_idx++) {
            raw_ddm_buffer[lag_idx*raw_ddm_size+iq] -=
                2*data_in[lag_shifts[lag_idx]][0]/(Nt-lag_shifts[lag_idx]);
            raw_ddm_buffer[lag_idx*raw_ddm_size+iq] *= mean_weight;
        }

        fftwf_free(data_in);
        fftwf_free(data_out);
    }

    fftwf_destroy_plan(forward_plan);
    fftwf_destroy_plan(backward_plan);
    fftwf_free(in);
    fftwf_free(out);
}