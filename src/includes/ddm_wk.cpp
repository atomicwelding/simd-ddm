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

    // Size of zero-padded temporal vector used for temporal FFT calculation
    // It is a multiple of 2 for efficiency, and it allows an unbiased calculation
    // of the autocorrelation up to max_lag_shift.
    const size_t N2 = std::pow(2.,std::ceil(std::log2(Nt+max_lag_shift)));

    fftwf_complex* in = utils::allocate_complex_float_array(N2);
    fftwf_complex* out = utils::allocate_complex_float_array(N2);
	fftwf_plan forward_plan = fftwf_plan_dft_1d(N2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_plan backward_plan = fftwf_plan_dft_1d(N2, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    #pragma omp parallel for schedule(nonmonotonic:dynamic) 
    for(int iq=0; iq<raw_ddm_size; iq++) {
        // Allocation of local buffers for this thread
        fftwf_complex* data_in = utils::allocate_complex_float_array(N2);
        fftwf_complex* data_out = utils::allocate_complex_float_array(N2);

        // Copy data and pad with zeros
        int it;
        for(it=0; it<Nt; it++) {
            data_in[it][0] = stack_fft[raw_ddm_size*it+iq][0];
            data_in[it][1] = stack_fft[raw_ddm_size*it+iq][1];
        }
        for(it=Nt; it<N2; it++) {
            data_in[it][0] = 0;
            data_in[it][1] = 0;
        }

        // Calculation of accumulated sums for the first part of the ISF
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
            raw_ddm_buffer[raw_ddm_size*lag_idx+iq] = asums[lag_shifts[lag_idx]];

        // Calculation of autocorrelation function from Wiener-Khinchin theorem for
        // the second part of the ISF
        fftwf_execute_dft(forward_plan, data_in, data_out);
        for(it=0; it<N2; it++) {
            data_out[it][0] =
                (data_out[it][0]*data_out[it][0]+data_out[it][1]*data_out[it][1])/N2;
            data_out[it][1] = 0;
        }
        fftwf_execute_dft(backward_plan, data_out, data_in);

        // ISF is obtained by substracting 2x the real-part of the autocorrelation
        // function from the accumulated sums already calculated above, using normalization factors
        // mean_weight/(Nt-lag_shift).
        float mean_weight = 1. / ( 2*raw_ddm_size );
        for(int lag_idx=0; lag_idx<n_lags; lag_idx++) {
            raw_ddm_buffer[raw_ddm_size*lag_idx+iq] -= 2*data_in[lag_shifts[lag_idx]][0];
            raw_ddm_buffer[raw_ddm_size*lag_idx+iq] *= mean_weight/(Nt-lag_shifts[lag_idx]);
        }

        // destroy local thread buffer
        fftwf_free(data_in);
        fftwf_free(data_out);
    }

    fftwf_destroy_plan(forward_plan);
    fftwf_destroy_plan(backward_plan);
    fftwf_free(in);
    fftwf_free(out);
}