#ifndef DEV_DDM
#define DEV_DDM

#include <stdexcept>
#include <fftw3.h>
#include "utils.hpp"

namespace DDM {
    void ddm_loop_autovec(float* ddm, const fftwf_complex* stack_fft, const int fft_size, utils::Options &opt);
    void ddm_loop_avx(float* ddm, const fftwf_complex* stack_fft, const int fft_size, utils::Options &opt);

    void ddm_loop_log_autovec(float* ddm,
                              const fftwf_complex* stack_fft,
                              const int fft_size,
                              std::vector<int> delays,
                              utils::Options &opt);

    void ddm_loop_log_avx(float* ddm,
                              const fftwf_complex* stack_fft,
                              const int fft_size,
                              std::vector<int> delays,
                              utils::Options &opt);
};

#endif //DEV_DDM
