#ifndef DEV_DDM
#define DEV_DDM

#include <fftw3.h>
#include "utils.hpp"

namespace DDM {
    void ddm_loop_autovec(float* raw_ddm, const fftwf_complex* stack_fft, const int fft_size, utils::Options &opt);
    void ddm_loop_avx(float* raw_ddm, const fftwf_complex* stack_fft, const int fft_size, utils::Options &opt);

    void ddm_loop_log_autovec(float* raw_ddm,
                              const fftwf_complex* stack_fft,
                              const int fft_size,
                              std::vector<int> delays,
                              utils::Options &opt);

    void ddm_loop_log_avx(float* raw_ddm,
                          const fftwf_complex* stack_fft,
                          const int fft_size,
                          std::vector<int> delays,
                          utils::Options &opt);

    void ddmshift(float* raw_ddm,
                   float* ddm,
                   int width,
                   int height,
                   utils::Options &opt);
};

#endif //DEV_DDM
