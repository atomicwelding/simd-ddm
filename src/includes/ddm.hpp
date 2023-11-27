#ifndef DEV_DDM
#define DEV_DDM

#include <fftw3.h>
#include "stack.hpp"

#include "utils.hpp"

class DDM {
public:
    int ddm_width;
    int ddm_height;

    utils::Options& options;

    DDM(Stack* stack, utils::lspace<double> delays, utils::Options& options);

    utils::lspace<double> delays;

     int ddm_size() const {
        return this->ddm_width * this->ddm_height;
    }

    void save_ddm_buffer(std::string path);


private:
    int raw_ddm_width;
    int raw_ddm_height;

    int fft_size;



    float* raw_ddm_buffer;
    float* ddm_buffer;

    int raw_ddm_size() const  {
        return this->raw_ddm_width * this->raw_ddm_height;
    }

    fftwf_complex *stack_fft;

    void ddm_loop_avx_delays();
    void ddm_loop_autovec_delays();

    /*void ddm_loop_autovec(float* raw_ddm, const fftwf_complex* stack_fft,
                          const int fft_size, utils::Options &opt);

    void ddm_loop_avx(float* raw_ddm, const fftwf_complex* stack_fft,
                      const int fft_size, utils::Options &opt);

    void ddm_loop_log_autovec(float* raw_ddm,
                              const fftwf_complex* stack_fft,
                              const int fft_size,
                              std::vector<int> delays,
                              utils::Options &opt);

    void ddm_loop_log_avx(float* raw_ddm,
                          const fftwf_complex* stack_fft,
                          const int fft_size,
                          std::vector<int> delays,
                          utils::Options &opt);*/


    void ddmshift();

};
#endif //DEV_DDM
