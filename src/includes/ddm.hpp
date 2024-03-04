#ifndef DEV_DDM
#define DEV_DDM

#include <fftw3.h>

#include "stack.hpp"
#include "utils.hpp"
#include "delays.hpp"

template<typename T>
class DDM {
public:
    DDM(Stack &stack, Delays<T>& delays, utils::Options& options);
    ~DDM();

    Delays<T>& delays;

    void save();

    auto exposeDdmBuffer() const {
        return ddm_buffer;
    }

    const int raw_ddm_width, raw_ddm_height, raw_ddm_size;
    const int ddm_width, ddm_height, ddm_size;

private:
    Stack& stack;
    utils::Options options;

    float* raw_ddm_buffer;
    float* ddm_buffer;

    void computeFFT();
    void computeDDM();


    fftwf_complex *stack_fft;

    void ddm_loop_avx_delays();
    void ddm_loop_autovec_delays();



    void ddmshift();

};
#endif //DEV_DDM
