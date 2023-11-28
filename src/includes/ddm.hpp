#ifndef DEV_DDM
#define DEV_DDM

#include <fftw3.h>
#include "stack.hpp"

#include "utils.hpp"
#include "delays.hpp"


template<typename T>
class DDM {
public:
    int ddm_width;
    int ddm_height;


    DDM(Stack *stack, Delays<T>& delays, utils::Options& options);
    ~DDM();

    Delays<T> delays;

    void save();



private:
    Stack *stack;
    utils::Options options;

    int raw_ddm_width;
    int raw_ddm_height;

    int fft_size;

    float* raw_ddm_buffer;
    float* ddm_buffer;

    int raw_ddm_size() const  {
        return this->raw_ddm_width * this->raw_ddm_height;
    }

    void computeFFT();
    void computeDDM();


    fftwf_complex *stack_fft;

    void ddm_loop_avx_delays();
    void ddm_loop_autovec_delays();

    int ddm_size() const {
        return this->ddm_width * this->ddm_height;
    }

    void ddmshift();

};
#endif //DEV_DDM
