#ifndef APP_H
#define APP_H

#include <fftw3.h>
#include <string>

#include "stack.hpp"
#include "utils.hpp"



class App {
public:
    App(utils::Options &options);
    ~App();

    void run();

private:
	void ddm_loop_avx(float* ddm, const fftwf_complex* stack_fft, const int fft_size);
    void fit_routine(Stack<float>* stack, float* ddm,int tauMax, int fft_size);

    static double exp_to_fit(double x, double A, double B, double tau_k);

    utils::Options* options;
};

#endif // APP_H
