#ifndef APP_H
#define APP_H

#include <fftw3.h>
#include <string>

#include "stack.hpp"

struct Options {
    std::string path;
    int loadNframes;
    std::string encoding;
    bool do_normalize;
    int tauMax;
    std::string pathOutput;
    bool do_fit;
};

class App {
public:
    App(Options &options);
    ~App();

    void run();

private:
	void ddm_loop_autovec(float* ddm, const fftwf_complex* stack_fft, const int fft_size);
	void ddm_loop_avx(float* ddm, const fftwf_complex* stack_fft, const int fft_size);
    void fit_routine(Stack<float>* stack, float* ddm,int tauMax, int fft_size, int N);

    static double exp_to_fit(double x, double A, double B, double tau_k);

    Options* options;
};

#endif // APP_H
