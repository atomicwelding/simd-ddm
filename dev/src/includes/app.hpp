#ifndef APP_H
#define APP_H

#include <fftw3.h>
#include <string>

struct Options {
    std::string path;
    int loadNframes;
    std::string encoding;
    bool do_normalize;
    int tauMax;
    std::string pathOutput;
};

class App {
public:
    App(Options &options);
    ~App();

    void run();

private:
	void ddm_loop_autovec(float* ddm, const fftwf_complex* stack_fft, const int fft_size);
	void ddm_loop_avx(float* ddm, const fftwf_complex* stack_fft, const int fft_size);

    Options* options;
};

#endif // APP_H
