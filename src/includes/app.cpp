#include <vector>
#include <iostream>

#include <tinytiffwriter.h>
#include <omp.h>
#include <fftw3.h>
#include <cmath>
#include <ranges>
#include <fstream>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm.hpp"
#include "fitting.hpp"

App::App(utils::Options& options) : options(&options) {}
App::~App()= default;

void App::run() {

    if(utils::stoe(this->options->encoding) != Mono12Packed)
        throw std::runtime_error("Encoding not supported yet");

	Timer timer;

    std::cout << "* Loading images..." << std::flush;
	timer.start();
    auto* stack = new Stack(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->doNormalize,
                                           this->options->binFactor);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    // not the perfect place to do that
    int Nt = stack->times.size();
    double mean_sampling_time = (stack->times[Nt-1] - stack->times[0])/(Nt-1);
    std::vector<int> delays_filtered;
    if(this->options->doLogScale) {
        auto delays_time = utils::log_delays_in_time<float>(mean_sampling_time, this->options->delayMax, this->options->Ntau);
        auto delays = utils::log_delays_indexes(mean_sampling_time, this->options->delayMax, this->options->Ntau);


        if(delays_time[delays_time.size()-1] >= this->options->delayMax || delays[delays.size() - 1] >= this->options->loadNframes) {
            throw std::runtime_error("Error : the file you want to process is shorter than the max delay you're asking");
        }

        std::copy_if(delays.begin(), delays.end(), std::back_inserter(delays_filtered),
                     [&](int x) { return x < this->options->loadNframes; });
        std::cout << "/!\\ Ntau has changed, from " << this->options->Ntau << std::flush;
        this->options->Ntau = delays_filtered.size();
        std::cout << " to  " << this->options->Ntau << " due to log scaling !" << std::endl;


        std::ofstream delays_array_file;
        delays_array_file.open(this->options->pathOutput  + "_time_delays.dat");
        delays_array_file << "Shift" << "," << "Real_Time" << "\n";
        for(int i = 0; i < this->options->Ntau; i++) {
            delays_array_file << delays[i] << "," << delays_time[i] << "\n";
        }

        delays_array_file.close();
    }


    std::cout << "* Creating FFTW plan with " << omp_get_max_threads() << " threads..." << std::flush;

    int r = fftwf_init_threads();
    if(r == 0)
        throw std::runtime_error("Can't spawn threads");

	timer.start();
    fftwf_plan_with_nthreads(omp_get_max_threads());

    int rank = 2;
    int n_in[] = {stack->aoi_height, stack->aoi_width};
    int n_out[] = {stack->aoi_height, stack->aoi_width/2+1};
    int fft_size = n_out[0] * n_out[1];
    fftwf_complex* stack_fft  = fftwf_alloc_complex(fft_size * this->options->loadNframes);
    fftwf_plan plan = fftwf_plan_many_dft_r2c(rank, n_in, this->options->loadNframes,
                                              stack->images, n_in,
                                              1, stack->image_size,
                                              stack_fft, n_out,
                                              1, fft_size,
                                              FFTW_ESTIMATE);
	std::cout << "  " << timer.elapsedSec() << "s" << std::endl;


    std::cout << "* Performing DFT..." << std::flush;
    timer.start();
    fftwf_execute(plan);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;


    std::cout << "* Computing DDM differences..." << std::flush;
    timer.start();
    float* raw_ddm = fftwf_alloc_real(this->options->Ntau * fft_size);


    if(this->options->doLogScale) {
        #ifdef __AVX2__
                DDM::ddm_loop_log_avx(raw_ddm, stack_fft, fft_size, delays_filtered, *(this->options));
        #else
                DDM::ddm_loop_log_autovec(ddm, stack_fft, fft_size, delays_filtered, *(this->options));
        #endif
    } else {
        #ifdef __AVX2__
            DDM::ddm_loop_avx(raw_ddm, stack_fft, fft_size, *(this->options));
        #else
            DDM::ddm_loop_autovec(ddm, stack_fft, fft_size, *(this->options*));
        #endif
    }

	std::cout << "          " << timer.elapsedSec() << "s" << std::endl;

    std::cout << "* Mirroring DDM images ..." << std::flush;
    timer.start();

    int ddm_width = (n_out[1] * 2) - 1;
    int ddm_height = (n_out[0])+1;
    float* ddm = fftwf_alloc_real(this->options->Ntau * ddm_height * ddm_width);

    DDM::ddmshift(raw_ddm, n_out[1], n_out[0], ddm, *(this->options));
    std::cout << "              " << timer.elapsedSec() << "s" << std::endl;



    std::cout << "* Writing files ..." << std::flush;
    timer.start();
    TinyTIFFWriterFile* tif = TinyTIFFWriter_open(this->options->pathOutput.c_str(), 32, TinyTIFFWriter_Float,
                                                  1, ddm_width, ddm_height,
                                                  TinyTIFFWriter_Greyscale);
    if(!tif)
        throw std::runtime_error("Can't write files");

    for(int frame = 0; frame < this->options->Ntau; frame++) {
        float* data = &ddm[frame * ddm_height * ddm_width];
        TinyTIFFWriter_writeImage(tif, data);
    }
    timer.stop();
    std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;



    if(this->options->doFit) {
        std::cout << "* Fitting..." << std::flush;

        timer.start();
        auto exp_to_fit  = [](double tau, double A, double B, double f) -> double {
                    return A*(1-std::exp(-tau*f))+B;
        };

        // need to use it now in the fit routine
        auto ROI = fit::find_ROI(exp_to_fit, ddm, ddm_width, ddm_height, *this->options, mean_sampling_time);
        fit::fit_routine(exp_to_fit, stack, ddm, ddm_width, ddm_height, this->options->Ntau, fft_size);

        timer.stop();
        std::cout << "                           " << timer.elapsedSec() << "s" << std::endl;
    }


    std::cout << "* Cleaning ..." << std::endl;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    fftwf_free(stack_fft);
    delete stack;
    TinyTIFFWriter_close(tif);
};


