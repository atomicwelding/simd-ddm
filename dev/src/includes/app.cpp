#include <omp.h>
#include <fftw3.h>

#include <iostream>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"

App::App(Options& options) : options(&options) {}
App::~App(){}

void App::run() {
    if(utils::stoe(this->options->encoding) != Mono12Packed) {
        std::cout << "Encoding not supported yet" << std::endl;
        return;
    }

    std::cout << "* Loading images ..." << std::endl;
    Stack<float>* stack = new Stack<float>(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->do_normalize);

    // TODO:
    //     move code to stack?
    std::cout << "* Spawning " << omp_get_max_threads() << " threads ..." << std::endl;
    int r = fftwf_init_threads();
    if(r == 0) {
        std::cout << "err" << std::endl;
        return;
    }

    fftwf_plan_with_nthreads(omp_get_max_threads());
    fftwf_complex* out  = fftwf_alloc_complex(stack->image_size * this->options->loadNframes);
    fftwf_plan plan = fftwf_plan_dft_r2c_2d(stack->aoi_width,
                                            stack->aoi_height*this->options->loadNframes,
                                            stack->N_images_buffer,
                                            out,
                                            FFTW_ESTIMATE);

    std::cout << "* Performing DFT ..." << std::endl;
    fftwf_execute(plan);

    std::cout << "* Done !" << std::endl;

    std::cout << "* Cleaning ..." << std::endl;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    delete stack;
    fftw_free(out);
};
