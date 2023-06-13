#include <fftw3.h>
#include <iostream>

#include "app.hpp"
#include "utils/stack.hpp"

#define REAL 0
#define IMAG 1

App::App(const std::string& path) : path(path) {}
App::~App(){}

void App::run() {
    Stack* stack = new Stack(path);



    // |- shitty
    // v

    /*int width = stack->aoi_width + stack->stride;
    int height = stack->aoi_height;

    stack->load_next_frame();

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * stack->current_frame->im_bytes);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * stack->current_frame->im_bytes);

    fftw_plan plan = fftw_plan_dft_2d(width, height, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for(int i = 0; i < stack->current_frame->im_bytes; i++) {
        in[i][REAL] = stack->current_frame->im_data[i];
        in[i][IMAG] = 0;
    }

    fftw_execute(plan);
    for(int i = 0; i < stack->current_frame->im_bytes; i++) {
        std::cout << "real: " << out[i][REAL] << " imag: " << out[i][IMAG] << std::endl;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);*/
    delete stack;
};
