#include <fftw3.h>
#include <iostream>

#include "app.hpp"
#include "utils/stack.hpp"

#define REAL 0
#define IMAG 1

App::App(Options& options) : options(&options) {}
App::~App(){}

void App::run() {
    Stack* stack = new Stack(options->path);

    // struct passée par reference aussi ?
    for(int i = 1; i <= this->options->loadNframes; ++i) {
        std::cout << i << "/" << this->options->loadNframes << std::endl;
        // templating pas ouf
        stack->load_next_frame<float>();
    }

    delete stack;
};
