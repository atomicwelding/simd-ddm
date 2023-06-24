#include <fftw3.h>
#include <iostream>

#include "app.hpp"
#include "utils/stack.hpp"

#define REAL 0
#define IMAG 1

App::App(Options& options) : options(&options) {}
App::~App(){}

void App::run() {
    Stack* stack = new Stack(options->path, options->loadNframes, options->do_normalize);
    delete stack;
};
