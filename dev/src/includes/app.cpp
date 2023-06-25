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

    Stack<float>* stack = new Stack<float>(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->do_normalize);
    delete stack;
};
