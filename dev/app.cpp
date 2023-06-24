#include <fftw3.h>
#include <iostream>

#include "app.hpp"
#include "utils/stack.hpp"

#define REAL 0
#define IMAG 1

App::App(Options& options) : options(&options) {}
App::~App(){}

// à remplacer plus tard, dans une classe util fourre tout
int strentoint(std::string &s) {
    int id = -1;
    if(s == "Mono12Packed")
        id = 2;
    return id;
}

void App::run() {
    if(strentoint(this->options->encoding) == Mono12Packed) {
        Stack<float>* stack = new Stack<float>(this->options->path,
                                 strentoint(this->options->encoding),
                                 this->options->loadNframes,
                                 this->options->do_normalize);

        delete stack;
    }
    else
        std::cout << "osef" << std::endl;
};
