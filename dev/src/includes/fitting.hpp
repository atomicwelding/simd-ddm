#ifndef FITTING_HPP
#define FITTING_HPP

#include "stack.hpp"

namespace fit {
    void fit_routine(Stack<float>* stack, float* ddm, int tau_max, int fft_size);
    double exp_to_fit(double tau, double A, double B, double f);
}

#endif //FITTING_HPP
