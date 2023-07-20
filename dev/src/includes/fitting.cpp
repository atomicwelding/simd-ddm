#include "fitting.hpp"
#include "curve_fit.hpp"

#include <cmath>
#include <vector>
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <tinytiffwriter.h>


double fit::exp_to_fit(double tau, double A, double B, double f) {
    return A*(1-std::exp(-tau*f))+B;
}

void fit::fit_routine(Stack<float>* stack, float* ddm, int tau_max, int fft_size) {
    /**
     * Creates a 3-stacked TIFF images named "fit.tif",
     * containing values of parameters of the exponential fit ;
     * params to be fitted [A->f->B] :  A(1-exp[-tau*f])+B
     */
    int Nt = stack->times.size();
    int Ny = stack->aoi_height;
    int Nx = stack->aoi_width/2+1;

    double mean_sampling_time = (stack->times[Nt-1] - stack->times[0])/(Nt-1);

    std::vector<double> times;
    for(int tau = 1; tau <= tau_max; tau++)
        times.push_back(tau*mean_sampling_time);

    float *parameters = fftwf_alloc_real(fft_size*3);

    // shortcut ptrs
    float* As = parameters;
    float* Bs  = &parameters[fft_size];
    float* fs = &parameters[fft_size * 2];
    double A, B, f;

    std::vector<double> I_vals(tau_max);
    double I_min, I_max, I_thresh;
    int tau, iy, ix;

    for(int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            for(tau = 0; tau < tau_max; tau++)
                I_vals[tau] = ddm[ix + iy*Nx + tau * fft_size];

            // Estimation of the mode amplitude and noise
            A = I_vals[tau_max - 1];
            B = 0.;

            // Estimation of the relaxation frequency
            I_min = *std::min_element(I_vals.begin(), I_vals.end());
            I_max = *std::max_element(I_vals.begin(), I_vals.end());
            I_thresh = I_min+(I_max-I_min)*(1-std::exp(-1));
            for(tau = 0; tau < tau_max; tau++) {
                if(I_vals[tau]>I_thresh)
                    break;
            }
            f = 1./times[tau];

            // Curve fitting
            auto params_fitted = curve_fit(fit::exp_to_fit, {A, B, f}, times, I_vals);
            As[ix + iy*Nx] = params_fitted[0];
            Bs[ix + iy*Nx] = params_fitted[1];
            fs[ix + iy*Nx] = params_fitted[2];
        }
    }

    TinyTIFFWriterFile* fit_tiff = TinyTIFFWriter_open("fit.tif", 32, TinyTIFFWriter_Float,
                                                       1, Nx, Ny,
                                                       TinyTIFFWriter_Greyscale);
    if(!fit_tiff) {
        std::cout << "Can't write fitting parameters into tiff!" << std::endl;
    }

    for(int param = 0; param < 3; param++)
        TinyTIFFWriter_writeImage(fit_tiff, &parameters[param*fft_size]);
    TinyTIFFWriter_close(fit_tiff);
}