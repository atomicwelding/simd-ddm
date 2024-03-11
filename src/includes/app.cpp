#include <vector>
#include <iostream>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm.hpp"
#include "fit.hpp"

App::App(utils::Options &options) : options(options) {}
App::~App()= default;

void App::run() {

	Timer timer;

    std::cout << "* Loading images..." << std::flush;
	timer.start();
    Stack stack(options);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    DDM DDMStack(stack, options);
    DDMStack.save();

    QuadraticSmoothingFit myfit(DDMStack, options);
    myfit.process();

    // à revoir
    /*if(this->options->doFit) {
        std::cout << "* Fitting..." << std::flush;

        timer.start();
        auto expToFit  = [](double tau, double A, double B, double f) -> double {
                    return A*(1-std::exp(-tau*f))+B;
        };

        // need to use it now in the fit routine
        auto ROI = fit::find_ROI(expToFit, ddm, ddm_width, ddm_height, *this->options, mean_sampling_time);

        if(this->options->doLogScale)
            fit::EXPERIMENTAL_fit_routine_log(expToFit, stack, ddm, ddm_width, ddm_height, delays_time, ROI);
        else
            fit::fit_routine(expToFit, stack, ddm, ddm_width, ddm_height, this->options->Ntau, fft_size);

        timer.stop();
        std::cout << "                           " << timer.elapsedSec() << "s" << std::endl;
    }*/
};


