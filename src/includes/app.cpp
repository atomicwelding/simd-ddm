#include <vector>
#include <iostream>
#include <memory>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm_diff.hpp"
#include "ddm_wk.hpp"
#include "fit.hpp"

App::App(utils::Options &options) : options(options) {}
App::~App()= default;

void App::run() {

	Timer timer;

    std::cout << "* Loading images...               " << std::flush;
	timer.start();
    Stack stack(options);
	std::cout << timer.elapsedSec() << "s" << std::endl;

    std::shared_ptr<DDM> ddm;
    if(options.ddm_algo=="diff")
        ddm = std::make_shared<DDMDiff>(stack, options);
    else if(options.ddm_algo=="wk")
        ddm = std::make_shared<DDMWK>(stack, options);
    ddm->save();

    QuadraticSmoothingFit myfit(*ddm, options);
    myfit.process();

    // à revoir
    /*if(this->options->fit) {
        std::cout << "* Fitting..." << std::flush;

        timer.start();
        auto expToFit  = [](double tau, double A, double B, double f) -> double {
                    return A*(1-std::exp(-tau*f))+B;
        };

        // need to use it now in the fit routine
        auto ROI = fit::find_ROI(expToFit, ddm, ddm_width, ddm_height, *this->options, mean_sampling_time);

        if(this->options->log_lags)
            fit::EXPERIMENTAL_fit_routine_log(expToFit, stack, ddm, ddm_width, ddm_height, delays_time, ROI);
        else
            fit::fit_routine(expToFit, stack, ddm, ddm_width, ddm_height, this->options->N_lags, fft_size);

        timer.stop();
        std::cout << "                           " << timer.elapsedSec() << "s" << std::endl;
    }*/
};


