#include "app.hpp"


#include <vector>
#include <iostream>


#include "stack.hpp"
#include "delays.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm.hpp"
#include "fit.hpp"
//#include "fitting.hpp"

App::App(utils::Options &options) : options(options) {}
App::~App()= default;

void App::run() {

    if(utils::stoe(options.encoding) != Mono12Packed)
        throw std::runtime_error("Encoding not supported yet");

	Timer timer;

    std::cout << "* Loading images..." << std::flush;
	timer.start();
    Stack stack(options);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;


    // we dont use camera's times, as we are using the mean sampling time
    // can we show that they are the same or not
    float mean_sampling_time = (stack.times.back() - stack.times.front())/(stack.times.size() - 1);

    std::string mode = options.doLogScale? "logarithmic" : "linear";
    Delays<float> delays(mean_sampling_time, options, mode);
    if( delays.getIndex().back() > options.loadNframes )
        throw std::runtime_error("Error : the file you want to process is shorter than desired max delay");
    delays.save();

    DDM<float> DDMStack(stack, delays, options);
    DDMStack.save();

    typedef std::function<float(float,float,float,float)> fitfn;
    fitfn expToFit = [](float tau, float A, float B, float f) -> float {
        return A*(1-std::exp(-tau*f))+B;
    };

    Fit2D<float, fitfn> myfit(DDMStack, options, expToFit);
    myfit.fit();

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


