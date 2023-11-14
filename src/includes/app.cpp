#include <vector>
#include <iostream>


#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "ddm.hpp"
#include "fitting.hpp"

App::App(utils::Options& options) : options(&options) {}
App::~App()= default;

void App::run() {

    if(utils::stoe(this->options->encoding) != Mono12Packed)
        throw std::runtime_error("Encoding not supported yet");

	Timer timer;

    std::cout << "* Loading images..." << std::flush;
	timer.start();
    auto* stack = new Stack(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->doNormalize,
                                           this->options->binFactor);
	std::cout << "                     " << timer.elapsedSec() << "s" << std::endl;

    // not the perfect place to do that
    int Nt = stack->times.size();
    double mean_sampling_time = (stack->times[Nt-1] - stack->times[0])/(Nt-1);
    std::vector<int> delays_filtered;

    std::vector<double> delays_time;
    std::vector<int> delays;
    if(this->options->doLogScale) {
        delays_time = utils::log_delays_in_time<double>(mean_sampling_time, this->options->delayMax, this->options->Ntau);
        delays = utils::log_delays_indexes(mean_sampling_time, this->options->delayMax, this->options->Ntau);


        if(delays_time[delays_time.size()-1] >= this->options->delayMax || delays[delays.size() - 1] >= this->options->loadNframes) {
            throw std::runtime_error("Error : the file you want to process is shorter than the max delay you're asking");
        }

        std::copy_if(delays.begin(), delays.end(), std::back_inserter(delays_filtered),
                     [&](int x) { return x < this->options->loadNframes; });
        std::cout << "/!\\ Ntau has changed, from " << this->options->Ntau << std::flush;
        this->options->Ntau = delays_filtered.size();
        std::cout << " to  " << this->options->Ntau << " due to log scaling !" << std::endl;


        std::ofstream delays_array_file;
        delays_array_file.open(this->options->pathOutput  + "_time_delays.dat");
        delays_array_file << "Shift" << "," << "Real_Time" << "\n";
        for(int i = 0; i < this->options->Ntau; i++) {
            delays_array_file << delays[i] << "," << delays_time[i] << "\n";
        }

        delays_array_file.close();
    }

    auto* DDMStack = new DDM(stack, *(this->options));
    DDMStack->save_ddm_buffer(this->options->pathOutput);


    // à revoir
    /*if(this->options->doFit) {
        std::cout << "* Fitting..." << std::flush;

        timer.start();
        auto exp_to_fit  = [](double tau, double A, double B, double f) -> double {
                    return A*(1-std::exp(-tau*f))+B;
        };

        // need to use it now in the fit routine
        auto ROI = fit::find_ROI(exp_to_fit, ddm, ddm_width, ddm_height, *this->options, mean_sampling_time);

        if(this->options->doLogScale)
            fit::EXPERIMENTAL_fit_routine_log(exp_to_fit, stack, ddm, ddm_width, ddm_height, delays_time, ROI);
        else
            fit::fit_routine(exp_to_fit, stack, ddm, ddm_width, ddm_height, this->options->Ntau, fft_size);

        timer.stop();
        std::cout << "                           " << timer.elapsedSec() << "s" << std::endl;
    }*/


    std::cout << "* Cleaning ..." << std::endl;
    delete stack;
};


