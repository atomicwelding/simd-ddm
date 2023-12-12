#include "delays.hpp"

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>

/*mode might be better if declared as enum val*/

template<typename T>
Delays<T>::Delays(T mean_sampling_time, utils::Options options, std::string mode)
          : mean_sampling_time(mean_sampling_time), options(options), mode(mode) {
       /**
        * Comment generated by GPT-3
        *
        * Construct a Delays object with specified parameters and compute delays.
        *
        * This constructor initializes a Delays object with the provided mean sampling time, options, and mode.
        * It performs the following steps:
        *   1. Initializes mean sampling time, options, and mode based on the provided parameters.
        *   2. Computes delays using the computeDelays() function.
        *
        * @tparam T Template parameter specifying the data type.
        * @param mean_sampling_time The mean sampling time for delay computation.
        * @param options Pointer to the Options structure containing computation parameters.
        * @param mode The mode for computing delays (linear or logarithmic).
        */
    computeDelays();
}


template<typename T>
void Delays<T>::computeDelays() {
    /**
     * Comment generated by GPT-3
     *
     * Compute delays based on the specified mode and parameters.
     *
     * This function computes delays based on the specified mode and parameters, such as linear or logarithmic spacing.
     * It performs the following steps:
     *   1. Computes delays using either linear or logarithmic spacing.
     *   2. Adjusts delays based on the specified mode.
     *   3. Removes duplicate delays and adjusts the number of time delays if needed.
     *   4. Prints a warning if the number of time delays differs from the requested number.
     *
     * @tparam T Template parameter specifying the data type.
     */
    if (mode == "linear")
        time = utils::linspace<T>(1, options.delayMax, options.Ntau);
    else if (mode == "logarithmic")
        time = utils::logspace<T>(0, std::log10(options.delayMax / mean_sampling_time), options.Ntau);
    else
        throw std::runtime_error("Error : delays(...) : mode doesn't exists");



    if(mode == "logarithmic")
        std::transform(time.begin(), time.end(), time.begin(), [this](T x) -> T {
            return (T) std::floor(x)*mean_sampling_time;
        });

    auto ittime = std::unique(time.begin(), time.end());
    time.resize(std::distance(time.begin(), ittime));


    // serve as initialization before transforming it
    index.resize(time.size());
    std::transform(time.begin(), time.end(), index.begin(), [this](T x) -> int {
        return std::floor(x/this->mean_sampling_time);
    });

    auto itidx =  std::unique(index.begin(), index.end());
    index.resize(std::distance(index.begin(), itidx));

    std::cout   << "/!\\ Number of time delays might be different from what have been asked.\n"
                << "Asked: " << options.Ntau << " Given: " << index.size() << std::endl;

}

template<typename T>
void Delays<T>::switchMode(std::string mode) {
    /**
     * Comment generated by GPT-3
     *
     * Switch the mode and recompute delays.
     *
     * This function switches the mode to the specified one and recomputes delays based on the updated mode.
     *
     * @tparam T Template parameter specifying the data type.
     * @param mode The new mode to switch to.
     */
    this->mode = mode;
    computeDelays();
}

template<typename T>
void Delays<T>::save() {
    /**
     * Comment generated by GPT-3
     *
     * Save the computed delays to a file.
     *
     * This function saves the computed delays to a file in CSV format. It includes the shift index and real time for each delay.
     *
     * @tparam T Template parameter specifying the data type.
     */
    std::ofstream delays_file;
    delays_file.open(options.pathOutput  + "_delays.dat");
    delays_file << "Shift" << "," << "Real_Time" << "\n";
    for(int i = 0; i < time.size(); i++)
        delays_file << index[i] << "," << time[i] << "\n";
    delays_file.close();
}

template class Delays<float>;
template class Delays<double>;

