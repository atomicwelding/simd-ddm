#include "delays.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
// changer les erreurs de runtime à la compilation si possible
// pouvoir switch sur les deux représentations ? (log et lin)
// à voir, pas essentiel pour le moment

template<typename T>
Delays<T>::Delays(T mean_sampling_time, utils::Options *options, std::string mode)
          : mean_sampling_time(mean_sampling_time), options(options), mode(mode) {
    this->oldNtau = options->Ntau;
    computeDelays();
}


template<typename T>
void Delays<T>::computeDelays() {
    options->Ntau = oldNtau;

    if (mode == "linear")
        time = utils::linspace<T>(1, options->delayMax, options->Ntau);
    else if (mode == "logarithmic")
        time = utils::logspace<T>(0, std::log10(options->delayMax / mean_sampling_time), options->Ntau);
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
                << "Asked: " << oldNtau << " Given: " << index.size() << std::endl;
    options->Ntau = index.size();
}

template<typename T>
void Delays<T>::switchMode(std::string mode) {
    this->mode = mode;
    computeDelays();
}

template<typename T>
void Delays<T>::save() {
    std::ofstream delays_file;
    delays_file.open(options->pathOutput  + "_delays.dat");
    delays_file << "Shift" << "," << "Real_Time" << "\n";
    for(int i = 0; i < time.size(); i++)
        delays_file << index[i] << "," << time[i] << "\n";
    delays_file.close();
}

template class Delays<float>;
template class Delays<double>;

