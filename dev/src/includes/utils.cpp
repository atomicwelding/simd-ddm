#include "utils.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

int utils::stoe(std::string &s) {
    int id = -1;

    if(s == "Mono16")
        id = Mono16;
    if(s == "Mono12")
        id = Mono12;
    if(s == "Mono12Packed")
        id = Mono12Packed;
    if(s == "Mono32")
        id = Mono32;

    return id;
}

std::vector<float> utils::logspace(double start, double end, int num_points) {
    std::vector<float> result;
    float step = (end - start) / (num_points - 1);

    for (int i = 0; i < num_points; i++) {
        float value = start + i * step;
        result.push_back(std::pow(10, value));
    }

    return result;
}

std::vector<float> utils::log_delays_in_time(float sampling_time, float delay_max, int N_delay) {
    std::vector<float> delays = utils::logspace(0, std::log10(delay_max / sampling_time), N_delay);

    auto it = std::unique(delays.begin(), delays.end());
    delays.resize(std::distance(delays.begin(), it));

    for(int i = 0; i < N_delay; i++)
        delays[i] = std::floor(delays[i])*sampling_time;

    return delays;
}

std::vector<int> utils::log_delays_indexes(float sampling_time, float delay_max, int N_delay) {
    std::vector<float> lspc = utils::logspace(0, std::log10(delay_max / sampling_time), N_delay);
    std::vector<int> delays(N_delay);
    std::transform(lspc.begin(), lspc.end(), delays.begin(), [](float x) -> int {
        return std::floor(x);
    });
    auto it = std::unique(delays.begin(), delays.end());
    delays.resize(std::distance(delays.begin(), it));

    return delays;
}
