#include "utils.hpp"

#include <string>
#include <vector>
#include <cmath>

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


std::vector<int> utils::log_delays_indexes(float sampling_time, float delay_max, int N_delay) {
    std::vector<float> lspc = utils::logspace<float>(0, std::log10(delay_max / sampling_time), N_delay);
    std::vector<int> delays(N_delay);
    std::transform(lspc.begin(), lspc.end(), delays.begin(), [](float x) -> int {
        return std::floor(x);
    });
    auto it = std::unique(delays.begin(), delays.end());
    delays.resize(std::distance(delays.begin(), it));

    return delays;
}
