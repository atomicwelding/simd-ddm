#ifndef SIMD_DDM_DELAYS
#define SIMD_DDM_DELAYS

#include "utils.hpp"
#include <vector>
#include <string>

template <typename T>
class Delays {
public:
    Delays(T mean_sampling_time, utils::Options *options, std::string mode = "linear");

    int oldNtau;

    void switchMode(std::string mode);
    void computeDelays();
    void save();

    auto getTime() {
        return time;
    }

    auto getIndex() {
        return index;
    }

private:
    utils::Options *options;
    T mean_sampling_time;

    std::vector<T> time;
    std::vector<int> index;

    std::string mode;
};



#endif //SIMD_DDM_DELAYS
