#ifndef SIMD_DDM_DELAYS
#define SIMD_DDM_DELAYS

#include "utils.hpp"
#include <vector>
#include <string>

template <typename T>
class Delays {
public:

    Delays(T mean_sampling_time, utils::Options options, std::string mode);
    void switchMode(std::string mode);
    void computeDelays();
    void computeDelays2();
    void save();

    const std::vector<T>& getTime() const {
        return time;
    }

    auto getIndex() const {
        return index;
    }

    auto getSamplingTime() const {
        return mean_sampling_time;
    }

private:
    const utils::Options options;
    T mean_sampling_time;

    std::vector<T> time;
    std::vector<int> index;

    std::string mode;
};



#endif //SIMD_DDM_DELAYS
