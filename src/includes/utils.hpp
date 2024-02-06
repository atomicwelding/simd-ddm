#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>

// === ENCODING ===
#define Mono16 0
#define Mono12 1
#define Mono12Packed 2
#define Mono32 3

#define REAL 0
#define IMAG 1

namespace utils {
    int stoe(std::string &s);

	template <typename T>
	T sqr(T val) {
		return val*val;
	}

    struct Options {
        std::string path;
        std::string encoding;
        std::string pathOutput;

        int loadNframes;
        int Ntau;
        int binFactor;

        float delayMax;
        float frequencyThreshold;

        bool doNormalize;
        bool doFit;
        bool doLogScale;

    };

    template<typename InputIterator, typename ValueType>
    int closest_index(InputIterator first, InputIterator last, ValueType value)
    {
        auto closest_it = std::min_element(first, last, [&](ValueType x, ValueType y)
        {
            return std::abs(x - value) < std::abs(y - value);
        });

        // Calculate the index of the closest element
        return static_cast<int>(std::distance(first, closest_it));
    }

    template<typename T>
    std::vector<T> logspace(double start, double end, int num_points) {
        std::vector<T> result;
        float step = (end - start) / (num_points - 1);

        for (int i = 0; i < num_points; i++) {
            float value = start + i * step;
            result.push_back(std::pow(10, value));
        }

        return result;
    }

    std::vector<int> log_delays_indexes(float sampling_time, float delay_max, int N_delay);

    template<typename T>
    std::vector<T> log_delays_in_time(float sampling_time, float delay_max, int N_delay) {
        std::vector<T> delays = utils::logspace<T>(0, std::log10(delay_max / sampling_time), N_delay);

        auto it = std::unique(delays.begin(), delays.end());
        delays.resize(std::distance(delays.begin(), it));

        for(int i = 0; i < N_delay; i++)
            delays[i] = std::floor(delays[i])*sampling_time;

        return delays;
    }
}

#endif // UTILS_H
