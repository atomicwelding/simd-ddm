#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <tiffio.h>

// === ENCODING ===
#define Mono16 0
#define Mono12 1
#define Mono12Packed 2
#define Mono32 3

#define REAL 0
#define IMAG 1

namespace utils {

    int stoe(std::string& s);

    template <typename T>
    T sqr(T val) {
        return val*val;
    }
    struct Options {
        std::string path;
        std::string pathOutput;

        // change name of Ntau

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
    };

    template<typename T>
    std::vector<T> linspace(double start, double end, int num_points) {
        std::vector<T> result(num_points);
        double step = (end-start)/(num_points-1);
        double value;
        for(int i = 0; i < num_points; i++) {
            value = start + i*step;
            if constexpr(std::is_integral_v<T>)
                result[i] = std::round(value);
            else
                result[i] = value;
        }
        return result;
    };

    template<typename T>
    std::vector<T> logspace(double start, double end, int num_points) {
        std::vector<T> result(num_points);
        double pow_base = end/start;
        double pow_step = 1. / (num_points - 1);
        double value;
        for (int i = 0; i < num_points; i++) {
            value = start*std::pow(pow_base, i*pow_step);
            if constexpr(std::is_integral_v<T>)
                result[i] = std::round(value);
            else
                result[i] = value;
        }
        return result;
    };

	template <typename ImType>
	bool libTIFFWriter_writeImage(TIFF* tif, ImType* img, int width, int height);
}

#endif // UTILS_H
