#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <type_traits>

#include <tiffio.h>
#include <fftw3.h>

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
        std::string output_path;
        std::string path_movie;
        std::string ddm_algo;

        int bin_factor;
        int N_frames;

        int N_lags;
        float max_lag_time;
        bool log_lags;

        bool fit;
        float max_decay_freq;
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

    /**
     * Helper method to allocate an array with N aligned complex floats.
     * 
     * We use _mm_malloc if AVX512 is supported, else the routine from FFTW.
    */
    fftwf_complex* allocate_complex_float_array(size_t N);

    /**
     * Helper method to allocate an array with N aligned floats.
     * 
     * We use _mm_malloc if AVX512 is supported, else the routine from FFTW.
    */
    float* allocate_float_array(size_t N);

    /**
     * LibTiff wrapper allowing to efficiently append an image in a multiframe tiff file.
    */
	template <typename ImType>
	bool libTIFFWriter_writeImage(TIFF* tif, ImType* img, int width, int height);
}

#endif // UTILS_H