#ifndef UTILS_H
#define UTILS_H

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

    struct Options {
        std::string path;
        std::string encoding;
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

    template<typename T>
    struct lspace {
        std::vector<T> time;
        std::vector<int> index;
    };

    template<typename InputIterator, typename ValueType>
    int closest_index(InputIterator first, InputIterator last, ValueType value);

    template<typename T>
    std::vector<T> linspace(T start, T end, int num_points);

    template<typename T>
    std::vector<T> logspace(T start, T end, int num_points);

    template<typename T>
    lspace<T>* delays(T sampling_time, Options& opt, std::string mode = "linear");
}

#endif // UTILS_H
