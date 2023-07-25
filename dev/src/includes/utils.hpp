#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

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

        bool doNormalize;
        bool doFit;
        bool doLogScale;

    };

    std::vector<float> logspace(double start, double end, int num_points);

    std::vector<float> log_delays_in_time(float sampling_time, float delay_max, int N_delay);
    std::vector<int> log_delays_indexes(float sampling_time, float delay_max, int N_delay);
}

#endif // UTILS_H
