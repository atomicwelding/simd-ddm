#ifndef UTILS_H
#define UTILS_H

#include <string>

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
        int loadNframes;
        std::string encoding;
        bool doNormalize;
        int tauMax;
        std::string pathOutput;
        bool doFit;
        int binFactor;
    };

}

#endif // UTILS_H
