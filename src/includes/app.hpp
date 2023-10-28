#ifndef APP_H
#define APP_H

#include <fftw3.h>
#include <string>

#include "stack.hpp"
#include "utils.hpp"



class App {
public:
    App(utils::Options &options);
    ~App();

    void run();

private:
    utils::Options* options;
};

#endif // APP_H
