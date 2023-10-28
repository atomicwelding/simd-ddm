#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
public:
	Timer() {}

    void start();
    void stop();
    
    double elapsedMilliSec();
    double elapsedSec();

private:
    std::chrono::time_point<std::chrono::system_clock> start_time;
    std::chrono::time_point<std::chrono::system_clock> end_time;
    bool                                               running = false;
};

 #endif
