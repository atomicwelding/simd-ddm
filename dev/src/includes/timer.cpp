#include "timer.hpp"


void Timer::start() {
	start_time = std::chrono::system_clock::now();
	running = true;
}

void Timer::stop() {
	end_time = std::chrono::system_clock::now();
	running = false;
}

double Timer::elapsedMilliSec()
{
	if(running) {
		end_time = std::chrono::system_clock::now();
		running = false;
	}

	return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
}

double Timer::elapsedSec()
{
	return elapsedMilliSec() / 1000.0;
}
