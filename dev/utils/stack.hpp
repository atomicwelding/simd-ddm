#ifndef STACK_H
#define STACK_H

#include <fstream>
#include <stdexcept>
#include <string>

#include "frame.hpp"

#define Mono16 0
#define Mono12 1
#define Mono12Packed 2
#define Mono32 3


class Stack {
private:
    std::ifstream acq;

    int current_byte();
    std::runtime_error error_reading(const std::string& msg);

public:
    auto current_frame;
    void load_next_frame();

    uint16_t stride;
    uint8_t encoding;
    uint64_t clock_frequency;
    uint16_t aoi_width;
    uint16_t aoi_height;

    Stack(const std::string& path);
    ~Stack();
};

#endif // STACK_H
