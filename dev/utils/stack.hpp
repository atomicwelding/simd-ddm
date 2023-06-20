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

    void load_next_M12P_frame(int offset);
    void load_M12P_images(int N);

public:
    void* N_images_buffer;


    uint16_t stride;
    uint8_t encoding;
    uint64_t clock_frequency;
    uint16_t aoi_width;
    uint16_t aoi_height;

    Stack(const std::string& path, int N);
    ~Stack();
};

#endif // STACK_H
