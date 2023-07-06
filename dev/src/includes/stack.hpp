#ifndef STACK_H
#define STACK_H

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

template<typename T>
class Stack {
private:
    std::ifstream acq;

    int current_byte();
    std::runtime_error error_reading(const std::string& msg);

    void load_next_M12P_frame(int offset);
    void load_M12P_images(int N);
    void normalize();

public:
    T* N_images_buffer;

    uint16_t stride;
    uint8_t encoding;
    uint64_t clock_frequency;
    uint16_t aoi_width;
    uint16_t aoi_height;

    int image_size;
    int len_images_buffer;

    Stack(const std::string& path, int encoding, int N, bool do_normalize);
    ~Stack();
};

#endif // STACK_H
