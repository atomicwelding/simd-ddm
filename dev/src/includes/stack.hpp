#ifndef STACK_H
#define STACK_H

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

template<typename T>
class Stack {
public:
    Stack(const std::string& path, int encoding, int N, bool do_normalize);
    ~Stack();

    T* images;

    uint16_t stride;
    uint8_t encoding;
    uint64_t clock_frequency;
    uint16_t aoi_width;
    uint16_t aoi_height;

    int image_size;
    int len_images_buffer;

private:
    std::ifstream acq;

	uint32_t N_bytes_block;
	std::vector<char> block_buffer;

    int current_byte();
    std::runtime_error error_reading(const std::string& msg);

    void load_next_M12P_frame(int offset);
    void load_M12P_images(int N);
    void normalize();

};

#endif // STACK_H
