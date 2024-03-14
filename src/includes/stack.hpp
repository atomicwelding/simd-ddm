#ifndef STACK_H
#define STACK_H

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "utils.hpp"

class Stack {
public:
    /*
     * Constructor. We create a fstream to the file at path location and
     * retrieve some infos from the custom header.
     *
     * N corresponds to the number of frames to load into the image buffer.
     */
    Stack(utils::Options& options);
    ~Stack();

    double mean_sampling_time() {
        return (times.back()-times.front()) / (times.size()-1);
    }

    float* images;

    uint64_t clock_frequency;
    uint16_t aoi_width;
    uint16_t aoi_height;

    size_t image_size;
    size_t len_images_buffer;

    int bin_factor;

    std::vector<double> times;

private:
    int current_byte();
    std::runtime_error error_reading(const std::string& msg);

    void load_next_M12P_frame(size_t offset);
    void load_M12P_images(size_t N);
	void normalize();

    void binning(int N);

    std::ifstream acq;

	uint32_t N_bytes_block;
	std::vector<char> block_buffer;

    uint16_t stride;
    uint8_t encoding;
};

#endif // STACK_H
