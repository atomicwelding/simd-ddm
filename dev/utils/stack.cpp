#include <iostream>
#include <fftw3.h>
#include "stack.hpp"

// === UTILS === //
int Stack::current_byte() {
    return (int) this->acq.tellg();
}

std::runtime_error Stack::error_reading(const std::string& msg) {
    return std::runtime_error("[" + std::to_string(this->current_byte()) + "] " + msg);
}


// === "MAIN" METHODS === //
template <typename T>
void Stack::load_next_frame() {
    /*
     * This method intend to load the next frame by reading the so-called frame block
     * and frame block. The new frame is stored in an object and a ref to this object is
     * maintained in a member variable named current_frame.
     *
     * We do sanity checks every time by checking if we're looking at the right CIDs for
     * each blocks.
     *
     * Could be great to see when to delete the previous Frame objects. We do not handle
     * it yet, cause we will need previous frames to make the calculations.
     *
     * ONLY HANDLES MONO12PACKED
     */

    uint32_t im_cid, im_bytes;
    uint32_t tk_cid, tk_bytes;

    // CID = 0 for frame block
    this->acq.read(reinterpret_cast<char*>(&im_cid), 4);
    if(im_cid != 0)
        throw this->error_reading("Frame block malformed: cannot read image");

    // read length of frame block + 4 bytes for CID
    this->acq.read(reinterpret_cast<char*>(&im_bytes), 4);
    im_bytes -= 4;

    // load image in buffer, strip the padding
    // only supports Mono12Packed
    int width_in_bytes = (this->aoi_width/2)*3;
    int padding_size = (this->stride - width_in_bytes);

    int image_size = this->aoi_width * this->aoi_height;
    T* image = (T*) fftw_malloc(image_size*sizeof(T));

    char buf[3];
    int count = 0;
    for(int i = 0; i < image_size; i += 2) {
        this->acq.read(reinterpret_cast<char*>(&buf), 3);

        image[i] = (buf[0] << 4) + (buf[1] & 0xF);
        image[i+1] = (buf[2] << 4) + (buf[1] >> 4);

        count += 3;
        if(count >= width_in_bytes) {
            this->acq.seekg(padding_size, std::ios_base::cur);
            count = 0;
        }
    }

    // seek after the remaining additional padding bytes at the end of the image
    int additional_empty_bytes = im_bytes - (width_in_bytes + padding_size) * this->aoi_height;
    this->acq.seekg(additional_empty_bytes, std::ios_base::cur);

    // CID = 1 for tick block
    this->acq.read(reinterpret_cast<char*>(&tk_cid), 4);
    if(tk_cid != 1)
        throw this->error_reading("Tick block malformed: cannot read tick infos");

    // read length of tick block + 4 bytes for CID
    this->acq.read(reinterpret_cast<char*>(&tk_bytes), 4);
    tk_bytes -= 4;

    // load ticks
    uint64_t tk_data;
    this->acq.read(reinterpret_cast<char*>(&tk_data), tk_bytes);

    // craft the new frame and keep the reference to it
    if(this->encoding == Mono12Packed)
        this->current_frame = new Frame<T>(im_bytes, tk_bytes, image, tk_data);

}

Stack::Stack(const std::string& path) {
    /*
     * Constructor. We create an fstream to the file at path location and
     * retrieve some infos from the custom header.
     */

    this->acq.open(path, std::ios::binary);

    uint32_t h_cid = 0, h_length = 0;
    this->acq.read(reinterpret_cast<char*>(&h_cid), 4);
    this->acq.read(reinterpret_cast<char*>(&h_length), 4);

    // sanity check
    if(h_cid != 2 || h_length != 19)
        throw this->error_reading("Wrong file format");

    // get stride
    this->acq.read(reinterpret_cast<char*>(&this->stride), 2);

    // get encoding (types of encoding are defined as #define in header file)
    this->acq.read(reinterpret_cast<char*>(&this->encoding), 1);
    if(this->encoding != Mono12Packed)
        throw this->error_reading("Pixel encoding not implemented yet");

    // get clock frequency
    this->acq.read(reinterpret_cast<char*>(&this->clock_frequency), 8);

    // get aoi infos ; AOIWidth is measured in pixels
    this->acq.read(reinterpret_cast<char*>(&this->aoi_width), 2);
    this->acq.read(reinterpret_cast<char*>(&this->aoi_height), 2);

    // load the first frame in memory !
    if(this->encoding == Mono12Packed)
        this->load_next_frame<float>();
}

Stack::~Stack() {
    if(this->encoding == Mono12Packed)
        delete (Frame<float>*) this->current_frame;
    this->acq.close();
}
