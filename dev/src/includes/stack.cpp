#include <iostream>
#include <fftw3.h>

#include "stack.hpp"
#include "utils.hpp"

// === UTILS === //
template<typename T>
int Stack<T>::current_byte() {
    return (int) this->acq.tellg();
}

template<typename T>
std::runtime_error Stack<T>::error_reading(const std::string& msg) {
    return std::runtime_error("[" + std::to_string(this->current_byte()) + "] " + msg);
}

// === ENCODING SPECIFIC === // 
template<typename T>
void Stack<T>::load_M12P_images(int N) {
    for(int i = 0; i < N; ++i) {
        this->load_next_M12P_frame(i * this->image_size);
    }
}

template<typename T>
void Stack<T>::load_next_M12P_frame(int offset) {
    /*
     * Load next image into the buffer ;
     * Performs sanity check to verify we're reading the right bytes
     */
    uint32_t im_cid, im_bytes, tk_cid;

    // CID = 0 for frame block
    this->acq.read(reinterpret_cast<char*>(&im_cid), 4);
    if(im_cid != 0)
        throw this->error_reading("Frame block malformed: cannot read image");

    // read length of frame block + 4 bytes for CID
    this->acq.read(reinterpret_cast<char*>(&im_bytes), 4);
    im_bytes -= 4;

    // load image in buffer, strip the padding
    int width_in_bytes = (this->aoi_width/2)*3;
    int padding_size = (this->stride - width_in_bytes);

    int image_size = this->aoi_width * this->aoi_height;

    char buf[3];
    for(int i = 0, count = 3; i < image_size; i += 2, count += 3) {
        this->acq.read(buf, 3);

        this->N_images_buffer[offset + i] = (buf[0] << 4) + (buf[1] & 0xF);
        this->N_images_buffer[offset + i+1] = (buf[2] << 4) + (buf[1] >> 4);

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

    // skip tick infos
    this->acq.seekg(12, std::ios_base::cur);
}

// === GENERIC METHODS === //
template<typename T>
void Stack<T>::normalize() {
    /*
     * Normalize the signal in the buffer
     * by dividing it by the average value of pixels
     */
    float mean = 0.0;
    for(int i = 0; i < this->len_images_buffer; ++i)
        mean += this->N_images_buffer[i];
    mean /= this->len_images_buffer;

    for(int i = 0; i < len_images_buffer; ++i)
        this->N_images_buffer[i] /= mean;
}

template<typename T>
Stack<T>::Stack(const std::string& path, int encoding, int N, bool do_normalize) {
    /*
     * Constructor. We create an fstream to the file at path location and
     * retrieve some infos from the custom header.
     *
     * N corresponds to the number of frames to load into the image buffer.
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

    // checks
    if(this->encoding != encoding)
        throw this->error_reading("Wrong encoding. Please, check out the format of the signal.");
    if(this->encoding != Mono12Packed)
        throw this->error_reading("Pixel encoding not implemented yet");
    // get clock frequency

    this->acq.read(reinterpret_cast<char*>(&this->clock_frequency), 8);

    // get aoi infos ; AOIWidth is measured in pixels
    this->acq.read(reinterpret_cast<char*>(&this->aoi_width), 2);
    this->acq.read(reinterpret_cast<char*>(&this->aoi_height), 2);

    this->image_size = this->aoi_width * this->aoi_height;
    this->len_images_buffer = this->image_size * N;

    this->N_images_buffer = reinterpret_cast<T*>(fftw_malloc(this->len_images_buffer*sizeof(T)));
    if(this->encoding == Mono12Packed) {
            this->load_M12P_images(N);
    }

    if(do_normalize)
        this->normalize();
}

template<typename T>
Stack<T>::~Stack() {
    fftw_free(this->N_images_buffer);
    this->acq.close();
}


// === EXPLICITLY INSTANTIATING THE TEMPLATE === //
template class Stack<float>;
template class Stack<double>;
