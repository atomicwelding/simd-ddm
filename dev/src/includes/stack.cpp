#include <iostream>
#include <vector>
#include <fftw3.h>

#include "stack.hpp"
#include "utils.hpp"

// === CONSTRUCTORS === //
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

	// We peek the length of the frame+tick block before reading all images in order to read
	// the file by big chunks

    // CID = 0 for frame block
    uint32_t im_cid;
    this->acq.read(reinterpret_cast<char*>(&im_cid), 4);
    if(im_cid != 0)
        throw this->error_reading("Frame block malformed: cannot read image");

    // read length of frame block + 4 bytes for CID, and deduce the whole length of a
	// frame+tick block, including CID, length and data subblocks.
    this->acq.read(reinterpret_cast<char*>(&N_bytes_block), 4);
    N_bytes_block += 20;
	block_buffer.resize(N_bytes_block);

	// Go back to the start of the first image block, and start reading images by chunk
	this->acq.seekg(-8, std::ios_base::cur);
    this->images = reinterpret_cast<T*>(fftw_malloc(this->len_images_buffer*sizeof(T)));
    if(this->encoding == Mono12Packed) {
        this->load_M12P_images(N);
    }

    if(do_normalize) {
        this->acc_signal /= this->len_images_buffer;
        for(int i = 0; i < len_images_buffer; ++i)
            this->images[i] /= acc_signal;
    }
}

template<typename T>
Stack<T>::~Stack() {
    fftw_free(this->images);
    this->acq.close();
}


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
	this->acq.read(block_buffer.data(), N_bytes_block);

    // CID = 0 for frame block and CID = 1 for tick block
    uint32_t im_cid, im_bytes, tk_cid, tk_bytes;
	for(int i=0; i<4; i++) {
		reinterpret_cast<char*>(&im_cid)[i] = block_buffer[i];
		reinterpret_cast<char*>(&im_bytes)[i] = block_buffer[4+i];
		reinterpret_cast<char*>(&tk_cid)[i] = block_buffer[N_bytes_block-16+i];
		reinterpret_cast<char*>(&tk_bytes)[i] = block_buffer[N_bytes_block-12+i];
	}

    if(im_cid!=0 || im_bytes!=N_bytes_block-20)
        throw this->error_reading("Frame block malformed: cannot read image");
    if(tk_cid!=1 || tk_bytes!=12)
        throw this->error_reading("Tick block malformed: cannot read tick infos");

	// shortcut ptrs for each row of the image
	char* buf_row;
	T* im_row;

	int iy, ix, ib;
	for(iy=0; iy<aoi_height; iy++) {
		buf_row = &block_buffer[8+iy*stride];
		im_row = &this->images[offset+iy*this->aoi_width];
		for(ix=0, ib=0; ix<aoi_width; ix+=2, ib+=3) {
			im_row[ix] = (buf_row[ib] << 4) + (buf_row[ib+1] & 0xF);
			im_row[ix+1] = (buf_row[ib+2] << 4) + (buf_row[ib+1] >> 4);

            this->acc_signal += im_row[ix] + im_row[ix + 1];
		}
    }

}

// === EXPLICITLY INSTANTIATING THE TEMPLATE === //
template class Stack<float>;
template class Stack<double>;
