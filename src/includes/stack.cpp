#include <iostream>
#include <fftw3.h>
#include <cmath>

#include <vector>
#include <tinytiffwriter.h>

#include "stack.hpp"
#include "utils.hpp"

// === CONSTRUCTORS === //
Stack::Stack(const std::string& path, int encoding, int N, bool do_normalize, int bin_factor) {
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
    this->images = reinterpret_cast<float*>(fftw_malloc(this->len_images_buffer*sizeof(float)));
    if(this->encoding == Mono12Packed)
        this->load_M12P_images(N);

    if(do_normalize)
		normalize();

    this->bin_factor = bin_factor;
    if(this->bin_factor != 1)
        this->binning(N);

}

Stack::~Stack() {
    fftwf_free(this->images);
    this->acq.close();
}


// === UTILS === //
int Stack::current_byte() {
    return (int) this->acq.tellg();
}

std::runtime_error Stack::error_reading(const std::string& msg) {
    return std::runtime_error("[" + std::to_string(this->current_byte()) + "] " + msg);
}

// === ENCODING SPECIFIC === //
void Stack::load_M12P_images(int N) {
    for(int i = 0; i < N; ++i) {
        this->load_next_M12P_frame(i * this->image_size);
    }
}

void Stack::load_next_M12P_frame(int offset) {
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

    uint64_t tick_count;
    for(int i = 0; i < 8; i++)
        reinterpret_cast<char*>(&tick_count)[i] = block_buffer[N_bytes_block-8+i];
    this->times.push_back(tick_count*1./this->clock_frequency);

    if(im_cid!=0 || im_bytes!=N_bytes_block-20)
        throw this->error_reading("Frame block malformed: cannot read image");
    if(tk_cid!=1 || tk_bytes!=12)
        throw this->error_reading("Tick block malformed: cannot read tick infos");

	// shortcut ptrs for each row of the image
	char* buf_row;
	float* im_row;

	int iy, ix, ib;
	for(iy=0; iy<aoi_height; iy++) {
		buf_row = &block_buffer[8+iy*stride];
		im_row = &this->images[offset+iy*this->aoi_width];
		for(ix=0, ib=0; ix<aoi_width; ix+=2, ib+=3) {
			im_row[ix] = (buf_row[ib] << 4) + (buf_row[ib+1] & 0xF);
			im_row[ix+1] = (buf_row[ib+2] << 4) + (buf_row[ib+1] >> 4);
		}
    }
}


void Stack::normalize() {
    /*
     * Normalize the signal in the buffer
     * by dividing it by the average value of pixels
     */
    float mean = 0.0;
	#pragma omp parallel for reduction(+:mean)
	for(int i = 0; i < this->len_images_buffer; ++i)
		mean += this->images[i];
	mean /= this->len_images_buffer;

	for(int i = 0; i < len_images_buffer; ++i)
		this->images[i] /= mean;
}


void Stack::binning(int N) {
    /*
     * This implementation is very naive. It assumes :
     *      * Square images (width=height)
     *      * That the image size is a multiple of the binning factor
     * This has also dirty side effects.
     */

    if(this->aoi_width != this->aoi_height)
        throw std::runtime_error("Use only square images with binning");

    if(this->aoi_width%this->bin_factor != 0 || this->aoi_height%this->bin_factor != 0)
        throw std::runtime_error("Do not use binning with this image size (height or width isn't a multiple of the bin_factor");


    int bin_dim = (int) std::floor(this->aoi_width/this->bin_factor);
    int bin_image_size = bin_dim*bin_dim;
    int bin_factor_sqr = this->bin_factor * this->bin_factor;

    // creates a new stack
    float* bin_stack = reinterpret_cast<float*>(fftw_malloc(N * bin_image_size * sizeof(float)));

    // binning
    for(int n = 0; n < N; n++) {

        float* img_ptr = &this->images[n*this->image_size];
        for(int y = 0; y < bin_dim; y++) {
            for(int x = 0; x < bin_dim; x++) {

                int sum = 0;
                for(int by = 0; by < this->bin_factor; by++) {
                    for(int bx = 0; bx < this->bin_factor; bx++)
                        sum += img_ptr[(y * this->bin_factor + by)*this->aoi_width + (x * this->bin_factor + bx)];
                }


                bin_stack[n * bin_dim*bin_dim + y*bin_dim + x] = static_cast<float>(sum / bin_factor_sqr);
            }
        }
    }

    // free the previous stack
    fftwf_free(this->images);

    // set by side effects
    this->image_size = bin_dim*bin_dim;
    this->aoi_width = this->aoi_height = bin_dim;
    this->images = bin_stack;

/* DEBUG to see what binning does.
 * TinyTIFFWriterFile* btif = TinyTIFFWriter_open("bin.tif", 32, TinyTIFFWriter_Float,
                                                       1, bin_dim, bin_dim,
                                                       TinyTIFFWriter_Greyscale);

    if(!btif) std::cout << "error cant write bin.tif" << std::endl;

    TinyTIFFWriter_writeImage(btif, bin_stack);
    TinyTIFFWriter_close(btif);*/
}
