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


// template <typename T>
// void Stack::load_N_images_Mono12Packed(int N) {
//     uint32_t im_cid, im_bytes; // j'en fous quoi de ces trucs ?
//     uint32_t tk_cid, tk_bytes;

//         // let's load the images
//     int width_in_bytes = (this->aoi_width/2)*3;
//     int padding_size = (this->stride - width_in_bytes);

//     int image_size = this->aoi_width * this->aoi_height;
//     T* images_stack = (T*) fftw_malloc(image_size * sizeof(T) * N);

//     char buf[3];
//     int idx_im;
//     for(int n = 0; n < N; n++) {
//         // performs some sanity checks
//         this->acq.read(reinterpret_cast<char*>(&im_cid), 4);
//         if(im_cid != 0) throw this->error_reading("Frame block malformed : cannot read image");

//         // read typical length of an image
//         this->acq.read(reinterpret_cast<char*>(&im_bytes), 4);
//         im_bytes -= 4;

//         idx_im = n * image_size;
//         for(int i = 0, count = 3; i < image_size; i+=2, count +=3) {
//             this->acq.read(reinterpret_cast<char*>(&buf), 3);

//             images_stack[idx_im + i] = (buf[0] << 4) + (buf[1] & 0xF);
//             images_stack[idx_im + i+1] = (buf[2] << 4) + (buf[1] >> 4);

//         }
//     }
// }

void Stack::load_M12P_images(int N) {
    // not implemented yet
}

void Stack::load_next_M12P_frame(int offset) {
    // not implemented yet
    // load next frame into buffer
}

// === "MAIN" METHODS === //
// faut que je fasse une methode qui load_next_frame N fois ;
// simplement que je modifie la taille du buffer (et que load_next_frame prenne en argument un offset)
// et que je charge le nombre de frames (et donc alloue le buffer) directement dans le constructeur
// apres je vire l'objet frame
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
    int width_in_bytes = (this->aoi_width/2)*3;
    int padding_size = (this->stride - width_in_bytes);

    int image_size = this->aoi_width * this->aoi_height;
    T* image = (T*) fftw_malloc(image_size*sizeof(T));

    if(this->encoding == Mono12Packed) {
        char buf[3];
        for(int i = 0, count = 3; i < image_size; i += 2, count += 3) {
            this->acq.read(reinterpret_cast<char*>(&buf), 3);

            image[i] = (buf[0] << 4) + (buf[1] & 0xF);
            image[i+1] = (buf[2] << 4) + (buf[1] >> 4);

            if(count >= width_in_bytes) {
                this->acq.seekg(padding_size, std::ios_base::cur);
                count = 0;
            }
        }

        // seek after the remaining additional padding bytes at the end of the image
        int additional_empty_bytes = im_bytes - (width_in_bytes + padding_size) * this->aoi_height;
        this->acq.seekg(additional_empty_bytes, std::ios_base::cur);
    }

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

Stack::Stack(const std::string& path, int N) {
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
    if(this->encoding != Mono12Packed)
        throw this->error_reading("Pixel encoding not implemented yet");

    // get clock frequency
    this->acq.read(reinterpret_cast<char*>(&this->clock_frequency), 8);

    // get aoi infos ; AOIWidth is measured in pixels
    this->acq.read(reinterpret_cast<char*>(&this->aoi_width), 2);
    this->acq.read(reinterpret_cast<char*>(&this->aoi_height), 2);

    if(this->encoding == Mono12Packed)
        this->load_M12P_images(N);
}

Stack::~Stack() {
    if(this->encoding == Mono12Packed)
        fftw_free(this->N_images_buffer);
    this->acq.close();
}
