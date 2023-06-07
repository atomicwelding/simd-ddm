#include "stack.hpp"


// === UTILS === //
std::string Stack::current_byte() {
    return std::to_string(this->acq.tellg());
}

std::runtime_error Stack::error_reading(const std::string& msg) {
    return std::runtime_error("[" + this->current_byte() + "] " + msg);
}


// === "MAIN" METHODS === //
void Stack::load_next_frame() {
    /*
     * This method intend to load the next frame by reading the so-called frame block
     * and frame block. The new frame is stored in a object and a ref to this object is
     * maintained in a member variable named current_frame.
     *
     * We do sanity checks every time by checking if we're looking at the right CIDs for
     * each blocks.
     *
     * Could be great to see when to delete the previous Frame objects. We do not handle
     * it yet, cause we will need previous frames to make the calculations.
     */

    uint32_t im_cid, im_bytes;
    uint32_t tk_cid, tk_bytes;

    // CID = 0 for frame block
    this->acq.read(reinterpret_cast<char*>(&im_cid), 4);
    if(im_cid != 0)
        throw this->error_reading("Frame block malformed: cannot read image");

    // read length of frame block and retrieve the 4 bytes for CID
    this->acq.read(reinterpret_cast<char*>(&im_bytes), 4);
    im_bytes -= 4;

    // load image in a buffer
    char im_data[im_bytes];
    this->acq.read(im_data, im_bytes);

    // CID = 1 for tick block
    this->acq.read(reinterpret_cast<char*>(&tk_cid), 4);
    if(tk_cid != 1)
        throw this->error_reading("Tick block malformed: cannot read tick infos");

    // read length of tick block and retrieve the 4 bytes for CID
    this->acq.read(reinterpret_cast<char*>(&tk_bytes), 4);
    tk_bytes -= 4;

    // load ticks
    uint64_t tk_data;
    this->acq.read(reinterpret_cast<char*>(&tk_data), tk_bytes);

    // craft the new frame and keep the reference to it
    this->current_frame = new Frame(im_bytes, tk_bytes, im_data, tk_data);
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

    // get aoi infos
    this->acq.read(reinterpret_cast<char*>(&this->aoi_width), 2);
    this->acq.read(reinterpret_cast<char*>(&this->aoi_height), 2);

    // load the first frame in memory !
    this->load_next_frame();
}

Stack::~Stack() {
    delete this->current_frame;
    this->acq.close();
}
