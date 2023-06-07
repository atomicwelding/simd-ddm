#include "frame.hpp"

Frame::Frame(uint32_t im_bytes, uint32_t tk_bytes, char* im_data, uint64_t tk_data)
    : im_bytes(im_bytes), tk_bytes(tk_bytes), im_data(im_data), tk_data(tk_data) {}
