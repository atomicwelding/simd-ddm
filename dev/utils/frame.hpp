#ifndef FRAME_H
#define FRAME_H

#include <stdint.h>

class Frame {
public:
    uint32_t im_bytes;
    uint32_t tk_bytes;

    char* im_data;
    uint64_t tk_data;

    Frame(uint32_t im_bytes, uint32_t tk_bytes, char* im_data, uint64_t tk_data);
};


#endif // FRAME_H
