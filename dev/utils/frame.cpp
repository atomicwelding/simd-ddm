#include "frame.hpp"

template<class T>
Frame<T>::Frame(uint32_t im_bytes, uint32_t tk_bytes, T* im_data, uint64_t tk_data)
    : im_bytes(im_bytes), tk_bytes(tk_bytes), im_data(im_data), tk_data(tk_data) {}

template class Frame<float>;
