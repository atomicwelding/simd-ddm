#include "utils.hpp"

int utils::stoe(std::string& s) {
    int id = -1;

    if(s == "Mono16")
        id = Mono16;
    if(s == "Mono12")
        id = Mono12;
    if(s == "Mono12Packed")
        id = Mono12Packed;
    if(s == "Mono32")
        id = Mono32;

    return id;
};