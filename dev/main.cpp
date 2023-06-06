#include <iostream>
#include <stdexcept>

#include "utils/stack.hpp"

int main() {
    try {
        Stack* stack = new Stack("sample.dat");
        delete stack;
    } catch(std::exception& ex) {
        std::cout << "Exception caught: " << ex.what() << std::endl;
        return -1;
    }

    return 0;
}
