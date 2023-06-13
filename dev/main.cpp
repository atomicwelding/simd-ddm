#include <iostream>
#include <stdexcept>

#include "app.hpp"

int main() {
    try {
        App* app = new App("sample.dat");
        app->run();
        delete app;
    } catch(std::exception& ex) {
        std::cout << "Exception caught: " << ex.what() << std::endl;
        return -1;
    }

    return 0;
}
