#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>

#include "app.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv) {

    std::string path;
    int loadNframes;
    std::string encoding;
    bool do_normalize;

    po::options_description desc("Options");
    desc.add_options()
        ("path,p", po::value<std::string>(&path)->default_value("sample.dat"), "Path to the file to process")
        ("loadNframes,N", po::value<int>(&loadNframes)->default_value(1), "Number of frames to be loaded")
        ("encoding,e", po::value<std::string>(&encoding)->default_value("Mono12Packed"), "Encoding type")
        ("normalize, n", po::bool_switch(&do_normalize), "Normalize signal data");

    // parse
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // recap options
    std::cout << "File: " << path << std::endl;
    std::cout << "Encoding: " << encoding << std::endl;
    std::cout << "Frames to load: " << std::to_string(loadNframes) << std::endl;
    std::cout << "Normalize signal? " << (do_normalize ? "yes" : "no") << std::endl;

    std::cout << "Starting to process..." << std::endl;
    try {
        App* app = new App(path);
        app->run();
        delete app;
    } catch(std::exception& ex) {
        std::cout << "Exception caught: " << ex.what() << std::endl;
        return -1;
    }

    std::cout << "Processing done!" << std::endl;
    return 0;
}



