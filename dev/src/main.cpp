#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>

#include "includes/app.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv) {

    Options options;
    bool do_print_help = false;

    std::string help =
"This software intends to perform image processing. \
List of arguments that can be passed at the command-line:\n\
        --path, -p: Path to the file containing the signal to process. Defaults to sample.dat. \n\
        --encoding, -e: Specify the frame encoding. Defaults to Mono12Packed.\n\
        --normalize, -n: Normalize the signal by averaging all the pixels. Defaults to false.\n\
        --loadNFrames, -N: Load a fixed number N of frames.\n\
        --help, -h: Print this help message.\n\
        --tau, -t: Choose the number of dt you're interested in to compute the differences. Example :\n\
                 If -t 3 is passed, it will compute differences between images separated by 1, 2 and 3 images.";

    po::options_description desc("Options");
    desc.add_options()
        ("path,p", po::value<std::string>(&options.path)->default_value("sample.dat"), "Path to the file containing the signal to process")
        ("loadNframes,N", po::value<int>(&options.loadNframes)->default_value(1), "Load a fixed number of N frames")
        ("encoding,e", po::value<std::string>(&options.encoding)->default_value("Mono12Packed"), "Specify the frame encoding")
        ("normalize,n", po::bool_switch(&options.do_normalize), "Normalize the signal data by averaging all the pixels")
        ("help,h", po::bool_switch(&do_print_help), "Print out help message")
        ("tau,t", po::value<int>(&options.tau)->default_value(1), "Choose the number of dt you're interested in to compute the differences");

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    if(do_print_help) {
        std::cout << help << std::endl;
        return 0;
    }

    // recap options
    std::cout << "File: " << options.path << std::endl;
    std::cout << "Encoding: " << options.encoding << std::endl;
    std::cout << "Frames to load: " << std::to_string(options.loadNframes) << std::endl;
    std::cout << "Differences between images: 1 -> " << std::to_string(options.tau) << std::endl;
    std::cout << "Normalize signal? " << (options.do_normalize ? "yes" : "no") << std::endl;

    std::cout << "Starting to process..." << std::endl;
    try {
        App* app = new App(options);
        app->run();
        delete app;
    } catch(std::exception& ex) {
        std::cout << "Exception caught: " << ex.what() << std::endl;
        return -1;
    }

    std::cout << "Processing done!" << std::endl;
    return 0;
}



