#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "includes/app.hpp"
#include "includes/utils.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv) {

    utils::Options options;
    bool do_print_help = false;

    std::string help =
"This software performs DDM image processing. \n\
List of arguments that can be passed at the command-line:\n\
        --help, -h: Print this help message.\n\
        --output, -o: output folder, default to '.'.\n\
        --path, -p: Path to the file containing the movie to process.\n\
                    Defaults to sample.dat.\n\
        --bin, -b: Binning factor for the input movie.\n\
                   Defaults to 1.\n\
        --Nframes, -N: Load a fixed number N of frames.\n\
        --Nlags, -M: Number of lag times for the DDM calculation. The lag times\n\
                   are varying from the acquisition sampling time up to maxLag,\n\
                   using a linear scale by default or a log scale if\n\
                   --logScale,-l is used.\n\
        --maxLagTime, -t: Max lag time for the DDM. Defaults to 2 seconds.\n\
        --DDMalgo, -a: Algorithm used for DDM calculation. Two choices are possible:\n\
                    diff (usual DDM calculation based on image differences), or\n\
                    wk (temporal FFT-based calculation using the Wiener-Khinchin\n\
                    theorem). Defaults to wk.\n\
        --logLags, -l: Use a logarithmic scale to distribute the lag times.\n\
        --fit, -f: Fit the DDM data along the lag time axis.\n\
        --maxDecayFreq, -r: Set decay frequency threshold to automatically\n\
                            determine the fit ROI. Defaults to 100 s^(-1).";

    po::options_description desc("Options");
    desc.add_options()
        ("help,h", po::bool_switch(&do_print_help), "Print out help message")
        ("output,o", po::value<std::string>(&options.output_path)->default_value("."),
            "Set the output folder")
        ("path,p", po::value<std::string>(&options.path_movie)->default_value("sample.dat"),
            "Path to the file containing the movie to process")
        ("bin,b", po::value<int>(&options.bin_factor)->default_value(1),
            "Binning factor for the input movie")
        ("Nframes,N", po::value<int>(&options.N_frames)->default_value(1),
            "Load a fixed number N of frames")
        ("Nlags,M", po::value<int>(&options.N_lags)->default_value(1),
            "Choose the number of lag time for the DDM calculation")
        ("maxLagTime,t", po::value<float>(&options.max_lag_time)->default_value(2.),
            "Max lag time for the DDM")
        ("logScale,l", po::bool_switch(&options.log_lags)->default_value(false),
            "Use a logarithmic scale to distribute the lag times")
        ("DDMalgo,a", po::value<std::string>(&options.ddm_algo)->default_value("wk"),
            "Algorithm used for DDM calculation")
        ("fit,f", po::bool_switch(&options.fit)->default_value(false),
            "Fit the DDM data along the lag time axis")
        ("maxDecayFreq,r", po::value<float>(&options.max_decay_freq)->default_value(100.),
            "Set decay frequency threshold to automatically determine the fit ROI");

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

    if(options.N_lags>=options.N_frames) {
        std::cout << "Set the number of tau lesser than number of images to load" << std::endl;
        return -1;
    }

    // recap options
    std::cout << "Output folder: " << options.output_path << std::endl;
    std::cout << "File: " << options.path_movie << std::endl;
    std::cout << "Binning factor: " << options.bin_factor << std::endl;
    std::cout << "N frames: " << options.N_frames << std::endl;
    std::cout << "N lags: " << options.N_lags << std::endl;
    std::cout << "Max lag time: " << options.max_lag_time << " s" << std::endl;
    std::cout << "Log scaling for lag times? "
              << (options.log_lags? "yes" : "no") << std::endl;
    std::cout << "Fit DDM? " << (options.fit ? "yes" : "no") << std::endl;
    std::cout << "Max decay frequency: " << options.max_decay_freq << " s^(-1)" << std::endl;

    // Append input basename to the output path for consistent naming of output files
	auto basename = boost::filesystem::path(options.path_movie).stem().string();
	options.output_path = (boost::filesystem::path(options.output_path) / basename).string();

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