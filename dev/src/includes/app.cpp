#include <omp.h>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include <iostream>

#include "app.hpp"
#include "stack.hpp"
#include "utils.hpp"

App::App(Options& options) : options(&options) {}
App::~App()= default;

void App::run() {
    if(utils::stoe(this->options->encoding) != Mono12Packed) {
        std::cout << "Encoding not supported yet" << std::endl;
        return;
    }

    std::cout << "* Loading images ..." << std::endl;
    auto* stack = new Stack<float>(this->options->path,
                                           utils::stoe(this->options->encoding),
                                           this->options->loadNframes,
                                           this->options->do_normalize);

    /*
     * TODO:
     *      add inplace fourier transform
     *      refactor code, move it to stack
     *      padding à gérer de la fftw dans l'alloc ?
     */

    std::cout << "* Spawning " << omp_get_max_threads() << " threads ..." << std::endl;
    int r = fftwf_init_threads();
    if(r == 0) {
        std::cout << "err" << std::endl;
        return;
    }

    fftwf_plan_with_nthreads(omp_get_max_threads());
    fftwf_complex* out  = fftwf_alloc_complex(stack->image_size * this->options->loadNframes);
    fftwf_plan plan = fftwf_plan_dft_r2c_2d(stack->aoi_width,
                                            stack->aoi_height*this->options->loadNframes,
                                            stack->N_images_buffer,
                                            out,
                                            FFTW_ESTIMATE);

    std::cout << "* Performing DFT ..." << std::endl;
    fftwf_execute(plan);

    std::cout << "* Computing differences ..." << std::endl;

    // on commence singlethreaded, pour les i tau
    // on garde en tableau les references vers les differents buffers;
    std::complex<float>* references[this->options->tau];

    // pour chaque tau
    for(int i = 1; i <= this->options->tau; i++){

        // on cree un nouveau buffer et on l'ajoute à la liste des references
        auto* buff = new std::complex<float>[stack->image_size];
        references[i-1] = buff;

        // nos indices aplatis
        int idx, idx_k = 0, idx_kp = 0;

        // pour chaque image
        for(int t = i; t < this->options->loadNframes; t++) { // ajouter le compteur dans la stack

            // pour chaque ligne
            for(int y = 0; y < stack->aoi_height; y++) {

                // pour chaque colonne
                for(int x = 0; x < stack->aoi_width; x++) {
                    // on calcule le nouvel indice
                    idx = x + y*stack->aoi_width;
                    idx_k =  idx + t*stack->image_size;
                    idx_kp = idx + (t-i)*stack->image_size;

                    // on calcule la difference
                    auto diff = std::complex<float>((out[idx_k][REAL] - out[idx_k][REAL]),
                                                    (out[idx_k][IMAG] - out[idx_kp][IMAG]));
                    // si on est à la premiere image, moyenne = juste la premiere difference
                    if(t == 0)
                        buff[idx] = std::pow(diff, 2);
                    else { // sinon, on met à jour la moyenne avec l'algo online
                        float scFactor = static_cast<float>((t-1)/t);
                        buff[idx] = scFactor * buff[idx] + std::pow(diff/(float) t, 2);
                    }
                }
            }
        }
    }


    std::cout << "* Cleaning ..." << std::endl;
    fftwf_cleanup_threads();
    fftwf_destroy_plan(plan);
    delete stack;
    fftw_free(out);
};
