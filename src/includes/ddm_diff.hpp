#ifndef DDM_DIFF_H
#define DDM_DIFF_H

#include "ddm.hpp"

class DDMDiff : public DDM {
public:
    /**
     * Construct a DDMDiff object, which uses the "difference" algorithm for calculating
     * the image structure function.
    */
    DDMDiff(Stack &stack, utils::Options& options);

private:
    virtual void compute_DDM();

    /**
     * Comment generated by GPT-3
     *
     * Perform the DDM main loop computation using AVX512 instructions.
     *
     * This function computes the DDM differences based on the FFT results, utilizing AVX512
     * instructions for vectorization. It performs the following steps:
     *   1. Initializes the raw DDM buffer to zeros.
     *   2. Utilizes AVX512 instructions for parallelized DDM computation.
     *   3. Computes squared differences and accumulates them in the raw DDM buffer.
     *   4. Normalizes the raw DDM buffer by the mean weight.
     */
    void ddm_loop_avx512();
    /**
    * Same as ddm_loop_avx512 with avx2 instruction sets
     */
    void ddm_loop_avx2();
    /**
     * Same as ddm_loop_avx512 with compiler autovecterization
     */
    void ddm_loop_autovec();
};

#endif