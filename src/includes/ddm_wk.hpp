#ifndef DDM_WK_H
#define DDM_WK_H

#include "ddm.hpp"

class DDMWK : public DDM {
public:
    /**
     * Construct a DDMDiff object, which uses the "Wiener-Khinchin" (WK) algorithm for
     * calculating the image structure function.
    */
    DDMWK(Stack &stack, utils::Options& options);
    virtual ~DDMWK() {}

private:
    virtual void compute_DDM();

    void ddm_loop_autovec();
};

#endif