#ifndef FIT_MODELS_H
#define FIT_MODELS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <vector>

// Purely static class for a fit model with single exp relaxation + noise
class SingleExpFitModel {
public:
    SingleExpFitModel() = delete;

    struct FitData {
        std::vector<double> tau_vals; // Time delays
        std::vector<double> ISF_vals; // Image structure function values as function of tau
    };

    static int fit_res(
        const gsl_vector* fit_params, void* fit_data, gsl_vector* res_vals);
    static int fit_jac(
        const gsl_vector* fit_params, void* fit_data, gsl_matrix* jac_vals);
    static int fit_hessian_curv(
        const gsl_vector* fit_params, const gsl_vector* delta_fit_params,
        void* fit_data, gsl_vector* curv_vals);

    static void estimate_initial_params(const FitData &fit_data, gsl_vector* initial_params);

    static const size_t n_params = 3;    
};

template <typename FitModelType>
class FitSolver {
public:
    /**
     * Default constructor for a fit with a given number of data points
    */
    FitSolver(const size_t &n_vals);
    /**
     * Copy constructor useful for private ompenmp construct
    */
   FitSolver(const FitSolver &fit_solver);
   /**
    * Destructor releasing GSL objects
   */
    ~FitSolver();

    /**
     * Run the fit optimization using the data points currently in the
     * public attribute fit_data.
    */
    void optimize();
    /**
     * Get the idx-th fit parameter
    */
    double get_fit_param(unsigned int idx);

    typename FitModelType::FitData fit_data;

private:
    const gsl_multifit_nlinear_type *solver_type = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *workspace;
    gsl_multifit_nlinear_parameters solver_params;
    gsl_multifit_nlinear_fdf fit_model;
    gsl_vector* initial_params;
};

#endif