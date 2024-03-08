#include "fit_models.hpp"

#include <cmath>
#include <iostream>

int SingleExpFitModel::fit_res(
        const gsl_vector* fit_params, void* fit_data, gsl_vector* res_vals) {
    
    FitData* d = (FitData*)fit_data;
    double A = gsl_vector_get(fit_params, 0);
    double B = gsl_vector_get(fit_params, 1);
    double freq = gsl_vector_get(fit_params, 2);
    for(size_t i=0; i<d->tau_vals.size(); ++i)
        gsl_vector_set(res_vals, i, A*(1-std::exp(-freq*d->tau_vals[i]))+B-d->ISF_vals[i]);    
    return GSL_SUCCESS;
}

int SingleExpFitModel::fit_jac(
        const gsl_vector* fit_params, void* fit_data, gsl_matrix* jac_vals) {
    
    FitData* d = (FitData*)fit_data;
    double A = gsl_vector_get(fit_params, 0);
    double freq = gsl_vector_get(fit_params, 2);
    double tau_i, exp_freq_tau_i;
    for(size_t i=0; i<d->tau_vals.size(); ++i) {
        tau_i = d->tau_vals[i];
        exp_freq_tau_i = std::exp(-freq*tau_i);
        gsl_matrix_set(jac_vals, i, 0, 1-exp_freq_tau_i);
        gsl_matrix_set(jac_vals, i, 1, 1);
        gsl_matrix_set(jac_vals, i, 2, A*tau_i*exp_freq_tau_i);
    }
    return GSL_SUCCESS;
}

int SingleExpFitModel::fit_hessian_curv(
        const gsl_vector* fit_params, const gsl_vector* delta_fit_params,
        void* fit_data, gsl_vector* curv_vals) {
    
    FitData* d = (FitData*)fit_data;
    double A = gsl_vector_get(fit_params, 0);
    double delta_A = gsl_vector_get(delta_fit_params, 0);
    double freq = gsl_vector_get(fit_params, 2);
    double delta_freq = gsl_vector_get(delta_fit_params, 2);    
    double tau_i;
    for(size_t i=0; i<d->tau_vals.size(); ++i) {
        tau_i = d->tau_vals[i];
        gsl_vector_set(
            curv_vals, i,
            tau_i*std::exp(-freq*tau_i)*delta_freq*(2*delta_A-A*tau_i*delta_freq));
    }
    return GSL_SUCCESS;
}

void SingleExpFitModel::estimate_initial_params(
        const FitData &fit_data, gsl_vector* initial_params) {

    auto &tau_vals = fit_data.tau_vals;
    auto &ISF_vals = fit_data.ISF_vals;

    gsl_vector_set(initial_params, 0, ISF_vals.back());
    gsl_vector_set(initial_params, 1, 0.);

    // freq estimation
    const double &ISF_min = *std::min_element(ISF_vals.begin(), ISF_vals.end());
    const double &ISF_max = *std::max_element(ISF_vals.begin(), ISF_vals.end());
    double ISF_thresh = ISF_min + (ISF_max - ISF_min) * (1 - std::exp(-1));
    for(int it=0; it<tau_vals.size(); it++) {
        if(ISF_vals[it]>ISF_thresh) {
            gsl_vector_set(initial_params, 2, 1./tau_vals[it]);
            break;
        }
    }
}

template <typename FitModelType>
FitSolver<FitModelType>::FitSolver(const size_t &n_vals) :
    solver_params(gsl_multifit_nlinear_default_parameters()),
    fit_model(gsl_multifit_nlinear_fdf()),
    initial_params(gsl_vector_alloc(FitModelType::n_params)) {

    fit_data.tau_vals.resize(n_vals);
    fit_data.ISF_vals.resize(n_vals);

    fit_model.f = &FitModelType::fit_res;
    fit_model.df = &FitModelType::fit_jac;
    fit_model.fvv = &FitModelType::fit_hessian_curv;
    fit_model.p = FitModelType::n_params;
    fit_model.n = n_vals;
    fit_model.params = &fit_data;
    
    solver_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    workspace = gsl_multifit_nlinear_alloc(
        solver_type, &solver_params, fit_model.n, fit_model.p);
}

template <typename FitModelType>
FitSolver<FitModelType>::FitSolver(const FitSolver &fit_solver) :
    FitSolver(fit_solver.fit_model.n) {

    fit_data = fit_solver.fit_data;
}

template <typename FitModelType>
FitSolver<FitModelType>::~FitSolver() {
    gsl_vector_free(initial_params);
    gsl_multifit_nlinear_free(workspace);
}

template <typename FitModelType>
void FitSolver<FitModelType>::optimize() {

    // Estimate initial fit parameters if none are provided
    FitModelType::estimate_initial_params(fit_data, initial_params);

    // Initialize solver and iterate until converged
    gsl_multifit_nlinear_init(initial_params, &fit_model, workspace);

    const size_t max_iter = 200;
    const double xtol = 1.0e-8;
    const double gtol = 1.0e-8;
    const double ftol = 1.0e-8;
    int info;
    gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, workspace);
}

template <typename FitModelType>
double FitSolver<FitModelType>::get_fit_param(unsigned int idx) {
    return gsl_vector_get(gsl_multifit_nlinear_position(workspace), idx);
}

template class FitSolver<SingleExpFitModel>;