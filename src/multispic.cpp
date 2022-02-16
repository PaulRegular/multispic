#include <TMB.hpp>
#include "utils.hpp" // Atomic functions positive function and for censored likelihoods

template<class Type>
Type objective_function<Type>::operator() ()
{

    // Data inputs
    DATA_VECTOR(L);
    DATA_IVECTOR(L_species);
    DATA_IVECTOR(L_year);
    DATA_VECTOR(I);
    DATA_IVECTOR(I_species);
    DATA_IVECTOR(I_survey);  // index vector for survey parameters
    DATA_IVECTOR(I_sy);      // species-year index that corresponds to L row number
    DATA_IVECTOR(keep);      // useful for cross-validation
    DATA_SCALAR(min_B);
    DATA_INTEGER(log_K_option);
    DATA_INTEGER(log_sd_B_option);
    DATA_INTEGER(log_B0_option);
    DATA_INTEGER(log_r_option);
    DATA_INTEGER(log_q_option);
    DATA_INTEGER(log_sd_I_option);
    DATA_INTEGER(logit_rho_option);
    DATA_INTEGER(logit_phi_option);
    DATA_INTEGER(K_betas_option);
    DATA_INTEGER(pe_betas_option);
    DATA_SCALAR(lower_log_K);
    DATA_SCALAR(upper_log_K);
    DATA_SCALAR(lower_log_sd_B);
    DATA_SCALAR(upper_log_sd_B);
    DATA_SCALAR(lower_log_B0);
    DATA_SCALAR(upper_log_B0);
    DATA_SCALAR(lower_log_r);
    DATA_SCALAR(upper_log_r);
    DATA_SCALAR(lower_log_q);
    DATA_SCALAR(upper_log_q);
    DATA_SCALAR(lower_log_sd_I);
    DATA_SCALAR(upper_log_sd_I);
    DATA_SCALAR(lower_logit_rho);
    DATA_SCALAR(upper_logit_rho);
    DATA_SCALAR(mean_logit_phi);
    DATA_SCALAR(sd_logit_phi);
    DATA_SCALAR(lower_logit_phi);
    DATA_SCALAR(upper_logit_phi);
    DATA_SCALAR(mean_K_betas);
    DATA_SCALAR(sd_K_betas);
    DATA_SCALAR(lower_K_betas);
    DATA_SCALAR(upper_K_betas);
    DATA_SCALAR(mean_pe_betas);
    DATA_SCALAR(sd_pe_betas);
    DATA_SCALAR(lower_pe_betas);
    DATA_SCALAR(upper_pe_betas);
    DATA_SCALAR(dmuniform_sd);
    DATA_MATRIX(survey_covariates);
    DATA_MATRIX(pe_covariates);
    DATA_MATRIX(K_covariates);
    DATA_IVECTOR(K_map);
    DATA_MATRIX(B_groups);

    // Parameters
    PARAMETER_MATRIX(log_B);
    PARAMETER(mean_log_sd_B);
    PARAMETER(log_sd_log_sd_B);
    PARAMETER_VECTOR(log_sd_B);
    PARAMETER(mean_logit_rho);
    PARAMETER(log_sd_logit_rho);
    PARAMETER_VECTOR(logit_rho);
    PARAMETER(logit_phi);
    PARAMETER_VECTOR(log_K);
    PARAMETER(mean_log_K);
    PARAMETER(log_sd_log_K);
    PARAMETER(mean_log_B0);
    PARAMETER(log_sd_log_B0);
    PARAMETER_VECTOR(log_B0);
    PARAMETER(mean_log_r);
    PARAMETER(log_sd_log_r);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER(mean_log_q);
    PARAMETER(log_sd_log_q);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_q_betas);
    PARAMETER(mean_log_sd_I);
    PARAMETER(log_sd_log_sd_I);
    PARAMETER_VECTOR(log_sd_I);
    PARAMETER_VECTOR(log_sd_I_betas);
    PARAMETER_VECTOR(K_betas);
    PARAMETER_VECTOR(pe_betas);

    // Dim
    int nK = log_K.size();         // number of Ks estimated
    int nY = log_B.rows();         // number of years
    int nS = log_B.cols();         // number of species
    int nL = L.size();
    int nI = I.size();
    int n_surveys = survey_covariates.rows();

    // Containers
    matrix<Type> pred_B(nY, nS);
    matrix<Type> log_pred_B(nY, nS);
    array<Type> delta(nY, nS); // AR1 function expects an array
    matrix<Type> tot_B_mat(nY, nS);
    matrix<Type> tot_B(nY, nK);
    matrix<Type> log_tot_B(nY, nK);
    matrix<Type> log_grouped_K(nY, nK);
    vector<Type> B_vec(nL);
    vector<Type> log_B_vec(nL);
    vector<Type> log_res_pe(nL);
    vector<Type> log_std_res_pe(nL);
    vector<Type> log_pe(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_I_res(nI);
    vector<Type> log_I_std_res(nI);
    matrix<Type> L_mat(nY, nS);
    vector<Type> F(nL);
    vector<Type> log_F(nL);
    matrix<Type> K_mat(nY, nS);
    vector<Type> K_vec(nL);
    vector<Type> log_K_vec(nL);
    vector<Type> pred_log_q(n_surveys);
    vector<Type> pred_log_sd_I(n_surveys);

    // Transformations
    matrix<Type> B = exp(log_B.array());
    vector<Type> B0 = exp(log_B0);
    vector<Type> log_I = log(I);
    vector<Type> sd_B = exp(log_sd_B);
    vector<Type> rho = 2.0 / (1.0 + exp(-logit_rho)) - 1.0; // want cor to be between -1 and 1
    Type phi = 1.0 / (1.0 + exp(-logit_phi)); // want temporal cor to be between 0 and 1
    vector<Type> K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    Type sd_log_sd_B = exp(log_sd_log_sd_B);
    Type sd_log_K = exp(log_sd_log_K);
    Type sd_log_B0 = exp(log_sd_log_B0);
    Type sd_log_r = exp(log_sd_log_r);
    Type sd_log_q = exp(log_sd_log_q);
    Type sd_log_sd_I = exp(log_sd_log_sd_I);
    Type sd_logit_rho = exp(log_sd_logit_rho);


    // Set-up a vector of B, landings matrix, and pe covariate effects
    vector<Type> pe_covar_vec = pe_covariates * pe_betas;
    vector<Type> K_covar_vec = K_covariates * K_betas;
    matrix<Type> pe_covar_mat(nY, nS);
    matrix<Type> K_covar_mat(nY, nS);
    for (int i = 0; i < nL; i++) {
        L_mat(L_year(i), L_species(i)) = L(i);
        log_B_vec(i) = log_B(L_year(i), L_species(i));
        B_vec(i) = exp(log_B_vec(i));
        pe_covar_mat(L_year(i), L_species(i)) = pe_covar_vec(i);
        K_covar_mat(L_year(i), L_species(i)) = K_covar_vec(i);
    }

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Priors / random effects
    if (log_K_option > 1) {
        for(int i = 0; i < log_K.size(); i++) {
            if (log_K_option == 4) {
                nll += dmuniform(log_K(i), lower_log_K, upper_log_K, dmuniform_sd);
            } else {
                nll -= dnorm(log_K(i), mean_log_K, sd_log_K, true);
            }
        }
    }
    if (log_B0_option > 1) {
        for(int i = 0; i < log_B0.size(); i++) {
            if (log_B0_option == 4) {
                nll += dmuniform(log_B0(i), lower_log_B0, upper_log_B0, dmuniform_sd);
            } else {
                nll -= dnorm(log_B0(i), mean_log_B0, sd_log_B0, true);
            }
        }
    }
    if (log_sd_B_option > 1) {
        for(int i = 0; i < log_sd_B.size(); i++) {
            if (log_sd_B_option == 4) {
                nll += dmuniform(log_sd_B(i), lower_log_sd_B, upper_log_sd_B, dmuniform_sd);
            } else {
                nll -= dnorm(log_sd_B(i), mean_log_sd_B, sd_log_sd_B, true);
            }
        }
    }
    if (log_r_option > 1) {
        for(int i = 0; i < log_r.size(); i++) {
            if (log_r_option == 4) {
                nll += dmuniform(log_r(i), lower_log_r, upper_log_r, dmuniform_sd);
            } else {
                nll -= dnorm(log_r(i), mean_log_r, sd_log_r, true);
            }
        }
    }

    if (log_q_option == 2) {
        pred_log_q = log_q;
    } else {
        pred_log_q = survey_covariates * log_q_betas;
    }
    if (log_q_option > 1) {
        for(int i = 0; i < pred_log_q.size(); i++) {
            if (log_q_option == 4) {
                nll += dmuniform(pred_log_q(i), lower_log_q, upper_log_q, dmuniform_sd);
            } else {
                nll -= dnorm(pred_log_q(i), mean_log_q, sd_log_q, true);
            }
        }
    }

    if (log_sd_I_option == 2) {
        pred_log_sd_I = log_sd_I;
    } else {
        pred_log_sd_I = survey_covariates * log_sd_I_betas;
    }
    if (log_sd_I_option > 1) {
        for(int i = 0; i < pred_log_sd_I.size(); i++) {
            if (log_sd_I_option == 4) {
                nll += dmuniform(pred_log_sd_I(i), lower_log_sd_I, upper_log_sd_I, dmuniform_sd);
            } else {
                nll -= dnorm(pred_log_sd_I(i), mean_log_sd_I, sd_log_sd_I, true);
            }
        }
    }
    if (logit_rho_option > 1) {
        for(int i = 0; i < logit_rho.size(); i++) {
            if (logit_rho_option == 4) {
                nll += dmuniform(logit_rho(i), lower_logit_rho, upper_logit_rho, dmuniform_sd);
            } else {
                nll -= dnorm(logit_rho(i), mean_logit_rho, sd_logit_rho, true);
            }
        }
    }
    if (logit_phi_option > 1) {
        if (logit_phi_option == 4) {
            nll += dmuniform(logit_phi, lower_logit_phi, upper_logit_phi, dmuniform_sd);
        } else {
            nll -= dnorm(logit_phi, mean_logit_phi, sd_logit_phi, true);
        }
    }
    if (pe_betas_option > 1) {
        for(int i = 0; i < pe_betas.size(); i++) {
            if (pe_betas_option == 4) {
                nll += dmuniform(pe_betas(i), lower_pe_betas, upper_pe_betas, dmuniform_sd);
            } else {
                nll -= dnorm(pe_betas(i), mean_pe_betas, sd_pe_betas, true);
            }
        }
    }
    if (K_betas_option > 1) {
        for(int i = 0; i < K_betas.size(); i++) {
            if (K_betas_option == 4) {
                nll += dmuniform(K_betas(i), lower_K_betas, upper_K_betas, dmuniform_sd);
            } else {
                nll -= dnorm(K_betas(i), mean_K_betas, sd_K_betas, true);
            }
        }
    }


    // Process equation
    using namespace density;

    vector<Type> grouped_B(nS);
    vector<Type> B_row(nS);
    vector<Type> B_groups_row(nS);
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            K_mat(i, j) = K(K_map(j)) * exp(K_covar_mat(i, j));
            if (i == 0) {
                pred_B(i, j) = B0(j) * exp(pe_covar_mat(i, j));
            } else {
                B_row = B.row(i - 1);
                B_groups_row = B_groups.row(j);
                tot_B_mat(i - 1, j) = (B_row * B_groups_row).sum(); // use 0s and 1s to group B sum
                pred_B(i, j) = (B(i - 1, j) + (r(j) / (m(j) - 1.0)) * B(i - 1, j) *
                    (1.0 - pow((tot_B_mat(i - 1, j) / K_mat(i, j)), m(j) - 1.0)) -
                    L_mat(i - 1, j)) * exp(pe_covar_mat(i, j));
            }
            pred_B(i, j) = pos_fun(pred_B(i, j), min_B, pen);
            log_pred_B(i, j) = log(pred_B(i, j));
            delta(i, j) = log_B(i, j) - log_pred_B(i, j);
        }
    }
    if (nS == 1) { // there were issues with the matrix math when attempting to fit to one species using the 'else' code below
        Type sd_B_scalar = sqrt(((sd_B(0) * sd_B(0)) / (1 - (phi * phi))));
        vector<Type> delta_vec = delta.col(0);
        nll += SCALE(AR1(phi), sd_B_scalar)(delta_vec);
    } else {
        array<Type> delta_transpose = delta.transpose();
        vector<Type> sd_B_vec = sqrt(((sd_B * sd_B) / (1 - (phi * phi))));
        nll += AR1(phi, VECSCALE(UNSTRUCTURED_CORR(rho), sd_B_vec))(delta_transpose);
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = pred_log_q(I_survey(i)) + log_B_vec(I_sy(i));
        if (keep(i) == 1) {
            nll -= dnorm(log_I(i), log_pred_I(i), exp(pred_log_sd_I(I_survey(i))), true);
        }
        log_I_res(i) = log_I(i) - log_pred_I(i);
        log_I_std_res(i) = log_I_res(i) / exp(pred_log_sd_I(I_survey(i)));
    }


    // Calculate process residuals (aka process error, aka delta) and F
    for (int i = 0; i < nL; i++) {
        log_pe(i) = pe_covar_mat(L_year(i), L_species(i)) + delta(L_year(i), L_species(i)); // process error = covariate effect + residual process error
        log_res_pe(i) = delta(L_year(i), L_species(i));                                     // residual process error not explained by the covariate
        log_std_res_pe(i) = log_res_pe(i) / sd_B(L_species(i));
        F(i) = L(i) / B_vec(i);
        log_F(i) = log(F(i));
        K_vec(i) = K_mat(L_year(i), L_species(i));
        log_K_vec(i) = log(K_vec(i));
    }

    // Total biomass by group
    tot_B.setZero();
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            tot_B(i, K_map(j)) += B(i, j);
        }
    }
    log_tot_B = log(tot_B.array());

    // K by group
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            log_grouped_K(i, K_map(j)) = log(K_mat(i, j));
        }
    }

    // AD report values
    ADREPORT(log_B_vec);
    ADREPORT(log_pred_I);
    ADREPORT(log_F);
    ADREPORT(log_K_vec);
    ADREPORT(log_tot_B);
    ADREPORT(pred_log_q);
    ADREPORT(pred_log_sd_I);
    ADREPORT(log_grouped_K);

    REPORT(B);
    REPORT(B_vec);
    REPORT(F);
    REPORT(K_vec);
    REPORT(tot_B);
    REPORT(log_res_pe);
    REPORT(log_std_res_pe);
    REPORT(log_pe);
    REPORT(log_I_res);
    REPORT(log_I_std_res);
    REPORT(log_pred_I);
    REPORT(K_mat);
    REPORT(tot_B_mat);

    REPORT(pen);

    nll += pen;
    return nll;

}

