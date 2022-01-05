#include <TMB.hpp>

// Function for keeping values positive
template<class Type>
Type pos_fun(Type x, Type eps, Type &pen){
    pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x - eps, 2), Type(0));
    return CppAD::CondExpGe(x, eps, x, eps / (Type(2) - x / eps));
}

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
    DATA_INTEGER(log_sd_B_option);
    DATA_INTEGER(log_B0_option);
    DATA_INTEGER(log_r_option);
    DATA_INTEGER(log_q_option);
    DATA_INTEGER(log_sd_I_option);
    DATA_INTEGER(logit_cor_option);
    DATA_MATRIX(pe_covariates);
    DATA_MATRIX(K_covariates);

    // Parameters
    PARAMETER_MATRIX(delta); // process error
    PARAMETER(mean_log_sd_B);
    PARAMETER(log_sd_log_sd_B);
    PARAMETER_VECTOR(log_sd_B);
    PARAMETER(mean_logit_cor);
    PARAMETER(log_sd_logit_cor);
    PARAMETER_VECTOR(logit_cor);
    PARAMETER(mean_log_B0);
    PARAMETER(log_sd_log_B0);
    PARAMETER_VECTOR(log_B0);
    PARAMETER(log_K);
    PARAMETER(mean_log_r);
    PARAMETER(log_sd_log_r);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER(mean_log_q);
    PARAMETER(log_sd_log_q);
    PARAMETER_VECTOR(log_q);
    PARAMETER(mean_log_sd_I);
    PARAMETER(log_sd_log_sd_I);
    PARAMETER_VECTOR(log_sd_I);
    PARAMETER_VECTOR(pe_betas);
    PARAMETER_VECTOR(K_betas);

    // Dim
    int nY = delta.rows();         // number of years
    int nS = delta.cols();         // number of species
    int nL = L.size();
    int nI = I.size();

    // Containers
    matrix<Type> B(nY, nS);
    matrix<Type> log_B(nY, nS);
    vector<Type> tot_B(nY);
    vector<Type> log_tot_B(nY);
    vector<Type> B_vec(nL);
    vector<Type> log_B_vec(nL);
    vector<Type> log_B_res(nL);
    vector<Type> log_B_std_res(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_I_res(nI);
    vector<Type> log_I_std_res(nI);
    matrix<Type> L_mat(nY, nS);
    vector<Type> F(nL);
    vector<Type> log_F(nL);
    vector<Type> K_vec(nY);
    vector<Type> log_K_vec(nY);

    // Transformations
    vector<Type> B0 = exp(log_B0);
    vector<Type> log_I = log(I);
    vector<Type> sd_B = exp(log_sd_B);
    vector<Type> cor = 2.0 / (1.0 + exp(-logit_cor)) - 1.0; // want cor to be between -1 and 1
    Type K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);
    Type sd_log_sd_B = exp(log_sd_log_sd_B);
    Type sd_log_B0 = exp(log_sd_log_B0);
    Type sd_log_r = exp(log_sd_log_r);
    Type sd_log_q = exp(log_sd_log_q);
    Type sd_log_sd_I = exp(log_sd_log_sd_I);
    Type sd_logit_cor = exp(log_sd_logit_cor);


    // Set-up a landings matrix
    for (int i = 0; i < nL; i++) {
        L_mat(L_year(i), L_species(i)) = L(i);
    }

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Priors / random effects
    if (log_B0_option > 1) {
        for(int i = 0; i < log_B0.size(); i++) {
            nll -= dnorm(log_B0(i), mean_log_B0, sd_log_B0, true);
        }
    }
    if (log_sd_B_option > 1) {
        for(int i = 0; i < log_sd_B.size(); i++) {
            nll -= dnorm(log_sd_B(i), mean_log_sd_B, sd_log_sd_B, true);
        }
    }
    if (log_r_option > 1) {
        for(int i = 0; i < log_r.size(); i++) {
            nll -= dnorm(log_r(i), mean_log_r, sd_log_r, true);
        }
    }
    if (log_q_option > 1) {
        for(int i = 0; i < log_q.size(); i++) {
            nll -= dnorm(log_q(i), mean_log_q, sd_log_q, true);
        }
    }
    if (log_sd_I_option > 1) {
        for(int i = 0; i < log_sd_I.size(); i++) {
            nll -= dnorm(log_sd_I(i), mean_log_sd_I, sd_log_sd_I, true);
        }
    }
    if (logit_cor_option > 1) {
        for(int i = 0; i < logit_cor.size(); i++) {
            nll -= dnorm(logit_cor(i), mean_logit_cor, sd_logit_cor, true);
        }
    }


    // Process equation
    using namespace density;
    vector<Type> pe_covar_effect = pe_covariates * pe_betas;
    vector<Type> K_covar_effect = K_covariates * K_betas;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            K_vec(i) = K + K_covar_effect(i);
            if (i == 0) {
                B(i, j) = B0(j) * exp(pe_covar_effect(i) + delta(i, j));
            } else {
                B(i, j) = (B(i - 1, j) + (r(j) / (m(j) - 1.0)) * B(i - 1, j) *
                    (1.0 - pow((B.row(i - 1).sum() / K_vec(i)), m(j) - 1.0)) -
                    L_mat(i - 1, j)) * exp(pe_covar_effect(i) + delta(i, j));
            }
            B(i, j) = pos_fun(B(i, j), min_B, pen);
            log_B(i, j) = log(B(i, j));
        }
        nll += VECSCALE(UNSTRUCTURED_CORR(cor), sd_B)(delta.row(i));
    }
    log_K_vec = log(K_vec);

    // Transform B matrix to vector
    for (int i = 0; i < nL; i++) {
        log_B_vec(i) = log_B(L_year(i), L_species(i));
        B_vec(i) = exp(log_B_vec(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_B_vec(I_sy(i));
        if (keep(i) == 1) {
            nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        }
        log_I_res(i) = log_I(i) - log_pred_I(i);
        log_I_std_res(i) = log_I_res(i) / sd_I(I_survey(i));
    }


    // Calculate standardized residuals (process errors) and F
    for (int i = 0; i < nL; i++) {
        log_B_res(i) = delta(L_year(i), L_species(i));
        log_B_std_res(i) = log_B_res(i) / sd_B(L_species(i));
        F(i) = L(i) / B_vec(i);
        log_F(i) = log(F(i));
    }

    // Total biomass
    tot_B = B.rowwise().sum();
    log_tot_B = log(tot_B);

    // AD report values
    ADREPORT(log_B_vec);
    ADREPORT(log_pred_I);
    ADREPORT(log_F);
    ADREPORT(log_K_vec);
    ADREPORT(log_tot_B);

    REPORT(log_B_res);
    REPORT(log_B_std_res);
    REPORT(log_I_res);
    REPORT(log_I_std_res);
    REPORT(log_pred_I);

    REPORT(pen);

    nll += pen;
    return nll;

}

