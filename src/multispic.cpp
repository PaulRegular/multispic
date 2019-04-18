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
    DATA_SCALAR(min_P);
    DATA_INTEGER(log_sd_P_option);
    DATA_INTEGER(log_P0_option);
    DATA_INTEGER(log_r_option);
    DATA_INTEGER(log_q_option);
    DATA_INTEGER(log_sd_I_option);
    DATA_MATRIX(covariates);

    // Parameters
    PARAMETER_MATRIX(log_P);
    PARAMETER(mean_log_sd_P);
    PARAMETER(log_sd_log_sd_P);
    PARAMETER_VECTOR(log_sd_P);
    PARAMETER_VECTOR(logit_cor);
    PARAMETER(mean_log_P0);
    PARAMETER(log_sd_log_P0);
    PARAMETER_VECTOR(log_P0);
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
    PARAMETER_VECTOR(betas);

    // Dim
    int nY = log_P.rows();         // number of years
    int nS = log_P.cols();         // number of species
    int nL = L.size();
    int nI = I.size();

    // Containers
    matrix<Type> P(nY, nS);
    matrix<Type> pred_P(nY, nS);
    matrix<Type> log_pred_P(nY, nS);
    vector<Type> P_vec(nL);
    vector<Type> log_P_vec(nL);
    vector<Type> B_vec(nL);
    vector<Type> log_B_vec(nL);
    vector<Type> log_P_res(nL);
    vector<Type> log_P_std_res(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_I_res(nI);
    vector<Type> log_I_std_res(nI);
    matrix<Type> L_mat(nY, nS);
    vector<Type> F(nL);
    vector<Type> log_F(nL);

    // Transformations
    vector<Type> P0 = exp(log_P0);
    vector<Type> log_I = log(I);
    vector<Type> sd_P = exp(log_sd_P);
    vector<Type> cor = 2.0 / (1.0 + exp(-logit_cor)) - 1.0; // want cor to be between -1 and 1
    Type K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);
    Type sd_log_sd_P = exp(log_sd_log_sd_P);
    Type sd_log_P0 = exp(log_sd_log_P0);
    Type sd_log_r = exp(log_sd_log_r);
    Type sd_log_q = exp(log_sd_log_q);
    Type sd_log_sd_I = exp(log_sd_log_sd_I);

    // Set-up a landings matrix and vector of P
    for (int i = 0; i < nL; i++) {
        L_mat(L_year(i), L_species(i)) = L(i);
        log_P_vec(i) = log_P(L_year(i), L_species(i));
        P_vec(i) = exp(log_P_vec(i));
    }

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Priors / random effects
    if (log_P0_option > 0) {
        for(int i = 0; i < log_P0.size(); i++) {
            nll -= dnorm(log_P0(i), mean_log_P0, sd_log_P0, true);
        }
    }
    if (log_sd_P_option > 0) {
        for(int i = 0; i < log_sd_P.size(); i++) {
            nll -= dnorm(log_sd_P(i), mean_log_sd_P, sd_log_sd_P, true);
        }
    }
    if (log_r_option > 0) {
        for(int i = 0; i < log_r.size(); i++) {
            nll -= dnorm(log_r(i), mean_log_r, sd_log_r, true);
        }
    }
    if (log_q_option > 0) {
        for(int i = 0; i < log_q.size(); i++) {
            nll -= dnorm(log_q(i), mean_log_q, sd_log_q, true);
        }
    }
    if (log_sd_I_option > 0) {
        for(int i = 0; i < log_sd_I.size(); i++) {
            nll -= dnorm(log_sd_I(i), mean_log_sd_I, sd_log_sd_I, true);
        }
    }

    // Process equation
    using namespace density;
    vector<Type> covar_effect = covariates * betas;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            if (i == 0) {
                pred_P(i, j) = P0(j);
            } else {
                P(i - 1, j) = exp(log_P(i - 1, j));
                pred_P(i, j) = P(i - 1, j) + (r(j) / (m(j) - 1.0)) * P(i - 1, j) *
                    (1.0 - pow(P.row(i - 1).sum(), m(j) - 1.0)) -
                    (L_mat(i - 1, j) / K) + covar_effect(i - 1);
            }
            pred_P(i, j) = pos_fun(pred_P(i, j), min_P, pen);
            log_pred_P(i, j) = log(pred_P(i, j));
        }
        nll += VECSCALE(UNSTRUCTURED_CORR(cor), sd_P)(log_P.row(i) - log_pred_P.row(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++) {
        log_pred_I(i) = log_q(I_survey(i)) + log_P_vec(I_sy(i)) + log_K;
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        log_I_res(i) = log_I(i) - log_pred_I(i);
        log_I_std_res(i) = log_I_res(i) / sd_I(I_survey(i));
    }


    // Calculate residuals and F
    for (int i = 0; i < nL; i++) {
        B_vec(i) = P_vec(i) * K;
        log_B_vec(i) = log(B_vec(i));
        log_P_res(i) = log_P(L_year(i), L_species(i)) - log_pred_P(L_year(i), L_species(i));
        log_P_std_res(i) = log_P_res(i) / sd_P(L_species(i));
        F(i) = L(i) / B_vec(i);
        log_F(i) = log(F(i));
    }

    // AD report values
    ADREPORT(log_P_vec);
    ADREPORT(log_B_vec);
    ADREPORT(log_pred_I);
    ADREPORT(log_F);

    REPORT(log_P_res);
    REPORT(log_P_std_res);
    REPORT(log_I_res);
    REPORT(log_I_std_res);

    REPORT(pen);

    nll += pen;
    return nll;

}

