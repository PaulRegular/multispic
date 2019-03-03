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
    DATA_SCALAR(min_B);
    DATA_INTEGER(log_q_option);
    DATA_INTEGER(log_sd_I_option);

    // Parameters
    PARAMETER_MATRIX(log_B);
    PARAMETER_VECTOR(log_sd_B);
    PARAMETER_VECTOR(logit_cor);
    PARAMETER(log_K);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER(mean_log_q);
    PARAMETER(log_sd_log_q);
    PARAMETER_VECTOR(log_q);
    PARAMETER(mean_log_sd_I);
    PARAMETER(log_sd_log_sd_I);
    PARAMETER_VECTOR(log_sd_I);

    // Dim
    int nY = log_B.rows();         // number of years
    int nS = log_B.cols();         // number of species
    int nL = L.size();
    int nI = I.size();

    // Containers
    matrix<Type> B(nY, nS);
    matrix<Type> pred_B(nY, nS);
    matrix<Type> log_pred_B(nY, nS);
    vector<Type> B_vec(nL);
    vector<Type> log_B_vec(nL);
    vector<Type> log_B_res(nL);
    vector<Type> log_B_std_res(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_I_res(nI);
    vector<Type> log_I_std_res(nI);
    matrix<Type> L_mat(nY, nS);

    // Transformations
    vector<Type> log_I = log(I);
    vector<Type> sd_B = exp(log_sd_B);
    vector<Type> cor = 2.0 / (1.0 + exp(-logit_cor)) - 1.0; // want cor to be between -1 and 1
    Type K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);
    Type sd_log_q = exp(log_sd_log_q);
    Type sd_log_sd_I = exp(log_sd_log_sd_I);

    // Set-up a landings matrix and vector of B
    for (int i = 0; i < nL; i++) {
        L_mat(L_year(i), L_species(i)) = L(i);
        log_B_vec(i) = log_B(L_year(i), L_species(i));
        B_vec(i) = exp(log_B_vec(i));
    }

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Priors / random effects
    if (log_q_option > 0) {
        for(int i = 0; i < log_q.size(); i++) {
            nll -= dnorm(log_q(i), mean_log_q, sd_log_q, true);
        }
    }
    if (log_sd_I_option > 0) {
        for(int i = 0; i < log_q.size(); i++) {
            nll -= dnorm(log_sd_I(i), mean_log_sd_I, sd_log_sd_I, true);
        }
    }

    // Process equation
    using namespace density;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            if (i > 0) {
                B(i - 1, j) = exp(log_B(i - 1, j));
                pred_B(i, j) = (r(j) / (m(j) - 1.0)) * B(i - 1, j) *
                    (1.0 - pow((B.row(i - 1).sum() / K), m(j) - 1.0)) -
                    L_mat(i - 1, j);
                pred_B(i, j) = pos_fun(pred_B(i, j), min_B, pen);
                log_pred_B(i, j) = log(pred_B(i, j));
            }
        }
        nll += VECSCALE(UNSTRUCTURED_CORR(cor), sd_B)(log_B.row(i) - log_pred_B.row(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_B_vec(I_sy(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        log_I_res(i) = log_I(i) - log_pred_I(i);
        log_I_std_res(i) = log_I_res(i) / sd_I(I_survey(i));
    }


    // Calculate residuals
    for (int i = 0; i < nL; i++) {
        log_B_res(i) = log_B(L_year(i), L_species(i)) - log_pred_B(L_year(i), L_species(i));
        log_B_std_res(i) = log_B_res(i) / sd_B(L_species(i));
    }

    // AD report values
    ADREPORT(log_B_vec);
    ADREPORT(log_pred_I);

    REPORT(log_B_res);
    REPORT(log_B_std_res);
    REPORT(log_I_res);
    REPORT(log_I_std_res);

    REPORT(pen);

    nll += pen;
    return nll;

}

