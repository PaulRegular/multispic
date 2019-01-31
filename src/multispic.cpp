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
    DATA_IVECTOR(I_survey);  // note: survey factor is species specific
    DATA_IVECTOR(I_sy);      // species-year index that corresponds to L row number
    DATA_SCALAR(min_P);

    // Parameters
    PARAMETER_MATRIX(log_P);
    PARAMETER_VECTOR(log_P0);
    PARAMETER_VECTOR(log_sd_P);
    PARAMETER_VECTOR(logit_cor);
    PARAMETER_VECTOR(log_K);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_sd_I);

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
    vector<Type> log_P_res(nL);
    vector<Type> B_vec(nL);
    vector<Type> log_B_vec(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_res_I(nI);
    vector<Type> std_res_I(nI);
    matrix<Type> L_mat(nY, nS);

    // Transformations
    vector<Type> P0 = exp(log_P0);
    vector<Type> log_I = log(I);
    vector<Type> sd_P = exp(log_sd_P);
    vector<Type> cor = 2.0 / (1.0 + exp(-logit_cor)) - 1.0; // want cor to be between -1 and 1
    vector<Type> K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);

    // Set-up a landings matrix and vector of P
    for (int i = 0; i < nL; i++) {
        L_mat(L_year(i), L_species(i)) = L(i);
        log_P_vec(i) = log_P(L_year(i), L_species(i));
        P_vec(i) = exp(log_P_vec(i));
        B_vec(i) = P_vec(i) * K(L_species(i));
        log_B_vec(i) = log(B_vec(i));
    }

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Process equation
    using namespace density;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nS; j++) {
            if (i == 0) {
                pred_P(i, j) = P0(j);
            } else {
                P(i - 1, j) = exp(log_P(i - 1, j));
                pred_P(i, j) = P(i - 1, j) +
                    (r(j) / (m(j) - 1.0)) *
                    (1.0 - pow(P(i - 1, j), m(j) - 1.0)) -
                    (L_mat(i - 1, j) / K(j));
                pred_P(i, j) = pos_fun(pred_P(i, j), min_P, pen);
                log_pred_P(i, j) = log(pred_P(i, j));
            }
        }
        nll += VECSCALE(UNSTRUCTURED_CORR(cor), sd_P)(log_P.row(i) - log_pred_P.row(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_P_vec(I_sy(i)) + log_K(I_species(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        log_res_I(i) = log_I(i) - log_pred_I(i);
        std_res_I(i) = log_res_I(i) / sd_I(I_survey(i));
    }


    // Calculate residuals
    for (int i = 0; i < nL; i++) {
        log_P_res(i) = log_P(L_year(i), L_species(i)) - log_pred_P(L_year(i), L_species(i));
    }

    // AD report values
    ADREPORT(log_P_vec);
    ADREPORT(log_P_res);
    ADREPORT(log_B_vec);
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);
    ADREPORT(std_res_I);

    REPORT(pen);

    nll += pen;
    return nll;

}

