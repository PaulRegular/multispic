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
    DATA_IVECTOR(I_sy);          // species-year index that corresponds to L row number
    DATA_IVECTOR(I_survey);      // note: survey factor is species specific
    DATA_SCALAR(min_P);

    // Parameters
    PARAMETER_MATRIX(pe);         // process error
    PARAMETER_VECTOR(logit_P0);   // pop size in first year
    PARAMETER_VECTOR(log_sd_pe);
    PARAMETER_VECTOR(logit_cor);
    PARAMETER_VECTOR(log_K);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_sd_I);

    // Transformations
    vector<Type> log_I = log(I);
    vector<Type> sd_pe = exp(log_sd_pe);
    vector<Type> P0 = 2.0 / (1.0 + exp(-logit_cor));  // want P0 to be between 0 and 2
    vector<Type> cor = 2.0 / (1.0 + exp(-logit_cor)) - 1.0; // want cor to be between -1 and 1
    vector<Type> K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);

    // Dim
    int nY = pe.rows();
    int nL = L.size();
    int nI = I.size();

    // Containers
    vector<Type> log_P(nL);
    vector<Type> pred_P(nL);
    vector<Type> log_pred_P(nL);
    vector<Type> log_res_P(nL);
    vector<Type> P(nL);
    vector<Type> B(nL);
    vector<Type> log_B(nL);
    vector<Type> log_pred_I(nI);
    vector<Type> log_res_I(nI);
    vector<Type> std_res_I(nI);

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Process equation
    for (int i = 0; i < nL; i++){
        if (L_year(i) == 0) {
            pred_P(i) = P0(L_species(i));
            P(i) = pred_P(i) * exp(pe(L_year(i), L_species(i)));
        } else {
            pred_P(i) = P(i - 1) +
                (r(L_species(i)) / (m(L_species(i)) - 1)) *
                (Type(1) - pow(P(i - 1), m(L_species(i)) - Type(1))) -
                (L(i - 1) / K(L_species(i)));
            pred_P(i) = pos_fun(pred_P(i), min_P, pen);
            P(i) = pred_P(i) * exp(pe(L_year(i), L_species(i)));
        }
        B(i) = P(i) * K(L_species(i));
        log_P(i) = log(P(i));
    }

    using namespace density;
    for (int i = 0; i < nY; i++) {
        nll += VECSCALE(UNSTRUCTURED_CORR(cor), sd_pe)(pe.row(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_P(I_sy(i)) + log_K(I_species(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        log_res_I(i) = log_I(i) - log_pred_I(i);
        std_res_I(i) = log_res_I(i) / sd_I(I_survey(i));
    }

    // More transformations
    log_B = log(B);
    log_pred_P = log(pred_P);
    log_res_P = log_P - log_pred_P;

    // AD report values
    ADREPORT(log_P);
    ADREPORT(log_pred_P);
    ADREPORT(log_res_P);
    ADREPORT(log_B);
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);
    ADREPORT(std_res_I);

    REPORT(pen);
    nll += pen;
    return nll;

}

