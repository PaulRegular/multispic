#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
    pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x - eps, 2), Type(0));
    return CppAD::CondExpGe(x, eps, x, eps / (Type(2) - x / eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{

    // Data inputs
    DATA_VECTOR(L);
    DATA_VECTOR(I);
    DATA_IVECTOR(I_survey);
    DATA_IVECTOR(I_year);
    DATA_SCALAR(min_P);

    // Parameters
    PARAMETER_VECTOR(log_P);
    PARAMETER(log_sd_P);
    PARAMETER(log_K);
    PARAMETER(log_r);
    PARAMETER(log_m);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_sd_I);

    // Transformations
    vector<Type> log_I = log(I);
    vector<Type> P = exp(log_P);
    Type sd_P = exp(log_sd_P);
    Type K = exp(log_K);
    Type r = exp(log_r);
    Type m = exp(log_m);
    vector<Type> sd_I = exp(log_sd_I);

    // Containers
    int n_years = L.size();
    vector<Type> pred_P(n_years);
    vector<Type> log_pred_P(n_years);
    vector<Type> log_res_P(n_years);
    vector<Type> B(n_years);
    vector<Type> log_B(n_years);
    int n_obs = I.size();
    vector<Type> log_pred_I(n_obs);
    vector<Type> log_res_I(n_obs);

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Process equation
    nll -= dnorm(log_P(0), Type(0), sd_P, true);
    for (int i = 1; i < n_years; i++){
        pred_P(i) = P(i - 1) + (r / (m - 1)) * (Type(1) - pow(P(i - 1), m - Type(1))) - (L(i - 1) / K);
        pred_P(i) = posfun(pred_P(i), min_P, pen);
        nll -= dnorm(log_P(i), log(pred_P(i)), sd_P, true);
        SIMULATE {
            log_P(i) = rnorm(log(pred_P(i)), sd_P);
        }
    }

    // Observation equations
    for (int i = 0; i < n_obs; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_P(I_year(i)) + log_K;
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        SIMULATE {
            log_I(i) = rnorm(log_pred_I(i), sd_I(I_survey(i)));
        }
    }

    // Derived quantities
    log_pred_P = log(pred_P);
    log_res_P = log_P - log_pred_P;
    B = P * K;
    log_B = log(B);
    log_res_I = log_I - log_pred_I;

    nll += pen;

    // Report values
    REPORT(pen);
    ADREPORT(log_P);
    ADREPORT(log_pred_P);
    ADREPORT(log_res_P);
    ADREPORT(log_B);
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);

    SIMULATE {
        REPORT(log_P);
        REPORT(log_I);
    }

    return nll;

}


