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

    // Parameters
    PARAMETER_VECTOR(log_B);
    PARAMETER(log_sd_B);
    PARAMETER(log_K);
    PARAMETER(log_r);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_sd_I);

    // Transformations
    vector<Type> log_I = log(I);
    vector<Type> B = exp(log_B);
    Type sd_B = exp(log_sd_B);
    Type K = exp(log_K);
    Type r = exp(log_r);
    vector<Type> sd_I = exp(log_sd_I);

    // Process equation
    Type pen = Type(0);
    Type nll = Type(0);

    int n_years = L.size();
    vector<Type> pred_L(n_years);
    vector<Type> pred_B(n_years);
    for (int i = 1; i < n_years; i++){
        pred_B(i) = B(i - 1) + r * B(i - 1) * (Type(1) - B(i - 1) / K) - L(i - 1);
        pred_B(i) = posfun(pred_B(i), Type(1.0e-6), pen);
        nll -= dnorm(log_B(i), log(pred_B(i)), sd_B, true);
    }
    vector<Type> log_pred_B = log(pred_B);
    vector<Type> log_res_B = log_B - log_pred_B;

    // Observation equations
    int n_obs = I.size();
    vector<Type> log_pred_I(n_obs);
    for (int i = 0; i < n_obs; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_B(I_year(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
    }
    vector<Type> log_res_I = log_I - log_pred_I;

    nll += pen;

    // Report values
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);
    ADREPORT(log_B);
    ADREPORT(log_pred_B);
    ADREPORT(log_res_B);
    REPORT(pen);

    return nll;

}


