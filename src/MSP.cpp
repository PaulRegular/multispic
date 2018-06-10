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
    PARAMETER(log_sd_L);

    // Transformations
    vector<Type> log_L = log(L);
    vector<Type> log_I = log(I);
    vector<Type> B = exp(log_B);
    Type sd_B = exp(log_sd_B);
    Type K = exp(log_K);
    Type r = exp(log_r);
    vector<Type> sd_I = exp(log_sd_I);
    Type sd_L = exp(log_sd_L);

    // Process equation
    Type pen = Type(0);
    Type nll = Type(0);

    int n_years = L.size();
    vector<Type> log_pred_L(n_years);
    vector<Type> pred_L = exp(log_pred_L);
    vector<Type> log_pred_B(n_years);
    vector<Type> pred_B = exp(log_pred_B);
    vector<Type> log_res_B = log_B - log_pred_B;
    for (int i = 1; i < n_years; i++){
        pred_B(i) = B(i - 1) + r * B(i - 1) * (Type(1) - B(i - 1) / K) - pred_L(i - 1);
        pred_B(i) = posfun(pred_B(i), Type(1.0e-6), pen);
        nll -= dnorm(log_B(i), log_pred_B(i), sd_B, true);
    }

    // Observation equations
    vector<Type> log_res_L = log_L - log_pred_L;
    for (int i = 1; i < n_years; i++){
        nll -= dnorm(log_L(i), log_pred_L(i), sd_L, true);

    int n_obs = I.size();
    vector<Type> log_pred_I(n_obs);
    vector<Type> log_res_I = log_I - log_pred_I;
    for (int i = 0; i < n_obs; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_B(I_year(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
    }

    nll += pen;

    // Report values
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);
    ADREPORT(log_pred_L);
    ADREPORT(log_res_L);
    ADREPORT(log_B);
    ADREPORT(log_pred_B);
    ADREPORT(log_res_B);
    REPORT(pen);

    return nll;

}
