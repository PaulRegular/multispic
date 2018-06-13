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
    PARAMETER_VECTOR(log_P);
    PARAMETER_VECTOR(log_sd_P);
    PARAMETER_VECTOR(log_K);
    PARAMETER_VECTOR(log_r);
    PARAMETER_VECTOR(log_m);
    PARAMETER_VECTOR(log_gamma); // year effect parameter
    PARAMETER(log_sd_gamma);
    PARAMETER_VECTOR(log_q);
    PARAMETER_VECTOR(log_sd_I);

    // Transformations
    vector<Type> log_I = log(I);
    vector<Type> P = exp(log_P);
    vector<Type> sd_P = exp(log_sd_P);
    vector<Type> K = exp(log_K);
    vector<Type> r = exp(log_r);
    vector<Type> m = exp(log_m);
    Type sd_gamma = exp(log_sd_gamma);
    vector<Type> sd_I = exp(log_sd_I);

    // Containers
    int nL = L.size();
    vector<Type> pred_P(nL);
    vector<Type> log_pred_P(nL);
    vector<Type> log_res_P(nL);
    vector<Type> B(nL);
    vector<Type> log_B(nL);
    int nI = I.size();
    vector<Type> log_pred_I(nI);
    vector<Type> log_res_I(nI);

    // Initalize nll
    Type pen = Type(0);
    Type nll = Type(0);

    // Process equation
    for (int i = 0; i < nL; i++){
        nll -= dnorm(log_gamma(i), Type(0), sd_gamma);
        if (L_year(i) == 0) {
            nll -= dnorm(log_P(i), Type(0), sd_P(L_species(i)), true);
            SIMULATE {
                log_P(i) = rnorm(Type(0), sd_P(L_species(i)));
            }
        } else {
            pred_P(i) = (P(i - 1) +
                (r(L_species(i)) / (m(L_species(i)) - 1)) *
                (Type(1) - pow(P(i - 1), m(L_species(i)) - Type(1))) -
                (L(i - 1) / K(L_species(i)))) *
                exp(log_gamma(i));
            pred_P(i) = pos_fun(pred_P(i), min_P, pen);
            nll -= dnorm(log_P(i), log(pred_P(i)), sd_P(L_species(i)), true);
            SIMULATE {
                log_P(i) = rnorm(log(pred_P(i)), sd_P(L_species(i)));
            }
        }
        B(i) = P(i) * K(L_species(i));
    }

    // Observation equations
    for (int i = 0; i < nI; i++){
        log_pred_I(i) = log_q(I_survey(i)) + log_P(I_sy(i)) + log_K(I_species(i));
        nll -= dnorm(log_I(i), log_pred_I(i), sd_I(I_survey(i)), true);
        SIMULATE {
            log_I(i) = rnorm(log_pred_I(i), sd_I(I_survey(i)));
        }
    }

    // More transformations
    log_B = log(B);
    log_pred_P = log(pred_P);

    // Residuals
    log_res_P = log_P - log_pred_P;
    log_res_I = log_I - log_pred_I;

    // AD report values
    ADREPORT(log_P);
    ADREPORT(log_pred_P);
    ADREPORT(log_res_P);
    ADREPORT(log_B);
    ADREPORT(log_gamma);
    ADREPORT(log_pred_I);
    ADREPORT(log_res_I);

    // Report simulated values
    SIMULATE {
        REPORT(log_P);
        REPORT(log_I);
    }

    REPORT(pen);
    nll += pen;
    return nll;

}

