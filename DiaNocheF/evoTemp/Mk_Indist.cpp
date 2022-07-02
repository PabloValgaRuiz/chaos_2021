#include "config.hpp"
#include "Mk_Indist.hpp"
#include "benchmark.hpp"


Mk_Indist::Mk_Indist(double _lambda, double _p, const MobMatrix& M)
:lambda{_lambda}, p{_p}{
    PROFILE_FUNCTION();

    n_eff.resize(M.N); I_eff.resize(M.N); probInf.resize(M.N);
    rhoE.resize(M.N); rhoA.resize(M.N); rhoP.resize(M.N); rhoI.resize(M.N); rhoD.resize(M.N);  rhoR.resize(M.N); PI.resize(M.N);
}


void Mk_Indist::calc_neff(const MobMatrix& M){


    for(int i = 0; i < M.N; i++)
        n_eff[i] = 0;

    for(int i = 0; i < M.N; i++){
        for(int j = 0; j < M.vecinos[i]; j++){
            n_eff[i] += (1 - p) * M.Mpesos[i][j];
            n_eff[M.Mvecinos[i][j]] += p * M.Mpesos[i][j];
        }
    }
}

void Mk_Indist::calc_BetaIeff(const MobMatrix& M){
    PROFILE_FUNCTION();

    for(int i = 0; i < M.N; i++)
        I_eff[i] = 0;
    for(int i = 0; i < M.N; i++){

        for(int j = 0; j < M.vecinos[i]; j++){

        #ifdef SIS                                                      //Delta de kronecker (delta_ij)
            I_eff[i] += lambda * rhoE[i] * (1 - p) * M.population[i] * (i == M.Mvecinos[i][j]);
            I_eff[M.Mvecinos[i][j]] += lambda * rhoE[i] * p * M.Mpesos[i][j];
        #endif

        #ifdef SIR
            I_eff[i] += lambda * rhoE[i] * (1 - p) * M.population[i] * (i == M.Mvecinos[i][j]);
            I_eff[M.Mvecinos[i][j]] += lambda * rhoE[i] * p * M.Mpesos[i][j];
        #endif

        #ifdef SEAPIDR
            I_eff[i] += (lambdaA * rhoA[i] + lambdaP * rhoP[i] + lambdaI * rhoI[i]) * (1 - p) * M.population[i] * (i == M.Mvecinos[i][j]);
            I_eff[M.Mvecinos[i][j]] += (lambdaA * rhoA[i] + lambdaP * rhoP[i] + lambdaI * rhoI[i]) * p * M.Mpesos[i][j];
        #endif
        }
    }
    
}