#include "config.hpp"
#include "Mk_IndistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <math.h>
#include <iostream>


Mk_IndistDiaNocheF::Mk_IndistDiaNocheF(double _lambda, double _p, const MobMatrix& M)
    :Mk_IndistDiaNocheF::Mk_Indist(_lambda, _p, M), probInfDia{probInf}{
    PROFILE_FUNCTION();

    I_effN.resize(M.N); probInfNoche.resize(M.N); fvector.resize(M.N);
    calculaLambda0(M);
    calc_neff(M);
    calculaZs(M);
    calculafvector(M);
}


void Mk_IndistDiaNocheF::inicializar(double _rhoInicial){
    PROFILE_FUNCTION();

    for(int i = 0; i < rhoE.size(); i++){
        rhoE[i] = _rhoInicial;
    }
}

void Mk_IndistDiaNocheF::inicializar(std::vector<double> _rhoInicial){
    PROFILE_FUNCTION();

    for(int i = 0; i < rhoE.size(); i++){
        rhoE[i] = _rhoInicial[i];
    }
}


void Mk_IndistDiaNocheF::update(const MobMatrix& T){
    PROFILE_FUNCTION();

    calc_BetaIeff(T);
    calc_BetaIeffN(T);
}


void Mk_IndistDiaNocheF::iteracion(const MobMatrix& T){
    PROFILE_FUNCTION();
    
    update(T);
    
    for(int i = 0; i < T.N; i++){
        if(n_eff[i] != 0)
            probInfDia[i] = 1 - pow((1 - I_eff[i]/n_eff[i]), zD * f(i, T));
        else probInfDia[i] = 0;

        if(T.population[i] != 0)
            probInfNoche[i] = 1 - pow((1 - I_effN[i]/T.population[i]), zN * sigma);
        else probInfNoche[i] = 0;
    }
    

    for(int i = 0; i < T.N; i++){
        PI[i] = (1 - p) * (probInfDia[i] + (1 - probInfDia[i]) * probInfNoche[i]);
    }
    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++){
            if(T.population[i] != 0)
                PI[i] += p * (T.Mpesos[i][j]/T.population[i]) * (probInfDia[T.Mvecinos[i][j]] + (1 - probInfDia[T.Mvecinos[i][j]]) * probInfNoche[i]);
        }
    }

    for(int i = 0; i < T.N; i++){

        double rhoEtemp = rhoE[i];
        double rhoAtemp = rhoA[i];
        double rhoPtemp = rhoP[i];
        double rhoItemp = rhoI[i];
        double rhoDtemp = rhoD[i];
        #ifdef SIS
        rhoE[i] = rhoEtemp * (1 - mu) + (1 - rhoEtemp) * PI[i]; //SIS O SIR
        #endif
        #ifdef SIR
        rhoE[i] = rhoEtemp * (1 - mu) + (1 - rhoEtemp - rhoR[i] ) * PI[i];
        rhoR[i] += rhoEtemp * mu;
        #endif
        #ifdef SEAPIDR
        rhoE[i] = rhoEtemp * (1 - nu) + (1 - rhoEtemp - rhoAtemp - rhoPtemp - rhoItemp - rhoDtemp - rhoR[i] ) * PI[i];
        rhoA[i] = rhoAtemp * (1 - muA) + rhoEtemp * (1 - x) * nu;
        rhoP[i] = rhoPtemp * (1 - alpha) + rhoEtemp * x * nu;
        rhoI[i] = rhoItemp * (1 - delta - muI) + rhoPtemp * alpha;
        rhoD[i] = rhoDtemp * (1 - gamma) + rhoItemp * delta;
        rhoR[i] += rhoAtemp * muA + rhoItemp * muI + rhoDtemp * gamma;
        #endif
            
    }
    
}


void Mk_IndistDiaNocheF::calc_BetaIeffN(const MobMatrix& T){
    PROFILE_FUNCTION();
    #ifdef SIS
    for(int i = 0; i < T.N; i++){
        I_effN[i] = lambda * T.population[i] * rhoE[i];
    }
    #endif
    #ifdef SIR
    for(int i = 0; i < T.N; i++){
        I_effN[i] = lambda * T.population[i] * rhoE[i];
    }
    #endif
    #ifdef SEAPIDR
    for(int i = 0; i < T.N; i++){
        I_effN[i] = T.population[i] * (lambdaA * rhoA[i] + lambdaP * rhoP[i] + lambdaI * rhoI[i]);
    }
    #endif
}

double Mk_IndistDiaNocheF::f(int i, const MobMatrix& T) const{

    if(T.area[i] != 0)
        return n_eff[i]/T.area[i];
    else return 0;
}
double Mk_IndistDiaNocheF::f(int i, double _p, const MobMatrix& T) const{

    if (T.area[i] != 0){
        if(_p == 0){
            return T.population[i]/T.area[i];
        }
        else{
            double tempNeff = 0;
            for(int j = 0; j < T.vecinos[i]; j++){
                tempNeff += T.Mpesos[i][j] * (1 - _p); 
            }
            for(int j = 0; j < T.vecinosT[i]; j++){
                tempNeff += T.MpesosT[i][j] * _p;
            }
            return tempNeff/T.area[i];
        }
    }
    else return 0;
}

void Mk_IndistDiaNocheF::calculafvector(const MobMatrix& T){
    PROFILE_FUNCTION();

    for(int i = 0; i < T.N; i++){
        fvector[i] = f(i, T);
    }
}

void Mk_IndistDiaNocheF::calculaZs(const MobMatrix& M){

    std::vector<double> n_eff_p1; n_eff_p1.resize(M.N);
    for(int i = 0; i < M.N; i++)
        for(int j = 0; j < M.vecinos[i]; j++)
            n_eff_p1[M.Mvecinos[i][j]] += M.Mpesos[i][j];

    double temp1 = 0;
    for(int i = 0; i < M.N; i++){
        temp1 += n_eff_p1[i] * f(i, 1, M);
    }
    zD = M.Pob * kD / temp1;

    temp1 = 0;
    for(int i = 0; i < M.N; i++){
        temp1 += M.population[i] * sigma;
    }
    zN = M.Pob * kN / temp1;
}


void Mk_IndistDiaNocheF::calculaLambda0(const MobMatrix& M){
    PROFILE_FUNCTION();

    double zD = 0, zN = 0;

    std::vector<double> n_eff_p1; n_eff_p1.resize(M.N);
    for(int i = 0; i < M.N; i++)
        for(int j = 0; j < M.vecinos[i]; j++)
            n_eff_p1[M.Mvecinos[i][j]] += M.Mpesos[i][j];

    double temp1 = 0;
    for(int i = 0; i < M.N; i++){
        temp1 += n_eff_p1[i] * f(i, 1, M);
    }
    zD = M.Pob * kD / temp1;

    temp1 = 0;
    for(int i = 0; i < M.N; i++){
        temp1 += M.population[i] * sigma;
    }
    zN = M.Pob * kN / temp1;
    
    //Calcular el maximo de zD*fi + zN*sigma -> maximo de fi si las sigmas son iguales
    double F = 0, SIGMA = 0;
    double temp2 = 0; temp1 = 0;
    for(int i = 0; i < M.N; i++){
        temp1 = f(i, 0, M);
        temp2 = sigma; //sigma[i] si fuese un vector
        if(zD*F + zN*SIGMA < zD*temp1 + zN*temp2) {F = temp1;  SIGMA = temp2;}
    }
#ifdef SIS
    lambda0 = mu / (zD * F + zN * SIGMA);
#endif

#ifdef SIR
    lambda0 = mu / (zD * F + zN * SIGMA);
#endif

    
#ifdef SEAPIDR
    lambda0 = 1 / ( (zD * F + zN * SIGMA) * ( (1 - x)/(2*muA) + x/(delta + muI) + x/alpha ) );
#endif
    
}

void Mk_IndistDiaNocheF::contarInfectados(const MobMatrix& M){
    PROFILE_FUNCTION();

    totalRhoS = totalRhoE = totalRhoA = totalRhoP = totalRhoI = totalRhoD = totalRhoR = 0;
    for(int i = 0; i < M.N; i++){
        totalRhoE += M.population[i] * rhoE[i];
        totalRhoA += M.population[i] * rhoA[i];
        totalRhoP += M.population[i] * rhoP[i];
        totalRhoI += M.population[i] * rhoI[i];
        totalRhoD += M.population[i] * rhoD[i];
        totalRhoR += M.population[i] * rhoR[i];
    }
    totalRhoS = M.Pob - totalRhoE - totalRhoA - totalRhoP - totalRhoI - totalRhoD - totalRhoR;
}