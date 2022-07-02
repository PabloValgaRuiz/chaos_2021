#include "config.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <math.h>
#include <iostream>


Mk_DistDiaNocheF::Mk_DistDiaNocheF(double _lambda, double _p, const MobMatrix& M)
    :Mk_DistDiaNocheF::Mk_Dist(_lambda, _p, M), probInfDia{probInf}{
    PROFILE_FUNCTION();

    I_effN.resize(M.N); probInfNoche.resize(M.N); fvector.resize(M.N);
    calculaLambda0(M);
    calc_neff(M);
    calculaZs(M);
    calculafvector(M);
}


void Mk_DistDiaNocheF::inicializar(double _rhoInicial){
    PROFILE_FUNCTION();

    for(int i = 0; i < rhoE.size(); i++){
        for(int j = 0; j < rhoE[i].size(); j++){
            rhoE[i][j] = _rhoInicial;
        }
    }
}


void Mk_DistDiaNocheF::update(const MobMatrix& T){
    PROFILE_FUNCTION();

    calc_BetaIeff(T);
    calc_BetaIeffN(T);
}


void Mk_DistDiaNocheF::iteracion(const MobMatrix& T){
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
        for(int j = 0; j < T.vecinos[i]; j++){
            PI[i][j] = 0;
        }
    }
    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++){
            PI[i][j] = (1 - p) * (probInfDia[i] + (1 - probInfDia[i]) * probInfNoche[i]) + p * (probInfDia[T.Mvecinos[i][j]] +
                                                                    (1 - probInfDia[T.Mvecinos[i][j]]) * probInfNoche[i]);
        }
    }

    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++){
            double rhoEtemp = rhoE[i][j];
            double rhoAtemp = rhoA[i][j];
            double rhoPtemp = rhoP[i][j];
            double rhoItemp = rhoI[i][j];
            double rhoDtemp = rhoD[i][j];
            #ifdef SIS
            rhoE[i][j] = rhoEtemp * (1 - mu) + (1 - rhoEtemp) * PI[i][j]; //SIS O SIR
            #endif
            #ifdef SIR
            rhoE[i][j] = rhoEtemp * (1 - mu) + (1 - rhoEtemp - rhoR[i][j] ) * PI[i][j];
            rhoR[i][j] += rhoEtemp * mu;
            #endif
            #ifdef SEAPIDR
            rhoE[i][j] = rhoEtemp * (1 - nu) + (1 - rhoEtemp - rhoAtemp - rhoPtemp - rhoItemp - rhoDtemp - rhoR[i][j] ) * PI[i][j];
            rhoA[i][j] = rhoAtemp * (1 - muA) + rhoEtemp * (1 - x) * nu;
            rhoP[i][j] = rhoPtemp * (1 - alpha) + rhoEtemp * x * nu;
            rhoI[i][j] = rhoItemp * (1 - delta - muI) + rhoPtemp * alpha;
            rhoD[i][j] = rhoDtemp * (1 - gamma) + rhoItemp * delta;
            rhoR[i][j] += rhoAtemp * muA + rhoItemp * muI + rhoDtemp * gamma;
            #endif
            
        }
    }
    
}


void Mk_DistDiaNocheF::calc_BetaIeffN(const MobMatrix& T){
    PROFILE_FUNCTION();
    #ifdef SIS
    for(int i = 0; i < T.N; i++){
        I_effN[i] = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            I_effN[i] += lambda * T.Mpesos[i][j] * rhoE[i][j];
        }
    }
    #endif
    #ifdef SIR
    for(int i = 0; i < T.N; i++){
        I_effN[i] = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            I_effN[i] += lambda * T.Mpesos[i][j] * rhoE[i][j];
        }
    }
    #endif
    #ifdef SEAPIDR
    for(int i = 0; i < T.N; i++){
        I_effN[i] = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            I_effN[i] += T.Mpesos[i][j] * (lambdaA * rhoA[i][j] + lambdaP * rhoP[i][j] + lambdaI * rhoI[i][j]);
        }
    }
    #endif
}

double Mk_DistDiaNocheF::f(int i, const MobMatrix& T) const{

    if(T.area[i] != 0)
        return n_eff[i]/T.area[i];
    else return 0;
}
double Mk_DistDiaNocheF::f(int i, double _p, const MobMatrix& T) const{

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

void Mk_DistDiaNocheF::calculafvector(const MobMatrix& T){
    PROFILE_FUNCTION();

    for(int i = 0; i < T.N; i++){
        fvector[i] = f(i, T);
    }
}

void Mk_DistDiaNocheF::calculaZs(const MobMatrix& M){

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


void Mk_DistDiaNocheF::calculaLambda0(const MobMatrix& M){
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

void Mk_DistDiaNocheF::contarInfectados(const MobMatrix& M){
    PROFILE_FUNCTION();

    totalRhoS = totalRhoE = totalRhoA = totalRhoP = totalRhoI = totalRhoD = totalRhoR = 0;
    for(int i = 0; i < M.N; i++){
        for(int j = 0; j < M.vecinos[i]; j++){
            totalRhoE += M.Mpesos[i][j] * rhoE[i][j];
            totalRhoA += M.Mpesos[i][j] * rhoA[i][j];
            totalRhoP += M.Mpesos[i][j] * rhoP[i][j];
            totalRhoI += M.Mpesos[i][j] * rhoI[i][j];
            totalRhoD += M.Mpesos[i][j] * rhoD[i][j];
            totalRhoR += M.Mpesos[i][j] * rhoR[i][j];
        }
    }
    totalRhoS = M.Pob - totalRhoE - totalRhoA - totalRhoP - totalRhoI - totalRhoD - totalRhoR;
}