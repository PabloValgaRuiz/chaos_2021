#pragma once
#include "config.hpp"
#include "MobMatrix.hpp"
#include <vector>

class Mk_Dist{

public:

    Mk_Dist(double _lambda, double _p, const MobMatrix& M);

    virtual void calculaLambda0(const MobMatrix& M) = 0;
    virtual void inicializar(double _rhoInicial) = 0;
    virtual void inicializar() = 0;
    virtual void iteracion(const MobMatrix& T) = 0;
    virtual void contarInfectados(const MobMatrix& T) = 0;

    virtual void calc_neff(const MobMatrix& M);
    virtual void calc_BetaIeff(const MobMatrix& M);
    void setLambda(double _lambda){

    #ifdef SIS
        lambda = _lambda * lambda0;
    #endif

    #ifdef SIR
        lambda = _lambda * lambda0;
    #endif

    #ifdef SEAPIDR
        lambdaP = _lambda * lambda0;
        lambdaI = _lambda * lambda0;
        lambdaA = _lambda * lambda0 / 2;
    #endif
    
    }


public:
    double lambda;  //Probabilidad de infeccion por contacto
    double lambdaP = 0.07, lambdaI = 0.07, lambdaA = 0.035;
    double p;		//Probabilidad de desplazamiento
    double lambda0, mu = 0.2;
    double nu = 1.0/2.6, alpha = 1.0/2.6, delta = 1.0/3.0, gamma = 1.0/14.0, muI = 1.0/4.2, muA = 1.0/6.8, x = 0.35;
    std::vector<std::vector<double>> rhoE, rhoA, rhoP, rhoI, rhoD, rhoR, PI;
    std::vector<double> n_eff, I_eff, probInf;

public:
    double getLambda0() const{return lambda0;}
    const std::vector<std::vector<double>>& getRhoR() const{return rhoR;}

};