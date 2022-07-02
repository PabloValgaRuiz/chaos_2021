#pragma once
#include "config.hpp"
#include "MobMatrix.hpp"
#include "Mk_Dist.hpp"

class Mk_DistDiaNocheF : public Mk_Dist{
public:
    Mk_DistDiaNocheF(double _lambda, double _p, const MobMatrix& M);
    virtual void calculaLambda0(const MobMatrix& M);
    void inicializar(double _rhoInicial);
    void inicializar(){
        rhoE[0][0] = 0.00001;
    }
    virtual void update(const MobMatrix& T);
    virtual void iteracion(const MobMatrix& T);
    
    double totalRhoS, totalRhoE, totalRhoA, totalRhoP, totalRhoI, totalRhoD, totalRhoR;
    virtual void contarInfectados(const MobMatrix& M);

public:
    double kD = 8, kN = 3, sigma = 5;
    double zD, zN;
    std::vector<double>& probInfDia;
    std::vector<double> probInfNoche;
    std::vector<double> I_effN;
    std::vector<double> fvector;
    double f(int i, double _p, const MobMatrix& T) const;
    double f(int i, const MobMatrix& T) const;
    void calc_BetaIeffN(const MobMatrix& T);
    void calculaZs(const MobMatrix& M);
    void calculafvector(const MobMatrix& T);
};