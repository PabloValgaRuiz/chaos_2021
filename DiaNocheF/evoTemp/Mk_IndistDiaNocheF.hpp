#pragma once
#include "config.hpp"
#include "MobMatrix.hpp"
#include "Mk_Indist.hpp"

class Mk_IndistDiaNocheF : public Mk_Indist{
public:
    Mk_IndistDiaNocheF(double _lambda, double _p, const MobMatrix& M);
    virtual void calculaLambda0(const MobMatrix& M);
    void inicializar(double _rhoInicial);
    void inicializar(std::vector<double> _rhoInicial);
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