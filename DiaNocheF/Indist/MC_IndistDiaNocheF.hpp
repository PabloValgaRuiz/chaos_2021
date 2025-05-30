#pragma once
#include "config.hpp"
#include "MobMatrix.hpp"
#include "MC_Indist.hpp"

class MC_IndistDiaNocheF : public MC_Indist{

public:
    virtual void iteracion(const MobMatrix& T);
    virtual void update(const MobMatrix& T);

    MC_IndistDiaNocheF(double _lambda, double _p, const MobMatrix& M);
    void calculaLambda0(const MobMatrix& T);
    void generafvector(const MobMatrix& T);

protected:
    std::vector<double> probInfDes, probInfOrg, fvector;
    double kD = 8, kN = 3, sigma = 5;
    double zD, zN;

    double f(int i, double _p, const MobMatrix &T) const;
    double f(int i, const MobMatrix &T) const;
    void calculateZs(const MobMatrix &T);
};