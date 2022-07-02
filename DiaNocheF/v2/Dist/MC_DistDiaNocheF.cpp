#include "config.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <chrono>
#include <random>
#include <math.h>

const auto time0 = std::chrono::steady_clock::now();

MC_DistDiaNocheF::MC_DistDiaNocheF(double _lambda, double _p, const MobMatrix& M)
: MC_DistDiaNocheF::MC_Dist{_lambda,_p,M}
{
    PROFILE_FUNCTION();
    
    calculaLambda0(M);
    probInfDes.resize(M.N); probInfOrg.resize(M.N); fvector.resize(M.N);
}

void MC_DistDiaNocheF::update(const MobMatrix& T){
    PROFILE_FUNCTION();

    this->desplaza();
    this->calc_neff();
    generafvector(T);
    this->calculateZs(T);
    this->calcInfectados();
}

void MC_DistDiaNocheF::iteracion(const MobMatrix& T){
    PROFILE_FUNCTION();

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    this->update(T);

    std::vector<int> contactosD(T.N);
    std::vector<int> contactosN(T.N);
    std::vector<double> probContactoD(T.N);
    std::vector<double> probContactoN(T.N);
    std::vector<std::vector<size_t>> positionsMC(T.N); for(int i = 0; i < T.N; i++) positionsMC[i].resize(n_eff[i]);

    for(int i = 0; i < T.N; i++){
        contactosD[i] = static_cast<int>(zD * fvector[i]);
        contactosN[i] = static_cast<int>(zN * sigma);
        probContactoD[i] = zD * fvector[i] - static_cast<double>(contactosD[i]);
        probContactoN[i] = zN * sigma - static_cast<double>(contactosN[i]);
    }

{   
    InstrumentationTimer timer("Iteracion sobre individuos");
    //Contagios por el dia
    for(int k = 0; k < Est.size(); k++){
        switch(Est[k]){
            case S:
                int contactos;
                if(dist(mt) < probContactoD[Desplazamiento[k]])
                    contactos = contactosD[Desplazamiento[k]] + 1;
                else contactos = contactosD[Desplazamiento[k]];
                int pop_eff = n_eff[Desplazamiento[k]] - 1;
                int infected = InfectadosDes[Desplazamiento[k]];

                if(contactos >= pop_eff){
                    if(dist(mt) < 1 - pow(1 - lambda, infected)){
                        Est[k] = E;
                        goto getOut;
                    }
                }
                else{
                    for(int i = 0; i < contactos; i++){
                        if(dist(mt) < (static_cast<double>(infected) / pop_eff)){
                            if(dist(mt) < lambda){
                                Est[k] = E;
                                goto getOut;
                            }
                            else{
                                infected--;
                                pop_eff--;
                            }
                        }
                        else{
                            pop_eff--;
                        }
                    }
                }
                if(dist(mt) < probContactoN[Org[k]])
                    contactos = contactosN[Org[k]] + 1;
                else contactos = contactosN[Org[k]];
                pop_eff = T.population[Org[k]] - 1;
                infected = InfectadosOrg[Org[k]];
                if(contactos >= pop_eff){
                    if(dist(mt) < 1 - pow(1 - lambda, infected)){
                        Est[k] = E;
                        goto getOut;
                    }
                }
                else{
                    for(int i = 0; i < contactos; i++){
                        if(dist(mt) < static_cast<double>(infected) / pop_eff){
                            if(dist(mt) < lambda){
                                Est[k] = E;
                                goto getOut;
                            }
                            else{
                                infected--;
                                pop_eff--;
                            }
                        }
                        else{
                            pop_eff--;
                        }
                    }
                }

                getOut:
                break;
            
            case E:
#ifdef SIS
                if(dist(mt) < mu)
                    Est[k] = S;
#endif
#ifdef SIR
                if(dist(mt) < mu)
                    Est[k] = R;
#endif
#ifdef SEAPIDR
                double temp = dist(mt);
                if(temp < (1-x)*nu)
                    Est[k] = A;
                else if(temp < nu)
                    Est[k] = P;
#endif
                break;
#ifdef SEAPIDR
            case A:
                if(dist(mt) < muA)
                    Est[k] = R;
                break;
            
            case P:
                if(dist(mt) < alpha)
                    Est[k] = I;
                break;
            
            case I:
                temp = dist(mt);
                if(temp < delta)
                    Est[k] = D;
                else if(temp < delta + muI)
                    Est[k] = R;
                break;
            case D:
                if(dist(mt) < gamma)
                    Est[k] = R;
#endif
        }
    }
}
    
}


void MC_DistDiaNocheF::calculaLambda0(const MobMatrix& T){
    PROFILE_FUNCTION();

    double zD = 0, zN = 0;

    std::vector<double> n_eff_p1; n_eff_p1.resize(T.N);
    for(int i = 0; i < T.N; i++)
        for(int j = 0; j < T.vecinos[i]; j++)
            n_eff_p1[T.Mvecinos[i][j]] += T.Mpesos[i][j];

    double temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 += n_eff_p1[i] * f(i, 1, T);
    }
    zD = T.Pob * kD / temp1;

    temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 += T.population[i] * sigma;
    }
    zN = T.Pob * kN / temp1;
    
    //Calcular el maximo de zD*fi + zN*sigma -> maximo de fi si las sigmas son iguales
    double F = 0, SIGMA = 0;
    double temp2 = 0; temp1 = 0;
    for(int i = 0; i < T.N; i++){
        temp1 = f(i, 0, T);
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

double MC_DistDiaNocheF::f(int i, double _p, const MobMatrix& T) const{

    double result = 0;  //Calcular la poblacion efectiva neff
    for(int j = 0; j < T.vecinosT[i]; j++){
        result += T.MpesosT[i][j] * _p;
    }
    for(int j = 0; j < T.vecinos[i]; j++){
        result += T.Mpesos[i][j] * (1 - _p);
    }
    if(T.area[i] != 0)
        return result/T.area[i];
    else return 0;
}

double MC_DistDiaNocheF::f(int i, const MobMatrix& T) const{
    return MC_DistDiaNocheF::f(i, this->p, T);
}

void MC_DistDiaNocheF::generafvector(const MobMatrix& T){
    PROFILE_FUNCTION();

    for(int i = 0; i < T.N; i++){
        fvector[i] = f(i, T);
    }
}

void MC_DistDiaNocheF::calculateZs(const MobMatrix&T){
    PROFILE_FUNCTION();

    //Calcular las z en funcion de p=1

    zD = zN = 0;
    double prov1 = 0, prov2 = 0;

    std::vector<double> n_eff_p1; n_eff_p1.resize(T.N);
    for(int i = 0; i < T.N; i++)
        for(int j = 0; j < T.vecinos[i]; j++)
            n_eff_p1[T.Mvecinos[i][j]] += T.Mpesos[i][j];


    for(int i = 0; i < T.N; i++){
        prov1 += n_eff_p1[i] * f(i, 1, T);
        prov2 += T.population[i] * sigma;
    }
    zD = T.Pob * kD / prov1;
    zN = T.Pob * kN / prov2;
}