#include "config.hpp"
#include "MobMatrix.hpp"
#include "MC_IndistDiaNocheF.hpp"
#include "Mk_IndistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <future>


static std::string name;
static std::string city = "0";


int simu(int argc, char* argv[]){
    double lambda = 4;
    double p = 0.2;
    name = "02";
    int tiempo = 1;//NO SE UTILIZA
    int hilos = 20;
    std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
    if(argc == 2){name = argv[1]; p = std::stoi(name)/10.0;}

    std::cout << p << " " << name << std::endl;
    std::cout << "cities2/oldnetwork"+city+".txt" << std::endl;
    MobMatrix T{"cities2/oldnetwork"+city+".txt", "areas/area"+city+".txt"};
    std::cout << T.Pob <<std::endl;

    MC_IndistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(lambda);
    montecarlo.inicializar(0.002);

    Mk_IndistDiaNocheF markov{0, p, T};
    markov.inicializar(0.002);
    markov.setLambda(lambda);


    //std::ofstream fileMC("out/MC_00_p" + name + ".txt");
    std::ofstream fileMC("out/simu"+city+"_p" + name + "(5).txt");
    int t = 0; double markovInfect, montecarloInfect;
    do{
        markov.iteracion(T);   
        markov.contarInfectados(T);

        montecarlo.iteracion(T);
        montecarlo.calcInfectados();
        markovInfect = markov.totalRhoE + markov.totalRhoA + markov.totalRhoP + markov.totalRhoI + markov.totalRhoD;
        montecarloInfect = montecarlo.totalE + montecarlo.totalA + montecarlo.totalP + montecarlo.totalI + montecarlo.totalD;
        fileMC << t++ << "\t" << markovInfect/T.Pob << "\t" << montecarloInfect/T.Pob << "\t" << markov.totalRhoR/T.Pob << "\t" << (double)montecarlo.totalR/T.Pob << std::endl;
        std::cout<<t<<"\t"<< markovInfect/T.Pob<<"\t" << montecarloInfect/T.Pob << "\t" << markov.totalRhoR/T.Pob << "\t" << (double)montecarlo.totalR/T.Pob << std::endl;
    }while( markovInfect > 1.0 && montecarloInfect > 0);
    
    fileMC.close();

    return 0;
}