#include "config.hpp"
#include "MobMatrix.hpp"
#include "MC_IndistDiaNocheF.hpp"
#include "Mk_IndistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <future>



static std::mutex s_Mutex;
static std::string name;
static std::string city = "0";
double montecarlo0, markov0;
double ceros = 0;


static void iteracionesMk(const MobMatrix& T, double lambda, double p, int tiempo){
    //std::ofstream fileMk("out/Mk_0_p" + name + ".txt");
    std::ofstream fileMk("out/Mk_"+city+"_p" + name + ".txt");
    for(double lambda = 0; lambda <= 5; lambda += 5.0/256){    // 0.02
        Mk_IndistDiaNocheF markov{0, p, T};
        markov.inicializar(0.002);
        markov.setLambda(lambda);
        do{
            markov.iteracion(T);
            markov.contarInfectados(T);
        }while((markov.totalRhoE + markov.totalRhoA + markov.totalRhoP + markov.totalRhoI + markov.totalRhoD)/T.Pob > 0.0000001);
        markov.contarInfectados(T);
        fileMk << lambda << " " << markov.totalRhoR/T.Pob << std::endl;
    }
    fileMk.close();
}

static void iteracionesMC(const MobMatrix& T, double &infectados, double lambda, double p, int tiempo, int i){
    MC_IndistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(lambda);
    montecarlo.inicializar(0.002);

    //NO VALIDO PARA SIS
    do{
        montecarlo.iteracion(T);
        montecarlo.calcInfectados();
        std::cout << montecarlo.totalE + montecarlo.totalA + montecarlo.totalP + montecarlo.totalI + montecarlo.totalD << std::endl;
    }while(montecarlo.totalE + montecarlo.totalA + montecarlo.totalP + montecarlo.totalI + montecarlo.totalD != 0);
    std::lock_guard<std::mutex> lock(s_Mutex);
    infectados += montecarlo.totalR;
}
int mkmc(int argc, char* argv[]){
    double lambda = 2;
    double p = 0.2;
    name = "02";
    int tiempo = 1;//NO SE UTILIZA
    int hilos = 20;

    if(argc == 2){name = argv[1]; p = std::stoi(name)/10.0;}

    std::cout << p << " " << name << std::endl;
    std::cout << "cities2/oldnetwork"+city+".txt" << std::endl;
    MobMatrix T{"cities2/oldnetwork"+city+".txt", "areas/area"+city+".txt"};
    std::cout << T.Pob <<std::endl;

    std::ofstream fileMC("out/MC_"+city+"_p" + name + ".txt");

    auto futureMk = std::async(std::launch::async, iteracionesMk, std::ref(T), lambda, p, tiempo);
    for(lambda = 2.5 + 1*5.0/16 ; lambda <= 5; lambda += 5.0/16){        
        double infectados = 0;
        std::vector<std::future<void>> futures;
        for(int i = 0; i < hilos; i++)
            futures.push_back(std::async(std::launch::async, iteracionesMC, std::ref(T), std::ref(infectados), lambda, p, tiempo, i));
        futures.clear();

        fileMC << lambda << " " << infectados/(hilos*T.Pob) << std::endl;
    }
    futureMk.wait();
    fileMC.close();

    return 0;
}