#include "config.hpp"
#include "MobMatrix.hpp"
#include "Mk_IndistDiaNocheF.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <iostream>
#include <fstream>
#include <math.h>


int main(){

    std::string city = "bogota";
    //MobMatrix T{"cities3/" + city + "/mobnetwork.txt", "cities3/" + city + "/Poparea.txt"};
    MobMatrix T{"bogota/mobnetwork.txt", "bogota/Poparea.txt"};
    Mk_IndistDiaNocheF indist{0, 1, T};
    Mk_DistDiaNocheF dist{0, 1, T};

    
    std::vector<double> semilla; semilla.resize(T.N);
    int i;
    std::cin >> i;
    semilla[i] = 0.001;
    int pasos = 100;

    indist.inicializar(semilla);
    dist.inicializar(semilla);

    //ny: 0.0466719
    indist.setLambda(1);
    dist.setLambda(1);

    std::ofstream file{"out/prueba.txt"};
    std::ofstream filetiempo{"out/pruebatiempo.txt"};
    std::ofstream filesimple{"out/pruebasimple.txt"};

    indist.contarInfectados(T);
    dist.contarInfectados(T);

    for(size_t i = 0; i < T.N; ++i){
        file << 0 << "\t" << i << "\t" << indist.rhoE[i] << "\t" << dist.patchRhoE[i] << "\n";
    }

    double minimum = 0.05;
    std::vector<bool> isCountedIndist; isCountedIndist.resize(T.N);
    std::vector<bool> isCountedDist; isCountedDist.resize(T.N);
    std::vector<int> timeIndist; timeIndist.resize(T.N);
    std::vector<int> timeDist; timeDist.resize(T.N);
    for(size_t t = 1; t <= pasos; t++){
        indist.iteracion(T);
        dist.iteracion(T);
        indist.contarInfectados(T);
        dist.contarInfectados(T);
        for(size_t i = 0; i < T.N; ++i){
            file << t << "\t" << i << "\t" << indist.rhoE[i] << "\t" << dist.patchRhoE[i] << "\n";

            if(indist.rhoE[i] > minimum && !(isCountedIndist[i]) ){
                timeIndist[i] = t;
                isCountedIndist[i] = true;
            }
            if(dist.patchRhoE[i] > minimum && !(isCountedDist[i]) ){
                timeDist[i] = t;
                isCountedDist[i] = true;
            }
        }
        filesimple << t << "\t" << indist.totalRhoE/T.Pob << "\t" << dist.totalRhoE/T.Pob << "\n";
    }
    for(size_t i = 0; i < T.N; ++i){
        if(timeIndist[i] == 0) timeIndist[i] = pasos;
        if(timeDist[i] == 0) timeDist[i] = pasos;
        filetiempo << i << "\t" << timeIndist[i] << "\t" << timeDist[i] << "\n";

    }
    // for(size_t i = 0; i < T.N; ++i){
    //  file << 0 << "\t" << i << "\t" << indist.rhoE[i] << "\t" << dist.patchRhoE[i] << "\n";
    // }
    // std::cout << 0 << "\t" << indist.totalRhoE << "\t" << dist.totalRhoE << std::endl;

    // for(size_t t = 1; t <= pasos; t++){
    //     indist.iteracion(T);
    //     indist.contarInfectados(T);
    //     dist.iteracion(T);
    //     dist.contarInfectados(T);
    //     for(size_t i = 0; i < T.N; ++i){
    //         file << t << "\t" << i << "\t" << indist.rhoE[i] << "\t" << dist.patchRhoE[i] << "\n";
    //     }
    //     std::cout << t << "\t" << indist.totalRhoE << "\t" << dist.totalRhoE << std::endl;
    // }

    file.close();
    filetiempo.close();
    filesimple.close();

    return 0;
}