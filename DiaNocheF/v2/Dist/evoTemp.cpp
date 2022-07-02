#include "evoTemp.hpp"

int evoTemp(){

    std::string city = "ny";

    MobMatrix T{"cities3/" + city + "/mobnetwork.txt", "cities3/" + city + "/Poparea.txt"};

    MC_DistDiaNocheF montecarlo{0, 1, T};
    Mk_DistDiaNocheF markov{0,1,T};

    montecarlo.inicializar(0.001);
    markov.inicializar(0.001);

    montecarlo.setLambda(4 * 0.0466719);
    markov.setLambda(4 * 0.0466719);

    montecarlo.calcInfectados();
    markov.contarInfectados(T);
    std::cout << 0 << "\t" << montecarlo.pobInf << "\t" << markov.totalRhoE << std::endl;

    for(size_t t = 1; t <= 100; t++){
        montecarlo.iteracion(T);
        montecarlo.calcInfectados();
        markov.iteracion(T);
        markov.contarInfectados(T);
        std::cout << t << "\t" << montecarlo.pobInf << "\t" << markov.totalRhoE << std::endl;
    }

    return 0;
}