#include "mkmc.hpp"

static std::mutex s_Mutex;
static std::string name;
static std::string city = "ny";
static double montecarlo0, markov0;
static int termalizacion = 150; //150

static void iteracionesMk(const MobMatrix& T, double p, int tiempo){
    //std::ofstream fileMk("out/Mk_0_p" + name + ".txt");
    std::ofstream fileMk("out/werror/"+city+"/Mk_p" + name + ".txt");
    for(double lambda = 0; lambda <= 3; lambda += 3.0/256){    // 0.01171875
        Mk_DistDiaNocheF markov{0, p, T};
        markov.inicializar(0.001);
        markov.setLambda(lambda);

        for(size_t t = 0; t < tiempo; t++){
            markov.iteracion(T);
        }
        markov.contarInfectados(T);
        fileMk << lambda << " " << markov.totalRhoE/T.Pob << std::endl;
    }
    fileMk.close();
}

static void iteracionesMC(const MobMatrix& T, double &infectados, double& infectadosCuadrado, double lambda, double p, int tiempo, int i){
    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(lambda);
    montecarlo.inicializar(0.002);
    montecarlo.calcInfectados();
    for(size_t t = 0; t < tiempo && montecarlo.pobInf; t++){
        montecarlo.iteracion(T);
    }
    std::lock_guard<std::mutex> lock(s_Mutex);
    infectados += static_cast<double>(montecarlo.pobInf) / T.Pob;//NO TOTALE EN SIS
    infectadosCuadrado += (static_cast<double>(montecarlo.pobInf) / T.Pob) * (static_cast<double>(montecarlo.pobInf) / T.Pob);
}
int mkmc(int argc, char* argv[]){


    double lambda{};
    double p = 1;
    name = "10";
    int tiempo = 200;
    int hilos = 24;

    std::cout << p << " " << city << std::endl;
    MobMatrix T{"cities3/" + city +"/mobnetwork.txt", "cities3/" + city +"/Poparea.txt"};

    std::cout << T.Pob << " cities3/" + city +"/mobnetwork.txt" <<std::endl;

    std::ofstream fileMC("out/werror/"+city+"/MC_p" + name + ".txt");

    auto futureMk = std::async(std::launch::async, iteracionesMk, std::ref(T), p, tiempo);
#if 1
    for(lambda = 0; lambda <= 3; lambda += 3.0/16){
        std::cout << "lambda = " << lambda << std::endl;
        double infectados = 0;
        double infectadosCuadrado  = 0;
        std::vector<std::future<void>> futures;
        for(int i = 0; i < hilos; i++)
            futures.push_back(std::async(std::launch::async, iteracionesMC, std::ref(T), std::ref(infectados), std::ref(infectadosCuadrado), lambda, p, tiempo, i));
        futures.clear();

        fileMC << lambda << " " << infectados/hilos << " " << sqrt((infectadosCuadrado - infectados*infectados/hilos)/(hilos * (hilos - 1))) << std::endl;
        std::cout << lambda << " " << infectados/hilos << " " << sqrt((infectadosCuadrado - infectados*infectados/hilos)/(hilos * (hilos - 1))) << std::endl;
    }
#endif
    futureMk.wait();
    fileMC.close();

    return 0;
}