#include <precompiled.hpp>
#include "config.hpp"
#include "MobMatrix.hpp"
#include "MergeMatrix.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <math.h>
#include <future>
#include <unistd.h>
#include <tuple>

#define zcta


static std::mutex s_Mutex;
static std::mutex threads_Mutex;
static int threads;
static std::string name;
static double LAMBDA0;

std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const Mk_DistDiaNocheF& markov);

static void iteracionesMk(const MobMatrix& T, double p, std::string name, std::string resol){
    //std::ofstream fileMk("out/Mk_0_p" + name + ".txt");
    std::ofstream fileMk("out/Mk_" + resol + name + ".txt");
    for(double lambda = 0; lambda <= 16; lambda += 16.0/128){    // 0.02
        Mk_DistDiaNocheF markov{0, p, T};
        markov.inicializar(0.002);
        double lambda0 = markov.getLambda0();
        markov.setLambda(lambda*LAMBDA0/lambda0);//0.00274911 es el lambda0 de zcta
        do{
            markov.iteracion(T);
            markov.contarInfectados(T);
        }while((markov.totalRhoE + markov.totalRhoA + markov.totalRhoP + markov.totalRhoI + markov.totalRhoD)/T.Pob > 0.0000001);
        markov.contarInfectados(T);
        fileMk << lambda << " " << markov.totalRhoR/T.Pob << std::endl;
    }
    fileMk.close();
}
static int diagramasfases(){
    double p;
    std::vector<std::future<void>> futures;
    MobMatrix T1{"cities/zctanetwork.txt", "cities/zcta_areapop.txt"};
    MobMatrix T2{"cities/subctynetwork.txt", "cities/Subcounties_areapop.txt"};
    MobMatrix T3{"cities/ctynetwork.txt", "cities/Counties_areapop.txt"};

    for(std::string name : {"00","02","04","06","08","10"}){
        p = std::stoi(name)/10.0;
        futures.push_back(std::async(std::launch::async, iteracionesMk, T1, p, name, "zcta"));
        futures.push_back(std::async(std::launch::async, iteracionesMk, T2, p, name, "subcty"));
        futures.push_back(std::async(std::launch::async, iteracionesMk, T3, p, name, "cty"));
    }
    futures.clear();
    return 0;
}

static std::vector<std::vector<double>> markovres(const MobMatrix& T, double p, double lambda){

    Mk_DistDiaNocheF markov{0, p, T};
    markov.inicializar(0.001);

    double lambda0 = markov.getLambda0();
    markov.setLambda(lambda*LAMBDA0/lambda0);
    do{
        markov.iteracion(T);
        markov.contarInfectados(T);
    }while((markov.totalRhoE + markov.totalRhoA + markov.totalRhoP + markov.totalRhoI + markov.totalRhoD)/T.Pob > 0.000001);
    markov.contarInfectados(T);
    return markov.getRhoR();
}

void iteracionHeatMap(const MobMatrix& T, double lambda, double p, int tiempo, std::vector<std::vector<double>>& heatMap, int i, int j){
{
    Mk_DistDiaNocheF markov{0, p, T};
    markov.inicializar(0.0000001);
    double lambda0 = markov.getLambda0();
    markov.setLambda(lambda*LAMBDA0/lambda0);
    do{
        markov.iteracion(T);
        markov.contarInfectados(T);
    }while((markov.totalRhoE + markov.totalRhoA + markov.totalRhoP + markov.totalRhoI + markov.totalRhoD)/T.Pob > 0.00000002);

    markov.contarInfectados(T);

    //Cada hilo cambia un par (i,j) distinto, no es necesario un mutex
    heatMap[i][j] = markov.totalRhoR/T.Pob;
}//Cerrar el scope antes de disminuir el valor de threads, para que se libere la memoria de las variables locales antes de crear el siguiente
    std::lock_guard<std::mutex> lock(threads_Mutex);
    threads--;
}
static int heatMap(){
    const double LAMBDA_MIN = 0.96875;
    const double LAMBDA_MAX = 1.03125;
    const double P_MIN = 0;
    const double P_MAX = 1;
    
    const int MAX_THREADS = 24;

    MobMatrix T{"cities/zctanetwork.txt", "cities/zcta_areapop.txt"};

    std::vector<std::vector<double>> heatMap;
    heatMap.resize(65); //Numero de p's //2^7 + 1 = 129
    for(auto& i : heatMap){i.resize(65);} //Numero de lambdas

    std::ofstream fileMk("out/heatmap_zcta65.txt");
    std::vector<std::future<void>> futures;
    double i = 0;
    for(double p = P_MIN; p <= P_MAX; p += (P_MAX - P_MIN)/(heatMap.size() - 1)){

        double j = 0;
        for(double lambda = LAMBDA_MIN; lambda <= LAMBDA_MAX; lambda += (LAMBDA_MAX - LAMBDA_MIN)/(heatMap[i].size() - 1)){
            std::cout << i  << "    " << j << std::endl;
            while(threads >= MAX_THREADS){usleep(1000);} //Esperar a que baje el numero de hilos (LINUX)

            threads_Mutex.lock();
            threads++;  //Sumar un hilo
            futures.push_back(std::async(std::launch::async, iteracionHeatMap, std::ref(T), lambda, p, 0, std::ref(heatMap), i, j));
            threads_Mutex.unlock();
            //La funcion al terminar resta el hilo
            j++;
        }
        futures.clear();

        j = 0;
        for(double lambda = LAMBDA_MIN; lambda <= LAMBDA_MAX; lambda += (LAMBDA_MAX - LAMBDA_MIN)/(heatMap[i].size() - 1)){
            fileMk << p << "    " << lambda << "    " << heatMap[i][j] << "\n";
            j++;
        }
        fileMk.flush();
        i++;
    }

    fileMk.close();
    return 0;
}


typedef std::vector<double> Vector;
typedef std::vector<Vector> Matriz;

static Vector reduce(MobMatrix T, Matriz rhoR){
    Vector rho_zcta;
    rho_zcta.resize(T.N);
    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++)
            rho_zcta[i] += rhoR[i][j] * T.Mpesos[i][j];
        if(T.population[i] != 0)
            rho_zcta[i] /= T.population[i];
    }
    return rho_zcta;
}
static Vector merge(MobMatrix T, MergeMatrix Fraction, Vector rho_zcta){

    Vector rho_county, pop_county, area_county;
    rho_county.resize(Fraction.N2);
    pop_county.resize(Fraction.N2);
    area_county.resize(Fraction.N2);

    for(int i = 0; i < T.N; i++){
        for(int k = 0; k < Fraction.vecinos[i]; k++){
            rho_county[Fraction.Mvecinos[i][k]] += T.population[i] * Fraction.Mpesos[i][k] * rho_zcta[i];
            pop_county[Fraction.Mvecinos[i][k]] += T.population[i] * Fraction.Mpesos[i][k];
            area_county[Fraction.Mvecinos[i][k]] += T.area[i] * Fraction.Mpesos[i][k];
        }
    }
    for(int k = 0; k < Fraction.N2; k++){
        if(pop_county[k] != 0) {rho_county[k] /= pop_county[k];}
        std::cout << k << "\t" << area_county[k] << "\t" << pop_county[k] << std::endl;
    }
    return rho_county;
}

int threshold(){
    
    MobMatrix T1{"cities/ctynetwork.txt", "cities/Counties_areapop.txt"};
    
    std::ofstream output{"out/threshold_cty.txt"};
    for(double p = 1.0/128; p <= 1; p += 1.0/128){
        Mk_DistDiaNocheF markov0{0, p, T1};
        markov0.setLambda(1);
        LAMBDA0 = markov0.getLambda0();
        output << p << "\t" << 1.0/iteracion(T1, markov0).second << std::endl;
        }
    output.close();

    return 0;
}


int main(int argc, char* argv[]){
    
    Instrumentor::Get().BeginSession("Session Name");        // Begin session 
{
    InstrumentationTimer timer("Program");
    double p = 1;
    double lambda = 1.4; //lambda0 = 0.00414403 zcta, 0.00760323 subcty, 0.00716934 cty
    name = "10";
    if(argc == 2){name = argv[1]; p = std::stoi(name)/10.0;}

    {
    MobMatrix T1{"cities/zctanetwork.txt", "cities/zcta_areapop.txt"};
    Mk_DistDiaNocheF markov0{0, 0, T1};
    LAMBDA0 = markov0.getLambda0();
    }
    
    return heatMap();
    // MobMatrix T{"cities/zctanetwork.txt", "cities/zcta_areapop.txt"};
    // Mk_DistDiaNocheF markov{0, 0, T};
    // std::cout << markov.zD << std::endl;
    // return 0;

    //return heatMap();
    
//zcta
{
    MobMatrix T{"cities/zctanetwork.txt", "cities/zcta_areapop.txt"};
    MergeMatrix Fraction{"cities/Mergelabelzctacounty.txt"};
    
    auto rhoR_zcta = reduce(T, markovres(T, p, lambda));
    auto rhoR_county = merge(T, Fraction, rhoR_zcta);

    std::ofstream output("out/result_zcta.txt");
    for(int i = 0; i < Fraction.N2; i++){
        output << i << "\t" << rhoR_county[i] << std::endl;
    }
    output.close();
}
//subcty
{
    MobMatrix T{"cities/subctynetwork.txt", "cities/Subcounties_areapop.txt"};
    MergeMatrix Fraction{"cities/Mergelabelsubcountycounty.txt"};

    auto rhoR_zcta = reduce(T, markovres(T, p, lambda));
    auto rhoR_county = merge(T, Fraction, rhoR_zcta);
    
    std::ofstream output("out/result_subcounty.txt");
    for(int i = 0; i < Fraction.N2; i++){
        output << i << "\t" << rhoR_county[i] << std::endl;
    }
    output.close();
}
//cty
{
    MobMatrix T{"cities/ctynetwork.txt", "cities/Counties_areapop.txt"};
    auto rhoR_zcta = reduce(T, markovres(T, p, lambda));

    std::ofstream output("out/result_county.txt");
    for(int i = 0; i < T.N; i++){
        output << i << "\t" << rhoR_zcta[i] << std::endl;
    }
    output.close();
}

}
    Instrumentor::Get().EndSession();       // End Session
    return 0;
}
