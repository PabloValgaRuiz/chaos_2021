#include <precompiled.hpp>
#include "config.hpp"
#include "MobMatrix.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <math.h>
#include <future>
#include <unistd.h>
#include "mkmc.hpp"
#include "evoTemp.hpp"

static std::mutex s_Mutex;
static std::mutex threads_Mutex;
static int threads;
void iteracionHeatMap(const MobMatrix& T, double lambda, double p, int tiempo, std::vector<std::vector<double>>& heatMap, int i, int j);
int main(int argc, char* argv[]){
    
    return mkmc(argc, argv);

    Instrumentor::Get().BeginSession("Session Name");        // Begin session 
{
    InstrumentationTimer timer("Program");

    std::cout << "Input city: ";
    std::string name{};
    std::getline(std::cin, name);
    
    const double LAMBDA_MIN = 0;
    const double LAMBDA_MAX = 3;
    const double P_MIN = 0;
    const double P_MAX = 1;
    
    int tiempo = 1600;
    const int MAX_THREADS = 24;

    MobMatrix T{"cities3/" + name +"/mobnetwork.txt", "cities3/" + name +"/Poparea.txt"};
    std::cout << T.Pob << " cities3/" + name +"/mobnetwork.txt" <<std::endl;

    std::vector<std::vector<double>> heatMap;
    heatMap.resize(129); //Numero de p's //2^7 + 1
    for(auto& i : heatMap){i.resize(129);} //Numero de lambdas

    std::ofstream fileMk("out/" + name + "/heatmapDist.txt");
    std::vector<std::future<void>> futures;
    double i = 0;
    threads = 0;
    for(double p = P_MIN; p <= P_MAX; p += (P_MAX - P_MIN)/(heatMap.size() - 1)){
        InstrumentationTimer timer("p valor");
        double j = 0;
        
        for(double lambda = LAMBDA_MIN; lambda <= LAMBDA_MAX; lambda += (LAMBDA_MAX - LAMBDA_MIN)/(heatMap[i].size() - 1)){
            std::cout << i  << "    " << j << std::endl;
            while(threads >= MAX_THREADS){usleep(1000);} //Esperar a que baje el numero de hilos (LINUX)

            threads_Mutex.lock();
            threads++;  //Sumar un hilo
            futures.push_back(std::async(std::launch::async, iteracionHeatMap, std::ref(T), lambda, p, tiempo, std::ref(heatMap), i, j));
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
}
    Instrumentor::Get().EndSession();       // End Session
    return 0;
}

void iteracionHeatMap(const MobMatrix& T, double lambda, double p, int tiempo, std::vector<std::vector<double>>& heatMap, int i, int j){
{
    PROFILE_FUNCTION();
    Mk_DistDiaNocheF markov{0, p, T};
    markov.inicializar(0.0000001);
    markov.setLambda(lambda);

    for(size_t t = 0; t < tiempo; t++){
        markov.iteracion(T);
    }
    markov.contarInfectados(T);

    //Cada hilo cambia un par (i,j) distinto, no es necesario un mutex
    heatMap[i][j] = markov.totalRhoE/T.Pob;
}//Cerrar el scope antes de disminuir el valor de threads, para que se libere la memoria de las variables locales antes de crear el siguiente
    std::lock_guard<std::mutex> lock(threads_Mutex);
    threads--;
}