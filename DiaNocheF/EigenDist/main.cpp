#include <precompiled.hpp>
#include "config.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <unistd.h>
#include <future>
#include "ThreadPool.hpp"

std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const Mk_DistDiaNocheF& markov);
std::pair<Eigen::ArrayXd, double> eigen_p_1(int N, Eigen::MatrixXi mPesosDense);

static std::mutex threads_Mutex;
static int threads;



void iteracionWrapper(const MobMatrix& T, const Mk_DistDiaNocheF markov, double& result){
    try{
    result = iteracion(T, markov).second;
    }catch(...){
        std::cout << "Excepcion" << std::endl;
    }
    std::lock_guard<std::mutex> lock(threads_Mutex);
    threads--;
}

int main(){
    Instrumentor::Get().BeginSession("Mi sesion");
{
PROFILE_SCOPE("Program");
    // char city_char[100];
    // std::cout << "Input city: ";
    // std::cin.getline(city_char, 100);
    // std::string city{city_char};

    // MobMatrix T{"cities3/" + city + "/mobnetwork.txt", "cities3/" + city + "/Poparea.txt"};
    // std::cout << T.Pob <<std::endl;
    
    // const double P_MIN = 0;
    // const double P_MAX = 1;

    // const int MAX_THREADS = 6;
    
    // std::vector<double> vectorP;
    // vectorP.resize(129); //Numero de p's //2^7 + 1

    // std::vector<std::future<void>> futures;
    // int i = 0;
    // for(double p = P_MIN; p <= P_MAX; p += (P_MAX - P_MIN)/(vectorP.size() - 1)){

    //     while(threads >= MAX_THREADS){usleep(1000);} //Esperar a que baje el numero de hilos (LINUX)

    //     Mk_DistDiaNocheF markov{0, p, T};
    //     markov.setLambda(1);
    //     threads_Mutex.lock();
    //     threads++;  //Sumar un hilo
    //     threads_Mutex.unlock();
    //     futures.push_back(std::async(std::launch::async, iteracionWrapper, std::ref(T), markov, std::ref(vectorP[i])));
    //     i++;
    // }
    // futures.clear();

    // std::ofstream file("out/" + city + "/thresholdDist.txt");
    // i = 0;
    // for(double p = P_MIN; p <= P_MAX; p += (P_MAX - P_MIN)/(vectorP.size() - 1)){
    //     file << p << "  " << 1.0/vectorP[i] << "\n";
    //     // Mk_DistDiaNocheF markov{0, p, T};
    //     // markov.setLambda(1);
    //     // file << p << "  " << 1.0/iteracion(T, markov).second << "\n";
    //     i++;
    // }
    // file.flush();
    // file.close();
    
    
   //OBTENER EL EIGENVECTOR DE TODAS LAS CIUDADES

    ThreadPool pool{6};

    std::array<std::string, 1> names = {"bogota"};//{"austin", "baltimore", "dc", "detroit", "los angeles", "ma", "miami", "ny", "seattle"};

    for(auto name : names){
        pool.enqueue([=]{

            MobMatrix T{"cities3/"+name+"/mobnetwork.txt", "cities3/"+name+"/Poparea.txt"};
            std::ofstream output("out/eigenvectors3/" + name + "_10.txt");
            Mk_DistDiaNocheF markov(0, 1, T);
            markov.setLambda(1);

            auto[eigenvector, eigenvalue] = iteracion(T, markov);
            for(int i = 0; i < T.N; i++){
                for(int j = 0; j < T.N; j++){
                    output << i << "\t" << j << "\t" << eigenvector[i*T.N+j] << std::endl;
                }
            }
            std::cout << 1 << "\t" << eigenvalue << std::endl;
            output.close();

        });
    }

}
    Instrumentor::Get().EndSession();
}