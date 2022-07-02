#include "mkmc.hpp"

static std::mutex s_Mutex;
static std::string name;
static double montecarlo0, markov0;
static int termalizacion = 20; //150

static void iteracionesMk(const MobMatrix& T, double lambda, double p, int tiempo){
    //std::ofstream fileMk("out/Mk_0_p" + name + ".txt");
    std::ofstream fileMk("out/benchmarking" + name + ".txt");
    for(double lambda = 0; lambda < 10.00000000001; lambda += 0.04){    // 0.02
        Mk_DistDiaNocheF markov{0, p, T};
        markov.inicializar(0.008);
        markov.setLambda(lambda);
        for(int i = 0; i < 20; i++){ //500 EN VEZ DE TERMALIZACION
            markov.iteracion(T);
            markov.contarInfectados(T);
        }
        markov.contarInfectados(T);
        fileMk << lambda << " " << markov.totalRhoE << std::endl;
    }
    fileMk.close();
}

static void iteracionesMC(const MobMatrix& T, double &infectados, double lambda, double p, int tiempo, int i){
    MC_DistDiaNocheF montecarlo{0, p, T};
    montecarlo.setLambda(lambda);
    montecarlo.inicializar(0.008);

    for(int k = 0; k < termalizacion; k++){ //Termalizacion
        montecarlo.iteracion(T);
    }
    for(int k = 0; k < tiempo; k++){
        montecarlo.iteracion(T);
        montecarlo.calcInfectados();
        std::lock_guard<std::mutex> lock(s_Mutex);
        infectados += montecarlo.totalE;

    }
}
int mkmc(int argc, char* argv[]){
    Instrumentor::Get().BeginSession("Session Name");        // Begin session 
{
    InstrumentationTimer timer("Program");

    double lambda = 2;
    double p = 0.2;
    name = "02";
    int tiempo = 2;
    int hilos = 2;

    if(argc == 2){name = argv[1]; p = std::stoi(name)/10.0;}

    std::cout << p << " " << name << std::endl;
    
    MobMatrix T{"cities2/oldnetwork14.txt", "areas/area14.txt"};
    std::cout << T.Pob <<std::endl;

    //std::ofstream fileMC("out/MC_00_p" + name + ".txt");
    std::ofstream fileMC("out/benchmarking2" + name + ".txt");

    auto futureMk = std::async(std::launch::async, iteracionesMk, std::ref(T), lambda, p, tiempo);
    for(lambda = 0; lambda < 10.000000000001; lambda += 2.5){
        InstrumentationTimer timer2{"Loop for lambda"};
        
        double infectados = 0;
        std::vector<std::future<void>> futures;
        for(int i = 0; i < hilos; i++)
            futures.push_back(std::async(std::launch::async, iteracionesMC, std::ref(T), std::ref(infectados), lambda, p, tiempo, i));
        futures.clear();

        fileMC << lambda << " " << infectados/(tiempo * hilos) << std::endl;
    }
    futureMk.wait();
    fileMC.close();    
    

    

}
    Instrumentor::Get().EndSession();       // End Session
    return 0;
}