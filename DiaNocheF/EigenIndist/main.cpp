#include "config.hpp"
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <iostream>
#include <tuple>
#include "Mk_IndistDiaNocheF.hpp"
#include "benchmark.hpp"

std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const Mk_IndistDiaNocheF& markov);
std::pair<Eigen::ArrayXd, double> eigen_p_1(int N, Eigen::MatrixXi mPesosDense);

int main(){
    Instrumentor::Get().BeginSession("Mi sesion");
{
PROFILE_SCOPE("Program");

    std::string city;
    std::cout << "Insert city: ";
    std::cin >> city;
    MobMatrix T{"cities3/" + city + "/mobnetwork.txt", "cities3/" + city + "/Poparea.txt"};
    std::cout << T.Pob <<std::endl;

    std::ofstream file{"out/" + city + "/thresholdIndist.txt"};
    for(double p = 0; p <= 1; p += 1.0/128){
        Mk_IndistDiaNocheF markov{0, p, T};
        markov.setLambda(1);
        auto pair = iteracion(T, markov);
        file << p << "  " << 1.0/pair.second << std::endl;
    }

    file.close();


}
    Instrumentor::Get().EndSession();
    return 0;
}