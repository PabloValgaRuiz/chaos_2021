#include "config.hpp"
#include "MobMatrix.hpp"
#include "benchmark.hpp"
#include <fstream>
#include <iostream>


MobMatrix::MobMatrix(const std::string& _mobility_network, const std::string& _pop_area){
    PROFILE_FUNCTION();


    this->mobility_network = _mobility_network;
    this->pop_area = _pop_area;
    readN();
    leer_vecinos();
    leer_vecinosT();
    leer_matrices();
    calculaTraspuesta(Mpesos, MpesosT);
    readPopArea();
}

void MobMatrix::readN(){

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    int a,b,c;
    int NN = 0,LL = 0;
    inFile >> N >> Links;
    while(inFile >> a >> b >> c){
        LL++;
        if(a >= NN) NN = a + 1;
    }

    std::cout << N << " " << Links << std::endl;
    inFile.close();

    if (NN != N || LL != Links) {
        std::cout << NN << "    " << N << "    " << LL << "    " << Links << std::endl;
        std::cout << "Discrepancia de N o Links";
        std::exit(1); // terminate with error
    }

}

void MobMatrix::readPopArea(){
    Pob = 0;
    this->population.resize(N);
    this->area.resize(N);

    std::ifstream inFile;
    inFile.open(this->pop_area);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }
    int i, k = 0;
    while(inFile >> i >> area[i]){
        if(i != k++) {std::cerr << "Fichero de áreas incompleto" << std::endl;}
    }
    inFile.close();
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < vecinos[i]; j++){
            population[i] += Mpesos[i][j];
        }
        //area[i] = population[i]/100.0;
        Pob += population[i];
    }
}

void MobMatrix::leer_vecinos()
{
    int i, j1, j2, trash1, trash2, trash;
    vecinos.resize(N);
    for ( i = 0 ; i < N ; i++)
        vecinos[i] = 0;
    
    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    inFile >> trash >> trash;
    while(inFile >> j1 >> trash1 >> trash){
        vecinos[j1]++;
    }
    inFile.close();

    // FILE *f;
    //f = fopen(this->mobility_network.c_str(),"rt");
    // fscanf(f,"%d %d", &trash, &trash);
    // while(feof(f) == 0)
    // {
    //     fscanf(f,"%d %d %d", &j1, &trash1, &trash);
    //     if(trash1 != trash2 || j1 != j2)
    //         vecinos[j1]++;
    //     j2 = j1;
    //     trash2 = trash1;
    // }
    // fclose(f);

    
}
void MobMatrix::leer_vecinosT()
{
    vecinosT.resize(N);
    int i, j1, j2, trash1, trash2, trash;
    for ( i = 0 ; i < N ; i++)
        vecinosT[i] = 0;

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    inFile >> trash >> trash;
    while(inFile >> trash1 >> j1 >> trash){
        vecinosT[j1]++;
    }
    inFile.close();

    // FILE *f;
    // f = fopen(this->mobility_network.c_str(),"rt");
    // fscanf(f,"%d %d", &trash, &trash);
    // while(feof(f) == 0)
    // {
    //     fscanf(f,"%d %d %d", &trash1, &j1, &trash);
    //     if(trash1 != trash2 || j1 != j2)
    //         vecinosT[j1]++;
    //     j2 = j1;
    //     trash2 = trash1;
    // }
    // fclose(f);

    
}

void MobMatrix::leer_matrices()
{
    Mvecinos.resize(N);
    Mpesos.resize(N);
    for(int i = 0; i < N; i++){
        Mvecinos[i].resize(vecinos[i]);
        Mpesos[i].resize(vecinos[i]);
    }

    int I,i,j,trash;

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    inFile >> trash >> trash;
    for(i = 0 ; i < N ; i++) {
        for ( j = 0 ; j < vecinos[i] ; j++) {
            inFile >> trash >> Mvecinos[i][j] >> Mpesos[i][j];
            Mpesos[i][j] = Ratio * Mpesos[i][j];
            //if(Mpesos[i][j] == 0) Mpesos[i][j] = 10;
            if(trash != i) { std::cout << "Error en la lectura de matrices"; exit(EXIT_FAILURE);}
        }
    }

    inFile.close();
}

void MobMatrix::calculaTraspuesta(const std::vector<std::vector<int>>& Matrix, std::vector<std::vector<int>>& MatrixT)
{  
    MvecinosT.resize(N);
    MatrixT.resize(N);
    for(int i = 0; i < N; i++){
        MvecinosT[i].resize(vecinosT[i]);
        MatrixT[i].resize(vecinosT[i]);
    }
	int k = 0;
    int i = 0, j = 0;
    int B;
	int* AUX = (int*)malloc(N * sizeof(int));
	for (i = 0; i < N; i++)
		AUX[i] = 0;

	for (i = 0; i < N; i++){
		for (j = 0; j < vecinos[i]; j++){
            B = Mvecinos[i][j];
            MvecinosT[B][AUX[B]] = i;
		    MatrixT[B][AUX[B]] = Matrix[i][j];
		    AUX[B] = AUX[B] + 1;
		}
	}

	free(AUX);
}
