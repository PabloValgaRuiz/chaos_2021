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
    readPopArea();
    leer_matrices();
    calculaTraspuesta(Mpesos, MpesosT);
}

void MobMatrix::readN(){

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    int a,b,c;
    std::string trash;
    N = Links = 0;
    inFile >> trash >> trash >> trash;
    while(inFile >> a >> b >> c){
        Links++;
        if(a >= N) N = a + 1;
    }

    std::cout << N << " " << Links << std::endl;
    inFile.close();
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
    
    std::string trash;
    inFile >> trash >> trash >> trash;
    double pop;
    while(inFile >> i >> area[i] >> pop){
        if(i != k++) {std::cerr << "Fichero de áreas incompleto" << std::endl;}

        population[i] = Ratio * pop;
        Pob += Ratio * pop;
    }
    inFile.close();
    
}

void MobMatrix::leer_vecinos()
{
    int i, j1, j2, trash1;
    vecinos.resize(N);
    for ( i = 0 ; i < N ; i++)
        vecinos[i] = 0;
    
    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    std::string trash;
    inFile >> trash >> trash >> trash;
    while(inFile >> j1 >> trash1 >> trash1){
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
    int i, j1, j2, trash1;
    for ( i = 0 ; i < N ; i++)
        vecinosT[i] = 0;

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }
    std::string trash;
    inFile >> trash >> trash >> trash;
    while(inFile >> trash1 >> j1 >> trash1){
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

    int I,i,j;

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }
    int trash1;
    std::string trash;
    inFile >> trash >> trash >> trash;

    std::vector<int> pesos; pesos.resize(N);
    for(i = 0 ; i < N ; i++) {
        for ( j = 0 ; j < vecinos[i] ; j++) {
            inFile >> trash1 >> Mvecinos[i][j] >> Mpesos[i][j];
            pesos[i] += Mpesos[i][j];

            //if(Mpesos[i][j] == 0) Mpesos[i][j] = 10;
            if(trash1 != i) { std::cout << "Error en la lectura de matrices"; exit(EXIT_FAILURE);}
        }
    }
    for(i = 0 ; i < N ; i++) {
        for ( j = 0 ; j < vecinos[i] ; j++) {
            Mpesos[i][j] = ((double)Mpesos[i][j] * (double)population[i])/pesos[i];
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
