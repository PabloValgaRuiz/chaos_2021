#include "config.hpp"
#include "MergeMatrix.hpp"
#include "benchmark.hpp"
#include <fstream>
#include <iostream>


MergeMatrix::MergeMatrix(const std::string& _mobility_network){
    PROFILE_FUNCTION();


    this->mobility_network = _mobility_network;
    readN();
    leer_vecinos();
    leer_matrices();
}

void MergeMatrix::readN(){

    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file" << std::endl;
        std::exit(1); // terminate with error
    }

    int a,b; double c;
    std::string trash;
    N1 = N2 = Links = 0;
    inFile >> trash >> trash >> trash;
    while(inFile >> a >> b >> c){
        Links++;
        if(a >= N1) N1 = a + 1;
        if(b >= N2) N2 = b + 1;
    }

    std::cout << N1 << "    " << N2 << " " << Links << std::endl;
    inFile.close();
}

void MergeMatrix::leer_vecinos()
{
    int i, j1, j2, trash1; double trash2;
    vecinos.resize(N1);
    for ( i = 0 ; i < N1 ; i++)
        vecinos[i] = 0;
    
    std::ifstream inFile;
    inFile.open(this->mobility_network);
    if (!inFile) {
        std::cout << "Unable to open file";
        std::exit(1); // terminate with error
    }

    std::string trash;
    inFile >> trash >> trash >> trash;
    while(inFile >> j1 >> trash1 >> trash2){
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

void MergeMatrix::leer_matrices()
{
    Mvecinos.resize(N1);
    Mpesos.resize(N1);
    for(int i = 0; i < N1; i++){
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
    for(i = 0 ; i < N1 ; i++) {
        for ( j = 0 ; j < vecinos[i] ; j++) {
            inFile >> trash1 >> Mvecinos[i][j] >> Mpesos[i][j];
            Mpesos[i][j] = Ratio * Mpesos[i][j];
            //if(Mpesos[i][j] == 0) Mpesos[i][j] = 10;
            if(trash1 != i) { std::cout << "Error en la lectura de matrices"; exit(EXIT_FAILURE);}
        }
    }

    inFile.close();
}
