#pragma once
#include "config.hpp"
#include <string>
#include <vector>


class MergeMatrix{//MODIFICADA PARA ESCALAS NUEVA YORK
protected:

    std::string mobility_network;

    void readN();
    void readCityPatch();
    void leer_vecinos(); //LOS FICHEROS NO PUEDEN ACABAR CON UNA LINEA EN BLANCO: DARÁ UN VECINO DE MÁS Y CRASHEARÁ
    void leer_matrices();


public:

    int N1 = 0, N2 = 0, Links = 0;
    double Ratio = 1;
    std::vector<int> population;
    std::vector<double> area;
    std::vector<int> vecinos;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<double>> Mpesos;

    MergeMatrix(const std::string& _mobility_network);

};