#pragma once
#include "config.hpp"
#include <string>
#include <vector>


class MobMatrix{//MODIFICADA PARA ESCALAS NUEVA YORK
protected:

    std::string mobility_network;
    std::string pop_area;

    void readN();
    void readCityPatch();
    void readPopArea();
    void leer_vecinos(); //LOS FICHEROS NO PUEDEN ACABAR CON UNA LINEA EN BLANCO: DARÁ UN VECINO DE MÁS Y CRASHEARÁ
    void leer_vecinosT();
    void leer_matrices();
    void calculaTraspuesta(const std::vector<std::vector<int> >& Matrix, std::vector<std::vector<int> >& MatrixT);


public:

    int N = 0, Links = 0, Pob = 0;
    double Ratio = 1;
    std::vector<int> population;
    std::vector<double> area;
    std::vector<int> vecinos;
    std::vector<int> vecinosT;
    std::vector<std::vector<int>> Mvecinos;
    std::vector<std::vector<int> > MvecinosT;
    std::vector<std::vector<int> > Mpesos;
    std::vector<std::vector<int> > MpesosT;

    MobMatrix(const std::string& _mobility_network, const std::string& _pop_area);

};