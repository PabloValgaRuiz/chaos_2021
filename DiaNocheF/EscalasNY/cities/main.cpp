#include <iostream>
#include <fstream>

int main(){

    std::ifstream input1{"Subcounties_areapop.txt"};
    std::string trash; input1 >> trash >> trash >> trash;
    int County,Pop;
    double Area;
    double temp = 0;
    int cty;
    while(input1 >> County >> Area >> Pop){
        if(Pop/Area > temp) {temp = Pop/Area; cty = County;}
    }
    std::cout << cty << "\t" << temp << std::endl;
    input1.close();
    return 0;





    // std::ifstream input("Subcounties_areapop.txt");

    // std::ofstream output("Subcounties_areapop(bis).txt");

    // std::string trash; input >> trash >> trash >> trash;
    // int County;
    // double Pop, Area;
    // output.precision(12);
    // output << "Subcounty\tArea\tPop" << std::endl; 
    // while(input >> County >> Pop >> Area){
    //     output << County << "\t" << Area << "\t" << (int)Pop << std::endl;
    // }
    // output.close();
}