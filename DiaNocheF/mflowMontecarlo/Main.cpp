#include <precompiled.hpp>
#include "config.hpp"
#include "MobMatrix.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <future>

char eigensfile[] = "eigenvectores/0_07.txt";
std::string carpeta = "out/montecarlo0.txt";	//Output

struct ChosenLink{
	int Nlc;
	std::vector<int> linksChosen;
	std::vector<int> popLinkChosen;
	std::vector<int> j_esimo;
	int popTotChosen;
};

int find_maximum(std::vector<double> a, int n);

void chooseLinks(const MobMatrix& T, ChosenLink& C); //Wrapper: pasa los argumentos agrupados en el struct ChosenLink
void chooseLinks(const MobMatrix& T, int NlinksChosen, std::vector<int>& linksChosen, std::vector<int>& popLinkChosen, std::vector<int>& j_esimo, int* popTotChosen);

int contar_Infectados_random (const MobMatrix& T, const MC_DistDiaNocheF& MC, size_t MUESTRA);

int contar_Infectados_chosen (const MobMatrix& T, const MC_DistDiaNocheF& MC, const ChosenLink& C, size_t MUESTRA);
int contar_Infectados_chosen (const MobMatrix& T, const MC_DistDiaNocheF& MC, int NlinksChosen, const std::vector<int>& linksChosen,
	const std::vector<int>& popLinkChosen, const std::vector<int>& j_esimo, int popTotChosen, size_t MUESTRA);

static double p = 0.7;				//Displacement probability
static int nThreads = 24; //About 8
static int nPasos = 100 ; //1000
static double beta = 0.4;

static std::mutex miMutex; //Llave (key) para bloquear la modificacion del heatmap

void iteraciones(const MobMatrix& T, const std::vector<ChosenLink>& vectorChosenLinks, std::vector<std::vector<double>>& heatMap){

	MC_DistDiaNocheF montecarlo{0, p, T};
	montecarlo.setLambda(beta);
	montecarlo.inicializar(0.0001);
	int PobInf;

	for (int t = 0; t < nPasos; t++){
		montecarlo.iteracion(T);
		std::cout << "Iteración: " << t << std::endl;
	}

	for(int i = 0; i < heatMap.size(); i++){//iteracion sobre links
		for(int j = 0; j < heatMap[i].size(); j++){//iteracion sobre tests

			int MUESTRA = ((j+1) * 10000) / heatMap[i].size();
			PobInf = contar_Infectados_chosen(T, montecarlo, vectorChosenLinks[i], MUESTRA);
			std::lock_guard<std::mutex> lock(miMutex);
			std::cout << vectorChosenLinks[i].Nlc << "\t" << ((j+1) * 10000) / heatMap[i].size() << "\t" << PobInf << std::endl;
			heatMap[i][j] += static_cast<double>(PobInf) / nThreads;
		}
	}
}

int main(int argc, char* argv[]){

	MobMatrix T{"cities2/oldnetwork0.txt", "areas/area0.txt"};
    std::cout << T.Pob <<std::endl;

#if 1
	//Heatmap dimensions
	size_t sizeLinks = 33;
	size_t sizeTests = 5;

	//Heatmap matrix

	std::vector<std::vector<double>> heatMap;
	heatMap.resize(sizeLinks);
	for(auto& v : heatMap)
		v.resize(sizeTests);

	//Different link collections to use
	std::vector<ChosenLink> vectorChosenLinks;
	vectorChosenLinks.resize(sizeLinks);

	// auto& C = vectorChosenLinks[0];
	// size_t NlcTemp = T.Links / (20 * sizeLinks + 1); //Care to not choose 0 or all the links (ex: if 1000 links, take from 166 to 866)
	// C.Nlc = NlcTemp;
	// C.linksChosen.resize(NlcTemp);
	// C.popLinkChosen.resize(NlcTemp);
	// C.j_esimo.resize(NlcTemp);
	// chooseLinks(T, C);
	// std::cout << 0 << "\t" << NlcTemp << std::endl;

	std::vector<std::future<void>> futures;
	for(size_t i = 0; i < sizeLinks; ++i){
		//futures.push_back(std::async(std::launch::async, [](const MobMatrix& T, ChosenLink& C, size_t i, size_t sizeLinks){
			auto& C = vectorChosenLinks[i];

			size_t NlcTemp = T.Links * (i+1) / (20 * sizeLinks + 1); //Care to not choose 0 or all the links (ex: if 1000 links, take from 166 to 866)
			C.Nlc = NlcTemp;
			C.linksChosen.resize(NlcTemp);
			C.popLinkChosen.resize(NlcTemp);
			C.j_esimo.resize(NlcTemp);

			//Choose the Nlc highest component links in the eigenvector
			chooseLinks(T, C);
			//This is the last time vectorChosenLinks will be modified

			std::cout << i << "\t" << NlcTemp << std::endl;

			//}, T, std::ref(vectorChosenLinks[i]), i, sizeLinks));
	}
	futures.clear();

	std::cout << "Links chosen" << std::endl;

	for(int l = 0; l < nThreads; l++){
		futures.push_back(std::async(std::launch::async, iteraciones, T, vectorChosenLinks, std::ref(heatMap)));
	}
	futures.clear();

	std::ofstream f(carpeta);
	for(int i = 0; i < heatMap.size(); i++)//iteracion sobre links
		for(int j = 0; j < heatMap[i].size(); j++)//iteracion sobre tests
			f << vectorChosenLinks[i].Nlc << "\t" << ((j+1) * 10000) / heatMap[i].size() << "\t" << heatMap[i][j] << "\n";
	f.close();
#endif
	return 0;
}
int contar_Infectados_random (const MobMatrix& T, const MC_DistDiaNocheF& MC, size_t MUESTRA){
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

	int PobInf = 0;
	int i,j,k,l;
	int persona, persona2;

	for(l = 0; l < MUESTRA; l++){
		persona = dist(mt) * T.Pob;
		if(MC.getEst()[persona] == E) PobInf++;
	}

	return PobInf;
}

int contar_Infectados_chosen (const MobMatrix& T, const MC_DistDiaNocheF& MC, const ChosenLink& C, size_t MUESTRA){
	return contar_Infectados_chosen(T, MC, C.Nlc, C.linksChosen, C.popLinkChosen, C.j_esimo, C.popTotChosen, MUESTRA);
}
int contar_Infectados_chosen(const MobMatrix& T, const MC_DistDiaNocheF& MC, int NlinksChosen, const std::vector<int>& linksChosen, 
	const std::vector<int>& popLinkChosen, const std::vector<int>& j_esimo, int popTotChosen, size_t MUESTRA)
{
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    int i,j,k;
    k=0;
    int PobInfChosen=0;
    //FILE *f;
    //f = fopen("Graficas/Infectados.txt","a");
	int persona, persona2, persona3;
	int l;
	int prov, prov2;
	for(l = 0; l < MUESTRA; l++){
		persona = dist(mt) * popTotChosen; //enésima persona del grupo de personas de los links escogidos
		prov = 0;
		for(k = 0; k < NlinksChosen; k++){
			//std::cout << "k = " << k << std::endl;
			//std::cout << linksChosen[k]/T.N << "\t" << j_esimo[k] << std::endl;
			prov += T.Mpesos[linksChosen[k]/T.N][j_esimo[k]];
			//std::cout << linksChosen[k]/T.N << "\t" << j_esimo[k] << std::endl;
			if(persona < prov){
				persona2 = persona - prov + T.Mpesos[linksChosen[k]/T.N][j_esimo[k]]; //enésima persona del link [linksChosen[k]/N][j_esimo[k]]
				prov2 = 0;
				for(i = 0; i < T.N; i++){
					//std::cout << "i = " << i << std::endl;
					for(j = 0; j < T.vecinos[i] && i < T.N; j++){
						//std::cout << "j = " << j << std::endl;
						if(i == linksChosen[k]/T.N && j == j_esimo[k]){
							persona3 = prov2; //enésima persona del total de la poblacion
							persona3 += persona2;
							
							if(MC.getEst()[persona3] == R) PobInfChosen++; //CONTAR A ESA PERSONA COMO INFECTADA

							i = j = T.N; //Salir de los loops
						}
						else prov2 += T.Mpesos[i][j];
					}
				}
				k = NlinksChosen; //Salir del loop
			}
		}

	}
	return PobInfChosen;
}

void chooseLinks(const MobMatrix& T, ChosenLink& C){
	chooseLinks(T, C.Nlc, C.linksChosen, C.popLinkChosen, C.j_esimo, &C.popTotChosen);
}

void chooseLinks(const MobMatrix& T, int NlinksChosen, std::vector<int>& linksChosen, std::vector<int>& popLinkChosen, std::vector<int>& j_esimo, int* popTotChosen){

	
	int i = 0, j = 0;
	int k;
	int a,b;
	double c;
	
	static bool control = 1; 
	static std::vector<double> eigensOriginal;
	if(control){ //Solo hacer esto la primera vez que se llama a la funcion
		eigensOriginal.resize(T.N * T.N);
		std::ifstream input(eigensfile);
		while(input >> a >> b >> eigensOriginal[i]){
			i++;
			if(i >= T.N*T.N) break;
		}
		input.close();
		control = 0;
		std::cout << "Done" << std::endl;
	}
	
	

	std::vector<double> eigens = eigensOriginal;	
	int check, prov;
	for(k = 0; k < NlinksChosen; k++){
		check = 1;
		do{
			prov = find_maximum(eigens, T.N*T.N);
			eigens[prov] = 0;
			if(T.vecinos[prov/T.N] == 0) continue;
			for(int j = 0; j < T.vecinos[prov/T.N]; j++){
				if(prov%T.N == T.Mvecinos[prov/T.N][j])
					check = 0;
			}
		}while(check);
		linksChosen[k] = prov;
		
		//std::cout << linksChosen[k] << std::endl;
		//printf("%d\t%d\n", linksChosen[k]/N, linksChosen[k]%N );
	}
	*popTotChosen = 0;

	for(k = 0; k < NlinksChosen; k++){
		int prov = 0;
		for(i = 0; i < T.N; i++){
			for(j = 0; j < T.vecinos[i]; j++){
				if(i == linksChosen[k]/T.N && T.Mvecinos[i][j] == linksChosen[k]%T.N){
					popLinkChosen[k] = prov;
					j_esimo[k] = j;
					*popTotChosen += T.Mpesos[i][j];
				}
				prov += T.Mpesos[i][j];
			}
		}
	}
}

int find_maximum(std::vector<double> a, int n) {
  int c, index = 0;
	for (c = 1; c < n; c++){
		if (a[c] > a[index]){
			index = c;
		}
	}
  return index;
}