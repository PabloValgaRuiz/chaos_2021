#include "config.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <iostream>
#include <tuple>
#include "MobMatrix.hpp"
#include "Mk_IndistDiaNocheF.hpp"
#include "benchmark.hpp"


std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const Mk_IndistDiaNocheF& markov)
{

    //==============================MATRIX CREATION====================================================//

InstrumentationTimer* timer = new InstrumentationTimer("Matrix construction");

	std::vector<Eigen::Triplet<double> > tripletList;
	
	tripletList.reserve(static_cast<int>(4 * T.N * T.N));
	Eigen::SparseMatrix<double, Eigen::ColMajor> MATRIX(T.N, T.N);	//SI SE UTILIZA ROWMAJOR CRASHEA EL GENEIGSSOLVER

	//TERMINOS CONSTANTES
	std::vector<double> Elem_dia; Elem_dia.resize(T.N);
	std::vector<double> Elem_noche; Elem_noche.resize(T.N);

	for(int i = 0; i < T.N; i++){
		if(markov.n_eff[i] == 0)
			Elem_dia[i] = 0;
		else
			Elem_dia[i] = markov.zD * markov.fvector[i] / markov.n_eff[i];

		if(T.population[i] == 0)
			Elem_noche[i] = 0;
		else
			Elem_noche[i] = markov.zN * markov.sigma / T.population[i];
	}

	int N = T.N;
	for(int i = 0; i < T.N; i++){

		tripletList.emplace_back(i, i, ((1-markov.p)*(1 - markov.p) * Elem_dia[i] + Elem_noche[i]) * T.population[i]);

		for(int j = 0; j < T.vecinos[i]; j++)
			if(T.population[i] != 0)
				tripletList.emplace_back(i, T.Mvecinos[i][j], markov.p*(1-markov.p) * Elem_dia[T.Mvecinos[i][j]] * T.population[T.Mvecinos[i][j]] * (T.Mpesos[i][j]/T.population[i]));

		for(int j = 0; j < T.vecinosT[i]; j++)
			tripletList.emplace_back(i, T.MvecinosT[i][j], markov.p*(1-markov.p) * Elem_dia[i] * T.MpesosT[i][j]);

		std::vector<double> temp; temp.resize(T.N);
		for(int l = 0; l < T.vecinos[i]; l++){
			for(int j = 0; j < T.vecinosT[T.Mvecinos[i][l]]; j++)	//Literalmente un infierno no se ni que he hecho
				if(T.population[i] != 0)
					temp[T.MvecinosT[T.Mvecinos[i][l]][j]] += (T.Mpesos[i][l]/T.population[i]) * T.MpesosT[T.Mvecinos[i][l]][j] * Elem_dia[T.Mvecinos[i][l]];
		}		
		for(int j = 0; j < N; j++){
			tripletList.emplace_back(i, j, temp[j] * markov.p * markov.p);
		}

		// for(int j = 0; j < T.N; j++){
		// 	double temp = 0;
		// 	for(int l = 0; l < T.N; l++)
		// 		temp += T.MpesosDense[i][l] * T.MpesosDense[j][l] * Elem_dia[l];
		// 	tripletList.emplace_back(i, j, temp * markov.p * markov.p * T.population[j]);
		// }
	}
	MATRIX.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();
delete timer;
	//==============================EIGEN CALCULATION===================================================//
InstrumentationTimer* timer2 = new InstrumentationTimer("Eigenvalues");
	// Construct matrix operation object using the wrapper class SparseGenMatProd
	Spectra::SparseGenMatProd<double> op(MATRIX);
	// Construct eigen solver object, requesting the largest eigenvalues
	Spectra::GenEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigs(&op, 1, 6);
	// Initialize and compute
	eigs.init();
	//std::cout << "Computing..." << std::endl;
	int nconv = eigs.compute();											//SI SE UTILIZA ROWMAJOR CRASHEA EL GENEIGSSOLVER

#ifdef SIS
	double R0 = (markov.lambda/markov.mu) * eigs.eigenvalues().real().coeffRef(0);
#endif

#ifdef SEAPIDR
	double R0 = (markov.lambdaA * (1 - markov.x) / markov.muA + markov.lambdaP * markov.x / markov.alpha + markov.lambdaI * markov.x / (markov.delta + markov.muI)) * 
		eigs.eigenvalues().real().coeffRef(0);
#endif

delete timer2;
	return {eigs.eigenvectors().real().array().abs(), R0};

}
std::pair<Eigen::ArrayXd, double> eigen_p_1(int N, Eigen::MatrixXi mPesosDense)
{
	int i = 0;
	double eigenvalue = static_cast<double>(mPesosDense.colwise().sum().maxCoeff(&i));
	mPesosDense.colwise().sum().maxCoeff(&i);
	Eigen::ArrayXd eigenvector(N*N); eigenvector.setZero();

	double valor = 1/pow(N, 0.5);
	for (int k = i; k < N*N; k += N)
		eigenvector(k) = valor;

	for(int i = 0; i < N; i++){	//Eliminar las componentes del autovector que no tengan poblacion
		for(int j = 0; j < N; j++){
			if(mPesosDense(i,j) == 0) eigenvector(i*N + j) = 0;
		}
	}

	return {eigenvector, eigenvalue};

}