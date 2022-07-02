#include <precompiled.hpp>
#include "config.hpp"
#include "MobMatrix.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"


std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const Mk_DistDiaNocheF& markov)
{

    //==============================MATRIX CREATION====================================================//

InstrumentationTimer* timer = new InstrumentationTimer("Matrix construction");

	std::vector<Eigen::Triplet<double> > tripletList;
	
	std::cout << "MegaBytes para almacenar:\t" << 4 * T.N * T.Links * sizeof(Eigen::Triplet<double>)/1000000 << "\n";
	tripletList.reserve(static_cast<int>(4 * T.N * T.Links));
	Eigen::SparseMatrix<double, Eigen::ColMajor> MATRIX(T.N * T.N, T.N * T.N);	//SI SE UTILIZA ROWMAJOR CRASHEA EL GENEIGSSOLVER

	//TERMINOS CONSTANTES
	std::vector<double> Elem_ik; Elem_ik.resize(T.N);
	std::vector<double> Elem_jk; Elem_jk.resize(T.N);
	std::vector<double> Elem_ki; Elem_ki.resize(T.N);
	std::vector<double> Elem_kj; Elem_kj.resize(T.N);
	for(int i = 0; i < T.N; i++){
		if(markov.n_eff[i] == 0){Elem_ik[i] = Elem_jk[i] = Elem_ki[i] = Elem_kj[i] = 0;
			if(T.population[i] != 0)
				Elem_ik[i] = markov.zN * markov.sigma / T.population[i];
		}
		else{
			if(T.population[i] == 0)
				Elem_ik[i] = (1 - markov.p)*(1 - markov.p) * markov.zD * markov.fvector[i] / markov.n_eff[i];
			else
				Elem_ik[i] = markov.zN * markov.sigma / T.population[i] + (1 - markov.p)*(1 - markov.p) * markov.zD * markov.fvector[i] / markov.n_eff[i];
			Elem_jk[i] = markov.p * (1 - markov.p) * markov.zD * markov.fvector[i] / markov.n_eff[i];
			Elem_ki[i] = markov.p * (1 - markov.p) * markov.zD * markov.fvector[i] / markov.n_eff[i];
			Elem_kj[i] = markov.p * markov.p * markov.zD * markov.fvector[i] / markov.n_eff[i];
		}
	}


	int N = T.N;
	//Version compacta: eliminando las columnas cuyas filas asociadas son cero
	for(int i = 0; i < T.N; i++){
		for(int j = 0; j < T.vecinos[i]; j++){
			for(int k = 0; k < T.vecinos[i]; k++)
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], i * N + T.Mvecinos[i][k], Elem_ik[i] * T.Mpesos[i][k] );
			for(int k = 0; k < T.vecinos[T.Mvecinos[i][j]]; k++)
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.Mvecinos[i][j] * N + T.Mvecinos[T.Mvecinos[i][j]][k], Elem_jk[T.Mvecinos[i][j]] * T.Mpesos[T.Mvecinos[i][j]][k]);
			for(int k = 0; k < T.vecinosT[i]; k++)
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[i][k] * N + i, Elem_ki[i] * T.MpesosT[i][k] );
			for(int k = 0; k < T.vecinosT[T.Mvecinos[i][j]]; k++)
					tripletList.emplace_back(i * N + T.Mvecinos[i][j], T.MvecinosT[T.Mvecinos[i][j]][k] * N + T.Mvecinos[i][j], Elem_kj[T.Mvecinos[i][j]] * T.MpesosT[T.Mvecinos[i][j]][k]);
		}
	}
	/*
	//Versión completa (mas sencilla de leer)
	for(int i = 0; i < T.N; i++){
		for(int j = 0; j < T.N; j++){
			for(int k = 0; k < T.vecinos[i]; k++)
			{
				//if(mPesosDense(i,j) != 0)//Eliminar las componentes sin poblacion
					tripletList.emplace_back(i * N + j, i * N + T.Mvecinos[i][k], Elem_ik[i] * T.Mpesos[i][k] );
			}
			for(int k = 0; k < T.vecinos[j]; k++){
				//if(mPesosDense(j,i) != 0)
					tripletList.emplace_back(i * N + j, j * N + T.Mvecinos[j][k], Elem_jk[j] * T.Mpesos[j][k]);
			}
			for(int k = 0; k < T.vecinosT[i]; k++)
			{
				//if(mPesosDense(i,j) != 0)//Eliminar las componentes sin poblacion
					tripletList.emplace_back(i * N + j, T.MvecinosT[i][k] * N + i, Elem_ki[i] * T.MpesosT[i][k] );
			}
			for(int k = 0; k < T.vecinosT[j]; k++){
				//if(mPesosDense(j,i) != 0)
					tripletList.emplace_back(i * N + j, T.MvecinosT[j][k] * N + j, Elem_kj[j] * T.MpesosT[j][k]);
			}
		}
	}
	*/
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