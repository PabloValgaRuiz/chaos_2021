#include "config.hpp"
#include "MobMatrix.hpp"
#include "MC_DistDiaNocheF.hpp"
#include "Mk_DistDiaNocheF.hpp"
#include "benchmark.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <future>

static void iteracionesMk(const MobMatrix& T, double lambda, double p, int tiempo);
static void iteracionesMC(const MobMatrix& T, double &infectados, double lambda, double p, int tiempo, int i);
int mkmc(int argc, char* argv[]);