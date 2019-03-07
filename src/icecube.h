#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <getopt.h>

#include "nuSQuIDS/nuSQuIDS.h"

#include "TStopwatch.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"

std::array < std::array < float,40 >, 150 > nuflux_init_pion, nuflux_init_kaon, nubarflux_init_pion, nubarflux_init_kaon;
std::array < std::array < float,21 >, 10 > data, mc;

int nmc = 8784618;
bool debug;
float mcscale;

squids::Const units;
float Etrue_min(1.e2*units.GeV), Etrue_max(1.e6*units.GeV), CZ_min(-1.), CZ_max(.24),Ereco_min(400.*units.GeV),Ereco_max(20000*units.GeV);
std::vector < int > vec_mc_pid;
std::vector < float > vec_mc_reco_energy, vec_mc_reco_cos_zenith, vec_mc_true_nu_energy, vec_mc_true_cos_zenith, vec_mc_weight;
std::vector < float > vec_mc_flux_pion_osc, vec_mc_flux_kaon_osc;

std::vector < float > vec_mc_flux_pion_nom, vec_mc_flux_kaon_nom;

float logfactorial(int i);
double ICMinimizer(const double * X);
