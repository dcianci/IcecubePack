#include <algorithm>
#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <getopt.h>

#include "nuSQuIDS/nuSQuIDS.h"

#include "TRandom3.h"

bool debug;

const squids::Const units;
const double Etrue_min(1.e2*units.GeV), Etrue_max(1.e6*units.GeV), CZ_min(-1.), CZ_max(.2),Ereco_min(400.*units.GeV),Ereco_max(20000*units.GeV);

// mc events
const int nMC(8784618);

// Bounds for nuisance params
const double x0_start(1.0),x1_start(1.0),x2_start(1.0),x3_start(.99);
const double x0_sigma(.1),x1_sigma(.05),x2_step(.2),x3_step(.02);
const double x0_lowbound(-3), x0_highbound(3);
const double x1_lowbound(-3), x1_highbound(3);
const double x2_lowbound(.6), x2_highbound(1.4);
const double x3_lowbound(.95),x3_highbound(1.089);
const double x2_sigma(.4);
TRandom3 RanGen(0);

// Bounds for grid scan of params
const double mnu_lowbound(.1), mnu_hibound(10.);
const double theta_lowbound(.01), theta_hibound(3.1415926/4.);
const double u_lowbound(.001), u_hibound(.6);
const int mgrd(50),ugrd(50),dgrd(1);

// SM mixing params
const float theta12(0.60214), theta13(0.147480), theta23(0.832522), dm21(7.55e-5), dm22(2.50e-3), delta13(1.32);
const int numneu = 4;

double ICMinimizer_2d(const std::array<double,3>&_X,const std::array<std::array<double,21>,10>&_data, const std::array<std::array<std::array<double,21>,10>,4>&_mc_fluxes, bool printout = false){
	// _X = {R_kpi, R_nunubar, Knorm}
	// _mc_fluxes = {nu_pion, nu_kaon, nubar_pion, nubar_kaon}
	double R_kpi(_X[0]), R_nunubar(_X[1]), Knorm(_X[2]);
	std::array < std::array < double,21 >,10 > mc;

	for(int emi = 0; emi < 10; emi++){
		for(int czi = 0; czi < 21; czi++){
			double fluxnu = Knorm * (_mc_fluxes[0][emi][czi] + _mc_fluxes[1][emi][czi] * R_kpi);
			double fluxnubar = Knorm * R_nunubar * (_mc_fluxes[2][emi][czi] + _mc_fluxes[3][emi][czi] * R_kpi);
			mc[emi][czi] = fluxnu + fluxnubar;
		}
	}

  double loglikelihood = 0;
	for(int cz = 0; cz < 21; cz++){
		for(int re = 0; re < 10; re++){
			if(mc[re][cz] > 0 and _data[re][cz] > 0){
				loglikelihood += (_data[re][cz] * log(mc[re][cz]) - mc[re][cz]);// - (_data[re][cz] * log(_data[re][cz]) - _data[re][cz]);
				std::cout << "BOOP: " << _data[re][cz] << " " << mc[re][cz] << " "  << _data[re][cz] * log(mc[re][cz]) - mc[re][cz] << " " << _data[re][cz] * log(_data[re][cz]) - _data[re][cz] << std::endl;
			}
		}
	}
	// Now add the penalty terms!
	double pull0,pull1,pull2;
	pull0 = .5 * pow(R_kpi-x0_start,2) / pow(x0_sigma,2);
  pull1 = .5 * pow(R_nunubar-x1_start,2) / pow(x1_sigma,2);
  pull2 = 0;//.5 * pow(Knorm-x2_start,2) / pow(x2_sigma,2);
	if(printout)
		std::cout << "PULLS: " << pull0 << " + " << pull1 << " + " << pull2 << " = " <<  pull1+pull2+pull0 << std::endl;

	return (loglikelihood + pull0 + pull1 + pull2 - 90056.10375420621);// - 600;
}
