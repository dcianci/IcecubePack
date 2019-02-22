#include "icecube.h"

using namespace nusquids;

int nuflux(){

	bool debug = false;

	std::cout << "Load up initial fluxes" << std::endl;
	
	std::array < float, 150 > vec_energy;
	std::array < float, 40 > vec_zenith;
	std::array < std::array < float,40 >, 150 > nuflux_init_pion, nuflux_init_kaon, nubarflux_init_pion, nubarflux_init_kaon;
	std::array < std::array < std::array < float,200 >,21 >, 10 > syst_resp_array_nu, syst_resp_array_nubar;
	std::array < std::array < float,21 >, 10 > nuflux_osc_pion, nuflux_osc_kaon, nubarflux_osc_pion, nubarflux_osc_kaon;

	////////////////////////////////////////////////////////
	// First, let's read everything in
	////////////////////////////////////////////////////////

	std::ifstream file;
	file.open("data/initial_pion_atmopheric_HondaGaisser.dat");
	for(int cosz = 0; cosz < 40; cosz++){
		for(int recoe = 0; recoe < 150; recoe++){
			file >> vec_zenith[cosz];
			file >> vec_energy[recoe];
			file >> nuflux_init_pion[recoe][cosz];
			file >> nubarflux_init_pion[recoe][cosz];
		}
	}
	file.close();
	file.open("data/initial_kaon_atmopheric_HondaGaisser.dat");
	for(int cosz = 0; cosz < 40; cosz++){
		for(int recoe = 0; recoe < 150; recoe++){
			file >> vec_zenith[cosz];
			file >> vec_energy[recoe];
			file >> nuflux_init_kaon[recoe][cosz];
			file >> nubarflux_init_kaon[recoe][cosz];
		}
	}
	file.close();

	if(debug)
		std::cout << "Hope none of these are zero! " << vec_zenith[0] << " " << vec_energy[0] << " " << nuflux_init_pion[0][0] << " " << nuflux_init_kaon[0][0] << " " << nubarflux_init_pion[0][0]  << " " << nubarflux_init_kaon[0][0] << std::endl;
	
	file.open("data/syst_resp_array_109_nu.txt");
	for(int pe = 0; pe < 10; pe ++){	// proxy (true) energy
		for(int cz = 0; cz < 21; cz ++){	// cos zenith angle
			for(int re = 0; re < 200; re ++){	// reco energy
				file >> syst_resp_array_nu_high[pe][cz][re];
			}
		}
	}
	file.close();

	file.open("data/syst_resp_array_109_nubar.txt");
	for(int pe = 0; pe < 10; pe ++){	// proxy (reco) energy
		for(int cz = 0; cz < 21; cz ++){	// cos zenith angle
			for(int re = 0; re < 200; re ++){	// true energy
				file >> syst_resp_array_nubar_high[pe][cz][re];
			}
		}
	}
	file.close();

	file.open("data/syst_resp_array_90_nu.txt");
	for(int pe = 0; pe < 10; pe ++){	// proxy (true) energy
		for(int cz = 0; cz < 21; cz ++){	// cos zenith angle
			for(int re = 0; re < 200; re ++){	// reco energy
				file >> syst_resp_array_nu_low[pe][cz][re];
			}
		}
	}
	file.close();

	file.open("data/syst_resp_array_90_nubar.txt");
	for(int pe = 0; pe < 10; pe ++){	// proxy (reco) energy
		for(int cz = 0; cz < 21; cz ++){	// cos zenith angle
			for(int re = 0; re < 200; re ++){	// true energy
				file >> syst_resp_array_nubar_low[pe][cz][re];
			}
		}
	}
	file.close();
	
	if(debug)
		std::cout << "Hope none of these are zero either! " << syst_resp_array_nu_high[0][0][0] << " " << syst_resp_array_nu_low[0][0][0] << " "  << syst_resp_array_nubar_low[0][0][0] << " " << syst_resp_array_nubar_high[0][0][0]  << std::endl;

	//////////////////////////////////////////////////////////////	
	// Propagate your fluxes through the earth!
	//////////////////////////////////////////////////////////////

	std::cout << "Create initial state" << std::endl;
	marray< double,4 > inistate_pions{40,150,2,4}, inistate_kaons{40,150,2,4};
	std::fill(inistate_pions.begin(),inistate_pions.end(),0);
	std::fill(inistate_kaons.begin(),inistate_kaons.end(),0);
	for(int ci = 0; ci < 40; ci++){
		for(int ei = 0; ei <  150; ei++){
			for(int flv = 0; flv < 4; flv++){
				inistate_pions[ci][ei][0][flv] = (flv == 1) ? nuflux_init_pion[ei][ci] : 0.0;	// set nonzero only for muons; add kaon and pion fluxes together
				inistate_pions[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_pion[ei][ci] : 0.0;	// same, but for nubars
				inistate_kaons[ci][ei][0][flv] = (flv == 1) ? nuflux_init_kaon[ei][ci] : 0.0;
				inistate_kaons[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_kaon[ei][ci] : 0.0;
			}
		}
	}
	
	std::cout << "Initialize nuSQuIDS" << std::endl;
	squids::Const units;

	double Emin = 2.e2*units.GeV;
	double Emax = 1.e6*units.GeV;
	double czmin = -1.;
	double czmax = .24;

	std::cout << "Construct nuSQuIDS-Atm object" << std::endl;
	nuSQUIDSAtm<> nus_atm(linspace(czmin,czmax,40),logspace(Emin,Emax,150),4,both,true);

	// Setup integration precision
	nus_atm.Set_rel_error(1.0e-6);
	nus_atm.Set_abs_error(1.0e-6);
	nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

	// Setup SM mixing params
	nus_atm.Set_MixingAngle(0,1,0.563942);
	nus_atm.Set_MixingAngle(0,2,0.154085);
	nus_atm.Set_MixingAngle(1,2,0.785398);
	nus_atm.Set_SquareMassDifference(1,7.65e-05);
	nus_atm.Set_SquareMassDifference(2,0.00247);
	nus_atm.Set_CPPhase(0,2,0);

	if(debug)
		nus_atm.Set_ProgressBar(true);
	else
		nus_atm.Set_ProgressBar(false);
	nus_atm.Set_IncludeOscillations(true);

	float dm241 = 1.2;
	float sin22th = .01;

	float theta24 = asin(sqrt(sin22th))/2;
	nus_atm.Set_SquareMassDifference(3,dm241);
	nus_atm.Set_MixingAngle(1,3,theta24);

	std::cout << "Setting initial pion state for nuSQuIDS" << std::endl;
	nus_atm.Set_initial_state(inistate_pions,flavor);

	std::cout << "Evolving states!" << std::endl;
	nus_atm.EvolveState();

	// Bin edges:
	int nZBins(21),nEBins(200);
	float centerE, centerCZ, binsE[nEBins+1],binsCosZ[nZBins+1];
	binsCosZ[0] = -1.f;
	for(int i = 0; i < nZBins+1; i++)
		binsCosZ[i] = std::max(-1.0,-1.02 + i*.06);
	for(int j = 0; j < nEBins+1; j++)
		binsE[j] = pow(10,log10(2e2) + j*log10(1e6/2e2)/nEBins);

	for(int cz = 0; cz < nZBins; cz++){
		for(int re = 0; re < nEBins; re++){
			centerCZ = binsCosZ[cz] + (binsCosZ[cz+1]-binsCosZ[cz])/2;
			centerE = pow(10,log10(binsE[re]) + log10(binsE[re+1]/binsE[re])/2.f);
			nuflux_osc_pion[cz][re] = nus_atm.EvalFlavor(1,centerCZ,centerE*units.GeV,0);
			nubarflux_osc_pion[cz][re] = nus_atm.EvalFlavor(1,centerCZ,centerE*units.GeV,1);
		}
	}

	std::cout << "Setting initial kaon state for nusquids" << std::endl;
	nus_atm.Set_initial_state(inistate_kaons,flavor);
	
	std::cout << "Evolving states!" << std::endl;
	nus_atm.EvolveState();
	
	for(int cz = 0; cz < nZBins; cz++){
		for(int re = 0; re < nEBins; re++){
			centerCZ = binsCosZ[cz] + (binsCosZ[cz+1]-binsCosZ[cz])/2;
			centerE = pow(10,log10(binsE[re]) + log10(binsE[re+1]/binsE[re])/2.f);
			nuflux_osc_kaon[re][cz] = nus_atm.EvalFlavor(1,centerCZ,centerE*units.GeV,0);
			nubarflux_osc_kaon[re][cz] = nus_atm.EvalFlavor(1,centerCZ,centerE*units.GeV,1);
		}
	}

	////////////////////////////////////////////////////////////////
	// Now, we can perform the fit!!!!
	////////////////////////////////////////////////////////////////


	std::array < std::array < float,21 >, 200 > nuflux_true, nubarflux_true;
	std::array < std::array < float,20 >, 10 > nuflux_reco, nubarflux_reco;	
	float etrue;

	// fit parameters. we'll just set 'em manually for now.
	float R_kpi(1), R_nunubar(1), Knorm(1), eff(.99), CRindex(0);
		
	// First, add our pion and kaon fluxes, with associated nuisance params
	for(int te = 0; te < nEBins; te++){
		etrue = pow(10,log10(binsE[te]) + log10(binsE[te+1]/binsE[te])/2.f);
		for(int cz = 0; cz < nZBins; cz++){
			nuflux_true[te][cz] = Knorm * ( nuflux_osc_pion[te][cz] + R_kpi * nuflux_osc_kaon[te][cz] ) * pow(etrue,-CRindex);
			nubarflux_true[te][cz] = Knorm * R_nunubar * ( nubarflux_osc_pion[te][cz] + R_kpi * nubarflux_osc_kaon[te][cz] ) * pow(etrue,-CRindex); 
		}
	}
	// next, convolve with systematic response matrices corresponding to our efficiency to get reco energy
	syst_resp_array_nu = interpSystRespArray_nu(eff);
	syst_resp_array_nubar = interpSystRespArray_nubar(eff);
	for(int cz = 0; cz < 21; cz++){
		for(int pe = 0; pe < 10; pe++){
			nuflux_reco[pe][cz] = 0;
			nubarflux_reco[pe][cz] = 0;
			for(int te = 0; te < 200; te++){
				nuflux_reco[pe][cz] += syst_resp_array_nu[pe][cz][te] * nuflux_true[te][cz];
				nubarflux_reco[pe][cz] += syst_resp_array_nubar[pe][cz][te] * nubarflux_true[te][cz];
			}
		}
	}
				
	// And now at long last, we can get the number of events per bin. whew.







	so, according to carlos' thesis,  the log likelihood is given by the 


}


std::array < std::array < std::array < float,200 >,21 >, 10 > interpSystRespArray_nu(float eff){
	std::array < std::array < std::array < float,200 >,21 >, 10 > resp;
	
	float higheff(1.089),loweff(.9);

	for(int pe = 0; pe < 10; pe++){
		for(int cz = 0; cz < 21; cz++){
			for(int te = 0; te < 200; te++){
				float slope = (syst_resp_array_nu_high[pe][cz][te]-syst_resp_array_nu_low[pe][cz][te])/(higheff-loweff);
				resp[pe][cz][te] = syst_resp_array_nu_low[pe][cz][te] + (eff - loweff) * slope;
			}
		}
	}
	
	return resp;
}

std::array < std::array < std::array < float,200 >,21 >, 10 > interpSystRespArray_nubar(float eff){
	std::array < std::array < std::array < float,200 >,21 >, 10 > resp;
	
	float higheff(1.089),loweff(.9);

	for(int pe = 0; pe < 10; pe++){
		for(int cz = 0; cz < 21; cz++){
			for(int te = 0; te < 200; te++){
				float slope = (syst_resp_array_nubar_high[pe][cz][te]-syst_resp_array_nubar_low[pe][cz][te])/(higheff-loweff);
				resp[pe][cz][te] = syst_resp_array_nubar_low[pe][cz][te] + (eff - loweff) * slope;
			}
		}
	}
	
	return resp;
}


int main(){

	nuflux();
	return 1;
}
