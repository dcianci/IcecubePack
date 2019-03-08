#include "icecube.h"

using namespace nusquids;

int nuflux(){

	debug = false;
	mcscale = .1;

	// to save time, we can reduce our mc by a fraction and then scale the weights accordingly
	nmc *= mcscale;

	vec_mc_pid.resize(nmc);
	vec_mc_reco_energy.resize(nmc);	vec_mc_reco_cos_zenith.resize(nmc);
	vec_mc_true_nu_energy.resize(nmc);	vec_mc_true_cos_zenith.resize(nmc);
	vec_mc_weight.resize(nmc);
	vec_mc_flux_kaon_osc.resize(nmc);
	vec_mc_flux_pion_osc.resize(nmc);

	vec_mc_flux_kaon_nom.resize(nmc);
	vec_mc_flux_pion_nom.resize(nmc);
	

////////////////////////////////////////////////////////
	// First, let's read everything in
	////////////////////////////////////////////////////////

	std::cout << "Load up initial fluxes" << std::endl;
	// Load in atmospheric fluxes
	float dummy;
	std::ifstream file;
	file.open("data/initial_pion_atmopheric_HondaGaisser.dat");
	for(int cz = 0; cz < 40; cz++){				// cos zenith angle
		for(int te = 0; te < 150; te++){		// true nu energy
			file >> dummy;
			file >> dummy;
			file >> nuflux_init_pion[te][cz];
			file >> nubarflux_init_pion[te][cz];
		}
	}
	file.close();
	file.open("data/initial_kaon_atmopheric_HondaGaisser.dat");
	for(int cz = 0; cz < 40; cz++){
		for(int te = 0; te < 150; te++){
			file >> dummy;
			file >> dummy;
			file >> nuflux_init_kaon[te][cz];
			file >> nubarflux_init_kaon[te][cz];
		}
	}
	file.close();

	if(debug)
		std::cout << "Hope none of these are zero! " << dummy << " " << nuflux_init_pion[0][0] << " " << nuflux_init_kaon[0][0] << " " << nubarflux_init_pion[0][0]  << " " << nubarflux_init_kaon[0][0] << std::endl;

	// Load in data!
	file.open("data/data_histo.txt");
	for(int cz = 0; cz < 21; cz ++){			// reco cos zenith
		for(int re = 0; re < 10; re ++){		// reco energy
			file >> data[re][cz];
		}
	}
	file.close();

	// load in mc!
	file.open("data/NuFSGenMC_nominal.dat");
	for(int i = 0; i < nmc; i++){
 		file >> vec_mc_pid[i];
		file >> vec_mc_reco_energy[i];
		file >> vec_mc_reco_cos_zenith[i];
		file >> vec_mc_true_nu_energy[i];
		file >> vec_mc_true_cos_zenith[i];
		file >> vec_mc_weight[i];
		file >> vec_mc_flux_pion_nom[i];
		file >> vec_mc_flux_kaon_nom[i];
	}
	file.close();

	if(debug)
		std::cout << vec_mc_pid[0] << " " << vec_mc_reco_energy[0] << " " << vec_mc_reco_cos_zenith[0] << " " << vec_mc_true_nu_energy[0] << " " << vec_mc_true_cos_zenith[0] << " " << vec_mc_weight[0] << std::endl;


	//////////////////////////////////////////////////////////////
	// Propagate your fluxes through the earth!
	//////////////////////////////////////////////////////////////

	std::cout << "Create initial states" << std::endl;
	marray< double,4 > inistate_pions{40,150,2,3}, inistate_kaons{40,150,2,3};
	std::fill(inistate_pions.begin(),inistate_pions.end(),0);
	std::fill(inistate_kaons.begin(),inistate_kaons.end(),0);
	for(int ci = 0; ci < 40; ci++){
		for(int ei = 0; ei <  150; ei++){
			for(int flv = 0; flv < 3; flv++){
				inistate_pions[ci][ei][0][flv] = (flv == 1) ? nuflux_init_pion[ei][ci] : 0.0;	// set nonzero only for muons
				inistate_pions[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_pion[ei][ci] : 0.0;	// same, but for nubars
				inistate_kaons[ci][ei][0][flv] = (flv == 1) ? nuflux_init_kaon[ei][ci] : 0.0;
				inistate_kaons[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_kaon[ei][ci] : 0.0;
			}
		}
	}

	std::cout << "Construct nuSQuIDS-Atm object" << std::endl;
	nuSQUIDSAtm<> nus_atm_pion(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),3,both,true);
	nuSQUIDSAtm<> nus_atm_kaon(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),3,both,true);

	// Setup integration precision
	nus_atm_pion.Set_rel_error(1.0e-6);								nus_atm_kaon.Set_rel_error(1.0e-6);
	nus_atm_pion.Set_abs_error(1.0e-6);								nus_atm_kaon.Set_abs_error(1.0e-6);
	nus_atm_pion.Set_GSL_step(gsl_odeiv2_step_rk4);		nus_atm_kaon.Set_GSL_step(gsl_odeiv2_step_rk4);

	// Setup SM mixing params
	nus_atm_pion.Set_MixingParametersToDefault();
	nus_atm_kaon.Set_MixingParametersToDefault();

	if(debug)
		nus_atm_pion.Set_ProgressBar(true);
	else
		nus_atm_pion.Set_ProgressBar(false);
	nus_atm_pion.Set_IncludeOscillations(true);
	nus_atm_kaon.Set_IncludeOscillations(true);

	std::cout << "Setting initial pion state for nuSQuIDS" << std::endl;
	nus_atm_pion.Set_initial_state(inistate_pions,flavor);
	nus_atm_kaon.Set_initial_state(inistate_kaons,flavor);

	std::cout << "Evolving PIONS" << std::endl;
	nus_atm_pion.EvolveState();
	std::cout << "Evolving KAONS" << std::endl;
	nus_atm_kaon.EvolveState();						

	// So we don't have to evaluate the flavors for the flux repeatedly for the fit, let's liberate them all Now
	for(int i = 0; i < nmc; i++){
		if(vec_mc_pid[i] > 0){
			vec_mc_flux_kaon_osc[i] = nus_atm_kaon.EvalFlavor(1,vec_mc_true_cos_zenith[i],vec_mc_true_nu_energy[i]*units.GeV,0);
			vec_mc_flux_pion_osc[i] = nus_atm_pion.EvalFlavor(1,vec_mc_true_cos_zenith[i],vec_mc_true_nu_energy[i]*units.GeV,0);
		}
		else{
			vec_mc_flux_kaon_osc[i] = nus_atm_kaon.EvalFlavor(1,vec_mc_true_cos_zenith[i],vec_mc_true_nu_energy[i]*units.GeV,1);
			vec_mc_flux_pion_osc[i] = nus_atm_pion.EvalFlavor(1,vec_mc_true_cos_zenith[i],vec_mc_true_nu_energy[i]*units.GeV,1);
		}
	}

	////////////////////////////////////////////////////////////////
	// Now, we can perform the fit!!!!
	////////////////////////////////////////////////////////////////

	std::cout << "now do the fit!" << std::endl;
	// fit parameters. we'll just set 'em manually for now.
	float R_kpi(1), R_nunubar(1), Knorm(1), eff(.99), CRindex(0);

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin","BFGS2");
	min->SetMaxIterations(1000);  // for GSL
	min->SetPrintLevel(1);
	min->SetPrecision(0.001);			//times 4 for normal
	min->SetTolerance(0.01);

	ROOT::Math::Functor f(&ICMinimizer, 5);
	min->SetFunction(f);

	// set vars
	// for now, just throw some fixed variables to make sure it's all good.
	min->SetVariable(0,"R_kpi",1., .01);
	min->SetVariableLimits(0,.7,1.3);
	min->SetVariable(1,"R_nunubar",1,.005);
	min->SetVariableLimits(1,.925,1.075);
	min->SetVariable(2,"K_norm",1,.01);
	min->SetVariableLimits(2,.9,1.1);
	min->SetFixedVariable(3,"eff",.99);
	min->SetVariable(4,"CRindex",0,.01);
	min->SetVariableLimits(4,0,.15);
	min->Minimize();

	const double *xs = min->X();
	std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "," << xs[3] << "," << xs[4] << "): "  << min->MinValue()  << std::endl;


	// Okay, let's see how things look after the fit?
	for(int cz = 0; cz < 21; cz++)
		for(int re = 0; re < 10; re++)
			mc[re][cz] = 0;

	float flux;
	int ei,czi;
	for(int i = 0; i < nmc; i++){
		// figure out bin coords
		ei = floor(10.f * log10(vec_mc_reco_energy[i]*units.GeV/Ereco_min) / log10(Ereco_max/Ereco_min));
		czi = floor((vec_mc_reco_cos_zenith[i] + 1.02)/.06);

		// get our flux
		if(vec_mc_pid[i] > 0)
			flux = xs[2] * (vec_mc_flux_pion_osc[i] + vec_mc_flux_kaon_osc[i] * xs[0] ) * pow(vec_mc_true_nu_energy[i]*units.GeV,-xs[4]);
		else
			flux = xs[2] * xs[1] * (vec_mc_flux_pion_osc[i] + vec_mc_flux_kaon_osc[i] * xs[0] ) * pow(vec_mc_true_nu_energy[i]*units.GeV,-xs[4]);
			
		// add the number of events to the bin!
		mc[ei][czi] += flux * vec_mc_weight[i] / mcscale;
	}
	for(int cz = 0; cz < 21; cz++)
		for(int re = 0; re < 10; re++)
			std::cout << "D: " << re << " " << cz << " " << data[re][cz] << " " <<  mc[re][cz] << std::endl;

	return 1;
}

double ICMinimizer(const double * X){

	// X = {R_kpi, R_nunubar, Knorm, eff, CRindex;}

	float flux;
	int ei, czi;

	// fill mc  with zeros
	for(int cz = 0; cz < 21; cz++)
		for(int re = 0; re < 10; re++)
			mc[re][cz] = 0;

	for(int i = 0; i < nmc; i++){
		// figure out bin coords
		ei = floor(10.f * log10(vec_mc_reco_energy[i]*units.GeV/Ereco_min) / log10(Ereco_max/Ereco_min));
		czi = floor((vec_mc_reco_cos_zenith[i] + 1.02)/.06);

		// get our flux
		if(vec_mc_pid[i] > 0)
			flux = X[2] * (vec_mc_flux_pion_osc[i] + vec_mc_flux_kaon_osc[i] * X[0] ) * pow(vec_mc_true_nu_energy[i]*units.GeV,-X[4]);
		else
			flux = X[2] * X[1] * (vec_mc_flux_pion_osc[i] + vec_mc_flux_kaon_osc[i] * X[0] ) * pow(vec_mc_true_nu_energy[i]*units.GeV,-X[4]);
			
		// add the number of events to the bin!
		mc[ei][czi] += flux * vec_mc_weight[i] / mcscale;
	}

	// next step is to apply dom efficiency, but i'll do  that laaater
	// Okay. Now, let's calculate the log likelihood.
	float loglikelihood = 0;
	for(int cz = 0; cz < 21; cz++){
		for(int re = 0; re < 10; re++){	
			std::cout << "B: " << re << " " << cz << " " << data[re][cz] << " " << mc[re][cz] << std::endl;
			float l1(0),l2(0),l3(0);
			if(mc[re][cz] > 0)
				l1 = data[re][cz] * log(mc[re][cz]);
			l2 = mc[re][cz];
			if(data[re][cz] > 0)
				l3 = logfactorial(data[re][cz]);
			
			std::cout << "l1: " << l1 << " l2: " << l2  << " l3: " <<  l3 << std::endl;
			loglikelihood += -l1 + l2 + l3;	
		}
	}
	// Now add the penalty terms!
	loglikelihood += .5 * pow(X[0]-1,2) / pow(.1,2);
	loglikelihood += .5 * pow(X[1]-1,2) / pow(.05,2);
	loglikelihood += .5 * pow(X[4]-0,2) / pow(.05,2);
	std::cout << "PENALTY: " << .5 * pow(X[0]-1,2) / pow(.1,2) << " " << .5 * pow(X[1]-1,2) / pow(.05,2) << " " << .5 * pow(X[4]-0,2) / pow(.05,2)  << std::endl;


	std::cout << "BOOP" << std::endl;
	std::cout << X[0] << " " << X[1]  << " " << X[2] << " " << X[3]  << " " << X[4] << std::endl;
	std::cout << loglikelihood << std::endl;
	return loglikelihood;
}

float logfactorial(int i){
	float val(0.f);
	for(int j = i; j > 0; j--){
		val += log(j);
	}
	return val;
}

int main(){

	nuflux();
	return 1;
}
