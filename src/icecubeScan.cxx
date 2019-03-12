#include "icecube.h"

using namespace nusquids;

int nuflux(int massIndex){
	
	bool debug = true;
	mcscale = .1;
	nmc *= mcscale;
	int mgrd(10),ugrd(10),dgrd(1);

	float dm241;
	dm241 = pow(10,(massIndex/float(mgrd)*log10(10./.1) + log10(.1)));

	vec_mc_pid.resize(nmc);
	vec_mc_reco_energy.resize(nmc);	vec_mc_reco_cos_zenith.resize(nmc);
	vec_mc_true_nu_energy.resize(nmc);	vec_mc_true_cos_zenith.resize(nmc);
	vec_mc_weight.resize(nmc);
	vec_mc_flux_kaon_osc.resize(nmc);
	vec_mc_flux_pion_osc.resize(nmc);

	////////////////////////////////////////////////////////
	// First, let's read everything in
	////////////////////////////////////////////////////////

	std::cout << "Load up initial fluxes" << std::endl;
	// Load in atmospheric fluxes
	std::ifstream file;
	float dummy;
	file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/initial_pion_atmopheric_HondaGaisser.dat");
	for(int cz = 0; cz < 40; cz++){				// cos zenith angle
		for(int te = 0; te < 150; te++){		// true nu energy
			file >> dummy;
			file >> dummy;
			file >> nuflux_init_pion[te][cz];
			file >> nubarflux_init_pion[te][cz];
		}
	}
	file.close();
	file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/initial_kaon_atmopheric_HondaGaisser.dat");
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
	file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/data_histo.txt");
	for(int cz = 0; cz < 21; cz ++){			// reco cos zenith
		for(int re = 0; re < 10; re ++){		// reco energy
			file >> data[re][cz];
		}
	}
	file.close();

	// load in mc!
	file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/NuFSGenMC_nominal.dat");
	for(int i = 0; i < nmc; i++){
		file >> vec_mc_pid[i];
		file >> vec_mc_reco_energy[i];
		file >> vec_mc_reco_cos_zenith[i];
		file >> vec_mc_true_nu_energy[i];
		file >> vec_mc_true_cos_zenith[i];
		file >> vec_mc_weight[i];
		file >> dummy;
		file >> dummy;
	}
	file.close();

	if(debug)
		std::cout << vec_mc_pid[0] << " " << vec_mc_reco_energy[0] << " " << vec_mc_reco_cos_zenith[0] << " " << vec_mc_true_nu_energy[0] << " " << vec_mc_true_cos_zenith[0] << " " << vec_mc_weight[0] << std::endl;

	std::cout << "Create output file: " << "/nevis/amsterdam/share/dcianci/out/icecube-likelihood_mi-"+std::to_string(massIndex)+".txt"  << std::endl;
	std::ofstream outfile;
	outfile.open("/nevis/amsterdam/share/dcianci/out/icecube-likelihood_mi-"+std::to_string(massIndex)+".txt");
	
	//////////////////////////////////////////////////////////////
	// Propagate your fluxes through the earth!
	//////////////////////////////////////////////////////////////

	std::cout << "Create initial states" << std::endl;
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

	std::cout << "Construct nuSQuIDS-Atm object" << std::endl;
	nuSQUIDSAtm<> nus_atm_pion(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),4,both,true);
	nuSQUIDSAtm<> nus_atm_kaon(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),4,both,true);

	// Setup integration precision
	nus_atm_pion.Set_rel_error(1.0e-6);								nus_atm_kaon.Set_rel_error(1.0e-6);
	nus_atm_pion.Set_abs_error(1.0e-6);								nus_atm_kaon.Set_abs_error(1.0e-6);
	nus_atm_pion.Set_GSL_step(gsl_odeiv2_step_rk4);		nus_atm_kaon.Set_GSL_step(gsl_odeiv2_step_rk4);

	// Setup SM mixing params
	nus_atm_pion.Set_MixingAngle(0,1,0.563942);				nus_atm_kaon.Set_MixingAngle(0,1,0.563942);
	nus_atm_pion.Set_MixingAngle(0,2,0.154085);				nus_atm_kaon.Set_MixingAngle(0,2,0.154085);
	nus_atm_pion.Set_MixingAngle(1,2,0.785398);				nus_atm_kaon.Set_MixingAngle(1,2,0.785398);
	nus_atm_pion.Set_SquareMassDifference(1,7.65e-05);nus_atm_kaon.Set_SquareMassDifference(1,7.65e-05);
	nus_atm_pion.Set_SquareMassDifference(2,0.00247);	nus_atm_kaon.Set_SquareMassDifference(2,0.00247);
	nus_atm_pion.Set_CPPhase(0,2,0);									nus_atm_kaon.Set_CPPhase(0,2,0);

	if(debug)
		nus_atm_pion.Set_ProgressBar(true);
	else
		nus_atm_pion.Set_ProgressBar(false);
	nus_atm_pion.Set_IncludeOscillations(true);
	nus_atm_kaon.Set_IncludeOscillations(true);



	// Here's where we'll separate out our oscillations..
	int ct = 1;
	for(int uei = 0; uei < ugrd; uei++){
		for(int umi = 0; umi < ugrd; umi++){
			for(int uti = 0; uti < 1; uti++){

				float u14, u24, u34;
				u14 = pow(10,(uei/float(ugrd)*log10(1./1e-3) + log10(1e-3)));
				u24 = pow(10,(umi/float(ugrd)*log10(1./1e-3) + log10(1e-3)));
				u34 = pow(10,(uti/float(ugrd)*log10(1./1e-3) + log10(1e-3)));
         
				float theta14 = asin(u14);
				float theta24 = asin(u24/cos(theta14));
				float theta34 = acos(u34/(cos(theta14)*cos(theta24)));

				for(int d24i = 0; d24i < dgrd; d24i++){
					for(int d34i = 0; d34i < dgrd; d34i++){
						std::cout << "Running model " << ct << " of " << ugrd * ugrd * ugrd * dgrd * dgrd  << std::endl;						
						

						float d24,d34;
						d24 = d24i/float(dgrd)*(2*3.1415926);
						d34 = d34i/float(dgrd)*(2*3.1415926);

						std::cout << "Setting new sterile parameters" << std::endl;
						nus_atm_pion.Set_SquareMassDifference(3,dm241);
						nus_atm_pion.Set_MixingAngle(0,3,theta14);
						nus_atm_pion.Set_MixingAngle(1,3,theta24);
						nus_atm_pion.Set_MixingAngle(2,3,theta34);
						nus_atm_pion.Set_CPPhase(1,3,d24);
						nus_atm_pion.Set_CPPhase(2,3,d34);

						nus_atm_kaon.Set_SquareMassDifference(3,dm241);
						nus_atm_kaon.Set_MixingAngle(0,3,theta14);
						nus_atm_kaon.Set_MixingAngle(1,3,theta24);
						nus_atm_kaon.Set_MixingAngle(2,3,theta34);
						nus_atm_kaon.Set_CPPhase(1,3,d24);
						nus_atm_kaon.Set_CPPhase(2,3,d34);


						std::cout << "Setting initial pion state for nuSQuIDS" << std::endl;
						nus_atm_pion.Set_initial_state(inistate_pions,flavor);
						nus_atm_kaon.Set_initial_state(inistate_kaons,flavor);

						std::cout << "Evolving states!" << std::endl;
						nus_atm_pion.EvolveState();
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

						ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Fumili2");
						min->SetMaxIterations(1000);  // for GSL
						min->SetPrintLevel(3);
						min->SetPrecision(0.001);			//times 4 for normal
						min->SetTolerance(.001);

						ROOT::Math::Functor f(&ICMinimizer, 5);
						min->SetFunction(f);

						// set vars
						// for now, just throw some fixed variables to make sure it's all good.
						min->SetVariable(0,"R_kpi",1.,.01);
						min->SetVariableLimits(0,.9,1.3);
						min->SetVariable(1,"R_nunubar",1,.01);
						min->SetVariableLimits(1,1,1.4);
						min->SetVariable(2,"K_norm",1,.05);
						min->SetVariableLimits(2,1,1.2);
						min->SetFixedVariable(3,"eff",.99);
						min->SetVariable(4,"CRindex",0,.01);
						min->SetVariableLimits(4,0,.15);
						min->Minimize();

						std::cout << "tol: " << min->Tolerance() << " prec: " << min->Precision() << std::endl;

						const double *xs = min->X();
						std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "," << xs[3] << "," << xs[4] << "): "  << min->MinValue()  << std::endl;

						// Write out our answer to a txt file!
						outfile << "logl: " << dm241 << " " << u14 << " " << u24 << " "  << u34 << " " << d24 << " " << d34  << " " << min->MinValue() << "\n";
						ct ++;
					}
				}
			}
		}
	}

	outfile.close();

	return  1;
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
			//std::cout << "B: " << re << " " << cz << " " << data[re][cz] << " " << mc[re][cz] << std::endl;
			float l1(0),l2(0),l3(0);
			if(mc[re][cz] > 0)
				l1 = data[re][cz] * log(mc[re][cz]);
			l2 = mc[re][cz];
			if(data[re][cz] > 0)
				l3 = logfactorial(data[re][cz]);
			
			//std::cout << "l1: " << l1 << " l2: " << l2  << " l3: " <<  l3 << std::endl;
			loglikelihood += -l1 + l2 + l3;	
		}
	}
	// Now add the penalty terms!
	loglikelihood += .5 * pow(X[0]-1,2) / pow(.1,2);
	loglikelihood += .5 * pow(X[1]-1,2) / pow(.05,2);
	loglikelihood += .5 * pow(X[4]-0,2) / pow(.05,2);
	//std::cout << "PENALTY: " << .5 * pow(X[0]-1,2) / pow(.1,2) << " " << .5 * pow(X[1]-1,2) / pow(.05,2) << " " << .5 * pow(X[4]-0,2) / pow(.05,2)  << std::endl;

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

int main(int argc, char* argv[]){

  int iarg = 0;
  opterr=1;
  int index = 0;

  const struct option longopts[] = {
    {"massindex", optional_argument,  0,  'm'},
	  {0,			no_argument, 		0,  0},
  };

  while(iarg != -1){
    iarg = getopt_long(argc,argv, "m:", longopts, &index);

    switch(iarg){
		  case 'm':
			  index = std::atoi(optarg);
			  break;
		}
  }

  nuflux(index);
  return 0;
}
