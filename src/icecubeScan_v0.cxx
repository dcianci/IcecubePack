// This version of the code is specifically for producing the ntuples needed for comparison to DUNE
// That is to say, these ntuples won't be very good for plotting (ie: you'll want to jigger with dimensionality for that)
// At this point, validation is done and we ain't changing the machinery

#include "icecube.h"

using namespace nusquids;

int nuflux(int stindex){
	
	debug = true;
	int numcores(160);

	const int numneu = 4;

	// 'final' distributions
	std::array < std::array < double, 21 >, 10> data_2d;
	std::array < std::array < double, 21 >, 10> flux_nu_kaon_2d, flux_nubar_kaon_2d, flux_nu_pion_2d,flux_nubar_pion_2d;

	// MC holders (these boys are too big for the stack so  we heapin' them)
	std::vector < int > mc_pid;
	std::vector < double > mc_emu_reco, mc_enu_true, mc_cosz_reco, mc_cosz_true, mc_weight;
	mc_pid.resize(nMC); 
	mc_emu_reco.resize(nMC);
  mc_enu_true.resize(nMC); 
  mc_cosz_reco.resize(nMC);
  mc_cosz_true.resize(nMC);
  mc_weight.resize(nMC);

	// flux init
	std::array < std::array < double, 40 >, 150 > nuflux_init_pion, nuflux_init_kaon, nubarflux_init_pion, nubarflux_init_kaon;

	////////////////////////////////////////////////////////
	// First, let's read everything in
	////////////////////////////////////////////////////////

	std::cout << "Load up initial fluxes" << std::endl;
	// Load in atmospheric fluxes
	std::ifstream file;
	float dummy;
	//file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/initial_pion_atmopheric_HondaGaisser.dat");
	file.open("initial_pion_atmopheric_HondaGaisser.dat");
	for(int cz = 0; cz < 40; cz++){				// cos zenith angle
		for(int te = 0; te < 150; te++){		// true nu energy
			file >> dummy;
			file >> dummy;
			file >> nuflux_init_pion[te][cz];
			file >> nubarflux_init_pion[te][cz];
		}
	}
	file.close();
	//file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/initial_kaon_atmopheric_HondaGaisser.dat");
	file.open("initial_kaon_atmopheric_HondaGaisser.dat");
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
	//file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/data_histo.txt");
	file.open("data_histo.txt");
	for(int cz = 0; cz < 21; cz ++){			// reco cos zenith
		for(int re = 0; re < 10; re ++){		// reco energy
			file >> data_2d[re][cz];
		}
	}
	file.close();

	// load in mc!
	std::cout << "Load the MC! (these are  chonky)" << std::endl;
	//file.open("/nevis/amsterdam/share/dcianci/IcecubePack/data/NuFSGenMC_nominal.dat");
	file.open("NuFSGenMC_nominal.dat");
	for(int i = 0; i < nMC; i++){
		file >> mc_pid[i];
		file >> mc_emu_reco[i];
		file >> mc_cosz_reco[i];
		file >> mc_enu_true[i];
		file >> mc_cosz_true[i];
		file >> mc_weight[i];
		file >> dummy;
		file >> dummy;
	}
	file.close();

	if(debug)
		std::cout << mc_pid[20] << std::endl;

	//////////////////////////////////////////////////////////////
	// Propagate your fluxes through the earth!
	//////////////////////////////////////////////////////////////

	std::cout << "Create initial states" << std::endl;
	marray< double, 4 > inistate_pions{40,150,2,numneu}, inistate_kaons{40,150,2,numneu};
	std::fill(inistate_pions.begin(),inistate_pions.end(),0);
	std::fill(inistate_kaons.begin(),inistate_kaons.end(),0);
	for(int ci = 0; ci < 40; ci++){
		for(int ei = 0; ei <  150; ei++){
			for(int flv = 0; flv < numneu; flv++){
				inistate_pions[ci][ei][0][flv] = (flv == 1) ? nuflux_init_pion[ei][ci] * 1.0e20 : 0.0;	// set nonzero only for muons; add kaon and pion fluxes together
				inistate_pions[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_pion[ei][ci] * 1.0e20 : 0.0;	// same, but for nubars
				inistate_kaons[ci][ei][0][flv] = (flv == 1) ? nuflux_init_kaon[ei][ci] * 1.0e20 : 0.0;
				inistate_kaons[ci][ei][1][flv] = (flv == 1) ? nubarflux_init_kaon[ei][ci] * 1.0e20 : 0.0;
			}
		}
	}

	std::cout << "Construct nuSQuIDS-Atm object" << std::endl;
	nuSQUIDSAtm<> nus_atm_pion(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),numneu,both,true);
	nuSQUIDSAtm<> nus_atm_kaon(linspace(CZ_min,CZ_max,40),logspace(Etrue_min,Etrue_max,150),numneu,both,true);
	
	// Setup integration precision
	nus_atm_pion.Set_rel_error(1.0e-3);								nus_atm_kaon.Set_rel_error(1.0e-3);
	nus_atm_pion.Set_abs_error(1.0e-3);								nus_atm_kaon.Set_abs_error(1.0e-3);
	nus_atm_pion.Set_GSL_step(gsl_odeiv2_step_rk4);		nus_atm_kaon.Set_GSL_step(gsl_odeiv2_step_rk4);

	// Setup SM mixing params
	float theta12(0.60214), theta13(0.147480), theta23(0.832522), dm21(7.55e-5), dm22(2.50e-3), delta13(1.32);
	nus_atm_pion.Set_MixingAngle(0,1,theta12);				nus_atm_kaon.Set_MixingAngle(0,1,theta12);
	nus_atm_pion.Set_MixingAngle(0,2,theta13);				nus_atm_kaon.Set_MixingAngle(0,2,theta13);
	nus_atm_pion.Set_MixingAngle(1,2,theta23);				nus_atm_kaon.Set_MixingAngle(1,2,theta23);
	nus_atm_pion.Set_SquareMassDifference(1,dm21);		nus_atm_kaon.Set_SquareMassDifference(1,dm21);
	nus_atm_pion.Set_SquareMassDifference(2,dm22);		nus_atm_kaon.Set_SquareMassDifference(2,dm22);
	nus_atm_pion.Set_CPPhase(0,2,theta13);						nus_atm_kaon.Set_CPPhase(0,2,delta13);

	if(debug)
		nus_atm_pion.Set_ProgressBar(true);
	else
		nus_atm_pion.Set_ProgressBar(false);
	nus_atm_pion.Set_IncludeOscillations(true);
	nus_atm_kaon.Set_IncludeOscillations(true);

	

	// Here's where we'll separate out our oscillations..
	// We've decided to stick with angles for the whole thing. Ok.
	// Welcome to dimension hell motherfucker
	int dm2grd(40), th14grd(15), th24grd(40), th34grd(10);
	//int npoints = dm2grd * th14grd * th24grd * th34grd;
	int npoints = dm2grd * th24grd;
	int stoffset = (npoints/numcores)*stindex;
	int ind, im, ith1, ith2, ith3;

	std::cout << "debug: dm2 theta14 theta24 theta34 delta delta logl_min" << std::endl;

	float m41,theta14,theta24,theta34,d24,d34;
	
	for(int univi = 0; univi < npoints/numcores; univi++){

		ind = univi+stoffset;
		im = ind % dm2grd;
		ith2 = ((ind-im) % (dm2grd*th24grd))/dm2grd;
		ith2 = ((ind-im-ith1) % (dm2grd*th14grd*th24grd))/(dm2grd*th14grd);
		ith3 = ((ind-im-ith1-ith2) % (dm2grd*th14grd*th24grd*th34grd))/(dm2grd*th14grd*th24grd);
		std::cout << "DIM: " << im << " " << ith1 << " " << ith2 << " " << ith3 << std::endl;
	
		m41 = pow(10,(im+.5)/float(dm2grd)*log10(mnu_hibound/mnu_lowbound) + log10(mnu_lowbound));
		theta14 = pow(10,(ith1-1+.5)/float(th14grd-2)*log10(theta_hibound/theta_lowbound) + log10(theta_lowbound));
		theta24 = pow(10,(ith2-1+.5)/float(th24grd-2)*log10(theta_hibound/theta_lowbound) + log10(theta_lowbound));
		theta34 = pow(10,(ith3+.5)/float(th34grd)*log10(theta_hibound/theta_lowbound) + log10(theta_lowbound));

	  	std::cout << "deborg: m41: " << m41 << " theta14: " << theta14 << "  theta24: " << theta24 << " theta34: " << theta34 << std::endl;

		d24 = 0;
		d34 = 0;

		std::cout << "Setting new sterile parameters" << std::endl;
		nus_atm_pion.Set_SquareMassDifference(3,pow(m41,2));
		nus_atm_pion.Set_MixingAngle(0,3,theta14);
		nus_atm_pion.Set_MixingAngle(1,3,theta24);
		nus_atm_pion.Set_MixingAngle(2,3,theta34);
		nus_atm_pion.Set_CPPhase(1,3,d24);
		nus_atm_pion.Set_CPPhase(2,3,d34);

		nus_atm_kaon.Set_SquareMassDifference(3,pow(m41,2));
		nus_atm_kaon.Set_MixingAngle(0,3,theta14);
		nus_atm_kaon.Set_MixingAngle(1,3,theta24);
		nus_atm_kaon.Set_MixingAngle(2,3,theta34);
		nus_atm_kaon.Set_CPPhase(1,3,d24);
		nus_atm_kaon.Set_CPPhase(2,3,d34);

		std::cout << "Setting initial pion state for nuSQuIDS" << std::endl;
		nus_atm_pion.Set_initial_state(inistate_pions,flavor);
		nus_atm_kaon.Set_initial_state(inistate_kaons,flavor);

		std::cout << "Evolving Pion flux!" << std::endl;
		nus_atm_pion.EvolveState();
		std::cout << "Evolving Kaon flux!" <<  std::endl;
		nus_atm_kaon.EvolveState();

		// Okay, now let's liberate the fluxes and spit them all into two arrays in reco energy x cosz
		std::cout << "Liberate the fluxes" << std::endl;
		for(int cz = 0; cz < 21; cz++) for(int em = 0; em < 10; em++){
			flux_nu_kaon_2d[em][cz] = 0.0;
			flux_nubar_kaon_2d[em][cz] = 0.0;
			flux_nu_pion_2d[em][cz] = 0.0;
			flux_nubar_pion_2d[em][cz] = 0.0;
		}		
		int emui, czrecoi;
		for(int j = 0; j < nMC; j++){
			// find indices for reco matrix
			emui = floor(10.f * log10(mc_emu_reco[j]*units.GeV/Ereco_min) / log10(Ereco_max/Ereco_min));
			czrecoi = floor((mc_cosz_reco[j] + 1.02)/.06);

			if(mc_pid[j] > 0){
				flux_nu_kaon_2d[emui][czrecoi] += std::max(nus_atm_kaon.EvalFlavor(1,mc_cosz_true[j],mc_enu_true[j]*units.GeV,0),0.0)*mc_weight[j]/1.0e20;
				flux_nu_pion_2d[emui][czrecoi] += std::max(nus_atm_pion.EvalFlavor(1,mc_cosz_true[j],mc_enu_true[j]*units.GeV,0),0.0)*mc_weight[j]/1.0e20;
			}
			else{
				flux_nubar_kaon_2d[emui][czrecoi] += std::max(nus_atm_kaon.EvalFlavor(1,mc_cosz_true[j],mc_enu_true[j]*units.GeV,1),0.0)*mc_weight[j]/1.0e20;
				flux_nubar_pion_2d[emui][czrecoi] += std::max(nus_atm_pion.EvalFlavor(1,mc_cosz_true[j],mc_enu_true[j]*units.GeV,1),0.0)*mc_weight[j]/1.0e20;
			}
		}
		std::array < std::array < std::array < double, 21 >, 10 >, 4 > mc_fluxes = {flux_nu_pion_2d, flux_nu_kaon_2d,flux_nubar_pion_2d,flux_nubar_kaon_2d};

		////////////////////////////////////////////////////////////////
		// Now, we can perform the fit!!!!
		////////////////////////////////////////////////////////////////

		std::cout << "now do the fit!" << std::endl;

		// fit parameters. we'll just set 'em manually for now.
		double R_kpi(x0_start), R_nunubar(x1_start), Knorm(x2_start), logl;
		double R_kpi_old, R_nunuubar_old, Knorm_old, logl_old, logl_min(9999999999.0);    
		std::array < double,3 > xpar, xpar_min, xpar_old = { R_kpi, R_nunubar, Knorm }; 

		// Okay, we're going to make our own minimizer because fuuuuuck this
		// Davio's MCMC
		logl_old = ICMinimizer_2d(xpar_old,data_2d,mc_fluxes,false);

		logl_min = logl_old;
		xpar = xpar_old;		

		int numit = 10000;
		int nchains = 40;
		double par0,par1,par2;
		// let's complicate things a little bit...
		for(int nch = 0; nch < nchains; nch++){

			// after first chain, we'll start in random place
			if(nch > 0){
				par0 = x0_start + (2*RanGen.Rndm()-1) * 2 * x0_sigma;
				par1 = x1_start + (2*RanGen.Rndm()-1) * 2 * x1_sigma;
				par2 = x2_start + (2*RanGen.Rndm()-1) * 2 * x2_sigma;
			}
		
			for(int nit = 0; nit < numit; nit++){
				par0 = xpar_old[0] + (2 * RanGen.Rndm() - 1) * x0_sigma;
				par1 = xpar_old[1] + (2 * RanGen.Rndm() - 1) * x1_sigma;
				par2 = xpar_old[2] + (2 * RanGen.Rndm() - 1) * x2_sigma;

				xpar[0] = std::min(par0,x0_highbound);  xpar[0] = std::max(xpar[0],x0_lowbound);	
				xpar[1] = std::min(par1,x1_highbound);  xpar[1] = std::max(xpar[1],x1_lowbound);	
				xpar[2] = std::min(par2,x2_highbound);  xpar[2] = std::max(xpar[2],x2_lowbound);	

				logl = ICMinimizer_2d(xpar,data_2d,mc_fluxes,false);
				if(logl < logl_min){
					logl_min = logl;
					xpar_min = xpar;
				}
				// Assume infinite temperature for now
				if(logl < logl_old){
					logl_old = logl;
					xpar_old = xpar;
				}

				//std::cout << nit << " logl: (" << xpar[0] << ", " <<  xpar[1] << ", " << xpar[2] << "): " << logl << std::endl;
			}
		
					
		std::cout << "logl: " << pow(m41,2) << " " << theta14 << " " << theta24 << " "  << theta34 << " " << d24 << " " << d34  << " " << logl_min << std::endl;
	}
	return  1;
}

int main(int argc, char* argv[]){

  int iarg = 0;
  opterr=1;
  int index = 0;

  const struct option longopts[] = {
    {"stindex", optional_argument,  0,  'm'},
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
