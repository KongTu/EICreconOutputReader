#include "pleaseIncludeMe.h"
int readSingleParticles(TString inname="input/input.root",TString outname="test"){

 	TString rec_file=inname;
	//ROOT::EnableImplicitMT(kNumThreads);
	ROOT::RDataFrame d("events", rec_file);

	auto d1 = d.Define("mult",getNtrk,{"ReconstructedChargedParticles"})
						 .Define("multMC",getNtrkMC,{"MCParticles"})
						 .Define("momentum",momenta_from_chargedparticles,{"ReconstructedChargedParticles"})
						 .Define("eta",getEta,{"momentum"})
						 .Define("pt",getPt,{"momentum"})
						 .Define("phi",getPhi,{"momentum"})
						 .Define("momentumMC",momenta_from_mcparticles,{"MCParticles"})
						 .Define("etaMC",getEta,{"momentumMC"})
						 .Define("ptMC",getPt,{"momentumMC"})
						 .Define("phiMC",getPhi,{"momentumMC"})
						 .Define("ptRes",pt_resolution,{"MCParticles","ReconstructedChargedParticles"})
						 ;

	TString output_name_dir = outname;
  	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

	auto h_mult_REC = d1.Histo1D({"h_mult_REC", "; N; counts", 10, -0.5, 9.5}, "mult");
	auto h_mult_MC = d1.Histo1D({"h_mult_MC", "; N; counts", 10, -0.5, 9.5}, "multMC");
	auto h_eta_REC = d1.Histo1D({"h_eta_REC", "; #eta; counts", 100, -5, 5}, "eta");
	auto h_pt_REC = d1.Histo1D({"h_pt_REC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "pt");
	auto h_phi_REC = d1.Histo1D({"h_phi_REC", "; #phi; counts", 100, -PI, PI}, "phi");
	auto h_eta_MC = d1.Histo1D({"h_eta_MC", "; #eta; counts", 100, -5, 5}, "etaMC");
	auto h_pt_MC = d1.Histo1D({"h_pt_MC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "ptMC");
	auto h_phi_MC = d1.Histo1D({"h_phi_MC", "; #phi; counts; counts", 100, -PI, PI}, "phiMC");
	auto h_pt_Res = d1.Histo1D({"h_pt_Res", "; Resolution; counts", 100, -1,1}, "ptRes");
	auto h_pt_Res2D = d1.Histo2D({"h_pt_Res2D", "; p_{T} (GeV/c); Resolution",100,0,5,100,-1,1},"pt","ptRes");

	//Q2,x_v cuts
	auto kineCut = [](const std::vector<double>& qsq, const std::vector<double>& y_elec) { 
		if(qsq.size()<1||y_elec.size()<1) return 0;
		if(qsq[0] > 1. && qsq[0] < 10. 
			&& y_elec[0] > 0.01 && y_elec[0] < 0.95) return 1;
		else return 0;
	};

	//MC & REC dis kinematics
	auto d2 = d.Define("scatMC",findScatElecMC,{"MCParticles"})
			   .Define("etaElecMC",getEta,{"scatMC"})
			   .Define("Q2elecMC",getQ2elec,{"scatMC"}).Define("YelecMC",getYelec,{"scatMC"}).Define("XelecMC",getXelec,{"scatMC"})
			   .Define("scatREC",findScatElecREC,{"EcalEndcapNClusters","ReconstructedChargedParticles"})
			   .Define("etaElecREC",getEta,{"scatREC"}).Define("EpzREC",getEpzREC,{"EcalEndcapNClusters","ReconstructedChargedParticles"})
			   .Define("Q2elecREC",getQ2elec,{"scatREC"}).Define("YelecREC",getYelec,{"scatREC"}).Define("XelecREC",getXelec,{"scatREC"})
			   .Filter(kineCut,{"Q2elecMC","YelecMC"}); 
			   ;

	auto h_Eta_Elect_MC = d2.Histo1D({"h_Eta_Elect_MC", "; #eta; counts", 150, -5, 10}, "etaElecMC");
	auto h_Q2elec_MC = d2.Histo1D({"h_Q2elec_MC", "; Q^{2}_{e}; counts", 100, 0,100}, "Q2elecMC");
	auto h_Yelec_MC = d2.Histo1D({"h_Yelec_MC", "; y_{e}; counts", 100, 0,1}, "YelecMC");
	auto h_Xelec_MC = d2.Histo1D({"h_Xelec_MC", "; x_{e}; counts", 1000, 0,1}, "XelecMC");
	auto h_Eta_Elect_REC = d2.Histo1D({"h_Eta_Elect_REC", "; #eta; counts", 150, -5, 10}, "etaElecREC");
	auto h_Q2elec_REC = d2.Histo1D({"h_Q2elec_REC", "; Q^{2}_{e}; counts", 100, 0,100}, "Q2elecREC");
	auto h_Yelec_REC = d2.Histo1D({"h_Yelec_REC", "; y_{e}; counts", 100, 0,1}, "YelecREC");
	auto h_Xelec_REC = d2.Histo1D({"h_Xelec_REC", "; x_{e}; counts", 1000, 0,1}, "XelecREC");
	auto h_Epz_REC = d2.Histo1D({"h_Epz_REC", "; E - P_{z} (GeV); counts", 200, 0,70}, "EpzREC");
	
	/*
	TODO: 
	1. NEED to add an association between rec and mc scat elec and plot background;
	2. NEED to reject those background based on pfRICH parametrization.
	*/

	auto d3 = d.Define("scatMC",findScatElecMC,{"MCParticles"})
			   .Define("etaElecMC",getEta,{"scatMC"})
			   .Define("Q2elecMC",getQ2elec,{"scatMC"}).Define("YelecMC",getYelec,{"scatMC"}).Define("XelecMC",getXelec,{"scatMC"})
			   .Define("scatRECbkg",findScatElecRECBkg,{"MCParticles","ReconstructedChargedParticles","ReconstructedChargedParticlesAssociations"})
			   .Define("scatRECbkg3Vect",tlorentzvector_to_tvector3,{"scatRECbkg"})
			   .Define("etaElecRECbkg",getEta,{"scatRECbkg3Vect"})
			   .Define("bkgProb",getPIDprob_pfRICH,{"scatRECbkg"})
			   .Define("bkgMass",getMass,{"scatRECbkg"})
			   .Filter(kineCut,{"Q2elecMC","YelecMC"}); 
			   ;

	auto h_Eta_Elect_REC_bkg = d3.Histo1D({"h_Eta_Elect_REC_bkg", "; #eta; counts", 150, -5, 10}, "etaElecRECbkg");
	auto h_Mass_Elect_REC_bkg = d3.Histo1D({"h_Mass_Elect_REC_bkg", "; mass; counts", 2000, -1, 19}, "bkgMass");
	auto h_etaVsPIDprob_bkg = d3.Histo2D({"h_etaVsPIDprob_bkg", "; #eta; PID probability",100,-5, 5,100,0,1},"etaElecRECbkg","bkgProb");

	//MC
	h_mult_MC->Write();
	h_eta_MC->Write();
	h_pt_MC->Write();
	h_phi_MC->Write();
	//REC
	h_mult_REC->Write();
	h_eta_REC->Write();
	h_pt_REC->Write();
	h_phi_REC->Write();
	//pt resolution
	h_pt_Res->Write();
	h_pt_Res2D->Write();

	//DIS
	h_Eta_Elect_MC->Write();
	h_Q2elec_MC->Write();
	h_Yelec_MC->Write();
	h_Xelec_MC->Write();
	h_Eta_Elect_REC->Write();
	h_Q2elec_REC->Write();
	h_Yelec_REC->Write();
	h_Xelec_REC->Write();
	h_Epz_REC->Write();

	//bkg
	h_Eta_Elect_REC_bkg->Write();
	h_Mass_Elect_REC_bkg->Write();
	h_etaVsPIDprob_bkg->Write();

	output->Write();
  	output->Close();

	return 0;
}
