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
						 .Define("pidProb",getPIDprob,{"momentum"})
						 .Define("momentumMC",momenta_from_mcparticles,{"MCParticles"})
						 .Define("etaMC",getEta,{"momentumMC"})
						 .Define("ptMC",getPt,{"momentumMC"})
						 .Define("phiMC",getPhi,{"momentumMC"})
						 .Define("pidProbMC",getPIDprobMC,{"momentumMC","MCParticles"})
						 .Define("ptRes",pt_resolution,{"MCParticles","ReconstructedChargedParticles"})
						 ;

	TString output_name_dir = outname;
  	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

	auto h_mult_REC = d1.Histo1D({"h_mult_REC", "; N; counts", 10, -0.5, 9.5}, "mult");
	auto h_mult_MC = d1.Histo1D({"h_mult_MC", "; N; counts", 10, -0.5, 9.5}, "multMC");
	auto h_eta_REC = d1.Histo1D({"h_eta_REC", "; #eta; counts", 100, -5, 5}, "eta");
	auto h_pt_REC = d1.Histo1D({"h_pt_REC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "pt");
	auto h_phi_REC = d1.Histo1D({"h_phi_REC", "; #phi; counts", 100, -PI, PI}, "phi");
	auto h_pidProb_REC = d1.Histo1D({"h_pidProb_REC", "; PID probability; counts", 100, 0, 1}, "pidProb");
	auto h_eta_MC = d1.Histo1D({"h_eta_MC", "; #eta; counts", 100, -5, 5}, "etaMC");
	auto h_pt_MC = d1.Histo1D({"h_pt_MC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "ptMC");
	auto h_phi_MC = d1.Histo1D({"h_phi_MC", "; #phi; counts; counts", 100, -PI, PI}, "phiMC");
	auto h_pt_Res = d1.Histo1D({"h_pt_Res", "; Resolution; counts", 100, -1,1}, "ptRes");
	auto h_pt_Res2D = d1.Histo2D({"h_pt_Res2D", "; p_{T} (GeV/c); Resolution",100,0,5,100,-1,1},"pt","ptRes");
	auto h_etaVsPIDprob = d1.Histo2D({"h_etaVsPIDprob", "; #eta; PID probability",100,-5, 5,100,0,1},"eta","pidProb");
	auto h_etaVsPIDprobMC = d1.Histo2D({"h_etaVsPIDprobMC", "; #eta; PID probability",100,-5, 5,100,0,1},"etaMC","pidProbMC");

	//MC dis kinematics
	auto d2 = d.Define("scatMC",findScatElecMC,{"MCParticles"})
			   .Define("etaElecMC",getEta,{"scatMC"})
			   .Define("Q2elecMC",getQ2elec,{"scatMC"}).Define("YelecMC",getYelec,{"scatMC"}).Define("XelecMC",getXelec,{"scatMC"})
			   ;

	auto h_Eta_Elect_MC = d2.Histo1D({"h_Eta_Elect_MC", "; #eta; counts", 100, -5, 5}, "etaElecMC");
	auto h_Q2elec_MC = d2.Histo1D({"h_Q2elec_MC", "; Q^{2}_{e}; counts", 1000, 0,1000}, "Q2elecMC");
	auto h_Yelec_MC = d2.Histo1D({"h_Yelec_MC", "; y_{e}; counts", 1000, 0,1}, "YelecMC");
	auto h_Xelec_MC = d2.Histo1D({"h_Xelec_MC", "; x_{e}; counts", 1000, 0,1}, "XelecMC");

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
	//pid
	h_pidProb_REC->Write();
	h_etaVsPIDprob->Write();
	h_etaVsPIDprobMC->Write();

	//DIS
	h_Eta_Elect_MC->Write();
	h_Q2elec_MC->Write();
	h_Yelec_MC->Write();
	h_Xelec_MC->Write();

	output->Write();
  	output->Close();

	return 0;
}
