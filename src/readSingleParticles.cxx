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
						 .Define("pidProb",getPIDprob,{"momentum",0})
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
	auto h_pidProb_REC = d1.Histo1D({"h_pidProb_REC", "; PID probability; counts", 100, 0, 1}, "pidProb");
	auto h_eta_MC = d1.Histo1D({"h_eta_MC", "; #eta; counts", 100, -5, 5}, "etaMC");
	auto h_pt_MC = d1.Histo1D({"h_pt_MC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "ptMC");
	auto h_phi_MC = d1.Histo1D({"h_phi_MC", "; #phi; counts; counts", 100, -PI, PI}, "phiMC");
	auto h_pt_Res = d1.Histo1D({"h_pt_Res", "; Resolution; counts", 100, -1,1}, "ptRes");
	auto h_pt_Res2D = d1.Histo2D({"h_pt_Res2D", "; p_{T} (GeV/c); Resolution",100,0,5,100,-1,1},"pt","ptRes");

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
	h_pidProb_REC->Write();
	//pt resolution
	h_pt_Res->Write();
	h_pt_Res2D->Write();

	output->Write();
  output->Close();

	return 0;
}
