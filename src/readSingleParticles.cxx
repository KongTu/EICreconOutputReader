#include "pleaseIncludeMe.h"
int readSingleParticles(TString inname="input/input.root",TString outname="test"){
	
 	TString rec_file=inname;
	//ROOT::EnableImplicitMT(kNumThreads);
	ROOT::RDataFrame d("events", rec_file);

	auto d1 = d.Define("mult",getNtrk,{"ReconstructedChargedParticles"})
						 .Define("momentum",momenta_from_chargedparticles,{"ReconstructedChargedParticles"})
						 .Define("eta",getEta,{"momentum"})
						 .Define("pt",getPt,{"momentum"})
						 .Define("momentumMC",momenta_from_mcparticles,{"MCParticles"})
						 .Define("etaMC",getEta,{"momentumMC"})
						 .Define("ptMC",getPt,{"momentumMC"})
						 .Define("pRes",momentum_resolution,{"MCParticles","ReconstructedChargedParticles"})
						 ;

	TString output_name_dir = outname;
  TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

	auto h_mult_REC = d1.Histo1D({"h_mult_REC", "; N; counts", 10, -0.5, 9.5}, "mult");
	auto h_eta_REC = d1.Histo1D({"h_eta_REC", "; #eta; counts", 100, -5, 5}, "eta");
	auto h_pt_REC = d1.Histo1D({"h_pt_REC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "pt");
	auto h_eta_MC = d1.Histo1D({"h_eta_MC", "; #eta; counts", 100, -5, 5}, "etaMC");
	auto h_pt_MC = d1.Histo1D({"h_pt_MC", "; p_{T} (GeV/c); counts", 100, 0, 5}, "ptMC");
	auto h_p_Res = d1.Histo1D({"h_p_Res", "; Resolution; counts", 100, -1,1}, "pRes");
	
	h_mult_REC->Write();
	h_eta_REC->Write();
	h_pt_REC->Write();
	h_eta_MC->Write();
	h_pt_MC->Write();
	h_p_Res->Write();

	output->Write();
  output->Close();

	return 0;
}
