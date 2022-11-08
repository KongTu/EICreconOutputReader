#include "pleaseIncludeMe.h"
int readSingleParticles(TString inname="input/input.root",TString outname="test"){
	
 	TString rec_file=inname;
	//ROOT::EnableImplicitMT(kNumThreads);
	ROOT::RDataFrame d("events", rec_file);

	auto d1 = d.Define("mult",getNtrk,{"ReconstructedChargedParticles"})
						 .Define("momentum",momenta_from_particles,{"ReconstructedChargedParticles"})
						 .Define("eta","momentum.Eta")
						 .Define("pt","momentum.Perp")
						 ;
	auto h_mult_REC = d1.Histo1D({"h_mult_REC", "; N; counts", 10, -0.5, 9.5}, "mult");
	auto h_eta_REC = d1.Histo1D({"h_eta_REC", "; N; counts", 100, -5, 5}, "eta");
	auto h_pt_REC = d1.Histo1D({"h_pt_REC", "; N; counts", 100, 0, 5}, "pt");

	TString output_name_dir = outname;
  TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");
	h_mult_REC->Write();
	h_eta_REC->Write();
	h_pt_REC->Write();
	
	output->Write();
  output->Close();

	return 0;
}
