#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TVector2.h"

//#include "nlohmann/json.hpp"
//#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"

auto getNtrk(const std::vector<edm4eic::ReconstructedParticleData>& parts)
{
  std::vector<int> mult;
  int n=0;
  for(auto& i1 : parts){
    if(i1.charge!=0) n++;
  }
  mult.push_back( n );
  return mult;
}

int readSingleParticles(TString inname="input/pi-_2GeV_130to177deg.0001.eicrecon.tree.edm4eic.root",TString outname="test"){
	
	
 	TString rec_file=inname;
	//ROOT::EnableImplicitMT(kNumThreads);
	ROOT::RDataFrame d("events", rec_file);

	auto d1 = d.Define("mult",getNtrk,{"ReconstructedChargedParticles"});
	auto h_mult_REC = d1.Histo1D({"h_mult_REC", "; N; counts", 10, -0.5, 9.5}, "mult");

	TString output_name_dir = outname;
    TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");
	h_mult_REC->Write();
	output->Write();
    output->Close();

	return 0;
}
