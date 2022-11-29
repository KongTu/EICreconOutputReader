#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include <Math/Vector4D.h>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "TLorentzRotation.h"
#include "TVector2.h"
#include "TVector3.h"

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4eic/ClusterData.h"
#include "edm4eic/MCRecoParticleAssociationData.h"
#include "edm4hep/MCParticleData.h"

auto findScatElecRECAssociation(const std::vector<edm4hep::MCParticleData>& mcs,
									const std::vector<edm4eic::ReconstructedParticleData>& parts,
										const std::vector<edm4eic::MCRecoParticleAssociationData>& assocs) 
{	
std::vector<TLorentzVector> momenta;
//simplest algo to find REC scat e'
double maxMom=0.;
TVector3 maxtrk(-1E10,-1E10,-1E10);
int elec_index=-1;
int index=-1;
for(auto& i2 : parts){
	index++;
	TVector3 trk(i2.momentum.x,i2.momentum.y,i2.momentum.z);
	if(i2.charge>0) continue;
	if(trk.Mag()>maxMom){
		maxMom=trk.Mag();
		maxtrk=trk;
		elec_index=index;
	}
}

//finding what sim ID this track corresponds to
int mc_elect_index=-1;
// std::cout << "List of association of rec and sim id: "<< std::endl;
for(auto& i3 : assocs){
	int rec_id = i3.recID;
		// std::cout << "rec id = " << rec_id << std::endl;
	int sim_id = i3.simID;
		// std::cout << "sim id = " << sim_id << std::endl;
	if (rec_id == elec_index) mc_elect_index=sim_id;
}

//finding what truth particle PID
index=-1;
int PDG=-99;
double mass=-0.1;
for(auto& i1 : mcs){
	index++;
	if(index == mc_elect_index){
		PDG=i1.PDG;
		mass=i1.mass;
	}
}

if(mass<0) {
	std::cout << "****** A new event ****** " << std::endl;
	std::cout << "reco level electron candidate index = " << elec_index << std::endl;
	std::cout << "No rec association found. " << std::endl;
}

//assigning the true MC mass; 
TLorentzVector maxtrk_4vect;
maxtrk_4vect.SetVectM(maxtrk,mass);

//only to save those electron candidate that doesn't have association.
if(mass<0.){
	momenta.push_back(maxtrk_4vect);
}
  return momenta;
}
auto getEta(const std::vector<TLorentzVector>& tracks)
{
	std::vector<double> eta;
	for(auto& i1 : tracks){eta.push_back(i1.Eta());}
	return eta;
}

int testAssociationExample(TString inname="input/input.root",TString outname="test"){

 	TString rec_file=inname;
	ROOT::RDataFrame d("events", rec_file);

	auto d1 = d.Define("failToAssociate",findScatElecRECAssociation,{"MCParticles","ReconstructedChargedParticles","ReconstructedChargedParticlesAssociations"})
			   .Define("etaOfFailure",getEta,{"failToAssociate"})
			   ;

	auto h_etaOfFailure = d1.Histo1D({"h_etaOfFailure", "; #eta; counts", 150, -5, 10}, "etaOfFailure");

	TString output_name_dir = outname;
  	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

	//write histo
	h_etaOfFailure->Write();

	output->Write();
  	output->Close();

	return 0;
}
