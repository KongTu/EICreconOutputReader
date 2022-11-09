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
#include "TVector3.h"

#include "fmt/color.h"
#include "fmt/core.h"

// #include "nlohmann/json.hpp"
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
// #include "edm4eic/ReconstructedParticleCollection.h"

//for pf-RICH only
#include "/tmp/EPIC-Kong/irt/delphes/include/DelphesConfig.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

auto ff = new TFile("pfRICH-configs/pfRICH-default-Nov8.root");
auto dconfig = dynamic_cast<DelphesConfig*>(ff->Get("DelphesConfigRICH"));

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
auto getNtrkMC(const std::vector<edm4hep::MCParticleData>& parts)
{
  std::vector<int> mult;
  int n=0;
  for(auto& i1 : parts){
    if(i1.charge!=0 && i1.generatorStatus==1) n++;
  }
  mult.push_back( n );
  return mult;
}
auto momenta_from_chargedparticles(const std::vector<edm4eic::ReconstructedParticleData>& parts) {
  std::vector<TVector3> momenta;
  for(auto& i1 : parts){
		TVector3 trk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
		if(i1.charge!=0) momenta.push_back(trk);
	}
  return momenta;
}
auto momenta_from_mcparticles(const std::vector<edm4hep::MCParticleData>& parts) {
  std::vector<TVector3> momenta;
  for(auto& i1 : parts){
		TVector3 trk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
		if(i1.charge!=0 && i1.generatorStatus==1) momenta.push_back(trk);
	}
  return momenta;
}
auto pt_resolution(const std::vector<edm4hep::MCParticleData>& mcs,
							const std::vector<edm4eic::ReconstructedParticleData>& recos){

	std::vector<double> resolution;
	for(auto& i1: recos){
		TVector3 trk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
		if(i1.charge==0) continue;
		double minR=99;
		TVector3 matchMCtrk(-99,-99,-99);
		for(auto& i2 : mcs){
			TVector3 trkMC(i2.momentum.x,i2.momentum.y,i2.momentum.z);
			if(i2.charge!=0 && i2.generatorStatus==1){
				if(trk.DeltaR(trkMC) < minR ){
					minR=trk.DeltaR(trkMC);
					matchMCtrk=trkMC;
				}
			}
		}
		double res= (matchMCtrk.Perp()-trk.Perp()) / matchMCtrk.Perp();
		resolution.push_back( res );
				
	}

	return resolution;
}
auto getPt(const std::vector<TVector3>& tracks)
{
	std::vector<double> Pt;
	for(auto& i1 : tracks){Pt.push_back(i1.Perp());}
	return Pt;
}
auto getEta(const std::vector<TVector3>& tracks)
{
	std::vector<double> eta;
	for(auto& i1 : tracks){eta.push_back(i1.Eta());}
	return eta;
}
auto getPhi(const std::vector<TVector3>& tracks)
{
	std::vector<double> Phi;
	for(auto& i1 : tracks){Phi.push_back(i1.Phi());}
	return Phi;
}
//pfRICH test
auto getPIDprob(const std::vector<TVector3>& tracks)
{	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	unsigned hdim = dconfig->GetMassHypothesisCount();
	std::vector<double> prob;
	for(auto& i1 : tracks){

		int hpid=0;
		double hmtx[hdim*hdim];
    	int ret = dconfig->GetSmearingMatrix(i1, hmtx);
    	if(ret!=0){prob.push_back(0.);}
    	else{
    		prob.push_back( hmtx[(hdim+1)*hpid] );
    	}	
	}
	return prob;
}
auto getPIDprobMC(const std::vector<TVector3>& tracks,
									const std::vector<edm4hep::MCParticleData>& mcs )
{	
	/*
	Need to write the association here. 
	Basically, one should match rec to mc, and use mc pdg 
	to determine which true PID is and set hpid; 
	-->	hpid==0,pion; hpid==1,kaon; hpod==2,proton
	*/

	unsigned hdim = dconfig->GetMassHypothesisCount();
	std::vector<double> prob;
	for(auto& i1 : tracks){
	//matching mc for now
		double minR=1e7;
		TVector3 matchMCtrk(0,0,0);
		int matchPID=-99;
		for(auto& i2 : mcs){
			TVector3 trkMC(i2.momentum.x,i2.momentum.y,i2.momentum.z);
			if(i2.charge!=0 && i2.generatorStatus==1){
				if(i1.DeltaR(trkMC) < minR ){
					minR=i1.DeltaR(trkMC);
					matchMCtrk=trkMC;
					matchPID=i2.PDG;
				}
			}
		}
		int hpid=-99;
		if(TMath::Abs(matchPID)==211) hpid=0;
		else if(TMath::Abs(matchPID)==321) hpid=1;
		else if(TMath::Abs(matchPID)==2212) hpid=2;
		else hpid=-99;
		
		double hmtx[hdim*hdim];
  	int ret = dconfig->GetSmearingMatrix(i1, hmtx);
  	if(ret!=0 || hpid==-99){
  		prob.push_back(0.);
  	}
  	else{
  		prob.push_back( hmtx[(hdim+1)*hpid] );
  	}	
	}
	return prob;
}