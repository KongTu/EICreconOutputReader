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
#include "edm4hep/MCParticleData.h"

//for pf-RICH only
#include "/gpfs02/eic/ztu/EPIC/software_tutorial/analysis/irt/delphes/include/DelphesConfig.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

auto ff = new TFile("pfRICH-configs/pfRICH-default-Nov8.root");
auto dconfig = dynamic_cast<DelphesConfig*>(ff->Get("DelphesConfigRICH"));

using namespace ROOT;

auto findScatElecMC(const std::vector<edm4hep::MCParticleData>& parts)
{
  std::vector<TVector3> momenta;
  double maxPt=0.;
  TVector3 leadingTrk(-1E10,-1E10,-1E10);
  //loop over all stable electrons;
  //find the leading pt one
  for(auto& i1 : parts){
  	TVector3 trk(-1E10,-1E10,-1E10);
    if(i1.generatorStatus==1&&i1.PDG==11) {
    	trk.SetXYZ(i1.momentum.x,i1.momentum.y,i1.momentum.z);
    	if(trk.Perp()>maxPt){
    		maxPt=trk.Perp();
    		leadingTrk=trk;
    	}
    }
  }
  momenta.push_back(leadingTrk);
  return momenta;
}
auto findScatElecREC(const std::vector<edm4eic::ClusterData>& clusters,
											const std::vector<edm4eic::ReconstructedParticleData>& parts) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  TLorentzVector escat(-1e10, -1e10, -1e10, -1e10);
  //EEMC
  double maxEnergy=0;
  for(auto& i1 : clusters){
    //need some projection, or matching the cluster here.
    auto energy=i1.energy;
    if(energy>maxEnergy){
    	maxEnergy=energy;
    }
  }
  double maxMom=0.;
  TVector3 maxtrk(-1E10,-1E10,-1E10);
  for(auto& i2 : parts){
  	TVector3 trk(i2.momentum.x,i2.momentum.y,i2.momentum.z);
  	if(trk.Mag()>maxMom){
  		maxMom=trk.Mag();
  		maxtrk=trk;
  	}
  }

  //electron hypothesis;
  double p = sqrt(maxEnergy*maxEnergy- MASS_ELECTRON*MASS_ELECTRON );
  double eta=maxtrk.Eta();
  double phi=maxtrk.Phi();
  double pt = TMath::Sin(maxtrk.Theta())*p;
  escat.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
  
  momenta.push_back(ROOT::Math::PxPyPzMVector{escat.Px(),escat.Py(),escat.Pz(),MASS_ELECTRON});
  return momenta;
}

auto getQ2elec(const std::vector<TVector3>& elec){

	std::vector<double> Q2elec;
	TLorentzVector ein(0,0,-18,18);//need to read from file;
	TLorentzVector scat;
	for(auto& i1 : elec){
		scat.SetPtEtaPhiM(i1.Pt(),i1.Eta(),i1.Phi(),MASS_ELECTRON);
		Q2elec.push_back(- (ein-scat).Mag2() );
	}
	return Q2elec;
}
auto getYelec(const std::vector<TVector3>& elec){

	std::vector<double> Yelec;
	TLorentzVector ein(0,0,-18,18);//need to read from file;
	TLorentzVector pin(-6.8742,0,274.9141,sqrt(275*275+MASS_PROTON*MASS_PROTON));
	TLorentzVector scat;
	TLorentzVector q;
	for(auto& i1 : elec){
		scat.SetPtEtaPhiM(i1.Pt(),i1.Eta(),i1.Phi(),MASS_ELECTRON);
		q=ein-scat;
		double Q2= -q.Mag2();
		double pq=pin.Dot(q);
		double y= pq/pin.Dot(ein);
		Yelec.push_back( y );
	}
	return Yelec;
}
auto getXelec(const std::vector<TVector3>& elec){

	std::vector<double> Xelec;
	TLorentzVector ein(0,0,-18,18);//need to read from file;
	TLorentzVector pin(-6.8742,0,274.9141,sqrt(275*275+MASS_PROTON*MASS_PROTON));
	TLorentzVector scat;
	TLorentzVector q;
	for(auto& i1 : elec){
		scat.SetPtEtaPhiM(i1.Pt(),i1.Eta(),i1.Phi(),MASS_ELECTRON);
		q=ein-scat;
		double Q2= -q.Mag2();
		double pq=pin.Dot(q);
		double x= Q2/(2.*pq);
		Xelec.push_back( x );
	}
	return Xelec;
}
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