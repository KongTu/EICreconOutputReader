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

auto findScatElecMC(const std::vector<edm4hep::MCParticleData>& parts)
{
	/*
	Comment - Things to think about:
		- how can we make sure the MC scat' e? 
		- we need to find a way to have this 
		information at the MC level directl
		from the generator, e.g., PYTHIA. 
	*/
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
	/*
	Comment - Things to think about:
		- We need to have some matching/projection
		code to match between tracks and clusters.
		Now, everything below is just for now. 
	*/
	std::vector<TVector3> momenta;
  TLorentzVector escat(-1E10, -1E10, -1E10, -1E10);
  //EEMC
  double maxEnergy=0;
  for(auto& i1 : clusters){
    auto energy=i1.energy;
    if(energy>maxEnergy){
    	maxEnergy=energy;
    }
  }
  double maxMom=0.;
  TVector3 maxtrk(-1E10,-1E10,-1E10);
  int elec_index=0;
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

  double Epz=0.;
  index=-1;
  for(auto& i2 : parts){
  	index++;
  	if(elec_index==index) continue;
  	TVector3 trk(i2.momentum.x,i2.momentum.y,i2.momentum.z);
  	Epz += sqrt(trk.Mag2()+MASS_PION*MASS_PION) - trk.Pz();
  }
  //3-second calibration.
  maxEnergy+=0.9;
  //electron hypothesis;
  double p = sqrt(maxEnergy*maxEnergy- MASS_ELECTRON*MASS_ELECTRON );
  double eta=maxtrk.Eta();
  double phi=maxtrk.Phi();
  double pt = TMath::Sin(maxtrk.Theta())*p;
  escat.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
  
  Epz += escat.E()-escat.Pz();
  if( escat.Eta()< 0. ) {
  	momenta.push_back(escat.Vect());
  }
  return momenta;
}
auto findScatElecRECBkg(const std::vector<edm4hep::MCParticleData>& mcs,
													const std::vector<edm4eic::ReconstructedParticleData>& parts,
															const std::vector<edm4eic::MCRecoParticleAssociationData>& assocs) 
{	
	std::vector<TLorentzVector> momenta;
	//finding REC scat e'
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

  //finding what truth particle ID
  int mc_elect_index=-1;
  for(auto& i3 : assocs){
  	int rec_id = i3.recID;
  	int sim_id = i3.simID;
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
  if(mass<0){
  	//let's do kinematic match for now.
  	double minR=999.;
  	for(auto& i1 : mcs){
  		TVector3 mcTrk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
  		if(mcTrk.DeltaR(maxtrk)<minR){
  			minR=mcTrk.DeltaR(maxtrk);
  			PDG=i1.PDG;
  			mass=i1.mass;
  		}
  	}
  	//finding the smallest deltaR as the MC particle
  }

  //assigning the true MC mass;
  TLorentzVector maxtrk_4vect;
  maxtrk_4vect.SetVectM(maxtrk,mass);
  
  if(PDG!=11 && maxtrk.Eta()<0){
  	momenta.push_back(maxtrk_4vect);
  }

  return momenta;
}
auto getEpzREC(const std::vector<edm4eic::ClusterData>& clusters,
											const std::vector<edm4eic::ReconstructedParticleData>& parts) 
{
	std::vector<double> EmPz;
  TLorentzVector escat(-1E10, -1E10, -1E10, -1E10);
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
  int elec_index=0;
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

  double Epz=0.;
  index=-1;
  for(auto& i2 : parts){
  	index++;
  	if(elec_index==index) continue;
  	TVector3 trk(i2.momentum.x,i2.momentum.y,i2.momentum.z);
  	Epz += sqrt(trk.Mag2()+MASS_PION*MASS_PION) - trk.Pz();
  }
  //3-second calibration.
  maxEnergy+=0.9;
  //electron hypothesis;
  double p = sqrt(maxEnergy*maxEnergy- MASS_ELECTRON*MASS_ELECTRON );
  double eta=maxtrk.Eta();
  double phi=maxtrk.Phi();
  double pt = TMath::Sin(maxtrk.Theta())*p;
  escat.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
  
  Epz += escat.E()-escat.Pz();
  EmPz.push_back(Epz);
  return EmPz;
}
auto getQ2elec(const std::vector<TVector3>& elec){

	std::vector<double> Q2elec;
	TLorentzVector ein(0,0,-18,18);//need to read from file;
	TLorentzVector scat;
	for(auto& i1 : elec){
		scat.SetPtEtaPhiM(i1.Pt(),i1.Eta(),i1.Phi(),MASS_ELECTRON);
		Q2elec.push_back(-(ein-scat).Mag2());
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
										const std::vector<edm4eic::ReconstructedParticleData>& recos,
											const std::vector<edm4eic::MCRecoParticleAssociationData>& assocs){

	std::vector<double> resolution;
	int index=-1;
	for(auto& i1: recos){
		index++;//start counting;
		if(i1.charge==0) continue;//only charged particles
		//find the matching sim_ID;
		int sim_id=-1
		for(auto& i3 : assocs){
			 if(index==i3.recID){
  		 		sim_id = i3.simID;
  		 }
		}
		//skip tracks now without matching; 
		//should be considered as fake rates;
		if(sim_id==-1) continue;
		
		//reco track 3 momentum
		TVector3 trk(i1.momentum.x,i1.momentum.y,i1.momentum.z);
		
		TVector3 matchMCtrk(-99,-99,-99);
		index=-1;
		for(auto& i2 : mcs){
			index++;
			if(index==sim_id
				&& i2.charge!=0 
					&& i2.generatorStatus==1)
			{
				matchMCtrk.SetXYZ(i2.momentum.x,i2.momentum.y,i2.momentum.z);
			}
		}
		if(matchMCtrk.X()!=-99){
			double res= (matchMCtrk.Perp()-trk.Perp()) / matchMCtrk.Perp();
			resolution.push_back( res );
		}
	}

	return resolution;
}
//if use this function to convert to TVector3 for other basic functions, 
//e.g., getPt,getEta,getPhi;
auto tlorentzvector_to_tvector3(const std::vector<TLorentzVector>& tracks)
{
	std::vector<TVector3> momenta;
	for(auto& i1 : tracks){
		TVector3 tmp = i1.Vect();
		momenta.push_back( tmp );
	}
	return momenta;
}
auto getMass(const std::vector<TLorentzVector>& tracks)
{
	std::vector<double> Mass;
	for(auto& i1 : tracks){
		Mass.push_back( i1.M() );
	}
	return Mass;
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
auto getPIDprob_pfRICH(const std::vector<TLorentzVector>& tracks)
{	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	unsigned hdim = dconfig->GetMassHypothesisCount();
	std::vector<double> prob;
	for(auto& i1 : tracks){
		int hpid=-1;
		if( fabs(i1.M()-MASS_PION)<1e-5) hpid=0;
		if( fabs(i1.M()-MASS_KAON)<1e-5) hpid=1;
		if( fabs(i1.M()-MASS_PROTON)<1e-5) hpid=2;

		double hmtx[hdim*hdim];
		TVector3 track = i1.Vect();
  	int ret = dconfig->GetSmearingMatrix(track, hmtx);
  	if(ret!=0 || hpid==-1){prob.push_back(0.);}
  	else{
  		prob.push_back( hmtx[(hdim+1)*hpid] );
  	}	
	}
	return prob;
}