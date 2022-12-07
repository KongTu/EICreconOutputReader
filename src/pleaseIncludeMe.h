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

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

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
//___________________________________________________________________________

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
//___________________________________________________________________________

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
//___________________________________________________________________________

auto getPt(const std::vector<TVector3>& tracks)
{
	std::vector<double> Pt;
	for(auto& i1 : tracks){Pt.push_back(i1.Perp());}
	return Pt;
}
//___________________________________________________________________________

auto getEta(const std::vector<TVector3>& tracks)
{
	std::vector<double> eta;
	for(auto& i1 : tracks){eta.push_back(i1.Eta());}
	return eta;
}
//___________________________________________________________________________

auto getPhi(const std::vector<TVector3>& tracks)
{
	std::vector<double> Phi;
	for(auto& i1 : tracks){Phi.push_back(i1.Phi());}
	return Phi;
}
//___________________________________________________________________________

auto getQ2elec(float electronEnergyInit, const TVector3 elec)
{
	TLorentzVector ein(0,0,-electronEnergyInit, electronEnergyInit);//need to read from file;
	
	TLorentzVector scat;
	scat.SetPtEtaPhiM(elec.Pt(),elec.Eta(),elec.Phi(),MASS_ELECTRON);
	
	double Q2elec = -(ein-scat).Mag2();
	
	return Q2elec;
}
//___________________________________________________________________________

auto getInelParamElectron_1(float protonEnergyInit, float electronEnergyInit, const TVector3 elec)
{
  //reverse fourmomentum of proton - for boost of electron to rest frame of proton
  TLorentzVector pFourMom_reverse(0, 0, -sqrt(electronEnergyInit*electronEnergyInit-MASS_PROTON*MASS_PROTON), electronEnergyInit);
  //pFourMom_reverse.SetXYZM(6.8742,0,-274.9141,sqrt(275*275+MASS_PROTON*MASS_PROTON));
  
  TLorentzVector e_init(0,0,-electronEnergyInit, electronEnergyInit);//need to read from file;
	
	TLorentzVector e_scat;
	e_scat.SetPtEtaPhiM(elec.Pt(),elec.Eta(),elec.Phi(),MASS_ELECTRON);
	
	e_init.Boost(pFourMom_reverse.BoostVector());
	e_scat.Boost(pFourMom_reverse.BoostVector());

  return (e_init.E() - e_scat.E())/e_scat.E();
}

auto getInelParamElectron_2(float protonEnergyInit, float electronEnergyInit, const TVector3 elec)
{

	TLorentzVector ein(0,0,-electronEnergyInit,electronEnergyInit);//need to read from file;
	TLorentzVector pin(0, 0, sqrt(electronEnergyInit*electronEnergyInit-MASS_PROTON*MASS_PROTON), electronEnergyInit);
	//TLorentzVector pin(-6.8742,0,274.9141,sqrt(275*275+MASS_PROTON*MASS_PROTON));
	TLorentzVector scat;
	TLorentzVector q;
	
	scat.SetPtEtaPhiM(elec.Pt(),elec.Eta(),elec.Phi(),MASS_ELECTRON);
	
	q=ein-scat;
	double Q2= -q.Mag2();
	double pq=pin.Dot(q);
	double Yelec = pq/pin.Dot(ein);
		
	return Yelec;
}

