#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include <Math/Vector4D.h>

#include "ROOT/RDataFrame.hxx"
#include <TH1.h>
#include <TH2.h>
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

//temporary hack to resolve issues
#define PODIO_JSON_OUTPUT

#include "nlohmann/json.hpp"
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"

//for pf-RICH only
//#include "/gpfs02/eic/ztu/EPIC/software_tutorial/analysis/irt/delphes/include/DelphesConfig.h"
#include "/gpfs02/eic/janvanek/eic/tmp/irt/delphes/include/DelphesConfig.h"
//#include "/gpfs02/eic/janvanek/eic/tmp/irt/delphes/include/DelphesConfigRICH.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

//pfRICH info input files
auto ff = new TFile("./pfRICH-configs/pfRICH-default-Nov8.root");
auto dconfig = dynamic_cast<DelphesConfig*>(ff->Get("DelphesConfigRICH"));

auto ff_e_pi = new TFile("./pfRICH-configs/pfRICH_e_pi.root");
auto dconfig_e_pi = dynamic_cast<DelphesConfig*>(ff_e_pi->Get("DelphesConfigRICH"));


//pfRICH TOF files
//"./pfRICH-configs/ETOF.e-pi.t000ps-pf10ps.root");
//"./pfRICH-configs/ETOF.e-pi.t025ps-pf10ps.root");
//"./pfRICH-configs/ETOF.e-pi.t025ps-pf20ps.root");

auto ff_e_pi_TOF = new TFile("./pfRICH-configs/ETOF.e-pi.t025ps-pf20ps.root");
auto dconfig_e_pi_TOF = dynamic_cast<DelphesConfig*>(ff_e_pi_TOF->Get("DelphesConfigTOF"));
//_________________________________________________________________________


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
//___________________________________________________________________________

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
//___________________________________________________________________________

//pfRICH test
auto getPIDprob_pfRICH(const std::vector<TLorentzVector>& tracks)
{	
  //hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	unsigned hdim = dconfig->GetMassHypothesisCount();
	std::vector<double> prob;
	
	for(auto& i1 : tracks)
	{
		int hpid=-1;
		if( fabs(i1.M()-MASS_PION)<1e-5) hpid=0;
		if( fabs(i1.M()-MASS_KAON)<1e-5) hpid=1;
		if( fabs(i1.M()-MASS_PROTON)<1e-5) hpid=2;

		double hmtx[hdim*hdim];
		TVector3 track = i1.Vect();
  	int ret = dconfig->GetSmearingMatrix(track, hmtx);
  	
  	if(ret!=0 || hpid==-1)
  	{
  	  prob.push_back(0.);
  	}
  	else
  	{
  		prob.push_back( hmtx[(hdim+1)*hpid] );
  	}	
	}
	
	return prob;
}
//___________________________________________________________________________

auto getPIDprob_pfRICH_single(TLorentzVector track)
{
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	unsigned hdim = dconfig->GetMassHypothesisCount();
	double prob;
	
	int hpid=-1;
	if( fabs(track.M()-MASS_PION)<1e-5) hpid=0;
	if( fabs(track.M()-MASS_KAON)<1e-5) hpid=1;
	if( fabs(track.M()-MASS_PROTON)<1e-5) hpid=2;

	double hmtx[hdim*hdim];

	int ret = dconfig->GetSmearingMatrix(track.Vect(), hmtx);
	
	if(ret!=0 || hpid==-1)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  prob = 0.;
	}
	else
	{
		prob = hmtx[(hdim+1)*hpid];
	}	
	
	return prob;
}
//___________________________________________________________________________

auto getPIDprob_pfRICH_MC(TLorentzVector track, int hpid, int config)
{
  //config == 0 - pi/K/p PID table
  //config == 1 - e/pi PID table
  
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	unsigned hdim;
	
	if(config == 0) 
	{
  	hdim = dconfig->GetMassHypothesisCount();
	}	
	
	if(config == 1) 
	{
	  hdim = dconfig_e_pi->GetMassHypothesisCount();
	}
	
	
	double hmtx[hdim*hdim];	
	int ret;
	
	if(config == 0) 
	{
   	ret = dconfig->GetSmearingMatrix(track.Vect(), hmtx);	
	}	
	
	if(config == 1) 
	{
	  ret = dconfig_e_pi->GetSmearingMatrix(track.Vect(), hmtx);	  
	}
	
	
	double prob;
	
	if(ret!=0 || hpid < 0 || hpid > 2)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  prob = 0.;
	}
	else
	{
		prob = hmtx[(hdim+1)*hpid];
	}	
	
	return prob;
}
//___________________________________________________________________________

auto getPIDprob_pfRICH_MC(TVector3 track, int hpid, int config)
{
  //config == 0 - pi/K/p PID table
  //config == 1 - e/pi PID table

	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	unsigned hdim;
	
	if(config == 0) 
	{
  	hdim = dconfig->GetMassHypothesisCount();
	}	
	
	if(config == 1) 
	{
	  hdim = dconfig_e_pi->GetMassHypothesisCount();
	}
	
	
	double hmtx[hdim*hdim];	
	int ret;
	
	if(config == 0) 
	{
   	ret = dconfig->GetSmearingMatrix(track, hmtx);	
	}	
	
	if(config == 1) 
	{
	  ret = dconfig_e_pi->GetSmearingMatrix(track, hmtx);	  
	}
	
	
	double prob;
	
	
	if(ret!=0 || hpid < 0 || hpid > 2)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  prob = 0.;
	}
	else
	{
		prob = hmtx[(hdim+1)*hpid];
	}	
	
	return prob;
}
//___________________________________________________________________________

int getPIDprob_pfRICH_mtx(TLorentzVector track, double (&PID_mtx)[9], int config)
{
  //config == 0 - pi/K/p PID table
  //config == 1 - e/pi PID table
  
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	unsigned hdim;
	
	if(config == 0) 
	{
  	hdim = dconfig->GetMassHypothesisCount();
	}	
	
	if(config == 1) 
	{
	  hdim = dconfig_e_pi->GetMassHypothesisCount();
	}
	
	
	double hmtx[hdim*hdim];	
	int ret;
	
	if(config == 0) 
	{
   	ret = dconfig->GetSmearingMatrix(track.Vect(), hmtx);	
	}	
	
	if(config == 1) 
	{
	  ret = dconfig_e_pi->GetSmearingMatrix(track.Vect(), hmtx);	  
	}
	
	
	if(ret!=0)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  return 0;
	}
	else
	{
  	for(unsigned int index = 0; index < hdim*hdim; index++)
    {
      PID_mtx[index] = hmtx[index]	;
    }
	
		return 1;
	}	
	
}
//___________________________________________________________________________

int getPIDprob_pfRICH_mtx(TVector3 track, double (&PID_mtx)[9], int config)
{
  //config == 0 - pi/K/p PID table
  //config == 1 - e/pi PID table
  
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	unsigned hdim;
	
	if(config == 0) 
	{
  	hdim = dconfig->GetMassHypothesisCount();
	}	
	
	if(config == 1) 
	{
	  hdim = dconfig_e_pi->GetMassHypothesisCount();
	}
	
	
	double hmtx[hdim*hdim];	
	int ret;
	
	if(config == 0) 
	{
   	ret = dconfig->GetSmearingMatrix(track, hmtx);	
	}	
	
	if(config == 1) 
	{
	  ret = dconfig_e_pi->GetSmearingMatrix(track, hmtx);	  
	}
	
	
	
	if(ret!=0)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  return 0;
	}
	else
	{
	  for(unsigned int index = 0; index < hdim*hdim; index++)
	  {
	    PID_mtx[index] = hmtx[index]	;
	  }
	
		return 1;
	}	
	
}
//___________________________________________________________________________

auto getPIDprob_TOF_MC(TVector3 track, int hpid, int config)
{
  //config == 0 - pi/K/p PID table
  
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	
	//config == 1 - e/pi PID table
	
	//hpid==0,electron
	//hpid==1,pion

	
	unsigned hdim = dconfig_e_pi_TOF->GetMassHypothesisCount();
	
	double hmtx[hdim*hdim];	
	int ret = dconfig_e_pi_TOF->GetSmearingMatrix(track, hmtx);;	
	
	double prob;
	
	if(ret!=0 || hpid < 0 || hpid > 2)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  prob = 0.;
	}
	else
	{
		prob = hmtx[(hdim+1)*hpid];
	}	
	
	return prob;
}
//___________________________________________________________________________

int getPIDprob_TOF_mtx(TLorentzVector track, double (&PID_mtx)[9], int config)
{
  //config == 0 - pi/K/p PID table
  
	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	
	//config == 1 - e/pi PID table
	
	//hpid==0,electron
	//hpid==1,pion
	
	unsigned hdim = dconfig_e_pi_TOF->GetMassHypothesisCount();
	
	double hmtx[hdim*hdim];	
	int ret = dconfig_e_pi_TOF->GetSmearingMatrix(track.Vect(), hmtx);
	
	
	if(ret!=0)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  return 0;
	}
	else
	{
  	for(unsigned int index = 0; index < hdim*hdim; index++)
    {
      PID_mtx[index] = hmtx[index]	;
    }
	
		return 1;
	}	
	
}
//___________________________________________________________________________

int getPIDprob_TOF_mtx(TVector3 track, double (&PID_mtx)[9], int config)
{
  //config == 0 - pi/K/p PID table

	//hpid==0,pion
	//hpid==1,kaon
	//hpod==2,proton
	
	
	//config == 1 - e/pi PID table
	
	//hpid==0,electron
	//hpid==1,pion
	
  unsigned hdim = dconfig_e_pi_TOF->GetMassHypothesisCount();
	
	double hmtx[hdim*hdim];	
	int ret = dconfig_e_pi_TOF->GetSmearingMatrix(track, hmtx);
	
	
	if(ret!=0)
	{
	  //cout<<"ret = "<<ret<<", hpid = "<<hpid<<endl;
	  return 0;
	}
	else
	{
	  for(unsigned int index = 0; index < hdim*hdim; index++)
	  {
	    PID_mtx[index] = hmtx[index]	;
	  }
	
		return 1;
	}	
	
}

