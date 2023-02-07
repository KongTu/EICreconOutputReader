#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include "pleaseIncludeMe.h"

using namespace std;

//int electronPionSeparation(TString inname="input/input.root",TString outname="test")
int kaonPurity(TString inname="./fileLists/flieList.list", TString outname="test", float e_energy = 18, float p_energy = 275)
{
  if( !(e_energy == 18 && p_energy == 275) && !(e_energy == 10 && p_energy == 100) && !(e_energy == 5 && p_energy == 41))
  {
    cout<<"Invalid beam energies."<<endl;

    return 0;
  }

  cout<<"Start"<<endl;
/*
  //old binnig for first testing
  const int nQ2bins = 2;
  float const Q2_bins[nQ2bins+1] = { 1., 2., 4. };

  const int nyInelParBins = 2;
  float const y_bins[nQ2bins+1] = { 0., 0.001, 0.002 };
*/

  const int nQ2bins = 4;
  float const Q2_bins[nQ2bins+1] = { 1,3,5,10,20 };

  const int nyInelParBins = 4;
  float const y_bins[nyInelParBins+1] = { 0.01,0.05,0.1,0.5,0.95 };

  const int nMomBins = 11;
  float const mom_bins[nMomBins+1] = { 0,0.5,1,1.5,2,3,4,5,6,7,10, 18 };
  
  TRandom3 *myRandom = new TRandom3();
  
  
  //pfRICH eta acceptance mc_mom.Eta() > -3.8 && mc_mom.Eta() < -1.5
  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -3.8, -3, -2.5, -2, -1.5};

  //____________________________________________________
  //pi eCAL info
  //values from Dimitry Kalinkin

  double array_mom_bins[8] = {0.1, 0.2, 0.5, 1., 2., 5., 10., 20.};
  double array_pi_false_rate_85[8] = {0.6332832672698294, 0.7495818158985306, 0.00930384575910461, 0.001692827694491846, 0.0001898238241173789, 0.00020018016214593134, 0.000536412900269677, 0.0006430092230696459};
  double array_pi_false_rate_95[8] = {0.8838860654090215, 0.9228366502089709, 0.02228665375203912, 0.0036237053948322794, 0.00048096291971113834, 0.0006523112180272588, 0.0022770026949322946, 0.0018829746368912276};


  TGraph *g_pi_false_rate_85 = new TGraph(8, array_mom_bins, array_pi_false_rate_85);

  TGraph *g_pi_false_rate_95 = new TGraph(8, array_mom_bins, array_pi_false_rate_95);
  //___________________________________________________


  //load files to TChain
  ifstream fileList;
  fileList.open(inname.Data());

  string fileFromList;

  auto myChain = new TChain("events");

  while(getline(fileList, fileFromList))
  {
    myChain->Add(fileFromList.c_str());
  }


  TTreeReader tree_reader(myChain);       // !the tree reader

  // MC particles
  TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
  TTreeReaderArray<int> mc_generatorStatus_array = {tree_reader, "MCParticles.generatorStatus"};


  // Reconstructed particles
  TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> reco_cahrge = {tree_reader, "ReconstructedChargedParticles.charge"};
  TTreeReaderArray<int> reco_PDG = {tree_reader, "ReconstructedChargedParticles.PDG"};

  TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
  TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

  

  //defining output file and histos.
  TString output_name_dir = outname;
	TFile* output = new TFile("output/"+output_name_dir+"-output_K.root","RECREATE");

  //MC histograms
  TH1F *h_eta_K_plus_MC[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC[nEtaBins+1];
  TH1F *h_p_K_minus_MC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_MC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_MC_pfRICH[nEtaBins+1];
  
  
/*  
  TH1F *h_eta_K_plus_MC_RC[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_RC[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_RC[nEtaBins+1];
  TH1F *h_p_K_minus_MC_RC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_MC_RC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_RC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_RC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_MC_RC_pfRICH[nEtaBins+1];
  
  
  
  TH1F *h_eta_K_plus_RC[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_RC[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_RC[nEtaBins+1];
  TH1F *h_p_K_minus_RC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_RC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_RC_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_RC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_RC_pfRICH[nEtaBins+1]; 
*/
  //____________________________________________________________


  TH1F *h_K_minus_purity_pfRICH_MC[nMomBins][nQ2bins][nyInelParBins];
  //TH1F *h_K_minus_purity_pfRICH_MC_RC[nMomBins][nQ2bins][nyInelParBins];
  //TH1F *h_K_minus_purity_pfRICH_RC[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_K_minus_purity_p_eta_pfRICH_MC[nMomBins][nEtaBins+1];
  //TH1F *h_K_minus_purity_p_eta_pfRICH_MC_RC[nMomBins][nEtaBins+1];
  //TH1F *h_K_minus_purity_p_eta_pfRICH_RC[nMomBins][nEtaBins+1];
  
  
  TH1F *h_K_plus_purity_pfRICH_MC[nMomBins][nQ2bins][nyInelParBins];
  //TH1F *h_K_plus_purity_pfRICH_MC_RC[nMomBins][nQ2bins][nyInelParBins];
  //TH1F *h_K_plus_purity_pfRICH_RC[nMomBins][nQ2bins][nyInelParBins];
  
  TH1F *h_K_plus_purity_p_eta_pfRICH_MC[nMomBins][nEtaBins+1];
  //TH1F *h_K_plus_purity_p_eta_pfRICH_MC_RC[nMomBins][nEtaBins+1];
  //TH1F *h_K_plus_purity_p_eta_pfRICH_RC[nMomBins][nEtaBins+1];


  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {
        //MC histograms
        h_eta_K_plus_MC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);       
        h_eta_K_minus_MC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);

        h_eta_K_plus_MC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_MC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_K_minus_MC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_MC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
/*       
        h_eta_K_plus_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);       
        h_eta_K_minus_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);

        h_eta_K_plus_RC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_K_minus_RC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        
        h_eta_K_plus_MC_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);       
        h_eta_K_minus_MC_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);

        h_eta_K_plus_MC_RC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_MC_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_K_minus_MC_RC_pfRICH[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_MC_RC_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
*/      

        //pfRICH
        h_K_minus_purity_pfRICH_MC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        //h_K_minus_purity_pfRICH_MC_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        //h_K_minus_purity_pfRICH_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        
        h_K_plus_purity_pfRICH_MC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_MC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        //h_K_plus_purity_pfRICH_MC_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_MC_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        //h_K_plus_purity_pfRICH_RC[mom_bin][Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);

      }

    }
  }
  
  //momentum histograms in eta bins
  for(unsigned int eta_bin = 0; eta_bin < nEtaBins+1; eta_bin++) 
  {
    //MC histograms
    h_p_K_plus_MC[eta_bin] = new TH1F(Form("h_p_K_plus_MC_eta_%i", eta_bin), Form("h_p_K_plus_MC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC[eta_bin] = new TH1F(Form("h_p_K_minus_MC_eta_%i", eta_bin), Form("h_p_K_minus_MC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_MC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_MC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_MC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_MC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_MC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
/*    
    h_p_K_plus_RC[eta_bin] = new TH1F(Form("h_p_K_plus_RC_eta_%i", eta_bin), Form("h_p_K_plus_RC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_RC[eta_bin] = new TH1F(Form("h_p_K_minus_RC_eta_%i", eta_bin), Form("h_p_K_minus_RC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    
    h_p_K_plus_MC_RC[eta_bin] = new TH1F(Form("h_p_K_plus_MC_RC_eta_%i", eta_bin), Form("h_p_K_plus_MC_RC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_RC[eta_bin] = new TH1F(Form("h_p_K_minus_MC_RC_eta_%i", eta_bin), Form("h_p_K_minus_MC_RC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_MC_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_MC_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_MC_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_MC_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_MC_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
*/  
        
    
    for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
    {
      h_K_minus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_MC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_MC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
      //h_K_minus_purity_p_eta_pfRICH_RC[mom_bin][eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_RC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_RC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
      //h_K_minus_purity_p_eta_pfRICH_MC_RC[mom_bin][eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_MC_RC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_MC_RC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
      
      h_K_plus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_MC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_MC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
      //h_K_plus_purity_p_eta_pfRICH_RC[mom_bin][eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_RC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_RC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
      //h_K_plus_purity_p_eta_pfRICH_MC_RC[mom_bin][eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_MC_RC_mom_%i_eta_%i" , mom_bin, eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_MC_RC_mom_%i_eta_%i" , mom_bin, eta_bin), 2, 0, 2);
    }
  
  }
  //______________________________________________________________________________________________________________________________________________________________________________________
  

	tree_reader.SetEntriesRange(0, myChain->GetEntries());

  while (tree_reader.Next())
  {
  	//MC analysis

    //finding the scattered electron
  	TLorentzVector scatMC(0,0,0,0);
  	int mc_elect_index = -1;
  	double maxP = -99.;
  	float scat_e_eta_MC = -99.;


  	for(int imc = 0; imc < mc_px_array.GetSize(); imc++)
  	{
  	
  	  if( mc_generatorStatus_array[imc] != 1 ) continue;

  		TVector3 mctrk(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

  		if(mc_pdg_array[imc] == 11 && mctrk.Mag() > maxP)
  		//if(mc_pdg_array[imc] == 11 && mctrk.Perp() > maxP)
  		{
  			maxP = mctrk.Mag();
  			scat_e_eta_MC = mctrk.Eta();
  			//maxP = mctrk.Perp();
  			mc_elect_index = imc;
  			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
  		}
   
  	} //end MC particles for loop
  	
  	//____________________________________________________________________________________________
    

  	double y_inelPar_e = getInelParamElectron_2(p_energy, e_energy, scatMC.Vect() );

    double Q2_electron = getQ2elec( e_energy, scatMC.Vect());


    //flag for wide MC Q2 and y cut
    int has_good_Q2_y_MC = 0;
    
    //if( Q2_electron > 1. && Q2_electron < 20. && y_inelPar_e > 0.01 && y_inelPar_e < 0.95) has_good_Q2_y_MC = 1;
    if( Q2_electron > 1. && y_inelPar_e < 0.95) has_good_Q2_y_MC = 1;

    
    //find bins
    int Q2_bin = -1;

    for(int j = 0; j < nQ2bins; j++)
    {
      if(Q2_electron > Q2_bins[j] && Q2_electron <= Q2_bins[j+1])
      {
        Q2_bin = j;
      }
     }


    int y_bin = -1;

    for(int j = 0; j < nyInelParBins; j++)
    {
      if(y_inelPar_e > y_bins[j] && y_inelPar_e <= y_bins[j+1])
      {
        y_bin = j;
      }
    }


    int mom_bin_scat_e = -1;

    for(int j = 0; j < nMomBins; j++) //loop over pT bins
    {
      if(scatMC.P() > mom_bins[j] && scatMC.P() <= mom_bins[j+1])
      {
        mom_bin_scat_e = j;
      }
    }
    
    
    int eta_bin_scat_e = -1;

    for(int j = 0; j < nEtaBins; j++) //loop over eta bins
    {
      if(scatMC.Eta() > eta_bins[j] && scatMC.Eta() <= eta_bins[j+1])
      {
        eta_bin_scat_e = j;
      }
    }

    //__________________________________________________________________________________

   
    //fill MC eta distributions
    
    //vectors for MC -> RC matchig
    //vector<int> K_minus_simID_vect;
    //vector<float> K_minus_eta_MC; //MC eta of pion
    
    //vector<int> K_plus_simID_vect;
    //vector<float> K_plus_eta_MC; //MC eta of pion

    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{
  	  if(has_good_Q2_y_MC != 1) continue;

  	  if( mc_generatorStatus_array[imc] != 1 ) continue;

  		TVector3 mc_mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
  		
      //only par in pfRICH acceptance
  		if(mc_mom.Eta() < -3.8 && mc_mom.Eta() > -1.5) continue;
	
  		//determine momentum and eta bin of particles in event
  		int mom_bin = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(mc_mom.Mag() > mom_bins[j] && mc_mom.Mag() <= mom_bins[j+1])
        {
          mom_bin = j;
        }

      }
      
      
      int eta_bin = -1;

      for(int j = 0; j < nEtaBins; j++) //loop over pT bins
      {
        if(mc_mom.Eta() > eta_bins[j] && mc_mom.Eta() <= eta_bins[j+1])
        {
          eta_bin = j;
        }

      }
      
      //_________________________________________________________________________
            
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(mc_mom, pfRICH_mtx, 0);      
      
      //generate random number for PID in pfRICH
      double rndm_PID_MC_loop = myRandom->Rndm();
      
      int is_pfRICH_kaon_MC = 0;
 
      if(good_pfRICH_mtx == 1)
      {
        //find matchig MC track
        //pi
        if( fabs(mc_pdg_array[imc]) == 211 )
        {          
          //pi mis-identified as K          
          if( pfRICH_mtx[1] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;
        }

        //K
        if( fabs(mc_pdg_array[imc]) == 321 )
        {          
          //K identified as K          
          if( pfRICH_mtx[4] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;
        }

        //p
        if( fabs(mc_pdg_array[imc]) == 2212 )
        {          
          //K identified as K          
          if( pfRICH_mtx[7] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;          
        }
      }
      //_________________________________________________________________________      

      //fill MC p distributions      
      //cout<<"p start"<<endl;
      
      if( mom_bin != -1 && eta_bin != -1)
      {
        //pure MC
        //K-
        if(mc_pdg_array[imc] == -321) 
        {
          h_p_K_minus_MC[eta_bin]->Fill(mc_mom.Mag());                 
          h_p_K_minus_MC[nEtaBins]->Fill(mc_mom.Mag()); //baseline in pfRICH eta acceptance
          
          //K_minus_simID_vect.push_back(imc);
          //K_minus_eta_MC.push_back(mc_mom.Eta());
                 
        }
        
        //K+
        if(mc_pdg_array[imc] == 321) 
        {
          h_p_K_plus_MC[eta_bin]->Fill(mc_mom.Mag());
          h_p_K_plus_MC[nEtaBins]->Fill(mc_mom.Mag()); //baseline in pfRICH eta acceptance
          
          //K_plus_simID_vect.push_back(imc);
          //K_plus_eta_MC.push_back(mc_mom.Eta());
                  
        }
        
        //pfRICH with MC                     
                            
    
    
        //______________________________________________________________________________________________________________________________________________________________________________________
        
    
        //K purity in p and eta bins
        //Q2 and y bin given by MC scattered e
      
        if( is_pfRICH_kaon_MC == 1 )
        {
          //K- identified using pfRICH (with mis-PID)
          if( mc_pdg_array[imc] < 0 )
          {
          
            h_p_K_minus_MC_pfRICH[eta_bin]->Fill(mc_mom.Mag());
            h_p_K_minus_MC_pfRICH[nEtaBins]->Fill(mc_mom.Mag());
            
            if( mc_pdg_array[imc] == -321 ) 
            {       
              h_K_minus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin]->Fill(1.5);
              h_K_minus_purity_p_eta_pfRICH_MC[mom_bin][nEtaBins]->Fill(1.5);           
            }      
            else //pi- and p-bar mis-identified as K-
            {
              h_K_minus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin]->Fill(0.5);
              h_K_minus_purity_p_eta_pfRICH_MC[mom_bin][nEtaBins]->Fill(0.5);
            }
          
          
          }
          
          //K+ identified using pfRICH (with mis-PID)
          if( mc_pdg_array[imc] > 0 )
          {
          
            h_p_K_plus_MC_pfRICH[eta_bin]->Fill(mc_mom.Mag());
            h_p_K_plus_MC_pfRICH[nEtaBins]->Fill(mc_mom.Mag());          
          
            if( mc_pdg_array[imc] == 321 ) 
            {       
              h_K_plus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin]->Fill(1.5);
              h_K_plus_purity_p_eta_pfRICH_MC[mom_bin][nEtaBins]->Fill(1.5);           
            }      
            else if( mc_pdg_array[imc] > 0 )//pi+ and p mis-identified as K+
            {
              h_K_plus_purity_p_eta_pfRICH_MC[mom_bin][eta_bin]->Fill(0.5);
              h_K_plus_purity_p_eta_pfRICH_MC[mom_bin][nEtaBins]->Fill(0.5);
            }
            
          }
          
        } //end if PID and acceptance
      
      } //end if bins
      //cout<<"p end"<<endl;
      
      //___________________________________________________________________________________________________________________________________________________________
      
      
      //fill eta distributions for MC particles
      if(Q2_bin < 0 || y_bin < 0 || mom_bin < 0) continue;
      
  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 &&  mc_mom.Pt() < maxP)
  		

  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		  h_eta_K_plus_MC[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());	   
  		}
  		//K-
  		if(mc_pdg_array[imc] == -321)
  		{
  		  h_eta_K_minus_MC[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
  		}
  		//_____________________________________________________________________________
  		
  		
  		if( is_pfRICH_kaon_MC == 1 )
      {
        //K- identified using pfRICH (with mis-PID)
        if( mc_pdg_array[imc] < 0 )
        {
          h_eta_K_minus_MC_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
          
          if( mc_pdg_array[imc] == -321 ) 
          {       
            h_K_minus_purity_pfRICH_MC[mom_bin][Q2_bin][y_bin]->Fill(1.5);
                      
          }      
          else //pi- and p-bar mis-identified as K-
          {
            h_K_minus_purity_pfRICH_MC[mom_bin][Q2_bin][y_bin]->Fill(0.5);
          }
        
        }
        
        //K+ identified using pfRICH (with mis-PID)
        if( mc_pdg_array[imc] > 0 )
        {
          h_eta_K_plus_MC_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
          
          if( mc_pdg_array[imc] == 321 ) 
          {       
            h_K_plus_purity_pfRICH_MC[mom_bin][Q2_bin][y_bin]->Fill(1.5);
                      
          }      
          else //pi+ and p mis-identified as K+
          {
            h_K_plus_purity_pfRICH_MC[mom_bin][Q2_bin][y_bin]->Fill(0.5);
          }
        
        }
        
      
      }
      

  		


  	}//end second MC particle loop
  	//___________________________________________________________________________________________________________________________________________________________

  	
  	//reconstructed charged particles analysis
  	//need to update
  	
/*
    //cout<<"RC - MC match start"<<endl;
    //MC -> RC matching
  	//find corresponding recID for scattered e and pi-
  	//int scat_e_recID = -1;
  	vector<int> K_minus_recID_vect;
  	
  	for(unsigned int ch_trk_assoc_i = 0; ch_trk_assoc_i < sim_id.GetSize(); ch_trk_assoc_i++)
  	{
  	  //if( mc_elect_index == sim_id[ch_trk_assoc_i] ) scat_e_recID = rec_id[ch_trk_assoc_i];
  	  
  	  for(unsigned int mc_K_index = 0; mc_K_index < K_minus_simID_vect.size(); mc_K_index++)
  	  {
  	    if( K_minus_simID_vect.at(mc_K_index) == sim_id[ch_trk_assoc_i] ) 
  	    {
  	      K_minus_recID_vect.push_back(rec_id[ch_trk_assoc_i]);
  	      
  	      break; //stop when matching RC track is found
  	    }      	    
  	  
  	  } 	  
  	    	
  	}

    
    //PID using pfRICH and MC -> RC tracks
    
    TVector3 scat_e_mom_RC_pfRICH(0,0,0);
    double maxP_pfRICH = -1;//to store leading momentum

    int scat_e_simID_pfRICH = -1; //simID of scattered e

    for(unsigned int iChTrack = 0; iChTrack < reco_px_array.GetSize(); iChTrack++)
    {
      //this momentum is used for PID efficiency histograms too
      TVector3 rc_mom( reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack] );    
   
    
      int eta_bin_RC = -1;
      
      for(int j = 0; j < nEtaBins; j++) //loop over pT bins
      {
        if(rc_mom.Eta() > eta_bins[j] && rc_mom.Eta() <= eta_bins[j+1])
        {
          eta_bin_RC = j;
        }

      }      
      
      //find momentum bin of K candidate
      int mom_bin_K_RC = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(rc_mom.Mag() > mom_bins[j] && rc_mom.Mag() <= mom_bins[j+1])
        {
          mom_bin_K_RC = j;
        }

      } 
      
          
      //MC -> RC matching
      if(eta_bin_RC != -1)
      {
                
        for( unsigned int K_minus_recID_i = 0; K_minus_recID_i < K_minus_recID_vect.size(); K_minus_recID_i++ )
        {
          if( K_minus_recID_vect.at(K_minus_recID_i) == iChTrack )
          {
            h_p_K_minus_MC_RC[eta_bin_RC]->Fill(rc_mom.Mag());
            h_p_K_minus_MC_RC[nEtaBins]->Fill(rc_mom.Mag());
           
            //change PID to MC
            double pfRICH_K_prob;            
            
            pfRICH_K_prob = getPIDprob_pfRICH_MC(rc_mom, 1, 0);            
            
            double noPID_K_prob = 1. - pfRICH_K_prob;

            if(  pfRICH_K_prob < 0.9999 && pfRICH_K_prob > 0)
            {
              //cout<<pfRICH_pi_prob<<endl;

              h_p_pi_minus_MC_RC_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag(), noPID_K_prob);

              
              //pfRICH eta acceptance
              if(rc_mom.Eta() > -3.8 && rc_mom.Eta() < -1.5)
              {
                h_p_pi_minus_MC_RC_pfRICH[nEtaBins]->Fill(rc_mom.Mag(), noPID_K_prob);              
              }
            }  
            
            break; //stop when matchig RC pion is found
            
          }   

        }      
      
      }
      
      
      //____________________________________________________________________________________________
      
      //pfRICH reconstruction
      //first, find leading momentum charged track      
      
      
      //find corresponding simID of charged track
      int simID = -1.;
      
      for(unsigned int ids = 0; ids < rec_id.GetSize(); ids++)
      {
        if( rec_id[ids] == iChTrack ) simID = sim_id[ids];      
      }
      
      if(simID < 0 )
      {
        //add some way to track bad simIDs
        continue;
      }
      
      
      TVector3 mc_mom;
      
      if(mc_generatorStatus_array[simID] == 1) mc_mom.SetXYZ(mc_px_array[simID], mc_py_array[simID], mc_pz_array[simID]); //matching track MC momentum
      

      if( rc_mom.Mag() > maxP_pfRICH && reco_cahrge[iChTrack] < 0) //store negative charged particle only
      {
        maxP_pfRICH = rc_mom.Mag();
        
        scat_e_mom_RC_pfRICH.SetXYZ( reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack] );

        scat_e_simID_pfRICH = simID; //save simID of leading momenum track, simID found earlier
      }
      //___________________________________________________________________________________________________________________________________________________________
      
      

      //fill pfRICH PID efficiency histograms for pi/K/p
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(rc_mom, pfRICH_mtx, 0);
      
      int is_pfRICH_kaon = 0; //for K purity study

      if(good_pfRICH_mtx == 1)
      {
        
        
        //generate random number for PID in pfRICH
        double rndm_PID_pfRICH_loop = myRandom->Rndm();
        
        //find matchig MC track
        //pi-
        if( mc_pdg_array[simID] == -211 )
        {                    
          //pi mis-identified as K          
          if( pfRICH_mtx[1] > rndm_PID_pfRICH_loop) is_pfRICH_kaon = 1;         
        }

        //K-
        if( mc_pdg_array[simID] == -321 )
        {         
          
          //K identified as K      
          if( pfRICH_mtx[4] > rndm_PID_pfRICH_loop) is_pfRICH_kaon = 1;
        }

        //p-bar
        if( mc_pdg_array[simID] == -2212 )
        { 
          //p-bar mis-identified as K      
          if( pfRICH_mtx[6] > rndm_PID_pfRICH_loop) is_pfRICH_kaon = 1;      
        }
      } 

      //K purity
      //Q2 and y bin given by MC scattered e
      if( mom_bin_K_RC != -1 )
      {
        if( is_pfRICH_kaon == 1 &&  rc_mom.Eta() > -3.8 && rc_mom.Eta() < -1.5 )
        {
          if( mc_pdg_array[simID] == -321 ) 
          {          
            h_K_minus_purity_pfRICH_RC[mom_bin_K_RC][Q2_bin][y_bin]->Fill(1.5);
            
            if(eta_bin_RC != -1 )
            {
              h_K_minus_purity_p_eta_pfRICH_RC[mom_bin_K_RC][eta_bin_RC]->Fill(1.5);
              h_K_minus_purity_p_eta_pfRICH_RC[mom_bin_K_RC][nEtaBins]->Fill(1.5);           
            }          
          }      
          else 
          {
            h_K_minus_purity_pfRICH_RC[mom_bin_K_RC][Q2_bin][y_bin]->Fill(0.5);
            
            if(eta_bin_RC != -1 )
            {
              h_K_minus_purity_p_eta_pfRICH_RC[mom_bin_K_RC][eta_bin_RC]->Fill(0.5);
              h_K_minus_purity_p_eta_pfRICH_RC[mom_bin_K_RC][nEtaBins]->Fill(0.5);            
            }
             
          }
          
        }
      
      } 
      
    }//end for over RC particles
    
    //___________________________________________________________________________________________________________________________________________________________
*/

  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}

