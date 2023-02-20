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
  TH1F *h_eta_K_plus_MC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC[nEtaBins+1];
  TH1F *h_p_K_minus_MC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_MC_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_pfRICH[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_MC_pfRICH[nEtaBins+1];
  
  

  TH1F *h_eta_K_plus_MC_RC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_RC[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_RC[nEtaBins+1];
  TH1F *h_p_K_minus_MC_RC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_MC_RC_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_MC_RC_pfRICH[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_MC_RC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_MC_RC_pfRICH[nEtaBins+1];
  
/*  
  
  TH1F *h_eta_K_plus_RC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_RC[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_RC[nEtaBins+1];
  TH1F *h_p_K_minus_RC[nEtaBins+1];    

  
  TH1F *h_eta_K_plus_RC_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_RC_pfRICH[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_plus_RC_pfRICH[nEtaBins+1];
  TH1F *h_p_K_minus_RC_pfRICH[nEtaBins+1]; 
*/
  //____________________________________________________________


  TH1F *h_K_minus_purity_pfRICH_MC[nQ2bins][nyInelParBins];
  TH1F *h_K_minus_purity_pfRICH_MC_RC[nQ2bins][nyInelParBins];
  //TH1F *h_K_minus_purity_pfRICH_RC[nQ2bins][nyInelParBins];
  
  TH1F *h_K_minus_purity_p_eta_pfRICH_MC[nEtaBins+1];
  TH1F *h_K_minus_purity_p_eta_pfRICH_MC_RC[nEtaBins+1];
  //TH1F *h_K_minus_purity_p_eta_pfRICH_RC[nEtaBins+1];
  
  
  TH1F *h_K_plus_purity_pfRICH_MC[nQ2bins][nyInelParBins];
  TH1F *h_K_plus_purity_pfRICH_MC_RC[nQ2bins][nyInelParBins];
  //TH1F *h_K_plus_purity_pfRICH_RC[nQ2bins][nyInelParBins];
  
  TH1F *h_K_plus_purity_p_eta_pfRICH_MC[nEtaBins+1];
  TH1F *h_K_plus_purity_p_eta_pfRICH_MC_RC[nEtaBins+1];
  //TH1F *h_K_plus_purity_p_eta_pfRICH_RC[nEtaBins+1];


  
  for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
  {
    for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
    {
      //K eta histograms
      h_eta_K_plus_MC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);       
      h_eta_K_minus_MC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);

      h_eta_K_plus_MC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      h_eta_K_minus_MC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      
/*       
      h_eta_K_plus_RC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_RC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);       
      h_eta_K_minus_RC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_RC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);

      h_eta_K_plus_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      h_eta_K_minus_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
*/      
      
      h_eta_K_plus_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);       
      h_eta_K_minus_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);

      h_eta_K_plus_MC_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      h_eta_K_minus_MC_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
  

      //K eta purity histograms
      h_K_minus_purity_pfRICH_MC[Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
      h_K_minus_purity_pfRICH_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
      //h_K_minus_purity_pfRICH_RC[Q2bin][y_bin] = new TH1F(Form("h_K_minus_purity_pfRICH_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_minus_purity_pfRICH_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
      
      h_K_plus_purity_pfRICH_MC[Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_MC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
      h_K_plus_purity_pfRICH_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
      //h_K_plus_purity_pfRICH_RC[Q2bin][y_bin] = new TH1F(Form("h_K_plus_purity_pfRICH_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_plus_purity_pfRICH_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);

    }

  }
  
  
  //momentum histograms in eta bins
  for(unsigned int eta_bin = 0; eta_bin < nEtaBins+1; eta_bin++) 
  {
    //K p histograms
    h_p_K_plus_MC[eta_bin] = new TH1F(Form("h_p_K_plus_MC_eta_%i", eta_bin), Form("h_p_K_plus_MC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC[eta_bin] = new TH1F(Form("h_p_K_minus_MC_eta_%i", eta_bin), Form("h_p_K_minus_MC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_MC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_MC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_MC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_MC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_MC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
/*    
    h_p_K_plus_RC[eta_bin] = new TH1F(Form("h_p_K_plus_RC_eta_%i", eta_bin), Form("h_p_K_plus_RC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_RC[eta_bin] = new TH1F(Form("h_p_K_minus_RC_eta_%i", eta_bin), Form("h_p_K_minus_RC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
*/    
    
    h_p_K_plus_MC_RC[eta_bin] = new TH1F(Form("h_p_K_plus_MC_RC_eta_%i", eta_bin), Form("h_p_K_plus_MC_RC_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_RC[eta_bin] = new TH1F(Form("h_p_K_minus_MC_RC_eta_%i", eta_bin), Form("h_p_K_minus_MC_RC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_K_plus_MC_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_plus_MC_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_plus_MC_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_K_minus_MC_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_K_minus_MC_RC_pfRICH_eta_%i", eta_bin), Form("h_p_K_minus_MC_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
 
        
    
    //K p purity histograms
    h_K_minus_purity_p_eta_pfRICH_MC[eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin), 2, 0, 2);
    //h_K_minus_purity_p_eta_pfRICH_RC[eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_RC_eta_%i" , eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_RC_eta_%i" , eta_bin), 2, 0, 2);
    h_K_minus_purity_p_eta_pfRICH_MC_RC[eta_bin] = new TH1F(Form("h_K_minus_purity_p_eta_pfRICH_MC_RC_eta_%i" , eta_bin), Form("h_K_minus_purity_p_eta_pfRICH_MC_RC_eta_%i" , eta_bin), 2, 0, 2);
    
    h_K_plus_purity_p_eta_pfRICH_MC[eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_MC_eta_%i" , eta_bin), 2, 0, 2);
    //h_K_plus_purity_p_eta_pfRICH_RC[eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_RC_eta_%i" , eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_RC_eta_%i" , eta_bin), 2, 0, 2);
    h_K_plus_purity_p_eta_pfRICH_MC_RC[eta_bin] = new TH1F(Form("h_K_plus_purity_p_eta_pfRICH_MC_RC_eta_%i" , eta_bin), Form("h_K_plus_purity_p_eta_pfRICH_MC_RC_eta_%i" , eta_bin), 2, 0, 2);
    
  
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
      
      if( eta_bin != -1)
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
              h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->Fill(1.5);
              h_K_minus_purity_p_eta_pfRICH_MC[nEtaBins]->Fill(1.5);           
            }      
            else //pi- and p-bar mis-identified as K-
            {
              h_K_minus_purity_p_eta_pfRICH_MC[eta_bin]->Fill(0.5);
              h_K_minus_purity_p_eta_pfRICH_MC[nEtaBins]->Fill(0.5);
            }
          
          
          }
          
          //K+ identified using pfRICH (with mis-PID)
          if( mc_pdg_array[imc] > 0 )
          {
          
            h_p_K_plus_MC_pfRICH[eta_bin]->Fill(mc_mom.Mag());
            h_p_K_plus_MC_pfRICH[nEtaBins]->Fill(mc_mom.Mag());          
          
            if( mc_pdg_array[imc] == 321 ) 
            {       
              h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->Fill(1.5);
              h_K_plus_purity_p_eta_pfRICH_MC[nEtaBins]->Fill(1.5);           
            }      
            else if( mc_pdg_array[imc] > 0 )//pi+ and p mis-identified as K+
            {
              h_K_plus_purity_p_eta_pfRICH_MC[eta_bin]->Fill(0.5);
              h_K_plus_purity_p_eta_pfRICH_MC[nEtaBins]->Fill(0.5);
            }
            
          }
          
        } //end if PID and acceptance
      
      } //end if bins
      //cout<<"p end"<<endl;
      
      //___________________________________________________________________________________________________________________________________________________________
      
      
      //fill eta distributions for MC particles
      if(Q2_bin < 0 || y_bin < 0 ) continue;
      
  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 &&  mc_mom.Pt() < maxP)
  		

  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		  h_eta_K_plus_MC[Q2_bin][y_bin]->Fill(mc_mom.Eta());	   
  		}
  		//K-
  		if(mc_pdg_array[imc] == -321)
  		{
  		  h_eta_K_minus_MC[Q2_bin][y_bin]->Fill(mc_mom.Eta());
  		}
  		//_____________________________________________________________________________
  		
  		
  		if( is_pfRICH_kaon_MC == 1 )
      {
        //K- identified using pfRICH (with mis-PID)
        if( mc_pdg_array[imc] < 0 )
        {
          h_eta_K_minus_MC_pfRICH[Q2_bin][y_bin]->Fill(mc_mom.Eta());
          
          if( mc_pdg_array[imc] == -321 ) 
          {       
            h_K_minus_purity_pfRICH_MC[Q2_bin][y_bin]->Fill(1.5);
                      
          }      
          else //pi- and p-bar mis-identified as K-
          {
            h_K_minus_purity_pfRICH_MC[Q2_bin][y_bin]->Fill(0.5);
          }
        
        }
        
        //K+ identified using pfRICH (with mis-PID)
        if( mc_pdg_array[imc] > 0 )
        {
          h_eta_K_plus_MC_pfRICH[Q2_bin][y_bin]->Fill(mc_mom.Eta());
          
          if( mc_pdg_array[imc] == 321 ) 
          {       
            h_K_plus_purity_pfRICH_MC[Q2_bin][y_bin]->Fill(1.5);
                      
          }      
          else //pi+ and p mis-identified as K+
          {
            h_K_plus_purity_pfRICH_MC[Q2_bin][y_bin]->Fill(0.5);
          }
        
        }       
      
      }   

  	}//end second MC particle loop
  	//___________________________________________________________________________________________________________________________________________________________

  	
  	//reconstructed charged particles analysis
  	//need to update
  	

    //cout<<"RC - MC match start"<<endl;
    //MC -> RC matching
  	//find corresponding recID for scattered e and K
  	int scat_e_recID = -1;
  	
  	for(unsigned int ch_trk_assoc_i = 0; ch_trk_assoc_i < sim_id.GetSize(); ch_trk_assoc_i++)
  	{
  	  if( mc_elect_index == sim_id[ch_trk_assoc_i] ) scat_e_recID = rec_id[ch_trk_assoc_i];	
  	}
  	
  	
  	//find RC scattered electron
  	TVector3 scatRC(0,0,0);
  	
    for(unsigned int iChTrack = 0; iChTrack < reco_px_array.GetSize(); iChTrack++)
    {
      //MC -> RC matched scattered electron
      if( scat_e_recID == iChTrack )
      {
        scatRC.SetXYZ(reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack]);
        
        break; //stop after matching RC track is found      
      }      
    
    }
    
    double y_inelPar_e_RC = getInelParamElectron_2(p_energy, e_energy, scatRC );

    double Q2_electron_RC = getQ2elec( e_energy, scatRC);

    //flag for wide MC Q2 and y cut
    int has_good_Q2_y_RC = 0;
    
    if( Q2_electron_RC > 1. && y_inelPar_e_RC < 0.95) has_good_Q2_y_RC = 1;

    
    //find bins
    int Q2_bin_RC = -1;

    for(int j = 0; j < nQ2bins; j++)
    {
      if(Q2_electron_RC > Q2_bins[j] && Q2_electron_RC <= Q2_bins[j+1])
      {
        Q2_bin_RC = j;
      }
     }


    int y_bin_RC = -1;

    for(int j = 0; j < nyInelParBins; j++)
    {
      if(y_inelPar_e_RC > y_bins[j] && y_inelPar_e_RC <= y_bins[j+1])
      {
        y_bin_RC = j;
      }
    }

    //__________________________________________________________________________________
    
    //PID using pfRICH and MC -> RC tracks
    
    //TVector3 scat_e_mom_RC_pfRICH(0,0,0);
    //double maxP_pfRICH = -1;//to store leading momentum


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
      
                          
 
      if(has_good_Q2_y_RC != 1) continue; 		
      
  			            
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(rc_mom, pfRICH_mtx, 0);  //pi/K/p matrix    
      
      //generate random number for PID in pfRICH
      double rndm_PID_RC_loop = myRandom->Rndm();
      
      int is_pfRICH_kaon_RC = 0;
 
      if(good_pfRICH_mtx == 1)
      {
        //find matchig MC track
        //pi
        if( fabs(reco_PDG[iChTrack]) == 211 )
        {          
          //pi mis-identified as K          
          if( pfRICH_mtx[1] > rndm_PID_RC_loop) is_pfRICH_kaon_RC = 1;
        }

        //K
        if( fabs(reco_PDG[iChTrack]) == 321 )
        {          
          //K identified as K          
          if( pfRICH_mtx[4] > rndm_PID_RC_loop) is_pfRICH_kaon_RC = 1;
        }

        //p
        if( fabs(reco_PDG[iChTrack]) == 2212 )
        {          
          //K identified as K          
          if( pfRICH_mtx[7] > rndm_PID_RC_loop) is_pfRICH_kaon_RC = 1;          
        }
      }
      //_________________________________________________________________________      


      //fill RC p distributions      
      //cout<<"p start"<<endl;
      
      if( eta_bin_RC != -1)
      {
        //pure MC
        //K-
        if(reco_PDG[iChTrack] == -321) 
        {
          h_p_K_minus_MC_RC[eta_bin_RC]->Fill(rc_mom.Mag());                 
          h_p_K_minus_MC_RC[nEtaBins]->Fill(rc_mom.Mag()); //baseline in pfRICH eta acceptance

                 
        }
        
        //K+
        if(reco_PDG[iChTrack] == 321) 
        {
          h_p_K_plus_MC_RC[eta_bin_RC]->Fill(rc_mom.Mag());
          h_p_K_plus_MC_RC[nEtaBins]->Fill(rc_mom.Mag()); //baseline in pfRICH eta acceptance

                  
        }       
           
        //K purity in p and eta bins
        //Q2 and y bin given by RC scattered e
      
        if( is_pfRICH_kaon_RC == 1 )
        {
          //K- identified using pfRICH (with mis-PID)
          if( reco_PDG[iChTrack] < 0 )
          {
          
            h_p_K_minus_MC_RC_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag());
            h_p_K_minus_MC_RC_pfRICH[nEtaBins]->Fill(rc_mom.Mag());
            
            if( reco_PDG[iChTrack] == -321 ) 
            {       
              h_K_minus_purity_p_eta_pfRICH_MC_RC[eta_bin_RC]->Fill(1.5);
              h_K_minus_purity_p_eta_pfRICH_MC_RC[nEtaBins]->Fill(1.5);           
            }      
            else //pi- and p-bar mis-identified as K-
            {
              h_K_minus_purity_p_eta_pfRICH_MC_RC[eta_bin_RC]->Fill(0.5);
              h_K_minus_purity_p_eta_pfRICH_MC_RC[nEtaBins]->Fill(0.5);
            }
          
          
          }
                    
          //K+ identified using pfRICH (with mis-PID)
          if( reco_PDG[iChTrack] > 0 )
          {
          
            h_p_K_plus_MC_RC_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag());
            h_p_K_plus_MC_RC_pfRICH[nEtaBins]->Fill(rc_mom.Mag());          
          
            if( reco_PDG[iChTrack] == 321 ) 
            {       
              h_K_plus_purity_p_eta_pfRICH_MC_RC[eta_bin_RC]->Fill(1.5);
              h_K_plus_purity_p_eta_pfRICH_MC_RC[nEtaBins]->Fill(1.5);           
            }      
            else //pi+ and p mis-identified as K+
            {
              h_K_plus_purity_p_eta_pfRICH_MC_RC[eta_bin_RC]->Fill(0.5);
              h_K_plus_purity_p_eta_pfRICH_MC_RC[nEtaBins]->Fill(0.5);
            }
            
          }
          
        } //end if PID and acceptance
      
      } //end if bins
      //cout<<"p end"<<endl;
      
      //___________________________________________________________________________________________________________________________________________________________
      
      
      //fill eta distributions for MC particles
      if(Q2_bin_RC < 0 || y_bin_RC < 0 ) continue;
      
  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 &&  mc_mom.Pt() < maxP)
  		

  		//K+
  		if(reco_PDG[iChTrack] == 321)
  		{
  		  h_eta_K_plus_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());	   
  		}
  		//K-
  		if(reco_PDG[iChTrack] == -321)
  		{
  		  h_eta_K_minus_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());
  		}
  		//_____________________________________________________________________________
  		
  		
  		if( is_pfRICH_kaon_RC == 1 )
      {
        //K- identified using pfRICH (with mis-PID)
        if( reco_PDG[iChTrack] < 0 )
        {
          h_eta_K_minus_MC_RC_pfRICH[Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());
          
          if( reco_PDG[iChTrack] == -321 ) 
          {       
            h_K_minus_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(1.5);
                      
          }      
          else //pi- and p-bar mis-identified as K-
          {
            h_K_minus_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(0.5);
          }
        
        }
        
        //K+ identified using pfRICH (with mis-PID)
        if( reco_PDG[iChTrack] > 0 )
        {
          h_eta_K_plus_MC_RC_pfRICH[Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());
          
          if( reco_PDG[iChTrack] == 321 ) 
          {       
            h_K_plus_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(1.5);
                      
          }      
          else //pi+ and p mis-identified as K+
          {
            h_K_plus_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(0.5);
          }
        
        }
        
      
      } //end if pfRICH PID            
          
    }//end loop over RC tracks

  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}

