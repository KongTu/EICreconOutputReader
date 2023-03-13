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

  const int nMomBins = 6;
  float const mom_bins[nMomBins+1] = { 1.,2.,4.,6.,8.,10.,12. };
  
  TRandom3 *myRandom = new TRandom3();
  
  
  //pfRICH eta acceptance mc_4mom.Eta() > -3.8 && mc_4mom.Eta() < -1.5
  const int nEtaBins = 4;
  float const eta_bins[nEtaBins+1] = { -3.8, -3, -2.5, -2, -1.5};

  //____________________________________________________
  //pi eCAL info
  //values from Dimitry Kalinkin

  double array_mom_bins[8] = {0.1, 0.2, 0.5, 1., 2., 5., 10., 20.};
  double array_pi_false_rate_85[8] = {0.6332832672698294, 0.7495818158985306, 0.00930384575910461, 0.001692827694491846, 0.0001898238241173789, 0.00020018016214593134, 0.000536412900269677, 0.0006430092230696459};
  double array_pi_false_rate_95[8] = {0.8838860654090215, 0.9228366502089709, 0.02228665375203912, 0.0036237053948322794, 0.00048096291971113834, 0.0006523112180272588, 0.0022770026949322946, 0.0018829746368912276};


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

  TH1F* h_hadron_zh = new TH1F("h_hadron_zh", "h_hadron_zh", 100, 0, 1);

  //MC histograms
  //p bins + p integrated
  //TH1F *h_eta_K_plus_lead_MC[nQ2bins][nyInelParBins];
  TH1F *h_eta_K_minus_lead_MC[nQ2bins][nyInelParBins];
  
  //TH1F *h_p_K_plus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_K_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  
  TH1F *h_p_pi_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  TH1F *h_p_p_minus_lead_MC_Q2_y[nQ2bins][nyInelParBins];
  
  
  TH1F *h_p_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  TH1F *h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  
  
  //for efficiency
  TH1F *h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins];
  
  TH1F *h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[nQ2bins][nyInelParBins];
  
  
 
  
  //MC -> RC histograms   
  //TH1F *h_eta_K_plus_lead_MC_RC_pfRICH[nQ2bins][nyInelParBins];
  //TH1F *h_eta_K_minus_lead_MC_RC_pfRICH[nQ2bins][nyInelParBins];
    
  //TH1F *h_p_K_plus_lead_MC_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  //TH1F *h_p_K_minus_lead_MC_RC_Q2_y_pfRICH[nQ2bins][nyInelParBins];
  
  //____________________________________________________________

  //K- 
  //TH1F *h_K_minus_lead_purity_pfRICH_MC_tot[nQ2bins][nyInelParBins];
  
  
  //TH1F *h_K_minus_lead_purity_pfRICH_MC_RC[nQ2bins][nyInelParBins];
  
/*  
  //K+
  TH1F *h_K_plus_lead_purity_pfRICH_MC_RC[nQ2bins][nyInelParBins]; 
*/


  for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
  {
    for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
    {
      //K eta histograms
      //MC
      //h_eta_K_plus_lead_MC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);       
      h_eta_K_minus_lead_MC[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
            
      //h_p_K_plus_lead_MC_Q2_y[Q2bin][y_bin] = new TH1F(Form("h_p_K_plus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_plus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      h_p_K_minus_lead_MC_Q2_y[Q2bin][y_bin] = new TH1F(Form("h_p_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      
      h_p_pi_minus_lead_MC_Q2_y[Q2bin][y_bin] = new TH1F(Form("h_p_pi_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_pi_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      h_p_p_minus_lead_MC_Q2_y[Q2bin][y_bin] = new TH1F(Form("h_p_p_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_p_minus_lead_MC_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      
      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_K_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin] = new TH1F(Form("h_p_K_as_K_minus_lead_MC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_as_K_minus_lead_MC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin), nMomBins, mom_bins);
         
      h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_pi_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_pi_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_p_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_p_as_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      
      h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_minus_lead_MC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
      h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2bin][y_bin] = new TH1F(Form("h_p_K_minus_lead_MC_pfRICH_eff_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_minus_lead_MC_pfRICH_Q2_eff_%i_y_%i" , Q2bin, y_bin), nMomBins, mom_bins);
      
      
      //h_K_minus_lead_purity_pfRICH_MC_tot[Q2bin][y_bin] = new TH1F(Form("h_K_minus_lead_purity_pfRICH_MC_tot_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_minus_lead_purity_pfRICH_MC_tot_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
           
      
      //MC->RC
      //h_eta_K_plus_lead_MC_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_plus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_plus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      //h_eta_K_minus_lead_MC_RC_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_eta_K_minus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_eta_K_minus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, -4, 0);
      
        
      //h_p_K_plus_lead_MC_RC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_K_plus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_plus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);    
      //h_p_K_minus_lead_MC_RC_Q2_y_pfRICH[Q2bin][y_bin] = new TH1F(Form("h_p_K_minus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_p_K_minus_lead_MC_RC_pfRICH_Q2_%i_y_%i" , Q2bin, y_bin), 100, 0, 20);
  

      //K eta purity histograms      
      //h_K_minus_lead_purity_pfRICH_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_K_minus_lead_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_minus_lead_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);
            
/*
      h_K_plus_lead_purity_pfRICH_MC_RC[Q2bin][y_bin] = new TH1F(Form("h_K_plus_lead_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), Form("h_K_plus_lead_purity_pfRICH_MC_RC_Q2_%i_y_%i" , Q2bin, y_bin), 2, 0, 2);      
*/
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

    //__________________________________________________________________________________

   
    //fill MC eta distributions
    
    //for leading K+/- (in pfRICH acceptance)
  	//TVector3 lead_K_minus_mom_MC(0,0,0);
  	//int mc_lead_K_PDG = -999;
  	//double maxP_lead_K_minus_MC = -99.;
  	//int p_bin_lead_K_minus_MC = -1;
  	
  	

  	
/*  	
  	TVector3 lead_K_plus_mom_MC(0,0,0);
  	//int mc_lead_K_plus_index = -1;
  	double maxP_lead_K_plus_MC = -99.;
  	int eta_bin_lead_K_plus_MC = -1;
  	
*/
    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{
  	  if(has_good_Q2_y_MC != 1) continue;

  	  if( mc_generatorStatus_array[imc] != 1 ) continue;

  		TLorentzVector mc_4mom;
  		
  		mc_4mom.SetXYZM(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc], mc_mass_array[imc]);
  		 		
  		
      //only par in pfRICH acceptance
  		//if(mc_4mom.Eta() < -3.8 && mc_4mom.Eta() > -1.5) continue;
	
  		//determine momentum and eta bin of particles in event     
      
      int mom_bin = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(mc_4mom.Vect().Mag() > mom_bins[j] && mc_4mom.Vect().Mag() <= mom_bins[j+1])
        {
          mom_bin = j;
        }

      }
      //_________________________________________________________________________
      
      float z_h = get_z_h(e_energy, scatMC.Vect(), p_energy, mc_4mom);
      
      h_hadron_zh->Fill(z_h);
      
      //if(/*mc_pdg_array[imc] == -321 &&*/ z_h > 0.2 && z_h < 1.0)
      if(/*mc_pdg_array[imc] == -321 &&*/ z_h > 0.0 && z_h < 1.0)
  		{
  		  TVector3 lead_K_minus_mom_MC(0,0,0);
  		  lead_K_minus_mom_MC.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
  		  
  		  //mc_lead_K_minus_index = imc;
  		  
  		  //p_bin_lead_K_minus_MC = mom_bin;
  		  
  		  //maxP_lead_K_minus_MC = mc_4mom.Mag(); 
  		  
  		  //mc_lead_K_PDG = mc_pdg_array[imc];
  		  
  		  
  		  //get pfRICH matrix
        double pfRICH_mtx[9];

        int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(mc_4mom.Vect(), pfRICH_mtx, 0);  //pi/K/p matrix    
        
        //generate random number for PID in pfRICH
        double rndm_PID_RC_loop = myRandom->Rndm();
        
        
        int is_pfRICH_kaon_MC = 0;
      	int is_pfRICH_kaon_MC_pi = 0; //mis-identified pi as K
      	int is_pfRICH_kaon_MC_p = 0; //mis-identified p as K
      	
      	double weight = 0;
        
        
   
        if(good_pfRICH_mtx == 1)
        {
          //find matchig MC track
          //pi
          if( fabs(mc_pdg_array[imc]) == 211 )
          {          
            //pi mis-identified as K          
            if( pfRICH_mtx[1] > rndm_PID_RC_loop) is_pfRICH_kaon_MC_pi = 1;
            
            weight = pfRICH_mtx[1];

          }

          //K
          if( fabs(mc_pdg_array[imc]) == 321 )
          {          
            //K identified as K          
            if( pfRICH_mtx[4] > rndm_PID_RC_loop) is_pfRICH_kaon_MC = 1;
            
            weight = pfRICH_mtx[4];

          }

          //p
          if( fabs(mc_pdg_array[imc]) == 2212 )
          {          
            //p identified as K          
            if( pfRICH_mtx[7] > rndm_PID_RC_loop) is_pfRICH_kaon_MC_p = 1;
            
            weight = pfRICH_mtx[7];
      
          }
                  
        }
        
        
        
        if(Q2_bin != -1 && y_bin != -1 )
      	{
      	  //true MC leading K- histos
      	  if( mc_pdg_array[imc] == -321)
      	  {
      	    
        	  h_eta_K_minus_lead_MC[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Eta()); //p integrated
        	  
        	  //h_eta_K_plus_lead_MC[Q2_bin][y_bin]->Fill(lead_K_plus_mom_MC.Eta()); 	  
        	  
        	  
        	  //fill p distributions in pfRICH eta acceptance
        	  if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    
        	    h_p_K_minus_lead_MC_Q2_y[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag());
        	    
        	    h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight); 
        	    h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight);
        	  
        	  }  	  
      	  
      	  }
      	  
      	  if( mc_pdg_array[imc] == -211)
      	  {
        	  //fill p distributions in pfRICH eta acceptance
        	  if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    
        	    h_p_pi_minus_lead_MC_Q2_y[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	    h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight); 
        	  
        	  }  	  
      	  
      	  }
      	  
      	  if( mc_pdg_array[imc] == -2212)
      	  {
      	    
        	  
        	  //fill p distributions in pfRICH eta acceptance
        	  if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {        	    
        	    h_p_p_minus_lead_MC_Q2_y[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	    
        	    h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight);        	  
        	  }  	  
      	  
      	  }
      	  
      	  if( mc_pdg_array[imc] == -211 || mc_pdg_array[imc] == -321 || mc_pdg_array[imc] == -2212)
      	  {
      	    if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  { 
        	    h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight); 
        	    h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag(), weight);
        	  }    	  
        	  
      	  }
/*	  
      	  //pi identified as K
      	  if( is_pfRICH_kaon_MC_pi )
      	  {
      	    if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    //negative particle
        	    if(mc_pdg_array[imc] < 0) h_p_pi_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	  
        	  }     	    
      	  }
      	  
      	  if( is_pfRICH_kaon_MC )
      	  {
      	    if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    
        	    //negative particle
        	    if(mc_pdg_array[imc] < 0) 
        	    {
        	      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	      h_p_K_as_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	    }
        	    
        	    
        	  
        	  }
      	  }
      	  
      	  if(is_pfRICH_kaon_MC_p)
      	  {
      	    if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    
        	    //negative particle
        	    if(mc_pdg_array[imc] < 0) h_p_p_as_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	  
        	  }
      	  }
      	  
      	  //all accepted hadrons as K
      	  if( is_pfRICH_kaon_MC_pi || is_pfRICH_kaon_MC || is_pfRICH_kaon_MC_p )
      	  {
      	  
      	    if( lead_K_minus_mom_MC.Eta() > -3.8 && lead_K_minus_mom_MC.Eta() < -1.5  )
        	  {
        	    
        	    //negative particle
        	    if(mc_pdg_array[imc] < 0)
        	      h_p_K_minus_lead_MC_Q2_y_pfRICH[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	    {
        	      h_p_K_minus_lead_MC_Q2_y_pfRICH_eff[Q2_bin][y_bin]->Fill(lead_K_minus_mom_MC.Mag()); 
        	    }
        	     
        	  
        	  }
        	  
        	  if(mc_pdg_array[imc] == -321)
      	    {  	      
      	      h_K_minus_lead_purity_pfRICH_MC_tot[Q2_bin][y_bin]->Fill(1.5);
      	    }
      	    else if( mc_pdg_array[imc] < 0 )
      	    {
      	      
      	      h_K_minus_lead_purity_pfRICH_MC_tot[Q2_bin][y_bin]->Fill(0.5); 
      	    }
      	  
      	  
      	  }  	  
      	     
        
*/
        }// end if bins
  		}//end if z_h

  	}//end second MC particle loop
  	//cout<<"MC loop end"<<endl;
  	
  	
  	
  	
  	
	
  	//___________________________________________________________________________________________________________________________________________________________

/*
    //MC -> RC matching
  	//find corresponding recID for scattered e
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
    
 	
  	TVector3 lead_K_minus_mom_RC_pfRICH(0,0,0);
  	double maxP_lead_K_minus_RC_pfRICH = -99.;
  	int mom_bin_lead_K_minus_RC_pfRICH = -1;
  	int lead_PDGid_muinus_RC = -99;
  	
  	
 	
  	TVector3 lead_K_plus_mom_RC_pfRICH(0,0,0);
  	double maxP_lead_K_plus_RC_pfRICH = -99.;
  	int eta_bin_lead_K_plus_RC_pfRICH = -1;


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
          //p identified as K          
          if( pfRICH_mtx[7] > rndm_PID_RC_loop) is_pfRICH_kaon_RC = 1;
    
        }
                
      }
      //_________________________________________________________________________
      
      
        
      
      if( is_pfRICH_kaon_RC == 1 )
      {
        //K- identified using pfRICH (with mis-PID)
        if( reco_PDG[iChTrack] < 0 )
        {
        
         
          //store RC leading K candidate
          if( rc_mom.Mag() > maxP_lead_K_minus_RC_pfRICH)
      		{
      		  lead_K_minus_mom_RC_pfRICH.SetXYZ(reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack]);
      		         		  
      		  mom_bin_lead_K_minus_RC_pfRICH = mom_bin_K_RC;
      		  
      		  maxP_lead_K_minus_RC_pfRICH = rc_mom.Mag();
      		  
      		  lead_PDGid_muinus_RC = reco_PDG[iChTrack];
      		  		  
      		}     
        }
                  
        //K+ identified using pfRICH (with mis-PID)
        if( reco_PDG[iChTrack] > 0 )
        {
        
                
          
          //store RC leading K candidate
          if( rc_mom.Mag() > maxP_lead_K_plus_RC_pfRICH)
      		{
      		  lead_K_plus_mom_RC_pfRICH.SetXYZ(reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack]);
      		         		  
      		  eta_bin_lead_K_plus_RC_pfRICH = eta_bin_RC;
      		  
      		  maxP_lead_K_plus_RC_pfRICH = rc_mom.Mag();
      		  
      		  lead_PDGid_plus_RC = reco_PDG[iChTrack];
      		  		  
      		}  
                  
        }
        
      } //end if PID and acceptance
        
          
    }//end loop over RC tracks
    
    //cout<<"RC loop end"<<endl;
    
    
    
    
    //fill RC leading K+- histos
  	if(Q2_bin_RC != -1 && y_bin_RC != -1)
  	{  
  	  if(lead_K_minus_mom_RC_pfRICH.Mag() > 0)
  	  {
  	    if(mom_bin_lead_K_minus_RC_pfRICH != -1 ) 
  	    {
  	      h_eta_K_minus_lead_MC_RC_pfRICH[mom_bin_lead_K_minus_RC_pfRICH][Q2_bin_RC][y_bin_RC]->Fill(lead_K_minus_mom_RC_pfRICH.Eta());
  	      h_p_K_minus_lead_MC_RC_Q2_y_pfRICH[mom_bin_lead_K_minus_RC_pfRICH][Q2_bin_RC][y_bin_RC]->Fill(lead_K_minus_mom_RC_pfRICH.Mag());
  	    }
  	    
  	    //p integrated
  	    h_eta_K_minus_lead_MC_RC_pfRICH[nMomBins][Q2_bin_RC][y_bin_RC]->Fill(lead_K_minus_mom_RC_pfRICH.Eta());  	    
  	    h_p_K_minus_lead_MC_RC_Q2_y_pfRICH[nMomBins][Q2_bin_RC][y_bin_RC]->Fill(lead_K_minus_mom_RC_pfRICH.Mag());
  	    
  	    if(lead_PDGid_muinus_RC == -321)
  	    {
  	      if(mom_bin_lead_K_minus_RC_pfRICH != -1 ) h_K_minus_lead_purity_pfRICH_MC_RC[mom_bin_lead_K_minus_RC_pfRICH][Q2_bin_RC][y_bin_RC]->Fill(1.5);
  	      h_K_minus_lead_purity_pfRICH_MC_RC[nMomBins][Q2_bin_RC][y_bin_RC]->Fill(1.5); //p integrated
  	    }
  	    else
  	    {
  	      if(mom_bin_lead_K_minus_RC_pfRICH != -1 ) h_K_minus_lead_purity_pfRICH_MC_RC[mom_bin_lead_K_minus_RC_pfRICH][Q2_bin_RC][y_bin_RC]->Fill(0.5);
  	      h_K_minus_lead_purity_pfRICH_MC_RC[nMomBins][Q2_bin_RC][y_bin_RC]->Fill(0.5); //p integrated
  	    }
  	  
  	  }
  	  
 	  
  	  if(lead_K_plus_mom_RC_pfRICH.Mag() > 0)
  	  {
  	    h_eta_K_plus_lead_MC_RC_pfRICH[Q2_bin_RC][y_bin_RC]->Fill(lead_K_plus_mom_RC_pfRICH.Eta());
  	    
  	    if(lead_PDGid_plus_RC == 321)
  	    {
  	      h_K_plus_lead_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(1.5);  	        	    
  	    }
  	    else
  	    {
  	      h_K_plus_lead_purity_pfRICH_MC_RC[Q2_bin_RC][y_bin_RC]->Fill(0.5);
  	    }
  	  
  	  }
 	  
  	}
  	
*/  	

  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}

