#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include "pleaseIncludeMe.h"

using namespace std;

int pfRICH_QA(TString inname="./fileLists/flieList.list", TString outname="test", float e_energy = 18, float p_energy = 275)
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

	//auto file = new TFile(inname);


  //TTreeReader tree_reader(tree);       // !the tree reader
  TTreeReader tree_reader(myChain);       // !the tree reader

  // MC particles
  TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
  //TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
  TTreeReaderArray<int> mc_generatorStatus_array = {tree_reader, "MCParticles.generatorStatus"};


  // Reconstructed particles
  TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  //TTreeReaderArray<float> reco_cahrge = {tree_reader, "ReconstructedChargedParticles.charge"};
  //TTreeReaderArray<int> reco_PDG = {tree_reader, "ReconstructedChargedParticles.PDG"};

  TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
  TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};



  //Reconstructed eCAL hits
  //for own clustering, does not have association info
 	//TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
  //TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
  //TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};



  //defining output file and histos.
  TString output_name_dir = outname;
	TFile* output = new TFile("output/"+output_name_dir+"-output_QA.root","RECREATE");

 
  

  //pi/K/p efficiency histogram
  TH1F *h_pi_pfRICH_PID_eff_MC[4]; //PID efficiency as a function of MC momentum - no PID baseline, good PID and mis-PID
  TH1F *h_K_pfRICH_PID_eff_MC[4];
  TH1F *h_p_pfRICH_PID_eff_MC[4];

  TH1F *h_pi_pfRICH_PID_eff_RC[4]; //PID efficiency as a function of RC momentum - no PID baseline, good PID and mis-PID
  TH1F *h_K_pfRICH_PID_eff_RC[4];
  TH1F *h_p_pfRICH_PID_eff_RC[4];

  TH1F *h_pi_pfRICH_PID_eff_MC_RC_match[4]; //PID efficiency as a function of MC momentum, from RC-MC matchig - no PID baseline, good PID and mis-PID
  TH1F *h_K_pfRICH_PID_eff_MC_RC_match[4];
  TH1F *h_p_pfRICH_PID_eff_MC_RC_match[4];
  
  TH2F *h_pi_pfRICH_MC_p_RC_p = new TH2F("h_pi_pfRICH_MC_p_RC_p", "h_pi_pfRICH_MC_p_RC_p", 100, 0, 20, 100, 0, 20);
  TH2F *h_K_pfRICH_MC_p_RC_p = new TH2F("h_K_pfRICH_MC_p_RC_p", "h_K_pfRICH_MC_p_RC_p", 100, 0, 20, 100, 0, 20);
  TH2F *h_p_pfRICH_MC_p_RC_p = new TH2F("h_p_pfRICH_MC_p_RC_p", "h_p_pfRICH_MC_p_RC_p", 100, 0, 20, 100, 0, 20);

  //e/pi PID efficiency histograms
  TH1F *h_e_pfRICH_PID_eff_MC[3];
  TH1F *h_e_pi_pfRICH_PID_eff_MC[3];

  TH1F *h_e_pfRICH_PID_eff_RC[3];
  TH1F *h_e_pi_pfRICH_PID_eff_RC[3];

  TH1F *h_e_pfRICH_PID_eff_MC_RC[3];
  TH1F *h_e_pi_pfRICH_PID_eff_MC_RC[3];
  
  
  //e/pi TOF PID efficiency histograms
  TH1F *h_e_pfRICH_TOF_PID_eff_MC[3];
  TH1F *h_e_pi_pfRICH_TOF_PID_eff_MC[3];

  //TH1F *h_e_pfRICH_TOF_PID_eff_RC[3];
  //TH1F *h_e_pi_pfRICH_TOF_PID_eff_RC[3];

  //TH1F *h_e_pfRICH_TOF_PID_eff_MC_RC[3];
  //TH1F *h_e_pi_pfRICH_TOF_PID_eff_MC_RC[3];

  for(unsigned int PID_bin = 0; PID_bin < 4; PID_bin++)
  {
    h_pi_pfRICH_PID_eff_MC[PID_bin] = new TH1F(Form("h_pi_pfRICH_PID_eff_MC_%i", PID_bin), Form("h_pi_pfRICH_PID_eff_MC_%i", PID_bin), 100, 0, 20); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_MC[PID_bin] = new TH1F(Form("h_K_pfRICH_PID_eff_MC_%i", PID_bin), Form("h_K_pfRICH_PID_eff_MC_%i", PID_bin), 100, 0, 20);
    h_p_pfRICH_PID_eff_MC[PID_bin] = new TH1F(Form("h_p_pfRICH_PID_eff_MC_%i", PID_bin), Form("h_p_pfRICH_PID_eff_MC_%i", PID_bin), 100, 0, 20);

    h_pi_pfRICH_PID_eff_RC[PID_bin] = new TH1F(Form("h_pi_pfRICH_PID_eff_RC_%i", PID_bin), Form("h_pi_pfRICH_PID_eff_RC_%i", PID_bin), 100, 0, 20); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_RC[PID_bin] = new TH1F(Form("h_K_pfRICH_PID_eff_RC_%i", PID_bin), Form("h_K_pfRICH_PID_eff_RC_%i", PID_bin), 100, 0, 20);
    h_p_pfRICH_PID_eff_RC[PID_bin] = new TH1F(Form("h_p_pfRICH_PID_eff_RC_%i", PID_bin), Form("h_p_pfRICH_PID_eff_RC_%i", PID_bin), 100, 0, 20);

    h_pi_pfRICH_PID_eff_MC_RC_match[PID_bin] = new TH1F(Form("h_pi_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), Form("h_pi_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), 100, 0, 20); //PID efficiency - good PID and mis-PID
    h_K_pfRICH_PID_eff_MC_RC_match[PID_bin] = new TH1F(Form("h_K_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), Form("h_K_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), 100, 0, 20);
    h_p_pfRICH_PID_eff_MC_RC_match[PID_bin] = new TH1F(Form("h_p_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), Form("h_p_pfRICH_PID_eff_MC_RC_match_%i", PID_bin), 100, 0, 20);

    if(PID_bin < 3)
    {
      //RICH PID
      h_e_pfRICH_PID_eff_MC[PID_bin] = new TH1F(Form("h_e_pfRICH_PID_eff_MC_%i", PID_bin), Form("h_e_pfRICH_PID_eff_MC_%i", PID_bin), 100, 0, 20);
      h_e_pi_pfRICH_PID_eff_MC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_PID_eff_MC_%i", PID_bin), Form("h_e_pi_pfRICH_PID_eff_MC_%i", PID_bin), 100, 0, 20);

      h_e_pfRICH_PID_eff_RC[PID_bin] = new TH1F(Form("h_e_pfRICH_PID_eff_RC_%i", PID_bin), Form("h_e_pfRICH_PID_eff_RC_%i", PID_bin), 100, 0, 20);
      h_e_pi_pfRICH_PID_eff_RC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_PID_eff_RC_%i", PID_bin), Form("h_e_pi_pfRICH_PID_eff_RC_%i", PID_bin), 100, 0, 20);

      h_e_pfRICH_PID_eff_MC_RC[PID_bin] = new TH1F(Form("h_e_pfRICH_PID_eff_MC_RC_%i", PID_bin), Form("h_e_pfRICH_PID_eff_MC_RC_%i", PID_bin), 100, 0, 20);
      h_e_pi_pfRICH_PID_eff_MC_RC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_PID_eff_MC_RC_%i", PID_bin), Form("h_e_pi_pfRICH_PID_eff_MC_RC_%i", PID_bin), 100, 0, 20);
      
      //TOF PID
      h_e_pfRICH_TOF_PID_eff_MC[PID_bin] = new TH1F(Form("h_e_pfRICH_TOF_PID_eff_MC_%i", PID_bin), Form("h_e_pfRICH_TOF_PID_eff_MC_%i", PID_bin), 100, 0, 20);
      h_e_pi_pfRICH_TOF_PID_eff_MC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_TOF_PID_eff_MC_%i", PID_bin), Form("h_e_pi_pfRICH_TOF_PID_eff_MC_%i", PID_bin), 100, 0, 20);

      //h_e_pfRICH_TOF_PID_eff_RC[PID_bin] = new TH1F(Form("h_e_pfRICH_TOF_PID_eff_RC_%i", PID_bin), Form("h_e_pfRICH_TOF_PID_eff_RC_%i", PID_bin), 100, 0, 20);
      //h_e_pi_pfRICH_TOF_PID_eff_RC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_TOF_PID_eff_RC_%i", PID_bin), Form("h_e_pi_pfRICH_TOF_PID_eff_RC_%i", PID_bin), 100, 0, 20);

      //h_e_pfRICH_TOF_PID_eff_MC_RC[PID_bin] = new TH1F(Form("h_e_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin), Form("h_e_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin), 100, 0, 20);
      //h_e_pi_pfRICH_TOF_PID_eff_MC_RC[PID_bin] = new TH1F(Form("h_e_pi_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin), Form("h_e_pi_pfRICH_TOF_PID_eff_MC_RC_%i", PID_bin), 100, 0, 20);
    }
  }
  //______________________________________________________________________________________________________________________________________________________________________________________
  
  cout<<"nEntries "<<myChain->GetEntries()<<endl;
	tree_reader.SetEntriesRange(0, myChain->GetEntries());
	//tree_reader.SetEntriesRange(0, 1000);
	

  //Long64_t iEvent = 0;

  while (tree_reader.Next())
  {
    //iEvent++;
  	//MC analysis

    //finding the scattering electron
  	TLorentzVector scatMC(0,0,0,0);

  	for(int imc = 0; imc < mc_px_array.GetSize(); imc++)
  	{
  	  if( mc_generatorStatus_array[imc] != 1 ) continue;

  		TVector3 mctrk(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

  		
  		//fill pfRICH PID efficiency histograms for pi/K/p for MC tracks
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(mctrk, pfRICH_mtx, 0);  
 
      if(good_pfRICH_mtx == 1)
      {
        //find matchig MC track
        //pi-
        if( mc_pdg_array[imc] == -211 )
        {
          h_pi_pfRICH_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_mtx[0] > 0 && pfRICH_mtx[0] < 1) h_pi_pfRICH_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_mtx[0]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[1] > 0 && pfRICH_mtx[1] < 1) h_pi_pfRICH_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_mtx[1]);//pi identified as K
          if(pfRICH_mtx[2] > 0 && pfRICH_mtx[2] < 1) h_pi_pfRICH_PID_eff_MC[3]->Fill(mctrk.Mag(), pfRICH_mtx[2]);//pi identified as p
          
          
        }

        //K-
        if( mc_pdg_array[imc] == -321 )
        {
          h_K_pfRICH_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_mtx[4] > 0 && pfRICH_mtx[4] < 1) h_K_pfRICH_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_mtx[4]);//K identified as K (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[3] > 0 && pfRICH_mtx[3] < 1) h_K_pfRICH_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_mtx[3]);//K identified as pi
          if(pfRICH_mtx[5] > 0 && pfRICH_mtx[5] < 1) h_K_pfRICH_PID_eff_MC[3]->Fill(mctrk.Mag(), pfRICH_mtx[5]);//K identified as p          

        }

        //p-bar
        if( mc_pdg_array[imc] == -2212 )
        {
          h_p_pfRICH_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_mtx[8] > 0 && pfRICH_mtx[8] < 1) h_p_pfRICH_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_mtx[8]);//p identified as p (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[6] > 0 && pfRICH_mtx[6] < 1) h_p_pfRICH_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_mtx[6]);//p identified as pi
          if(pfRICH_mtx[7] > 0 && pfRICH_mtx[7] < 1) h_p_pfRICH_PID_eff_MC[3]->Fill(mctrk.Mag(), pfRICH_mtx[7]);//p identified as K
          
        }
      }


       //e/pi pfRICH PID
      double pfRICH_mtx_e[9];

      int good_pfRICH_mtx_e = getPIDprob_pfRICH_mtx(mctrk, pfRICH_mtx_e, 1);

      if(good_pfRICH_mtx_e == 1)
      {
        //find matchig MC track
        //e
        if( mc_pdg_array[imc] == 11 )
        {
          h_e_pfRICH_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_mtx_e[0] > 0 && pfRICH_mtx_e[0] < 1) h_e_pfRICH_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_mtx_e[0]);//e identified as e (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[1] > 0 && pfRICH_mtx_e[1] < 1) h_e_pfRICH_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_mtx_e[1]);//e identified as pi
        }

        //pi-
        if( mc_pdg_array[imc] == -211 )
        {
          h_e_pi_pfRICH_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_mtx_e[3] > 0 && pfRICH_mtx_e[3] < 1) h_e_pi_pfRICH_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_mtx_e[3]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[2] > 0 && pfRICH_mtx_e[2] < 1) h_e_pi_pfRICH_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_mtx_e[2]);//pi identified as e
        }


      }
     

      //e/pi pfRICH TOF PID
      double pfRICH_TOF_mtx_e[9];

      int good_pfRICH_TOF_mtx_e = getPIDprob_TOF_mtx(mctrk, pfRICH_TOF_mtx_e, 0);
      
      //cout<<"test 2"<<endl;

      if(good_pfRICH_TOF_mtx_e == 1)
      {
        //find matchig MC track
        //e
        if( mc_pdg_array[imc] == 11 )
        {
          h_e_pfRICH_TOF_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_TOF_mtx_e[0] > 0 && pfRICH_TOF_mtx_e[0] < 1) h_e_pfRICH_TOF_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_TOF_mtx_e[0]);//e identified as e (diagonal term of pfRICH_mtx)
          if(pfRICH_TOF_mtx_e[1] > 0 && pfRICH_TOF_mtx_e[1] < 1) h_e_pfRICH_TOF_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_TOF_mtx_e[1]);//e identified as pi
        }
 
        //pi-
        if( mc_pdg_array[imc] == -211 )
        {
          h_e_pi_pfRICH_TOF_PID_eff_MC[0]->Fill(mctrk.Mag());//no PID baseline
          if(pfRICH_TOF_mtx_e[3] > 0 && pfRICH_TOF_mtx_e[3] < 1) h_e_pi_pfRICH_TOF_PID_eff_MC[1]->Fill(mctrk.Mag(), pfRICH_TOF_mtx_e[3]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_TOF_mtx_e[2] > 0 && pfRICH_TOF_mtx_e[2] < 1) h_e_pi_pfRICH_TOF_PID_eff_MC[2]->Fill(mctrk.Mag(), pfRICH_TOF_mtx_e[2]);//pi identified as e
        }

      }


  	} //end MC particles for loop
  	
  	//____________________________________________________________________________________________
    
    
    //PID using pfRICH and MC -> RC tracks


    int scat_e_simID_pfRICH = -1; //simID of scattered e
    //cout<<"RC track loop start"<<endl;
    for(unsigned int iChTrack = 0; iChTrack < reco_px_array.GetSize(); iChTrack++)
    {
      //this momentum is used for PID efficiency histograms too
      TVector3 rc_mom( reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack] );   
      
      
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
      //___________________________________________________________________________________________________________________________________________________________
      
      

      //fill pfRICH PID efficiency histograms for pi/K/p
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(rc_mom, pfRICH_mtx, 0);
      

      
      if( mc_pdg_array[simID] == -211 )
      {
        h_pi_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
      }

      if(good_pfRICH_mtx == 1)
      {
        
        
        //find matchig MC track
        //pi-
        if( mc_pdg_array[simID] == -211 )
        {
          //h_pi_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
          if(pfRICH_mtx[0] > 0 && pfRICH_mtx[0] < 1) h_pi_pfRICH_PID_eff_RC[1]->Fill(rc_mom.Mag(), pfRICH_mtx[0]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[1] > 0 && pfRICH_mtx[1] < 1) h_pi_pfRICH_PID_eff_RC[2]->Fill(rc_mom.Mag(), pfRICH_mtx[1]);//pi identified as K
          if(pfRICH_mtx[2] > 0 && pfRICH_mtx[2] < 1) h_pi_pfRICH_PID_eff_RC[3]->Fill(rc_mom.Mag(), pfRICH_mtx[2]);//pi identified as p
          
          if(mc_generatorStatus_array[simID] == 1)
          {
            h_pi_pfRICH_PID_eff_MC_RC_match[0]->Fill(mc_mom.Mag());//no PID baseline
            if(pfRICH_mtx[0] > 0 && pfRICH_mtx[0] < 1) h_pi_pfRICH_PID_eff_MC_RC_match[1]->Fill(mc_mom.Mag(), pfRICH_mtx[0]);//pi identified as pi (diagonal term of pfRICH_mtx)
            if(pfRICH_mtx[1] > 0 && pfRICH_mtx[1] < 1) h_pi_pfRICH_PID_eff_MC_RC_match[2]->Fill(mc_mom.Mag(), pfRICH_mtx[1]);//pi identified as K
            if(pfRICH_mtx[2] > 0 && pfRICH_mtx[2] < 1) h_pi_pfRICH_PID_eff_MC_RC_match[3]->Fill(mc_mom.Mag(), pfRICH_mtx[2]);//pi identified as p
            
            h_pi_pfRICH_MC_p_RC_p->Fill(mc_mom.Mag(), rc_mom.Mag());
          }

          

        }

        //K-
        if( mc_pdg_array[simID] == -321 )
        {
          h_K_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
          if(pfRICH_mtx[4] > 0 && pfRICH_mtx[4] < 1) h_K_pfRICH_PID_eff_RC[1]->Fill(rc_mom.Mag(), pfRICH_mtx[4]);//K identified as K (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[3] > 0 && pfRICH_mtx[3] < 1) h_K_pfRICH_PID_eff_RC[2]->Fill(rc_mom.Mag(), pfRICH_mtx[3]);//K identified as pi
          if(pfRICH_mtx[5] > 0 && pfRICH_mtx[5] < 1) h_K_pfRICH_PID_eff_RC[3]->Fill(rc_mom.Mag(), pfRICH_mtx[5]);//K identified as p

          if(mc_generatorStatus_array[simID] == 1)
          {
            h_K_pfRICH_PID_eff_MC_RC_match[0]->Fill(mc_mom.Mag());//no PID baseline
            if(pfRICH_mtx[4] > 0 && pfRICH_mtx[4] < 1) h_K_pfRICH_PID_eff_MC_RC_match[1]->Fill(mc_mom.Mag(), pfRICH_mtx[4]);//pi identified as pi (diagonal term of pfRICH_mtx)
            if(pfRICH_mtx[3] > 0 && pfRICH_mtx[3] < 1) h_K_pfRICH_PID_eff_MC_RC_match[2]->Fill(mc_mom.Mag(), pfRICH_mtx[3]);//pi identified as K
            if(pfRICH_mtx[5] > 0 && pfRICH_mtx[5] < 1) h_K_pfRICH_PID_eff_MC_RC_match[3]->Fill(mc_mom.Mag(), pfRICH_mtx[5]);//pi identified as p
            
            h_K_pfRICH_MC_p_RC_p->Fill(mc_mom.Mag(), rc_mom.Mag());
          }
          


        }

        //p-bar
        if( mc_pdg_array[simID] == -2212 )
        {
          h_p_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
          if(pfRICH_mtx[8] > 0 && pfRICH_mtx[8] < 1) h_p_pfRICH_PID_eff_RC[1]->Fill(rc_mom.Mag(), pfRICH_mtx[8]);//p identified as p (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx[6] > 0 && pfRICH_mtx[6] < 1) h_p_pfRICH_PID_eff_RC[2]->Fill(rc_mom.Mag(), pfRICH_mtx[6]);//p identified as pi
          if(pfRICH_mtx[7] > 0 && pfRICH_mtx[7] < 1) h_p_pfRICH_PID_eff_RC[3]->Fill(rc_mom.Mag(), pfRICH_mtx[7]);//p identified as K


          if(mc_generatorStatus_array[simID] == 1)
          {
            h_p_pfRICH_PID_eff_MC_RC_match[0]->Fill(mc_mom.Mag());//no PID baseline
            if(pfRICH_mtx[8] > 0 && pfRICH_mtx[8] < 1) h_p_pfRICH_PID_eff_MC_RC_match[1]->Fill(mc_mom.Mag(), pfRICH_mtx[8]);//pi identified as pi (diagonal term of pfRICH_mtx)
            if(pfRICH_mtx[6] > 0 && pfRICH_mtx[6] < 1) h_p_pfRICH_PID_eff_MC_RC_match[2]->Fill(mc_mom.Mag(), pfRICH_mtx[6]);//pi identified as K
            if(pfRICH_mtx[7] > 0 && pfRICH_mtx[7] < 1) h_p_pfRICH_PID_eff_MC_RC_match[3]->Fill(mc_mom.Mag(), pfRICH_mtx[7]);//pi identified as p
            
            h_p_pfRICH_MC_p_RC_p->Fill(mc_mom.Mag(), rc_mom.Mag());
          }

      
        }
      }


       //e/pi pfRICH PID
      double pfRICH_mtx_e[9];

      int good_pfRICH_mtx_e = getPIDprob_pfRICH_mtx(rc_mom, pfRICH_mtx_e, 1);
      
      
      if(good_pfRICH_mtx_e == 1)
      {

        //find matchig MC track
        //e
        if( mc_pdg_array[simID] == 11 )
        {
          h_e_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
          if(pfRICH_mtx_e[0] > 0 && pfRICH_mtx_e[0] < 1) h_e_pfRICH_PID_eff_RC[1]->Fill(rc_mom.Mag(), pfRICH_mtx_e[0]);//e identified as e (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[1] > 0 && pfRICH_mtx_e[1] < 1) h_e_pfRICH_PID_eff_RC[2]->Fill(rc_mom.Mag(), pfRICH_mtx_e[1]);//e identified as pi


          h_e_pfRICH_PID_eff_MC_RC[0]->Fill(mc_mom.Mag());//no PID baseline
          if(pfRICH_mtx_e[0] > 0 && pfRICH_mtx_e[0] < 1) h_e_pfRICH_PID_eff_MC_RC[1]->Fill(mc_mom.Mag(), pfRICH_mtx_e[0]);//e identified as e (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[1] > 0 && pfRICH_mtx_e[1] < 1) h_e_pfRICH_PID_eff_MC_RC[2]->Fill(mc_mom.Mag(), pfRICH_mtx_e[1]);//e identified as pi
        }

        //pi-
        if( mc_pdg_array[simID] == -211 )
        {
          h_e_pi_pfRICH_PID_eff_RC[0]->Fill(rc_mom.Mag());//no PID baseline
          if(pfRICH_mtx_e[3] > 0 && pfRICH_mtx_e[3] < 1) h_e_pi_pfRICH_PID_eff_RC[1]->Fill(rc_mom.Mag(), pfRICH_mtx_e[3]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[2] > 0 && pfRICH_mtx_e[2] < 1) h_e_pi_pfRICH_PID_eff_RC[2]->Fill(rc_mom.Mag(), pfRICH_mtx_e[2]);//pi identified as e


          h_e_pi_pfRICH_PID_eff_MC_RC[0]->Fill(mc_mom.Mag());//no PID baseline
          if(pfRICH_mtx_e[3] > 0 && pfRICH_mtx_e[3] < 1) h_e_pi_pfRICH_PID_eff_MC_RC[1]->Fill(mc_mom.Mag(), pfRICH_mtx_e[3]);//pi identified as pi (diagonal term of pfRICH_mtx)
          if(pfRICH_mtx_e[2] > 0 && pfRICH_mtx_e[2] < 1) h_e_pi_pfRICH_PID_eff_MC_RC[2]->Fill(mc_mom.Mag(), pfRICH_mtx_e[2]);//pi identified as e
        }


      }
      
      
            
    }//end for over RC particles
    
    //___________________________________________________________________________________________________________________________________________________________


  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}

