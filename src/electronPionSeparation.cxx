#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include "pleaseIncludeMe.h"

using namespace std;

//int electronPionSeparation(TString inname="input/input.root",TString outname="test")
int electronPionSeparation(TString inname="./fileLists/flieList.list", TString outname="test", float e_energy = 18, float p_energy = 275)
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

  //Reconstructed eCAL clusters
  TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
  TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
  TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};

  TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClustersAssociations.recID"};
  TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClustersAssociations.simID"};

  //Reconstructed eCAL hits
  //for own clustering, does not have association info
 	//TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
  //TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};
  //TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};



  //defining output file and histos.
  TString output_name_dir = outname;
	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

  //MC histograms
  TH1F* h_energy_MC = new TH1F("h_energy_MC","E_{MC} (GeV)",100,0,20);
  TH1F* h_energy_zoom_MC = new TH1F("h_energy_zoom_MC","E_{MC} (GeV)",20,e_energy-2,e_energy+2);

  TH1F* h_momentum_MC = new TH1F("h_momentum_MC","p_{MC} (GeV/c)",100,0,20);

  TH1F* h_Q2_MC = new TH1F("h_Q2_MC",";Q^{2}",100,0,20);
  TH1F* h_Q2_zoom_MC = new TH1F("h_Q2_zoom_MC",";Q^{2}",100,0,4);

  TH1F *h_y_inelPar_MC = new TH1F("h_y_inelPar_MC", "h_y_inelPar_MC", 100, 0, 1);
  TH1F *h_y_inelPar_zoom_MC = new TH1F("h_y_inelPar_zoom_MC", "h_y_inelPar_zoom_MC", 100, 0, 0.01);

  //eta distributions in multiple momentum bins
  TH1F *h_eta_scat_ele[nMomBins];

  TH1F *h_eta_ele[nMomBins];
  TH1F *h_eta_pi_plus[nMomBins];
  TH1F *h_eta_K_plus[nMomBins];
  TH1F *h_eta_proton[nMomBins];


  TH1F *h_eta_positron[nMomBins];
  TH1F *h_eta_pi_minus[nMomBins];
  TH1F *h_eta_K_minus[nMomBins];
  TH1F *h_eta_anti_proton[nMomBins];
  
  
  TH1F *h_eta_pi_minus_eCAL_85[nMomBins];
  TH1F *h_eta_pi_minus_eCAL_95[nMomBins];
    
  TH1F *h_eta_pi_minus_eCAL_85_pfRICH[nMomBins];
  TH1F *h_eta_pi_minus_eCAL_95_pfRICH[nMomBins];
    
  TH1F *h_eta_pi_minus_pfRICH[nMomBins];
  TH1F *h_eta_K_minus_pfRICH[nMomBins];
  TH1F *h_eta_anti_proton_pfRICH[nMomBins];

  //p distributions in multiple eta bins
  //use just wide Q^2 and y selection
  
  //eta bins defined higher + 1 eta bin for pfRICH acceptance window
  TH1F *h_p_scat_ele[nEtaBins+1];
  
  TH1F *h_p_pi_minus[nEtaBins+1];
  
  TH1F *h_p_pi_minus_eCAL_85[nEtaBins+1];
  TH1F *h_p_pi_minus_eCAL_95[nEtaBins+1];
  
  TH1F *h_p_pi_minus_eCAL_85_pfRICH[nEtaBins+1];
  TH1F *h_p_pi_minus_eCAL_95_pfRICH[nEtaBins+1];
  
  TH1F *h_p_pi_minus_pfRICH[nEtaBins+1];
  
  //MC tracks as a function of RC momentum
  //direct MC -> RC matching
  TH1F *h_p_scat_ele_MC_RC[nEtaBins+1];
  
  TH1F *h_p_pi_minus_MC_RC[nEtaBins+1];
  
  TH1F *h_p_pi_minus_MC_RC_eCAL_85[nEtaBins+1];
  TH1F *h_p_pi_minus_MC_RC_eCAL_95[nEtaBins+1];
  
  TH1F *h_p_pi_minus_MC_RC_eCAL_85_pfRICH[nEtaBins+1];
  TH1F *h_p_pi_minus_MC_RC_eCAL_95_pfRICH[nEtaBins+1];
  
  TH1F *h_p_pi_minus_MC_RC_pfRICH[nEtaBins+1];

  //____________________________________________________________

  //reco histograms with eCAL
  TH1F* h_energy_RC = new TH1F("h_energy_RC","E_{RC} (GeV)",100,0,20);

  TH1F* h_momentum_RC = new TH1F("h_momentum_RC","p_{RC} (GeV/c)",100,0,20);

  TH1F *h_E_over_p_RC = new TH1F("h_E_over_p_RC", "h_E_over_p_RC", 120, 0, 1.2);

  TH1F* h_Q2_RC = new TH1F("h_Q2_RC",";Q^{2}",100,0,20);

  TH1F *h_y_inelPar_RC = new TH1F("h_y_inelPar_RC", "h_y_inelPar_RC", 100, 0, 1);

  //reconstructed scattered electron with eCAL
  TH1F *h_eta_scat_ele_RC_eCAL[nMomBins];  
  TH1F *h_eta_scat_ele_RC_eCAL_E_over_p[nMomBins];
  
  //momentum distributions in eta bins
  TH1F *h_p_scat_ele_RC_eCAL[nEtaBins+1];  
  TH1F *h_p_scat_ele_RC_eCAL_E_over_p[nEtaBins+1];
  
  

  //background - negative charged particles
  //TH1F *h_eta_neg_ch_part_RC[nMomBins][nQ2bins][nyInelParBins];

  //scattered e purity
  TH1F *h_scat_ele_purity_eCAL[nMomBins];  
  TH1F *h_scat_ele_purity_eCAL_E_over_p[nMomBins];
  
  TH1F *h_scat_ele_purity_p_eCAL[nEtaBins+1];  
  TH1F *h_scat_ele_purity_p_eCAL_E_over_p[nEtaBins+1];
  
  //_______________________________________________________________

  //reco histograms with pfRICH
  TH1F* h_momentum_RC_pfRICH = new TH1F("h_momentum_RC_pfRICH","h_momentum_RC_pfRICH",100,0,20);

  TH1F* h_Q2_RC_pfRICH = new TH1F("h_Q2_RC_pfRICH","h_Q2_RC_pfRICH",100,0,20);

  TH1F *h_y_inelPar_RC_pfRICH = new TH1F("h_y_inelPar_RC_pfRICH", "h_momentum_RC_pfRICH", 100, 0, 1);

  TH1F *h_eta_scat_ele_RC_lead_p[nMomBins];
  TH1F *h_eta_scat_ele_RC_pfRICH[nMomBins];
  
  //momentum distributions in eta bins
  TH1F *h_p_scat_ele_RC_lead_p[nEtaBins+1];
  TH1F *h_p_scat_ele_RC_pfRICH[nEtaBins+1];  



  //for e purity with pfRICH selection
  TH1F *h_scat_ele_purity_lead_p[nMomBins];
  TH1F *h_scat_ele_purity_pfRICH[nMomBins];
  
  TH1F *h_scat_ele_purity_p_lead_p[nEtaBins+1];
  TH1F *h_scat_ele_purity_p_pfRICH[nEtaBins+1];



  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {   
    
    //MC histograms
    h_eta_scat_ele[mom_bin] = new TH1F(Form("h_eta_scat_ele_mom_%i" , mom_bin), Form("h_eta_scat_ele_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_ele[mom_bin] = new TH1F(Form("h_eta_ele_mom_%i" , mom_bin), Form("h_eta_ele_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_pi_plus[mom_bin] = new TH1F(Form("h_eta_pi_plus_mom_%i" , mom_bin), Form("h_eta_pi_plus_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_K_plus[mom_bin] = new TH1F(Form("h_eta_K_plus_mom_%i" , mom_bin), Form("h_eta_K_plus_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_proton[mom_bin] = new TH1F(Form("h_eta_proton_mom_%i" , mom_bin), Form("h_eta_proton_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_positron[mom_bin] = new TH1F(Form("h_eta_positron_mom_%i" , mom_bin), Form("h_eta_positron_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_pi_minus[mom_bin] = new TH1F(Form("h_eta_pi_minus_mom_%i" , mom_bin), Form("h_eta_pi_minus_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_pi_minus_eCAL_85[mom_bin] = new TH1F(Form("h_eta_pi_minus_eCAL_85_mom_%i" , mom_bin), Form("h_eta_pi_minus_eCAL_85_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_pi_minus_eCAL_95[mom_bin] = new TH1F(Form("h_eta_pi_minus_eCAL_95_mom_%i" , mom_bin), Form("h_eta_pi_minus_eCAL_95_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_pi_minus_pfRICH[mom_bin] = new TH1F(Form("h_eta_pi_minus_pfRICH_mom_%i" , mom_bin), Form("h_eta_pi_minus_pfRICH_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_pi_minus_eCAL_85_pfRICH[mom_bin] = new TH1F(Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i" , mom_bin), Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i" , mom_bin), 100, -4, 0);
    h_eta_pi_minus_eCAL_95_pfRICH[mom_bin] = new TH1F(Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i" , mom_bin), Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_K_minus[mom_bin] = new TH1F(Form("h_eta_K_minus_mom_%i" , mom_bin), Form("h_eta_K_minus_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_K_minus_pfRICH[mom_bin] = new TH1F(Form("h_eta_K_minus_pfRICH_mom_%i" , mom_bin), Form("h_eta_K_minus_pfRICH_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_anti_proton[mom_bin] = new TH1F(Form("h_eta_anti_proton_mom_%i" , mom_bin), Form("h_eta_anti_proton_mom_%i" , mom_bin), 100, -4, 0);

    h_eta_anti_proton_pfRICH[mom_bin] = new TH1F(Form("h_eta_anti_proton_pfRICH_mom_%i" , mom_bin), Form("h_eta_anti_proton_pfRICH_mom_%i" , mom_bin), 100, -4, 0);
    //___________________________________________________________________________________________________________________________________________________________

    //reconstructed particle histograms
    //eCAL
    h_eta_scat_ele_RC_eCAL[mom_bin] = new TH1F(Form("h_eta_scat_ele_RC_eCAL_mom_%i" , mom_bin), Form("h_eta_scat_ele_RC_eCAL_mom_%i" , mom_bin), 100, -4, 0);
    
    h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin] = new TH1F(Form("h_eta_scat_ele_RC_eCAL_E_over_p_mom_%i" , mom_bin), Form("h_eta_scat_ele_RC_eCAL_E_over_p_mom_%i" , mom_bin), 100, -4, 0);

    h_scat_ele_purity_eCAL[mom_bin] = new TH1F(Form("h_scat_ele_purity_eCAL_mom_%i" , mom_bin), Form("h_scat_ele_purity_eCAL_mom_%i" , mom_bin), 2, 0, 2);
    
    h_scat_ele_purity_eCAL_E_over_p[mom_bin] = new TH1F(Form("h_scat_ele_purity_eCAL_E_over_p_mom_%i" , mom_bin), Form("h_scat_ele_purity_eCAL_E_over_p_mom_%i" , mom_bin), 2, 0, 2);

    //pfRICH
    h_eta_scat_ele_RC_lead_p[mom_bin] = new TH1F(Form("h_eta_scat_ele_RC_lead_p_mom_%i" , mom_bin), Form("h_eta_scat_ele_RC_lead_p_mom_%i" , mom_bin), 100, -4, 0);
    
    h_eta_scat_ele_RC_pfRICH[mom_bin] = new TH1F(Form("h_eta_scat_ele_RC_pfRICH_mom_%i" , mom_bin), Form("h_eta_scat_ele_RC_pfRICH_mom_%i" , mom_bin), 100, -4, 0);

    
    h_scat_ele_purity_lead_p[mom_bin] = new TH1F(Form("h_scat_ele_purity_lead_p_mom_%i" , mom_bin), Form("h_scat_ele_purity_lead_p_mom_%i" , mom_bin), 2, 0, 2);        
    h_scat_ele_purity_pfRICH[mom_bin] = new TH1F(Form("h_scat_ele_purity_pfRICH_mom_%i" , mom_bin), Form("h_scat_ele_purity_pfRICH_mom_%i" , mom_bin), 2, 0, 2);



  }
  
  //momentum histograms in eta bins
  for(unsigned int eta_bin = 0; eta_bin < nEtaBins+1; eta_bin++) 
  {
    //MC histograms
    h_p_scat_ele[eta_bin] = new TH1F(Form("h_p_scat_ele_eta_%i", eta_bin), Form("h_p_scat_ele_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus[eta_bin] = new TH1F(Form("h_p_pi_minus_eta_%i", eta_bin), Form("h_p_pi_minus_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_eCAL_85[eta_bin] = new TH1F(Form("h_p_pi_minus_eCAL_85_eta_%i", eta_bin), Form("h_p_pi_minus_eCAL_85_eta_%i", eta_bin), 200, 0, 20);
    h_p_pi_minus_eCAL_95[eta_bin] = new TH1F(Form("h_p_pi_minus_eCAL_95_eta_%i", eta_bin), Form("h_p_pi_minus_eCAL_95_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_eCAL_85_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_eCAL_85_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_eCAL_85_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_pi_minus_eCAL_95_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_eCAL_95_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_eCAL_95_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    //MC tracks as a function of RC momentum
    //direct MC -> RC matching
    h_p_scat_ele_MC_RC[eta_bin] = new TH1F(Form("h_p_scat_ele_MC_RC_eta_%i", eta_bin), Form("h_p_scat_ele_MC_RC_eta_%i", eta_bin), 200, 0, 20);
        
    h_p_pi_minus_MC_RC[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_MC_RC_eCAL_85[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_eCAL_85_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_eCAL_85_eta_%i", eta_bin), 200, 0, 20);
    h_p_pi_minus_MC_RC_eCAL_95[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_eCAL_95_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_eCAL_95_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_MC_RC_eCAL_85_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_eCAL_85_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_eCAL_85_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    h_p_pi_minus_MC_RC_eCAL_95_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_eCAL_95_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_eCAL_95_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    h_p_pi_minus_MC_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_pi_minus_MC_RC_pfRICH_eta_%i", eta_bin), Form("h_p_pi_minus_MC_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    //eCAL histograms
    h_p_scat_ele_RC_eCAL[eta_bin] = new TH1F(Form("h_p_scat_ele_RC_eCAL_eta_%i", eta_bin), Form("h_p_scat_ele_RC_eCAL_eta_%i", eta_bin), 200, 0, 20);  
    h_p_scat_ele_RC_eCAL_E_over_p[eta_bin] = new TH1F(Form("h_p_scat_ele_RC_eCAL_E_over_p_eta_%i", eta_bin), Form("h_p_scat_ele_RC_eCAL_E_over_p_eta_%i", eta_bin), 200, 0, 20);
    
    h_scat_ele_purity_p_eCAL[eta_bin] = new TH1F(Form("h_scat_ele_purity_p_eCAL_eta_%i" , eta_bin), Form("h_scat_ele_purity_p_eCAL_eta_%i" , eta_bin), 2, 0, 2);        
    h_scat_ele_purity_p_eCAL_E_over_p[eta_bin] = new TH1F(Form("h_scat_ele_purity_p_eCAL_E_over_p_eta_%i" , eta_bin), Form("h_scat_ele_purity_p_eCAL_E_over_p_eta_%i" , eta_bin), 2, 0, 2);
    
    //pfRICH histograms
    h_p_scat_ele_RC_lead_p[eta_bin] = new TH1F(Form("h_p_scat_ele_RC_lead_p_eta_%i", eta_bin), Form("h_p_scat_ele_RC_lead_p_eta_%i", eta_bin), 200, 0, 20);
    h_p_scat_ele_RC_pfRICH[eta_bin] = new TH1F(Form("h_p_scat_ele_RC_pfRICH_eta_%i", eta_bin), Form("h_p_scat_ele_RC_pfRICH_eta_%i", eta_bin), 200, 0, 20);
    
    h_scat_ele_purity_p_lead_p[eta_bin] = new TH1F(Form("h_scat_ele_purity_p_lead_p_eta_%i" , eta_bin), Form("h_scat_ele_purity_p_lead_p_eta_%i" , eta_bin), 2, 0, 2);        
    h_scat_ele_purity_p_pfRICH[eta_bin] = new TH1F(Form("h_scat_ele_purity_p_pfRICH_eta_%i" , eta_bin), Form("h_scat_ele_purity_p_pfRICH_eta_%i" , eta_bin), 2, 0, 2);
    
   
  
  }
  //______________________________________________________________________________________________________________________________________________________________________________________
  

	tree_reader.SetEntriesRange(0, myChain->GetEntries());

  while (tree_reader.Next())
  {
  	//MC analysis

    //finding the scattering electron
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

  		//_______________________________________________________________

  		//fill pfRICH PID efficiency histograms for pi/K/p for MC tracks
      //get pfRICH matrix
      double pfRICH_mtx[9];

      int good_pfRICH_mtx = getPIDprob_pfRICH_mtx(mctrk, pfRICH_mtx, 0);
      
      
      //generate random number for PID in pfRICH
      TRandom3 *rndm_MC_loop = new TRandom3();
      double rndm_PID_MC_loop = rndm_MC_loop->Rndm();
      
      int is_pfRICH_kaon_MC = 0;
 
      if(good_pfRICH_mtx == 1)
      {
        //find matchig MC track
        //pi-
        if( mc_pdg_array[imc] == -211 )
        {
                  
          //pi mis-identified as K          
          if( pfRICH_mtx[1] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;
        }

        //K-
        if( mc_pdg_array[imc] == -321 )
        { 
          //K identified as K          
          if( pfRICH_mtx[4] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;

        }

        //p-bar
        if( mc_pdg_array[imc] == -2212 )
        {
          
          //K identified as K          
          if( pfRICH_mtx[7] > rndm_PID_MC_loop) is_pfRICH_kaon_MC = 1;
          
        }
       }



  	} //end MC particles for loop
  	
  	//____________________________________________________________________________________________
    

  	double y_inelPar_e = getInelParamElectron_2(p_energy, e_energy, scatMC.Vect() );

    double Q2_electron = getQ2elec( e_energy, scatMC.Vect());


    //flag for wide MC Q2 and y cut
    int has_good_Q2_y_MC = 0;
    
    //if( Q2_electron > 1. && Q2_electron < 20. && y_inelPar_e > 0.01 && y_inelPar_e < 0.95) has_good_Q2_y_MC = 1;
    if( Q2_electron > 1. && y_inelPar_e < 0.95) has_good_Q2_y_MC = 1;

    //fill truth scattered electron energy
    if(has_good_Q2_y_MC && Q2_electron < 10.) //additional Q2 cut to match Kong's selection
    {
      h_energy_MC->Fill(scatMC.E());
    	h_energy_zoom_MC->Fill(scatMC.E());

    	h_momentum_MC->Fill(scatMC.P());

      h_y_inelPar_MC->Fill(y_inelPar_e);
      h_y_inelPar_zoom_MC->Fill(y_inelPar_e);

      h_Q2_MC->Fill(Q2_electron);
      h_Q2_zoom_MC->Fill(Q2_electron);
    }


    //find bins

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

    //fill MC p distributions
    //scattered electrons
    if(has_good_Q2_y_MC == 1)
    {
      if( eta_bin_scat_e != -1 ) h_p_scat_ele[eta_bin_scat_e]->Fill(scatMC.P());
      if(scatMC.Eta() > -3.8 && scatMC.Eta() < -1.5) h_p_scat_ele[nEtaBins]->Fill(scatMC.P());
    }
    
    
    //fille MC eta distributions
    if( !( mom_bin_scat_e < 0) ) h_eta_scat_ele[mom_bin_scat_e]->Fill(scatMC.Eta());   

    
    vector<int> pi_simID_vect;
    vector<float> pi_eta_MC; //MC eta of pion

    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{

  	  if( mc_generatorStatus_array[imc] != 1 ) continue;

  		TVector3 mc_mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
  		TLorentzVector mc_4mom(0,0,0,0);  
  		mc_4mom.SetVectM(mc_mom, mc_mass_array[imc]);
  		

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
      
      //fill p distributions for MC particles
      if(has_good_Q2_y_MC == 1)
      {
        if(eta_bin != -1)
        {
          //pi-
          if(mc_pdg_array[imc] == -211) 
          {
            h_p_pi_minus[eta_bin]->Fill(mc_mom.Mag());          
            if(mc_mom.Eta() > -3.8 && mc_mom.Eta() < -1.5)
            {
              h_p_pi_minus[nEtaBins]->Fill(mc_mom.Mag()); //baseline in pfRICH eta acceptance            
              pi_simID_vect.push_back(imc);
              pi_eta_MC.push_back(mc_mom.Eta());
            }
            
            
            //suppression factors for eCAL
            //default value is 1 - i.e. no suppression
            double eCAL_suppress_85 = 1.;
            double eCAL_suppress_95 = 1.;

            //eta acceptance of eCAL
            if( mc_mom.Eta() > -3.14 && mc_mom.Eta() < -1.87 )
            {
              eCAL_suppress_85 = g_pi_false_rate_85->Eval(mc_mom.Mag());
              eCAL_suppress_95 = g_pi_false_rate_95->Eval(mc_mom.Mag());
              
              h_p_pi_minus_eCAL_85[eta_bin]->Fill(mc_mom.Mag(), eCAL_suppress_85);
              h_p_pi_minus_eCAL_95[eta_bin]->Fill(mc_mom.Mag(), eCAL_suppress_95);
              
              h_p_pi_minus_eCAL_85[nEtaBins]->Fill(mc_mom.Mag(), eCAL_suppress_85);
              h_p_pi_minus_eCAL_95[nEtaBins]->Fill(mc_mom.Mag(), eCAL_suppress_95);
            }
            
            
            //change PID to MC
            double pfRICH_pi_prob;
            
            //use e/pi PID table for p < 5 GeV/c
            if(mc_mom.Mag() < 5)
            {
              pfRICH_pi_prob = getPIDprob_pfRICH_MC(mc_4mom, 1, 1);        
            }
            else
            {
              pfRICH_pi_prob = getPIDprob_pfRICH_MC(mc_4mom, 0, 0);
            }
            
            double noPID_pi_prob = 1. - pfRICH_pi_prob;


            if(  pfRICH_pi_prob < 0.9999 && pfRICH_pi_prob > 0)
            {
              //cout<<pfRICH_pi_prob<<endl;

              h_p_pi_minus_pfRICH[eta_bin]->Fill(mc_mom.Mag(), noPID_pi_prob);

              h_p_pi_minus_eCAL_85_pfRICH[eta_bin]->Fill(mc_mom.Mag(), eCAL_suppress_85*noPID_pi_prob);
              h_p_pi_minus_eCAL_95_pfRICH[eta_bin]->Fill(mc_mom.Mag(), eCAL_suppress_95*noPID_pi_prob);
              
              //pfRICH eta acceptance
              if(mc_mom.Eta() > -3.8 && mc_mom.Eta() < -1.5)
              {
                h_p_pi_minus_pfRICH[nEtaBins]->Fill(mc_mom.Mag(), noPID_pi_prob);

                h_p_pi_minus_eCAL_85_pfRICH[nEtaBins]->Fill(mc_mom.Mag(), eCAL_suppress_85*noPID_pi_prob);
                h_p_pi_minus_eCAL_95_pfRICH[nEtaBins]->Fill(mc_mom.Mag(), eCAL_suppress_95*noPID_pi_prob);
              
              }
            }
          
          }
                  
        }                           
      }
      
      //______________________________________________________________________________________________________________________________________________________________________________________
      
      //fill eta distributions for MC particles
      if( mom_bin < 0) continue;      


  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 &&  mc_mom.Pt() < maxP)
  		if(mc_pdg_array[imc] == 11 && mc_mom.Mag() < maxP)
  		{
  			h_eta_ele[mom_bin]->Fill(mc_mom.Eta());
  		}
  		
  		//positrons
  		if(mc_pdg_array[imc] == -11 )
  		//if(mc_pdg_array[imc] == 11 && mc_mom.Mag() < maxP)
  		{
  			h_eta_positron[mom_bin]->Fill(mc_mom.Eta());
  		}
  		//_____________________________________________________________________________

  		//pi+
  		if(mc_pdg_array[imc] == 211)
  		{
  		   h_eta_pi_plus[mom_bin]->Fill(mc_mom.Eta());
  		}
  		//pi-
  		if(mc_pdg_array[imc] == -211)
  		{
        h_eta_pi_minus[mom_bin]->Fill(mc_mom.Eta());

        //suppression factors for eCAL
        //default value is 1 - i.e. no suppression
        double eCAL_suppress_85 = 1.;
        double eCAL_suppress_95 = 1.;

        //eta acceptance of eCAL
        if( mc_mom.Eta() > -3.14 && mc_mom.Eta() < -1.87 )
        {
          eCAL_suppress_85 = g_pi_false_rate_85->Eval(mc_mom.Mag());
          eCAL_suppress_95 = g_pi_false_rate_95->Eval(mc_mom.Mag());

          h_eta_pi_minus_eCAL_85[mom_bin]->Fill(mc_mom.Eta(), eCAL_suppress_85);
          h_eta_pi_minus_eCAL_95[mom_bin]->Fill(mc_mom.Eta(), eCAL_suppress_95);
        }


        //change PID to MC
        double pfRICH_pi_prob;
        
        //use e/pi PID table for p < % GeV/c
        if(mc_mom.Mag() < 5)
        {
          pfRICH_pi_prob = getPIDprob_pfRICH_MC(mc_4mom, 1, 1);        
        }
        else
        {
          pfRICH_pi_prob = getPIDprob_pfRICH_MC(mc_4mom, 0, 0);
        }
        
        double noPID_pi_prob = 1. - pfRICH_pi_prob;


        if(  pfRICH_pi_prob < 0.9999 && pfRICH_pi_prob > 0)
        {
          //cout<<pfRICH_pi_prob<<endl;

          h_eta_pi_minus_pfRICH[mom_bin]->Fill(mc_mom.Eta(), noPID_pi_prob);

          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin]->Fill(mc_mom.Eta(), eCAL_suppress_85*noPID_pi_prob);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin]->Fill(mc_mom.Eta(), eCAL_suppress_95*noPID_pi_prob);
        }



  		}
  		//_____________________________________________________________________________

  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		   h_eta_K_plus[mom_bin]->Fill(mc_mom.Eta());
  		}
  		//K-
  		if(mc_pdg_array[imc] == -321)
  		{
  		  h_eta_K_minus[mom_bin]->Fill(mc_mom.Eta());

  		  double pfRICH_K_prob = getPIDprob_pfRICH_MC(mc_4mom, 1, 0);
        double noPID_K_prob = 1 - pfRICH_K_prob;

        if( pfRICH_K_prob > 0 && pfRICH_K_prob < 0.9999)
        {
          h_eta_K_minus_pfRICH[mom_bin]->Fill(mc_mom.Eta(), noPID_K_prob);
        }

  		}
  		//_____________________________________________________________________________

  		//proton
  		if(mc_pdg_array[imc] == 2212)
  		{
  		   h_eta_proton[mom_bin]->Fill(mc_mom.Eta());
  		}
  		//anti-proton
  		if(mc_pdg_array[imc] == -2212)
  		{
  		  h_eta_anti_proton[mom_bin]->Fill(mc_mom.Eta());

  		  double pfRICH_anti_proton_prob = getPIDprob_pfRICH_MC(mc_4mom, 2, 0);
        double noPID_anti_proton_prob = 1 - pfRICH_anti_proton_prob;


        if( pfRICH_anti_proton_prob > 0 && pfRICH_anti_proton_prob < 0.9999)
        {

          h_eta_anti_proton_pfRICH[mom_bin]->Fill(mc_mom.Eta(), noPID_anti_proton_prob);
        }
  		}
  		//___________________________________________________________________________________________________________________________________________________________


  	}//end second MC particle loop
  	//___________________________________________________________________________________________________________________________________________________________
  	//cout<<"mc loop 2 end"<<endl;
  	
  	//reconstructed charged particles analysis

  	//using eCAL
		//find scattered electron candidate
		//look for cluster with highest energy

		//leading cluster
    double maxEnergy=-99.;
    int sim_id_scat_e_RC = -1;

    for(int iclus=0; iclus < em_energy_array.GetSize(); iclus++)
    {
  	  if(em_energy_array[iclus] > maxEnergy)
  	  {
  		  maxEnergy = em_energy_array[iclus];
  		  sim_id_scat_e_RC = em_sim_id_array[iclus];
  	  }
    }

    h_energy_RC->Fill(maxEnergy);

    //find charged track corresponding to leading culster
    //save scattered e candidate momentum

    TVector3 scat_e_mom_RC(0,0,0);
    int scat_e_iTrack = -1; //index of scattered e in charged track array

    for(unsigned int iChTrack = 0; iChTrack < sim_id.GetSize(); iChTrack++)
    {
      if( sim_id[iChTrack] == sim_id_scat_e_RC && reco_cahrge[iChTrack] < 0 )
      {
        scat_e_mom_RC.SetXYZ( reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack] );
        scat_e_iTrack = iChTrack;

        break; //stop when corresponding ch. track is found
      }

    }
    

    if(maxEnergy > 0)
    {

      h_momentum_RC->Fill(scat_e_mom_RC.Mag());

      h_E_over_p_RC->Fill( maxEnergy/scat_e_mom_RC.Mag() );

     
      
      double y_inelPar_e_RC = getInelParamElectron_2(p_energy, e_energy, scat_e_mom_RC);

      double Q2_electron_RC = getQ2elec( e_energy, scat_e_mom_RC);

      h_Q2_RC->Fill(Q2_electron_RC);
      h_y_inelPar_RC->Fill(y_inelPar_e_RC);


      //find bins for reconstructed scattered e and charged particles
      int mom_bin_scat_e_RC = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(scat_e_mom_RC.Mag() > mom_bins[j] && scat_e_mom_RC.Mag() <= mom_bins[j+1])
        {
          mom_bin_scat_e_RC = j;
        }

      }
      
      
      int eta_bin_scat_e_RC = -1;

      for(int j = 0; j < nEtaBins; j++) //loop over pT bins
      {
        if(scat_e_mom_RC.Eta() > eta_bins[j] && scat_e_mom_RC.Eta() <= eta_bins[j+1])
        {
          eta_bin_scat_e_RC = j;
        }

      }

      //__________________________________________________________________________________



      //reconstructed scattered electron
      if( !( mom_bin_scat_e_RC < 0) )
      {
        h_eta_scat_ele_RC_eCAL[mom_bin_scat_e_RC]->Fill(scat_e_mom_RC.Eta());

        if( mc_pdg_array[sim_id_scat_e_RC] == 11 ) h_scat_ele_purity_eCAL[mom_bin_scat_e_RC]->Fill(1.5); //RC scattered e candidate is MC electron
        else h_scat_ele_purity_eCAL[mom_bin_scat_e_RC]->Fill(0.5); //RC scattered e candidate is not MC electron
        
        //E/p cut
        if( fabs( 1. - maxEnergy/scat_e_mom_RC.Mag() ) < 0.1 )
        {
          h_eta_scat_ele_RC_eCAL_E_over_p[mom_bin_scat_e_RC]->Fill(scat_e_mom_RC.Eta());
          
          if( mc_pdg_array[sim_id_scat_e_RC] == 11 ) h_scat_ele_purity_eCAL_E_over_p[mom_bin_scat_e_RC]->Fill(1.5); //RC scattered e candidate is MC electron
          else h_scat_ele_purity_eCAL_E_over_p[mom_bin_scat_e_RC]->Fill(0.5); //RC scattered e candidate is not MC electron
        
        }

      }
      
      if( eta_bin_scat_e_RC != -1 )
      {      
        h_p_scat_ele_RC_eCAL[eta_bin_scat_e_RC]->Fill(scat_e_mom_RC.Mag());
        h_p_scat_ele_RC_eCAL[nEtaBins]->Fill(scat_e_mom_RC.Mag());
        
        if( mc_pdg_array[sim_id_scat_e_RC] == 11 ) 
        {
          h_scat_ele_purity_p_eCAL[eta_bin_scat_e_RC]->Fill(1.5); //RC scattered e candidate is MC electron
          h_scat_ele_purity_p_eCAL[nEtaBins]->Fill(1.5); //full pfRICH acceptance        
        
        }    
        else 
        {
          h_scat_ele_purity_p_eCAL[eta_bin_scat_e_RC]->Fill(0.5); //RC scattered e candidate is not MC electron
          h_scat_ele_purity_p_eCAL[nEtaBins]->Fill(0.5); //full pfRICH acceptance        
        }
        
        
        
        //E/p cut
        if( fabs( 1. - maxEnergy/scat_e_mom_RC.Mag() ) < 0.1 )
        {
          h_p_scat_ele_RC_eCAL_E_over_p[eta_bin_scat_e_RC]->Fill(scat_e_mom_RC.Mag());
          h_p_scat_ele_RC_eCAL_E_over_p[nEtaBins]->Fill(scat_e_mom_RC.Mag());
          
          if( mc_pdg_array[sim_id_scat_e_RC] == 11 ) 
          {
            h_scat_ele_purity_p_eCAL_E_over_p[eta_bin_scat_e_RC]->Fill(1.5); //RC scattered e candidate is MC electron
            h_scat_ele_purity_p_eCAL_E_over_p[nEtaBins]->Fill(1.5); //full pfRICH acceptance        
          
          }    
          else 
          {
            h_scat_ele_purity_p_eCAL_E_over_p[eta_bin_scat_e_RC]->Fill(0.5); //RC scattered e candidate is not MC electron
            h_scat_ele_purity_p_eCAL_E_over_p[nEtaBins]->Fill(0.5); //full pfRICH acceptance        
          } 
        }        
      
      }
      
    }//end if good leading cluster energy
    //___________________________________________________________________________________________________________________________________________________________

    //cout<<"RC - MC match start"<<endl;
    //MC -> RC matching
  	//find corresponding recID for scattered e and pi-
  	int scat_e_recID = -1;
  	vector<int> pi_recID_vect;
  	
  	for(unsigned int ch_trk_assoc_i = 0; ch_trk_assoc_i < sim_id.GetSize(); ch_trk_assoc_i++)
  	{
  	  if( mc_elect_index == sim_id[ch_trk_assoc_i] ) scat_e_recID = rec_id[ch_trk_assoc_i];
  	  
  	  for(unsigned int mc_pi_index = 0; mc_pi_index < pi_simID_vect.size(); mc_pi_index++)
  	  {
  	    if( pi_simID_vect.at(mc_pi_index) == sim_id[ch_trk_assoc_i] ) 
  	    {
  	      pi_recID_vect.push_back(rec_id[ch_trk_assoc_i]);
  	      
  	      break; //stop when matching RC track is found
  	    }      	    
  	  
  	  } 	  
  	    	
  	}

    
    //PID using pfRICH and MC -> RC tracks
    
    TVector3 scat_e_mom_RC_pfRICH(0,0,0);
    double maxP_pfRICH = -1;//to store leading momentum

    int scat_e_simID_pfRICH = -1; //simID of scattered e
    //cout<<"RC track loop start"<<endl;
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
      
          
      //MC -> RC matching
      if(eta_bin_RC != -1)
      {
        if( iChTrack == scat_e_recID ) 
        {
          h_p_scat_ele_MC_RC[eta_bin_RC]->Fill(rc_mom.Mag());
          h_p_scat_ele_MC_RC[nEtaBins]->Fill(rc_mom.Mag());      
        }
        
        for( unsigned int pi_recID_i = 0; pi_recID_i < pi_recID_vect.size(); pi_recID_i++ )
        {
          if( pi_recID_vect.at(pi_recID_i) == iChTrack )
          {
            h_p_pi_minus_MC_RC[eta_bin_RC]->Fill(rc_mom.Mag());
            h_p_pi_minus_MC_RC[nEtaBins]->Fill(rc_mom.Mag());

            //suppression factors for eCAL
            //default value is 1 - i.e. no suppression
            double eCAL_suppress_85 = 1.;
            double eCAL_suppress_95 = 1.;

            //eta acceptance of eCAL
            if( rc_mom.Eta() > -3.14 && rc_mom.Eta() < -1.87 )
            {
              eCAL_suppress_85 = g_pi_false_rate_85->Eval(rc_mom.Mag());
              eCAL_suppress_95 = g_pi_false_rate_95->Eval(rc_mom.Mag());
              
              h_p_pi_minus_MC_RC_eCAL_85[eta_bin_RC]->Fill(rc_mom.Mag(), eCAL_suppress_85);
              h_p_pi_minus_MC_RC_eCAL_95[eta_bin_RC]->Fill(rc_mom.Mag(), eCAL_suppress_95);
              
              h_p_pi_minus_MC_RC_eCAL_85[nEtaBins]->Fill(rc_mom.Mag(), eCAL_suppress_85);
              h_p_pi_minus_MC_RC_eCAL_95[nEtaBins]->Fill(rc_mom.Mag(), eCAL_suppress_95);
            }
            
            
            //change PID to MC
            double pfRICH_pi_prob;
            
            //use e/pi PID table for p < 5 GeV/c
            if(rc_mom.Mag() < 5)
            {
              pfRICH_pi_prob = getPIDprob_pfRICH_MC(rc_mom, 1, 1);        
            }
            else
            {
              pfRICH_pi_prob = getPIDprob_pfRICH_MC(rc_mom, 0, 0);
            }
            
            double noPID_pi_prob = 1. - pfRICH_pi_prob;


            if(  pfRICH_pi_prob < 0.9999 && pfRICH_pi_prob > 0)
            {
              //cout<<pfRICH_pi_prob<<endl;

              h_p_pi_minus_MC_RC_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag(), noPID_pi_prob);

              h_p_pi_minus_MC_RC_eCAL_85_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag(), eCAL_suppress_85*noPID_pi_prob);
              h_p_pi_minus_MC_RC_eCAL_95_pfRICH[eta_bin_RC]->Fill(rc_mom.Mag(), eCAL_suppress_95*noPID_pi_prob);
              
              //pfRICH eta acceptance
              if(rc_mom.Eta() > -3.8 && rc_mom.Eta() < -1.5)
              {
                h_p_pi_minus_MC_RC_pfRICH[nEtaBins]->Fill(rc_mom.Mag(), noPID_pi_prob);

                h_p_pi_minus_MC_RC_eCAL_85_pfRICH[nEtaBins]->Fill(rc_mom.Mag(), eCAL_suppress_85*noPID_pi_prob);
                h_p_pi_minus_MC_RC_eCAL_95_pfRICH[nEtaBins]->Fill(rc_mom.Mag(), eCAL_suppress_95*noPID_pi_prob);
              
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
      

      
    }//end for over RC particles
    
    //___________________________________________________________________________________________________________________________________________________________


    //analyze selected leading momentum tracks
    //find bins for reconstructed scattered e and charged particles
    double y_inelPar_e_RC_pfRICH = getInelParamElectron_2(p_energy, e_energy, scat_e_mom_RC_pfRICH);

    double Q2_electron_RC_pfRICH = getQ2elec( e_energy, scat_e_mom_RC_pfRICH);
    

    int mom_bin_scat_e_RC_pfRICH = -1;

    for(int j = 0; j < nMomBins; j++) //loop over pT bins
    {
      if(scat_e_mom_RC_pfRICH.Mag() > mom_bins[j] && scat_e_mom_RC_pfRICH.Mag() <= mom_bins[j+1])
      {
        mom_bin_scat_e_RC_pfRICH = j;
      }

    }
    
    
    int eta_bin_scat_e_RC_pfRICH = -1;

      for(int j = 0; j < nEtaBins; j++) //loop over pT bins
      {
        if(scat_e_mom_RC_pfRICH.Eta() > eta_bins[j] && scat_e_mom_RC_pfRICH.Eta() <= eta_bins[j+1])
        {
          eta_bin_scat_e_RC_pfRICH = j;
        }

      }
    
    
    double MC_PDG_ID_pfRICH = -99;

    if( scat_e_simID_pfRICH >= 0 && mc_generatorStatus_array[scat_e_simID_pfRICH] == 1 )
    {
      MC_PDG_ID_pfRICH = mc_pdg_array[scat_e_simID_pfRICH];
    
      if( !( mom_bin_scat_e_RC_pfRICH < 0) )
      {
        h_eta_scat_ele_RC_lead_p[mom_bin_scat_e_RC_pfRICH]->Fill(scat_e_mom_RC_pfRICH.Eta());
        
        //evaluate ekectron purity only in pfRICH acceptance
        if(scat_e_mom_RC_pfRICH.Eta() > -3.8 && scat_e_mom_RC_pfRICH.Eta() < -1.5)
        {
          if( MC_PDG_ID_pfRICH == 11 ) h_scat_ele_purity_lead_p[mom_bin_scat_e_RC_pfRICH]->Fill(1.5);
          else h_scat_ele_purity_lead_p[mom_bin_scat_e_RC_pfRICH]->Fill(0.5);      
        }
      }
      
      if( eta_bin_scat_e_RC_pfRICH != -1 )
      {        
        
        if( scat_e_mom_RC_pfRICH.Eta() > -3.8 && scat_e_mom_RC_pfRICH.Eta() < -1.5 )
        {
          h_p_scat_ele_RC_lead_p[eta_bin_scat_e_RC_pfRICH]->Fill(scat_e_mom_RC_pfRICH.Mag());
          h_p_scat_ele_RC_lead_p[nEtaBins]->Fill(scat_e_mom_RC_pfRICH.Mag());
          
          if( MC_PDG_ID_pfRICH == 11 ) 
          {
            h_scat_ele_purity_p_lead_p[eta_bin_scat_e_RC_pfRICH]->Fill(1.5);
            h_scat_ele_purity_p_lead_p[nEtaBins]->Fill(1.5);          
          
          }         
          else
          {
            h_scat_ele_purity_p_lead_p[eta_bin_scat_e_RC_pfRICH]->Fill(0.5);
            h_scat_ele_purity_p_lead_p[nEtaBins]->Fill(0.5);
          }
    
        }
      
      }
      
      
      
    }


    //generate random number for PID in pfRICH
    TRandom3 *rndm_pfRICH = new TRandom3();
    double rndm_PID_pfRICH = rndm_pfRICH->Rndm();

    int is_pfRICH_electron = 0; //flag if track is identified by pfRICH as electron
    int is_pfRICH_hadron = 0; //flag if track is identified by pfRICH as hadron
    

    //get pfRICH PID probablilities
    //also rejecting tracks with invalid pfRICH info


    //get pfRICH PID matrixes for mis-PID evaluation

    //pi/K/p matrix
    // 0 - pi identified as pi, 1 - pi identified as K, 2 - pi identified as p
    // 2 - K identified as pi, 3 - K identified as K, 4 - K identified as p
    // 5 - p identified as pi, 6 - p identified as K, 7 - p identified as p
    double PID_mtx[9];

    int pfRICH_good = getPIDprob_pfRICH_mtx(scat_e_mom_RC_pfRICH, PID_mtx, 0); // pi/K/p matrix


    //e/pi matrix
    //PID mtx has to have 9 elemets to be compatible with getPIDprob_pfRICH_mtx
    //for e/pi, only 4 elemets will be filled:
    // 0 - e identified as e, 1 - e identified as pi
    // 2 - pi identified as e, 3 - pi identified as pi
    double PID_mtx_e_pi[9];

    int pfRICH_good_e_pi = getPIDprob_pfRICH_mtx(scat_e_mom_RC_pfRICH, PID_mtx_e_pi, 1); // e/pi PID matrix



    //selection of electrons
    //direct for e
    //contamination by pi from e/pi PID table
    //veto on K and p

    //pi-
    if( MC_PDG_ID_pfRICH == -211 )
    {
      //PID of pi - not used now
      double pi_PID_prob_pfRICH = getPIDprob_pfRICH_MC(scat_e_mom_RC_pfRICH, 0, 0);

      if( pi_PID_prob_pfRICH > rndm_PID_pfRICH /*&& pi_PID_prob_pfRICH < 0.9999 */|| pi_PID_prob_pfRICH == 0) is_pfRICH_hadron = 1;

/*
      //pi mis-identified as e
      if(pfRICH_good_e_pi == 1)
      {
        //check if pi was mis-identified as e
        if( PID_mtx_e_pi[2] > rndm_PID_pfRICH) is_pfRICH_electron = 1;
      }
*/
    }

    //K-
    else if( MC_PDG_ID_pfRICH == -321 )
    {
      double K_PID_prob_pfRICH = getPIDprob_pfRICH_MC(scat_e_mom_RC_pfRICH, 1, 0);

      if( K_PID_prob_pfRICH > rndm_PID_pfRICH /*&& K_PID_prob_pfRICH < 0.9999 */ || K_PID_prob_pfRICH == 0)
      {
        is_pfRICH_hadron = 1; //for hadron veto in e selection
      }

    }

    //p-bar
    else if( MC_PDG_ID_pfRICH == -2212 )
    {
      double p_PID_prob_pfRICH = getPIDprob_pfRICH_MC(scat_e_mom_RC_pfRICH, 2, 0);

      if( p_PID_prob_pfRICH > rndm_PID_pfRICH /*&& p_PID_prob_pfRICH < 0.9999 */ || p_PID_prob_pfRICH == 0) is_pfRICH_hadron = 1;

    }

    //e-
    else if( MC_PDG_ID_pfRICH == 11 )
    {
      //if electron then accept

      double e_PID_prob_pfRICH = getPIDprob_pfRICH_MC(scat_e_mom_RC_pfRICH, 0, 1);

      if( e_PID_prob_pfRICH > rndm_PID_pfRICH )
      {
        is_pfRICH_electron = 1;
      }



    }
    else
    {
      //reject everything else
      //for testing
      is_pfRICH_hadron = 1;
    }
    
    

    //continue for good scattered electron candidates
    //accept only particles in pfRICH acceptance
    if( MC_PDG_ID_pfRICH != -99 && (is_pfRICH_electron == 1 || is_pfRICH_hadron != 1) &&  scat_e_mom_RC_pfRICH.Eta() > -3.8 && scat_e_mom_RC_pfRICH.Eta() < -1.5 )
    //if( MC_PDG_ID_pfRICH != -99 && is_pfRICH_electron == 1 &&  scat_e_mom_RC_pfRICH.Eta() > -3.8 && scat_e_mom_RC_pfRICH.Eta() < -1.5 )
    {      

      h_momentum_RC_pfRICH->Fill(scat_e_mom_RC_pfRICH.Mag());

      h_Q2_RC_pfRICH->Fill(Q2_electron_RC_pfRICH);
      h_y_inelPar_RC_pfRICH->Fill(y_inelPar_e_RC_pfRICH);
      

      if( !( mom_bin_scat_e_RC_pfRICH < 0) )
      {

        h_eta_scat_ele_RC_pfRICH[mom_bin_scat_e_RC_pfRICH]->Fill(scat_e_mom_RC_pfRICH.Eta());

        if( MC_PDG_ID_pfRICH == 11 ) h_scat_ele_purity_pfRICH[mom_bin_scat_e_RC_pfRICH]->Fill(1.5);
        else h_scat_ele_purity_pfRICH[mom_bin_scat_e_RC_pfRICH]->Fill(0.5);
      }
      
      
      if( eta_bin_scat_e_RC_pfRICH != -1 )
      {
        h_p_scat_ele_RC_pfRICH[eta_bin_scat_e_RC_pfRICH]->Fill(scat_e_mom_RC_pfRICH.Mag());
        h_p_scat_ele_RC_pfRICH[nEtaBins]->Fill(scat_e_mom_RC_pfRICH.Mag());
        
        if( MC_PDG_ID_pfRICH == 11 )
        {
          h_scat_ele_purity_p_pfRICH[eta_bin_scat_e_RC_pfRICH]->Fill(1.5);
          h_scat_ele_purity_p_pfRICH[nEtaBins]->Fill(1.5);        
        }    
        else
        {
          h_scat_ele_purity_p_pfRICH[eta_bin_scat_e_RC_pfRICH]->Fill(0.5);
          h_scat_ele_purity_p_pfRICH[nEtaBins]->Fill(0.5);
        }       
      
      }

    }//end if e selection



  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}

