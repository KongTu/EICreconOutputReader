#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include "pleaseIncludeMe.h"

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
  
  const int nMomBins = 5;
  float const mom_bins[nMomBins+1] = { 0,1,3,7,10, 18 };
  
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
  TTreeReaderArray<int> reco_cahrge = {tree_reader, "ReconstructedChargedParticles.charge"};
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
  TH1D* h_energy_MC = new TH1D("h_energy_MC","E_{MC} (GeV)",100,0,20);
  TH1D* h_energy_zoom_MC = new TH1D("h_energy_zoom_MC","E_{MC} (GeV)",20,e_energy-2,e_energy+2);
  
  TH1D* h_momentum_MC = new TH1D("h_momentum_MC","p_{MC} (GeV/c)",100,0,20);
  
  TH1D* h_Q2_MC = new TH1D("h_Q2_MC",";Q^{2}",100,0,20);
  TH1D* h_Q2_zoom_MC = new TH1D("h_Q2_zoom_MC",";Q^{2}",100,0,4);
  
  TH1D *h_y_inelPar_MC = new TH1D("h_y_inelPar_MC", "h_y_inelPar_MC", 100, 0, 1);
  TH1D *h_y_inelPar_zoom_MC = new TH1D("h_y_inelPar_zoom_MC", "h_y_inelPar_zoom_MC", 100, 0, 0.01);
  
  //eta distributions in multiple Q^2 and inelasticity bins
  TH1D *h_eta_scat_ele[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_ele[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_plus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_K_plus[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_proton[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_positron[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus_eCAL_85[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus_eCAL_85_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_K_minus[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_K_minus_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_anti_proton[nMomBins][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_anti_proton_pfRICH[nMomBins][nQ2bins][nyInelParBins];
  
  
  //reco histograms
  TH1D* h_energy_RC = new TH1D("h_energy_RC","E_{RC} (GeV)",100,0,20);
  
  TH1D* h_momentum_RC = new TH1D("h_momentum_RC","p_{RC} (GeV/c)",100,0,20);
  
  TH1D *h_E_over_p_RC = new TH1D("h_E_over_p_RC", "h_E_over_p_RC", 120, 0, 1.2);
  
  TH1D* h_Q2_RC = new TH1D("h_Q2_RC",";Q^{2}",100,0,20);
  
  TH1D *h_y_inelPar_RC = new TH1D("h_y_inelPar_RC", "h_y_inelPar_RC", 100, 0, 1);
  
  //reconstructed scattered electron
  TH1D *h_eta_scat_ele_RC[nMomBins][nQ2bins][nyInelParBins];
  
  //background - negative charged particles
  TH1D *h_eta_neg_ch_part_RC[nMomBins][nQ2bins][nyInelParBins];
  
  //negative charged particles after E/p cut
  //TH1D *h_eta_neg_ch_part_RC_eCAL_cut[nMomBins][nQ2bins][nyInelParBins];
  
  //negative charged par after pfRCIH veto
  //TH1D *h_eta_neg_ch_part_RC_pFRICH_cut[nMomBins][nQ2bins][nyInelParBins];
  
  //negative charged particles after eCAL+pfRCIH veto
  //TH1D *h_eta_neg_ch_part_RC_eCAL_pfRICH_cut[nMomBins][nQ2bins][nyInelParBins];
  
  //MC-RC matching
  TH1D *h_scat_ele_MC_RC_match[nMomBins][nQ2bins][nyInelParBins];

  //QA histograms
  TH1D *h_PID_pfRICH_mtx[nMomBins];
  TH1D *h_PID_pfRICH_mtx_nFill[nMomBins];
  
  TH1D *h_PID_pfRICH_pi_mtx[nMomBins];
  TH1D *h_PID_pfRICH_pi_mtx_nFill[nMomBins];
  
  TH1D *h_PID_pfRICH_K_mtx[nMomBins];
  TH1D *h_PID_pfRICH_K_mtx_nFill[nMomBins];
  
  TH1D *h_PID_pfRICH_p_mtx[nMomBins];
  TH1D *h_PID_pfRICH_p_mtx_nFill[nMomBins];
  
  
  for(unsigned int mom_bin = 0; mom_bin < nMomBins; mom_bin++)
  {
    //pfRCIH QA histograms
    h_PID_pfRICH_mtx[mom_bin] = new TH1D(Form("h_PID_pfRICH_mtx_mom_%i", mom_bin), Form("h_PID_pfRICH_mtx_mom_%i", mom_bin), 9, 0, 9);
    h_PID_pfRICH_mtx_nFill[mom_bin] = new TH1D(Form("h_PID_pfRICH_mtx_nFill_mom_%i", mom_bin), Form("h_PID_pfRICH_mtx_nFill_mom_%i", mom_bin), 9, 0, 9);
    
    h_PID_pfRICH_pi_mtx[mom_bin] = new TH1D(Form("h_PID_pfRICH_pi_mtx_mom_%i", mom_bin), Form("h_PID_pfRICH_pi_mtx_mom_%i", mom_bin), 9, 0, 9);
    h_PID_pfRICH_pi_mtx_nFill[mom_bin] = new TH1D(Form("h_PID_pfRICH_pi_mtx_nFill_mom_%i", mom_bin), Form("h_PID_pfRICH_pi_mtx_nFill_mom_%i", mom_bin), 9, 0, 9);
    
    h_PID_pfRICH_K_mtx[mom_bin] = new TH1D(Form("h_PID_pfRICH_K_mtx_mom_%i", mom_bin), Form("h_PID_pfRICH_K_mtx_mom_%i", mom_bin), 9, 0, 9);
    h_PID_pfRICH_K_mtx_nFill[mom_bin] = new TH1D(Form("h_PID_pfRICH_K_mtx_nFill_mom_%i", mom_bin), Form("h_PID_pfRICH_K_mtx_nFill_mom_%i", mom_bin), 9, 0, 9);
    
    h_PID_pfRICH_p_mtx[mom_bin] = new TH1D(Form("h_PID_pfRICH_p_mtx_mom_%i", mom_bin), Form("h_PID_pfRICH_p_mtx_mom_%i", mom_bin), 9, 0, 9);
    h_PID_pfRICH_p_mtx_nFill[mom_bin] = new TH1D(Form("h_PID_pfRICH_p_mtx_nFill_mom_%i", mom_bin), Form("h_PID_pfRICH_p_mtx_nFill_mom_%i", mom_bin), 9, 0, 9);
    //___________________________________________________________________________________________________________________________________________________________
    
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {
        //MC histograms
        h_eta_scat_ele[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_ele[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_K_plus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_K_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_proton[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_positron[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_positron_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_positron_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_pi_minus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_85_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_85_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_95_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_95_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_pi_minus_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_K_minus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_K_minus_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_K_minus_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
          
        h_eta_anti_proton[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_eta_anti_proton_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_anti_proton_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_anti_proton_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        //___________________________________________________________________________________________________________________________________________________________
        
        //reconstructed particle histograms
        h_eta_scat_ele_RC[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_scat_ele_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_scat_ele_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
  
        //background - negative charged particles
        h_eta_neg_ch_part_RC[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_neg_ch_part_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_neg_ch_part_RC_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        //negative charged particles after E/p cut
        //h_eta_neg_ch_part_RC_eCAL_cut[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_neg_ch_part_RC_eCAL_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_neg_ch_part_RC_eCAL_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        //negative charged par after pfRCIH veto
        //h_eta_neg_ch_part_RC_pFRICH_cut[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_neg_ch_part_RC_pFRICH_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_neg_ch_part_RC_pFRICH_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        //negative charged particles after eCAL+pfRCIH veto
        //h_eta_neg_ch_part_RC_eCAL_pfRICH_cut[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_neg_ch_part_RC_eCAL_pfRICH_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_neg_ch_part_RC_eCAL_pfRICH_cut_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 100, -4, 0);
        
        h_scat_ele_MC_RC_match[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_scat_ele_MC_RC_match_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_scat_ele_MC_RC_match_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 2, 0, 2);
        
      }
    
    }
  }
  
  //TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
  //TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);  
  //TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
  //TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC}",100,0,20,1000,-1,1);
  
  

	tree_reader.SetEntriesRange(0, myChain->GetEntries());
		
  while (tree_reader.Next())
  {

  	//MC analysis
  	
    //finding the scattering electron
  	TLorentzVector scatMC(0,0,0,0);
  	int mc_elect_index = -1;
  	double maxP = -99.;
  	
  	
  	for(int imc = 0; imc < mc_px_array.GetSize(); imc++)
  	{
  	  if( mc_generatorStatus_array[imc] != 1 ) continue;
  	
  		TVector3 mctrk(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);	
  		
  		if(mc_pdg_array[imc] == 11 && mctrk.Mag() > maxP)
  		//if(mc_pdg_array[imc] == 11 && mctrk.Perp() > maxP)
  		{
  			maxP = mctrk.Mag();
  			//maxP = mctrk.Perp();
  			mc_elect_index = imc;
  			scatMC.SetVectM(mctrk,mc_mass_array[imc]);  			
  		}
  	} 	
  	 	
  	//charge proton and electron energy loading - for now hard coded
    double y_inelPar_e = getInelParamElectron_2(p_energy, e_energy, scatMC.Vect() );
    
    double Q2_electron = getQ2elec( e_energy, scatMC.Vect());
    
    if(Q2_electron< 1. || Q2_electron > 20.) continue;
    if(y_inelPar_e < 0.01 ||y_inelPar_e > 0.95) continue;
    
    //fill truth scattered electron energy
    if(Q2_electron < 10.)
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
    
    //if(Q2_bin < 0 || y_bin < 0 || mom_bin_scat_e < 0) continue;
    //__________________________________________________________________________________
    
    if( !(Q2_bin < 0 || y_bin < 0 || mom_bin_scat_e < 0) ) h_eta_scat_ele[mom_bin_scat_e][Q2_bin][y_bin]->Fill(scatMC.Eta());
    
    
    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{
  	
  	  if( mc_generatorStatus_array[imc] != 1 ) continue; 
  	
  		TVector3 mc_mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
  		TLorentzVector mc_4mom(0,0,0,0);
  		mc_4mom.SetVectM(mc_mom, mc_mass_array[imc]);
  		
  		//determine momentum bin of particles in event
  		int mom_bin = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(mc_mom.Mag() > mom_bins[j] && mc_mom.Mag() <= mom_bins[j+1])
        {
          mom_bin = j;
        }
        
      }
      
      if(Q2_bin < 0 || y_bin < 0 || mom_bin < 0) continue;
      
      //_____________________________________________________________________________
      
  		
  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 &&  mc_mom.Pt() < maxP)
  		if(mc_pdg_array[imc] == 11 && mc_mom.Mag() < maxP)
  		{
  			h_eta_ele[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  			
  		}
  		//positrons
  		if(mc_pdg_array[imc] == -11 )
  		//if(mc_pdg_array[imc] == 11 && mc_mom.Mag() < maxP)
  		{
  			h_eta_positron[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  			
  		}
  		//_____________________________________________________________________________
  		
  		//pi+
  		if(mc_pdg_array[imc] == 211)
  		{
  		   h_eta_pi_plus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//pi-
  		if(mc_pdg_array[imc] == -211)
  		{
        h_eta_pi_minus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());

        //suppression factors for eCAL
        //default value is 1 - i.e. no suppression
        double eCAL_suppress_85 = 1.;
        double eCAL_suppress_95 = 1.;

        //eta acceptance of eCAL
        if( mc_mom.Eta() > -3.14 && mc_mom.Eta() < -1.87 )
        {
          eCAL_suppress_85 = g_pi_false_rate_85->Eval(mc_mom.Mag());
          eCAL_suppress_95 = g_pi_false_rate_95->Eval(mc_mom.Mag());

          h_eta_pi_minus_eCAL_85[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), eCAL_suppress_85);
          h_eta_pi_minus_eCAL_95[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), eCAL_suppress_95); 
        }
        
        
        //change PID to MC
        double pfRICH_pi_prob = getPIDprob_pfRICH_MC(mc_4mom, 0);    
        double noPID_pi_prob = 1. - pfRICH_pi_prob; 
        
        
        if(  pfRICH_pi_prob < 0.9999 && pfRICH_pi_prob > 0)
        {
          //cout<<pfRICH_pi_prob<<endl;
        
          h_eta_pi_minus_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), noPID_pi_prob);
        
          h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), eCAL_suppress_85*noPID_pi_prob);
          h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), eCAL_suppress_95*noPID_pi_prob);
        }
  		   
  		     

  		}
  		//_____________________________________________________________________________
  		
  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		   h_eta_K_plus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//K-
  		if(mc_pdg_array[imc] == -321)
  		{
  		  h_eta_K_minus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
  		  
  		  double pfRICH_K_prob = getPIDprob_pfRICH_MC(mc_4mom, 1);
        double noPID_K_prob = 1 - pfRICH_K_prob;
                       
        
        //if( mc_mom.Mag() > 4 && mc_mom.Mag() < 5 ) cout<<PID_mtx[4]<<endl;

        if( pfRICH_K_prob > 0 && pfRICH_K_prob < 0.9999)
        {
          //if(mom_bin == 2) cout<<PID_K_mtx[0]<<" "<<PID_K_mtx[4]<<" "<<PID_K_mtx[8]<<endl;
          //if(mom_bin == 2) cout<<PID_K_mtx[3]<<" "<<PID_K_mtx[4]<<" "<<PID_K_mtx[5]<<endl;
          //if( mc_mom.Mag() > 4 && mc_mom.Mag() < 5 ) cout<<pfRICH_K_prob<<endl;
                      
          h_eta_K_minus_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), noPID_K_prob);
        }
  		   		
  		}
  		//_____________________________________________________________________________
  		
  		//proton
  		if(mc_pdg_array[imc] == 2212)
  		{
  		   h_eta_proton[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//anti-proton
  		if(mc_pdg_array[imc] == -2212)
  		{
  		  h_eta_anti_proton[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
  		   
  		  double pfRICH_anti_proton_prob = getPIDprob_pfRICH_MC(mc_4mom, 2);
        double noPID_anti_proton_prob = 1 - pfRICH_anti_proton_prob;
        
     
        if( pfRICH_anti_proton_prob > 0 && pfRICH_anti_proton_prob < 0.9999)
        {
                 
          h_eta_anti_proton_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), noPID_anti_proton_prob);
        } 		
  		}
  		//___________________________________________________________________________________________________________________________________________________________
  		
  		
  		//general pfRICH info QA
  		double PID_mtx[9];
        
      //all particles
      int pfRICH_good = getPIDprob_pfRICH_QA(mc_4mom, PID_mtx);
      
      if( pfRICH_good == 0 ) continue;
      
      for(unsigned int mtx_bin = 0; mtx_bin < 9; mtx_bin++)
      {        
        if( PID_mtx[mtx_bin] > 0 && PID_mtx[mtx_bin] < 0.9999 )
        {
          h_PID_pfRICH_mtx[mom_bin]->Fill(mtx_bin+0.5, PID_mtx[mtx_bin]);
          h_PID_pfRICH_mtx_nFill[mom_bin]->Fill(mtx_bin+0.5);
          
          //pi-
          if(mc_pdg_array[imc] == -211)
          {
            h_PID_pfRICH_pi_mtx[mom_bin]->Fill(mtx_bin+0.5, PID_mtx[mtx_bin]);
            h_PID_pfRICH_pi_mtx_nFill[mom_bin]->Fill(mtx_bin+0.5);          
          }
          
          //K-
          if(mc_pdg_array[imc] == -321)
          {
            h_PID_pfRICH_K_mtx[mom_bin]->Fill(mtx_bin+0.5, PID_mtx[mtx_bin]);
            h_PID_pfRICH_K_mtx_nFill[mom_bin]->Fill(mtx_bin+0.5);          
          }
          
          //anti-proton
          if(mc_pdg_array[imc] == -2212)
          {
            h_PID_pfRICH_p_mtx[mom_bin]->Fill(mtx_bin+0.5, PID_mtx[mtx_bin]);
            h_PID_pfRICH_p_mtx_nFill[mom_bin]->Fill(mtx_bin+0.5);          
          }
                
        }
        
         
      }
      
  		
  	}//end secon MC particle loop
  	//___________________________________________________________________________________________________________________________________________________________
  	
  	
  	//reconstructed charged particles analysis
  		
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
      if( sim_id[iChTrack] == sim_id_scat_e_RC /*&& reco_cahrge[iChTrack] == -1 */)
      {
        scat_e_mom_RC.SetXYZ( reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack] );
        scat_e_iTrack = iChTrack;
        
        break; //stop when corresponding ch. track is found        
      }        
    
    }
    
    if(scat_e_mom_RC.Mag() == 0)
    {
      cout<<"No good RC scattered e found! Continue with next event..."<<endl;
      continue;    
    }
    
    h_momentum_RC->Fill(scat_e_mom_RC.Mag());
        
    h_E_over_p_RC->Fill( maxEnergy/scat_e_mom_RC.Mag() );
    
    //analyze RC tracks only when E/p cut is passed
    if( fabs( 1. - maxEnergy/scat_e_mom_RC.Mag() ) < 0.1 )
    {
      double y_inelPar_e_RC = getInelParamElectron_2(p_energy, e_energy, scat_e_mom_RC);
      
      double Q2_electron_RC = getQ2elec( e_energy, scat_e_mom_RC);
      
      h_Q2_RC->Fill(Q2_electron_RC);
      h_y_inelPar_RC->Fill(y_inelPar_e_RC);
      
      
      //find bins for reconstructed scattered e and charged particles
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
      
      
      int mom_bin_scat_e_RC = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(scat_e_mom_RC.Mag() > mom_bins[j] && scat_e_mom_RC.Mag() <= mom_bins[j+1])
        {
          mom_bin_scat_e_RC = j;
        }
        
      }
      
      //if(Q2_bin_RC < 0 || y_bin_RC < 0 || mom_bin_scat_e_RC < 0) continue;
      //__________________________________________________________________________________
      
      
      
      //reconstructed scattered electron
      if( !(Q2_bin_RC < 0 || y_bin_RC < 0 || mom_bin_scat_e_RC < 0) )  
      {
        h_eta_scat_ele_RC[mom_bin_scat_e_RC][Q2_bin_RC][y_bin_RC]->Fill(scat_e_mom_RC.Eta());
        
        if( mc_pdg_array[sim_id_scat_e_RC] == 11 ) h_scat_ele_MC_RC_match[mom_bin_scat_e_RC][Q2_bin_RC][y_bin_RC]->Fill(1.5); //RC scattered e candidate is MC electron
        else h_scat_ele_MC_RC_match[mom_bin_scat_e_RC][Q2_bin_RC][y_bin_RC]->Fill(0.5); //RC scattered e candidate is not MC electron
      
      }
      
      
      
      //loop over RC charged particles
      for(unsigned int iChTrack = 0; iChTrack < reco_px_array.GetSize(); iChTrack++)
      {
        if( reco_cahrge[iChTrack]  > 0 || iChTrack == scat_e_iTrack) continue; //analyze negative charge particles only, skip scattered e identified earlier
        
        TVector3 rc_mom(reco_px_array[iChTrack], reco_py_array[iChTrack], reco_pz_array[iChTrack]);
        
        //TLorentzVector rc_4mom(0,0,0,0);
		    //rc_4mom.SetVectM(rc_mom, mc_mass_array[imc]);
		    
		    //determine momentum bin of particles in event
		    int mom_bin_RC = -1;

        for(int j = 0; j < nMomBins; j++) //loop over pT bins
        {
          if(rc_mom.Mag() > mom_bins[j] && rc_mom.Mag() <= mom_bins[j+1])
          {
            mom_bin_RC = j;
          }
          
        }
        
        if(Q2_bin_RC < 0 || y_bin_RC < 0 || mom_bin_RC < 0) continue;
        
        //_____________________________________________________________________________
        
      
        //background - negative charged particles
        h_eta_neg_ch_part_RC[mom_bin_RC][Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());
        
        //negative charged particles after E/p cut
 /*       
        //first find corresponding cluster in eCAL
        float track_eCAL_energy = -99.;
        
        for(unsigned int iCluster = 0; iCluster < em_energy_array.GetSize(); iCluster++)
        {
          if( sim_id[iChTrack] == em_sim_id_array[iCluster] )
          {
            track_eCAL_energy = em_energy_array[iCluster];
            
            break; //stop when corresponding cluster is found
          
          }
        
        }
        
        //check that matchig cluser was found
        //apply same E/p cut as for scattered E - histogram will contain background for scattered e
        if(track_eCAL_energy > 0 && fabs( 1. - track_eCAL_energy/rc_mom.Mag() ) < 0.1 )
        {
            h_eta_neg_ch_part_RC_eCAL_cut[mom_bin_RC][Q2_bin_RC][y_bin_RC]->Fill(rc_mom.Eta());
        }
 */       
        
        //negative charged par after pfRCIH veto
        //h_eta_neg_ch_part_RC_pFRICH_cut[mom_bin_RC][Q2_bin_RC][y_bin_RC];
        
        //negative charged particles after eCAL+pfRCIH veto
        //h_eta_neg_ch_part_RC_eCAL_pfRICH_cut[mom_bin_RC][Q2_bin_RC][y_bin_RC];
      
      
      }    
      
    
    
    }//end if scattered e E/p cut
  	
   

  }//end while over TTree entries



	output->Write();
	output->Close();

	return 0;
}
