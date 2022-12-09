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
  //pi eCALion
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
  
  // MC particle pz array for each MC particle
  TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
  TTreeReaderArray<int> mc_generatorStatus_array = {tree_reader, "MCParticles.generatorStatus"};


  // Reconstructed particles pz array for each reconstructed particle
  TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<int> reco_PDG = {tree_reader, "ReconstructedChargedParticles.PDG"};

  TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
  TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};
  

  //defining output file and histos.
  TString output_name_dir = outname;
	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

  
  TH1D* h_energy_MC = new TH1D("h_energy_MC","E_{MC} (GeV)",100,0,20);
  TH1D* h_energy_zoom_MC = new TH1D("h_energy_zoom_MC","E_{MC} (GeV)",20,e_energy-2,e_energy+2);
  
  TH1D* h_momentum_MC = new TH1D("h_momentum_MC","p_{MC} (GeV/c)",100,0,20);
  
  TH1D* h_Q2_MC = new TH1D("h_Q2_MC",";Q^{2}",100,0,20);
  TH1D* h_Q2_zoom_MC = new TH1D("h_Q2_zoom_MC",";Q^{2}",100,0,4);
  
  TH1D *h_y_inelPar_MC = new TH1D("h_y_inelPar_MC", "h_y_inelPar_MC", 100, 0, 1);
  TH1D *h_y_inelPar_zoom_MC = new TH1D("h_y_inelPar_zoom_MC", "h_y_inelPar_zoom_MC", 100, 0, 0.01);
  
  //eta distributions in multiple Q^2 and inelasticity bins
  TH1D *h_eta_scat_ele[nMomBins+1][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_ele[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_plus[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_K_plus[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_proton[nMomBins+1][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_positron[nMomBins+1][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_85[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95[nMomBins+1][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_pi_minus_eCAL_85_pfRICH[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_pi_minus_eCAL_95_pfRICH[nMomBins+1][nQ2bins][nyInelParBins];
  
  TH1D *h_eta_K_minus[nMomBins+1][nQ2bins][nyInelParBins];
  TH1D *h_eta_anti_proton[nMomBins+1][nQ2bins][nyInelParBins];
  
  for(unsigned int mom_bin = 0; mom_bin < nMomBins+1; mom_bin++)
  {
    for(unsigned int Q2bin = 0; Q2bin < nQ2bins; Q2bin++)
    {
      for(unsigned int y_bin = 0; y_bin < nyInelParBins; y_bin++)
      {
        h_eta_scat_ele[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_scat_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_ele[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_ele_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_pi_plus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_K_plus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_K_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_plus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_proton[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        
        h_eta_positron[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_positron_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_positron_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        
        h_eta_pi_minus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_pi_minus_eCAL_85[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_85_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_85_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_pi_minus_eCAL_95[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_95_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_95_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        
        h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_85_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_pi_minus_eCAL_95_pfRICH_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
        
        h_eta_K_minus[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_K_minus_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);        
        h_eta_anti_proton[mom_bin][Q2bin][y_bin] = new TH1D(Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), Form("h_eta_anti_proton_mom_%i_Q2_%i_y_%i" , mom_bin, Q2bin, y_bin), 200, -4, 4);
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

  	//MCParticles
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
    
    if(Q2_bin < 0 || y_bin < 0 || mom_bin_scat_e < 0) continue;
    
    
    h_eta_scat_ele[mom_bin_scat_e][Q2_bin][y_bin]->Fill(scatMC.Eta());
    
     
    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{
  	  if( mc_generatorStatus_array[imc] != 1 ) continue; 
  	
  		TVector3 mc_mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
  		TLorentzVector mc_4mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc], mc_mass_array[imc]);
  		
  		//determine momentum bin of particles in event
  		int mom_bin = -1;

      for(int j = 0; j < nMomBins; j++) //loop over pT bins
      {
        if(mc_mom.Mag() > mom_bins[j] && mc_mom.Mag() <= mom_bins[j+1])
        {
          mom_bin = j;
        }
        
      }
      
      if(mom_bin < 0) continue;
  		
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
  		
  		//pi+
  		if(mc_pdg_array[imc] == 211)
  		{
  		   h_eta_pi_plus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//pi-
  		if(mc_pdg_array[imc] == -211)
  		{
  		   h_eta_pi_minus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());
  		   
  		   h_eta_pi_minus_eCAL_85[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), g_pi_false_rate_85->Eval(mc_mom.Mag()));
  		   h_eta_pi_minus_eCAL_95[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), g_pi_false_rate_95->Eval(mc_mom.Mag()));
  		   
  		   double pfRICH_pi_prob = getPIDprob_pfRICH_single(mc_4mom);
  		   double noPID_pi_prob = 1 - pfRICH_pi_prob;
  		   
  		   h_eta_pi_minus_eCAL_85_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), g_pi_false_rate_85->Eval(mc_mom.Mag())*noPID_pi_prob);
  		   h_eta_pi_minus_eCAL_95_pfRICH[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta(), g_pi_false_rate_95->Eval(mc_mom.Mag())*noPID_pi_prob);  

  		}
  		
  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		   h_eta_K_plus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//K-
  		if(mc_pdg_array[imc] == -321)
  		{
  		   h_eta_K_minus[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		
  		//proton
  		if(mc_pdg_array[imc] == 2212)
  		{
  		   h_eta_proton[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		//anti-proton
  		if(mc_pdg_array[imc] == -2212)
  		{
  		   h_eta_anti_proton[mom_bin][Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		
  	}//end secon particle loop
   

  }//end while over TTree entries

	output->Write();
	output->Close();

	return 0;
}
