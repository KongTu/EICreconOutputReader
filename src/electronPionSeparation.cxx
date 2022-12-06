#include "pleaseIncludeMe.h"

int electronPionSeparation(TString inname="input/input.root",TString outname="test")
{	
  const int nQ2bins = 2;
  float const Q2_bins[nQ2bins+1] = { 1., 2., 4. };
  
  const int nyInelParBins = 2;
  float const y_bins[nQ2bins+1] = { 0., 0.001, 0.002 };
	
	auto file = new TFile(inname);
	
	auto tree = (TTree *) file->Get("events");
  TTreeReader tree_reader(tree);       // !the tree reader
  
  // MC particle pz array for each MC particle
  TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};


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

  
  TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
  TH1D* h_energy_18_MC = new TH1D("h_energy_18_MC",";E_{MC} (GeV)",20,17,19);
  
  TH1D* h_Q2_MC = new TH1D("h_Q2_MC",";Q^{2}",100,0,20);
  TH1D* h_Q2_zoom_MC = new TH1D("h_Q2_zoom_MC",";Q^{2}",100,0,4);
  
  TH1D *h_y_inelPar_MC = new TH1D("h_y_inelPar_MC", "h_y_inelPar_MC", 100, 0, 1);
  TH1D *h_y_inelPar_zoom_MC = new TH1D("h_y_inelPar_zoom_MC", "h_y_inelPar_zoom_MC", 100, 0, 0.01);
  
  //eta distributions in multiple Q^2 and inelasticity bins
  TH1D *h_eta_scat_ele[nQ2bins+1][nyInelParBins+1];
  TH1D *h_eta_ele[nQ2bins+1][nyInelParBins+1];
  TH1D *h_eta_pi[nQ2bins+1][nyInelParBins+1];
  TH1D *h_eta_K[nQ2bins+1][nyInelParBins+1];
  TH1D *h_eta_proton[nQ2bins+1][nyInelParBins+1];
  
  for(unsigned int Q2bin = 0; Q2bin < nQ2bins+1; Q2bin++)
  {
    for(unsigned int y_bin = 0; y_bin < nyInelParBins+1; y_bin++)
    {
      h_eta_scat_ele[Q2bin][y_bin] = new TH1D(Form("h_eta_scat_ele_Q2_%i_y_%i", Q2bin, y_bin), Form("h_eta_scat_ele_Q2_%i_y_%i", Q2bin, y_bin), 200, -4, 4);
      h_eta_ele[Q2bin][y_bin] = new TH1D(Form("h_eta_ele_Q2_%i_y_%i", Q2bin, y_bin), Form("h_eta_ele_Q2_%i_y_%i", Q2bin, y_bin), 200, -4, 4);
      h_eta_pi[Q2bin][y_bin] = new TH1D(Form("h_eta_pi_Q2_%i_y_%i", Q2bin, y_bin), Form("h_eta_pi_Q2_%i_y_%i", Q2bin, y_bin), 200, -4, 4);
      h_eta_K[Q2bin][y_bin] = new TH1D(Form("h_eta_K_Q2_%i_y_%i", Q2bin, y_bin), Form("h_eta_K_Q2_%i_y_%i", Q2bin, y_bin), 200, -4, 4);
      h_eta_proton[Q2bin][y_bin] = new TH1D(Form("h_eta_proton_Q2_%i_y_%i", Q2bin, y_bin), Form("h_eta_proton_Q2_%i_y_%i", Q2bin, y_bin), 200, -4, 4);
    }
  
  }
  
  //TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
  //TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);  
  //TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
  //TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC}",100,0,20,1000,-1,1);
  
  

	tree_reader.SetEntriesRange(0, tree->GetEntries());
	
  while (tree_reader.Next())
  {

  	//MCParticles
    //finding the scattering electron
  	TLorentzVector scatMC(0,0,0,0);
  	int mc_elect_index = -1;
  	double maxPt = -99.;
  	
  	
  	for(int imc = 0; imc < mc_px_array.GetSize(); imc++)
  	{
  		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
  		
  		if(mc_pdg_array[imc] == 11 && mctrk.Mag() > maxPt)
  		//if(mc_pdg_array[imc] == 11 && mctrk.Perp() > maxPt)
  		{
  			maxPt = mctrk.Mag();
  			//maxPt = mctrk.Perp();
  			mc_elect_index = imc;
  			scatMC.SetVectM(mctrk,mc_mass_array[imc]);  			
  		}
  	} 	
  	
    //fill truth scattered electron energy
  	h_energy_MC->Fill(scatMC.E());
  	
    double y_inelPar_e = getInelParamElectron( 18, scatMC.E() );
    double Q2_electron = getQ2elec( 18, scatMC.Vect());
    
    h_y_inelPar_MC->Fill(y_inelPar_e);
    h_y_inelPar_zoom_MC->Fill(y_inelPar_e);
    
    h_Q2_MC->Fill(Q2_electron);
    h_Q2_zoom_MC->Fill(Q2_electron);
    
    //find bins
    int Q2_bin = -1;

    for(int j = 0; j < nQ2bins; j++) //loop over pT bins
    {
      if(Q2_electron > Q2_bins[j] && Q2_electron <= Q2_bins[j+1])
      {
        Q2_bin = j;
      }
      else if( Q2_electron > Q2_bins[nQ2bins] )
      {
        Q2_bin = nQ2bins;      
      }
    }
        
    
    int y_bin = -1;

    for(int j = 0; j < nyInelParBins; j++) //loop over pT bins
    {
      if(y_inelPar_e > y_bins[j] && y_inelPar_e <= y_bins[j+1])
      {
        y_bin = j;
      }
      else if( y_inelPar_e > y_bins[nyInelParBins] )
      {
        y_bin = nyInelParBins;      
      }
    }
    
    if(Q2_bin < 0 || y_bin < 0) continue;
    
    h_eta_scat_ele[Q2_bin][y_bin]->Fill(scatMC.Eta());
    
     
    //loop ove MC particles to fill distributions of produced particles
    for(int imc=0; imc < mc_px_array.GetSize(); imc++)
  	{
  		TVector3 mc_mom(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);	
  		
  		//all electrons except the scattered one (one with highest pT)
  		//if(mc_pdg_array[imc] == 11 && mc_mom.Perp() < maxPt)
  		if(mc_pdg_array[imc] == 11 && mc_mom.Mag() < maxPt)
  		{
  			h_eta_ele[Q2_bin][y_bin]->Fill(mc_mom.Eta());  			
  		}
  		
  		//pi+
  		if(mc_pdg_array[imc] == 211)
  		{
  		   h_eta_pi[Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		
  		//K+
  		if(mc_pdg_array[imc] == 321)
  		{
  		   h_eta_K[Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		
  		//proton
  		if(mc_pdg_array[imc] == 2212)
  		{
  		   h_eta_pi[Q2_bin][y_bin]->Fill(mc_mom.Eta());  		
  		}
  		
  	}//end secon particle loop
   

  }//end while over TTree entries

	output->Write();
	output->Close();

	return 0;
}
