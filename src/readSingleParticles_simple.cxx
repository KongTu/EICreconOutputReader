#include "pleaseIncludeMe.h"
int readSingleParticles_simple(TString inname="input/input.root",TString outname="test")
{	

	TString name_of_input = (TString) inname;	
	auto file=new TFile(name_of_input);
	auto tree = (TTree *) file->Get("events");
    TTreeReader tree_reader(tree);       // !the tree reader
    
    // MC particle pz array for each MC particle
    TTreeReaderArray<float> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<float> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<float> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};

    //Reconstructed EcalEndcapNClusters
    TTreeReaderArray<float> em_energy_array = {tree_reader, "EcalEndcapNClusters.energy"};
    TTreeReaderArray<float> em_x_array = {tree_reader, "EcalEndcapNClusters.position.x"};
    TTreeReaderArray<float> em_y_array = {tree_reader, "EcalEndcapNClusters.position.y"};

    //not sure how to use it (leave it for now)
    TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClustersAssociations.recID"};
    TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClustersAssociations.simID"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

    //defining output file and histos.
    TString output_name_dir = outname;
  	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

    TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
    TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
    TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
    TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC}",100,0,20,1000,-1,1);

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

    	//MCParticles
        //finding the scattering electron;
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
    		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
    		if(mc_pdg_array[imc]==11 	
    			&& mctrk.Perp()>maxPt){
    			maxPt=mctrk.Perp();
    			mc_elect_index=imc;
    			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
    		}
    	}
        //fill truth scattered electron energy
    	h_energy_MC->Fill(scatMC.E());

        //finding the leading energy cluster
    	double maxEnergy=-99.;
    	double xpos=-999.;
    	double ypos=-999.;
    	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++){
    		if(em_energy_array[iclus]>maxEnergy){
    			maxEnergy=em_energy_array[iclus];
    			xpos=em_x_array[iclus];
    			ypos=em_y_array[iclus];
    		}
    	}
    	
        //fill energy and x-y position;
		h_energy_REC->Fill(maxEnergy);
		h_emClus_position_REC->Fill(xpos,ypos);

        //calcuate the energy resolution;
		double res= (scatMC.E()-maxEnergy)/scatMC.E();
		h_energy_res->Fill(scatMC.E(), res);
	       
        //example loop over rec particles and fill eta 
    	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
    		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
    		h_eta->Fill(trk.Eta());
    	}
    }

	output->Write();
	output->Close();

	return 0;
}