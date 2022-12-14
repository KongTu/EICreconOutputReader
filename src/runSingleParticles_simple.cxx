#include "pleaseIncludeMe.h"
int runSingleParticles_simple(TString inname="input/input.root",TString outname="test")
{	

	TString name_of_input = (TString) inname;
	std::cout << "what is this rec_file = " << name_of_input << endl;
	
    auto tree = new TChain("events");
    tree->Add(name_of_input);
    TTreeReader tree_reader(tree);       // !the tree reader
    
    TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    
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
    TTreeReaderArray<float> emhits_x_array = {tree_reader, "EcalEndcapNRecHits.position.x"};
    TTreeReaderArray<float> emhits_y_array = {tree_reader, "EcalEndcapNRecHits.position.y"};

    TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader, "EcalEndcapNClustersAssociations.recID"};
    TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader, "EcalEndcapNClustersAssociations.simID"};

    // Reconstructed particles pz array for each reconstructed particle
    TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};

    TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticlesAssociations.recID"};
    TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticlesAssociations.simID"};

    TString output_name_dir = outname;
  	TFile* output = new TFile("output/"+output_name_dir+"-output.root","RECREATE");

    TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",100,0,20);
    TH1D* h_y_e = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
    TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);
    TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
    TH1D* h_energy_REC = new TH1D("h_energy_REC",";E_{REC} (GeV)",100,0,20);
    TH2D* h_emClus_position_REC = new TH2D("h_emClus_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
    TH2D* h_emHits_position_REC = new TH2D("h_emHits_position_REC",";x (cm);y (cm)",400,-800,800,400,-800,800);
    TH2D* h_energy_res = new TH2D("h_energy_res",";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} emcal",100,0,20,1000,-1,1);
    TH1D* h_energy_calibration_REC = new TH1D("h_energy_calibration_REC",";E (GeV)",200,0,2);

	tree_reader.SetEntriesRange(0, tree->GetEntries());
    while (tree_reader.Next()) {

        /*
        Beam particles
        */
        TLorentzVector ebeam(0,0,0,0);
        TLorentzVector pbeam(0,0,0,0);

    	//MCParticles
    	TLorentzVector scatMC(0,0,0,0);
    	int mc_elect_index=-1;
    	double maxPt=-99.;
    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
    		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
            if(mc_genStatus_array[imc]==4){
                if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
                if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON);
            }
            if(mc_genStatus_array[imc]!=1) continue;
            if(mc_pdg_array[imc]==11    
                && mctrk.Perp()>maxPt){
                maxPt=mctrk.Perp();
                mc_elect_index=imc;
                scatMC.SetVectM(mctrk,mc_mass_array[imc]);
            }
    	}
        TLorentzVector qbeam=ebeam-scatMC;
        double Q2=-(qbeam).Mag2();  
        double pq=pbeam.Dot(qbeam);
        double y= pq/pbeam.Dot(ebeam);
        
        //MC level phase space cut
        if(Q2<1.||Q2>10.) continue;
        if(y<0.01||y>0.95) continue;

        h_Q2_e->Fill(Q2);
        h_y_e->Fill(y);
        h_energy_MC->Fill(scatMC.E());

    	double maxEnergy=-99.;
    	double xpos=-999.;
    	double ypos=-999.;
        double xhitpos=-999.;
        double yhitpos=-999.;
    	for(int iclus=0;iclus<em_energy_array.GetSize();iclus++){
    		if(em_energy_array[iclus]>maxEnergy){
    			maxEnergy=em_energy_array[iclus];
    			xpos=em_x_array[iclus];
    			ypos=em_y_array[iclus];
                xhitpos=emhits_x_array[iclus];
                yhitpos=emhits_y_array[iclus];
    		}
    	}

    	maxEnergy *= 1.037;
		h_energy_calibration_REC->Fill( maxEnergy / scatMC.E() );
        // if(fabs(xpos)<100. && fabs(ypos)<100.) continue; //100mm square hole.

        double res= (scatMC.E()-maxEnergy)/scatMC.E();
        h_energy_res->Fill(scatMC.E(), res);
        h_energy_REC->Fill(maxEnergy);
        h_emClus_position_REC->Fill(xpos,ypos);
        h_emHits_position_REC->Fill(xhitpos,yhitpos);

    	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
    		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
    		h_eta->Fill(trk.Eta());
    	}
    }
	output->Write();
	output->Close();

	return 0;
}