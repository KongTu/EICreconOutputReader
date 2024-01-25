#include "pleaseIncludeMe.h"

auto giveme_t_method_L(TLorentzVector eIn, 
					   TLorentzVector eOut, 
					   TLorentzVector pIn, 
					   TLorentzVector vmOut)
{
	TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
	double method_L = 0;
	TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
	double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
	double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
	double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
	TLorentzVector a_beam_scattered_corr; 
	a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
	method_L = -(a_beam_scattered_corr-aInVec).Mag2();

	return method_L;
}
auto giveme_t_method_E(TLorentzVector eIn, 
					   TLorentzVector eOut, 
					   TLorentzVector pIn, 
					   TLorentzVector vmOut)
{
	double method_E;
	method_E = -(eIn-eOut-vmOut).Mag2();
	return method_E;
}
auto giveme_t_method_A(TLorentzVector eIn, 
					   TLorentzVector eOut, 
					   TLorentzVector pIn, 
					   TLorentzVector vmOut)
{
	double method_A;
	TVector2 sum_pt(vmOut.Px()+eOut.Px(), vmOut.Py()+eOut.Py());
	method_A = sum_pt.Mod2();
	return method_A;
}
auto giveme_pt2(TLorentzVector vmOut)
{
	double method_pt2;
	TVector2 sum_pt(vmOut.Px(), vmOut.Py());
	method_pt2 = sum_pt.Mod2();
	return method_pt2;
}

int photoproduction_vm_analysis(TString rec_file, TString outputfile, bool veto_)
{	
// read our configuration	
TString name_of_input = (TString) rec_file;
std::cout << "Input file = " << name_of_input << endl;
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
TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};

// Reconstructed particles pz array for each reconstructed particle
TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};

TTreeReaderArray<unsigned int> rec_id = {tree_reader, "ReconstructedChargedParticleAssociations.recID"};
TTreeReaderArray<unsigned int> sim_id = {tree_reader, "ReconstructedChargedParticleAssociations.simID"};

//FF system
TTreeReaderArray<float> zdc_ecal_cluster_x = {tree_reader, "EcalFarForwardZDCClusters.position.x"};
	
TTreeReaderArray<float> b0_hits_x = {tree_reader, "B0TrackerRecHits.position.x"};
TTreeReaderArray<float> b0_hits_z = {tree_reader, "B0TrackerRecHits.position.z"};
	
TTreeReaderArray<float> reco_RP_px = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
TTreeReaderArray<float> reco_RP_py = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
TTreeReaderArray<float> reco_RP_pz = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
	
TTreeReaderArray<float> reco_OM_px = {tree_reader, "ForwardOffMRecParticles.momentum.x"};
TTreeReaderArray<float> reco_OM_py = {tree_reader, "ForwardOffMRecParticles.momentum.y"};
TTreeReaderArray<float> reco_OM_pz = {tree_reader, "ForwardOffMRecParticles.momentum.z"};

TString output_name_dir = outputfile+"_output.root";
cout << "Output file = " << output_name_dir << endl;
TFile* output = new TFile(output_name_dir,"RECREATE");

//events
TH1D* h_Q2_e = new TH1D("h_Q2_e",";Q^{2}_{e,MC}",1000,0.001,0.1);
TH1D* h_y_e = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
TH1D* h_energy_MC = new TH1D("h_energy_MC",";E_{MC} (GeV)",100,0,20);
TH1D* h_t_MC = new TH1D("h_t_MC",";t_{MC}; counts",100,0,0.2);
TH1D* h_VM_mass_MC = new TH1D("h_VM_mass_MC",";mass (GeV)",200,0,4);

TH1D* h_Q2REC_e = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,20);
TH1D* h_yREC_e = new TH1D("h_yREC_e",";y_{e,REC}",100,0,1);

//track
TH1D* h_eta = new TH1D("h_eta",";#eta",100,-5,5);

//VM & t
TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC",";mass (GeV)",200,0,4);
TH1D* h_VM_pt_REC = new TH1D("h_VM_pt_REC",";p_{T} (GeV/c)",200,0,2);
TH2D* h_VM_res = new TH2D("h_VM_res",";p_{T,MC} (GeV); p_{T,MC}-E_{T,REC}/p_{T,MC}",100,0,2,1000,-1,1);
TH1D* h_t_REC = new TH1D("h_t_REC",";t_{REC} (GeV^{2}); counts",100,0,0.2);
TH2D* h_t_res = new TH2D("h_t_res",";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC}",100,0,0.2,1000,-10,10);
TH2D* h_t_2D = new TH2D("h_t_2D",";t_{MC} (GeV^{2}); t_{REC} (GeV^{2}) track-base",100,0,0.2,100,0,0.2);

tree_reader.SetEntriesRange(0, tree->GetEntries());
while (tree_reader.Next()) {

	/*
	Beam particles
	*/
	TLorentzVector ebeam(0,0,0,0);
	TLorentzVector pbeam(0,0,0,0);

	TLorentzVector vmMC(0,0,0,0);
	TLorentzVector kplusMC(0,0,0,0);
	TLorentzVector kminusMC(0,0,0,0);

	//MC level
	TLorentzVector scatMC(0,0,0,0);
	int mc_elect_index=-1;
	double maxPt=-99.;
	for(int imc=0;imc<mc_px_array.GetSize();imc++){
		TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);	
		if(mc_genStatus_array[imc]==4){//4 is Sartre.
			if(mc_pdg_array[imc]==11) ebeam.SetVectM(mctrk, MASS_ELECTRON);
				if(mc_pdg_array[imc]==1000791970) pbeam.SetVectM(mctrk, MASS_PROTON);
		}
		if(mc_genStatus_array[imc]!=1) continue;
		if(mc_pdg_array[imc]==11 	
			&& mctrk.Perp()>maxPt){
			maxPt=mctrk.Perp();
			mc_elect_index=imc;
			scatMC.SetVectM(mctrk,mc_mass_array[imc]);
		}
		if(mc_pdg_array[imc]==321
			&& mc_genStatus_array[imc]==1) kplusMC.SetVectM(mctrk,MASS_KAON);
		if(mc_pdg_array[imc]==-321
			&& mc_genStatus_array[imc]==1) kminusMC.SetVectM(mctrk,MASS_KAON);

	}
	vmMC=kplusMC+kminusMC;
	//protection.
	if(ebeam.E()==pbeam.E() && ebeam.E()==0) {
		std::cout << "problem with MC incoming beams" << std::endl;
		continue;
	}
	TLorentzVector qbeam=ebeam-scatMC;
	double Q2=-(qbeam).Mag2();  
	double pq=pbeam.Dot(qbeam);
	double y= pq/pbeam.Dot(ebeam);
	
	//MC level phase space cut
	if(Q2>1.) continue;
	// if(y<0.01||y>0.95) continue;

	h_Q2_e->Fill(Q2);
	h_y_e->Fill(y);
	h_energy_MC->Fill(scatMC.E());

	double t_MC=0.;
	if(vmMC.E()!=0 
		&& fabs(vmMC.Rapidity())<3.5)
	{
		double method_E = -(qbeam-vmMC).Mag2();
		t_MC=method_E;
		h_t_MC->Fill( method_E );
		h_VM_mass_MC->Fill( vmMC.M() );
	}

	//rec level
	//veto FFs
	if(veto_) 
	{
		//ZDC
		if(zdc_ecal_cluster_x.GetSize()>0) continue;
		//B0
		int b0hits[4]={0,0,0,0};
		for(int ihit=0;ihit<b0_hits_z.GetSize();ihit++){
			if(b0_hits_z[ihit]>5700 && b0_hits_z[ihit]<5900) b0hits[0]++;
			if(b0_hits_z[ihit]>6100 && b0_hits_z[ihit]<6200) b0hits[1]++;
			if(b0_hits_z[ihit]>6400 && b0_hits_z[ihit]<6500) b0hits[2]++;
			if(b0_hits_z[ihit]>6700 && b0_hits_z[ihit]<6750) b0hits[3]++;
		}
		if(b0hits[0]>0&&b0hits[1]>0&&b0hits[2]>0&&b0hits[3]>0) continue;
		//RP
		int num_RP_tracks=0;
		for(int irp=0;irp<reco_RP_px.GetSize();irp++){
			TVector3 prec_rp(reco_RP_px[irp], reco_RP_py[irp], reco_RP_pz[irp]);
			if(prec_rp.Mag()>55) num_RP_tracks++;
		}
		if(num_RP_tracks>0) continue;
		//OMD
		int num_OMD_tracks=0;
		for(int iomd=0;iomd<reco_OM_px.GetSize();iomd++){
			TVector3 prec_omd(reco_OM_px[iomd], reco_OM_py[iomd], reco_OM_pz[iomd]);
			if(prec_omd.Mag()>55) num_OMD_tracks++;
		}
		if(num_OMD_tracks>0) continue;
	}

    TLorentzVector hfs(0,0,0,0);
    TLorentzVector particle(0,0,0,0);
    TLorentzVector kplusREC(0,0,0,0);
    TLorentzVector kminusREC(0,0,0,0);
    TLorentzVector vmREC(0,0,0,0);

	//loop over track again;
	for(int itrk=0;itrk<reco_pz_array.GetSize();itrk++){
		TVector3 trk(reco_px_array[itrk],reco_py_array[itrk],reco_pz_array[itrk]);
		particle.SetVectM(trk,MASS_PION);//assume pions;

		hfs += particle; //hfs 4vector sum.
		//selecting phi->kk daughters;
		h_eta->Fill(trk.Eta());
		if(fabs(trk.Eta())<3.0){
			if(reco_charge_array[itrk]>0) kplusREC.SetVectM(trk,MASS_KAON);
			if(reco_charge_array[itrk]<0) kminusREC.SetVectM(trk,MASS_KAON);
		}
		
	}
	//4vector of VM;
	if(kplusREC.E()!=0. && kminusREC.E()!=0.){
		vmREC=kplusREC+kminusREC;
	}

	//cluster-base DIS kine;
	// TLorentzVector qbeamREC=ebeam-scatClusEREC;
	// double Q2REC=-(qbeamREC).Mag2();  
	// double pqREC=pbeam.Dot(qbeamREC);
	// double yREC= pqREC/pbeam.Dot(ebeam);
	// h_Q2REC_e->Fill(Q2REC);
	// h_yREC_e->Fill(yREC);

	// //Event selection:
	// if( EpzREC<27||EpzREC>40 ) continue;
	// if( EoverP<0.8||EoverP>1.18 ) continue;		

	// //REC level phase space cut
	// if(Q2REC>1) continue;
	// if(yREC<0.01||yREC>0.95) continue;

	//VM rec
	if(vmREC.E()==0) continue;
	double phi_mass = vmREC.M();
	h_VM_mass_REC->Fill(phi_mass);
	h_VM_pt_REC->Fill(vmREC.Pt());

	//select phi mass and rapidity window 
	if( fabs(phi_mass-1.02)<0.02
    		&& fabs(vmREC.Rapidity())<3.5 ){
    	//2 versions: track and energy cluster:
    	double t_REC = giveme_pt2(vmREC);
    	h_t_REC->Fill( t_REC );

		//t resolution;
		double res= (t_MC-t_REC)/t_MC;
		h_t_res->Fill(t_MC, res);
		
		//2D t
		h_t_2D->Fill(t_MC,t_REC);

    	//VM pt resolution;
    	res= (vmMC.Pt()-vmREC.Pt())/vmMC.Pt();
		h_VM_res->Fill(vmMC.Pt(), res);
	}

}
output->Write();
output->Close();

return 0;
}
