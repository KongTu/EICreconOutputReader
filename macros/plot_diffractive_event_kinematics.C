#include "RiceStyle.h"
using namespace std;
void plot_diffractive_event_kinematics(TString filename="phi", int doLog_=1){

	TFile* file = new TFile("../output/eicrecon-sartre_coherent_"+filename+"_output.root");
	//DIS kine
	TH1D* h_Q2_e = (TH1D*) file->Get("h_Q2_e");
    TH1D* h_y_e = (TH1D*) file->Get("h_y_e");  
 	TH1D* h_energy_MC = (TH1D*) file->Get("h_energy_MC");
    TH1D* h_Q2REC_e = (TH1D*) file->Get("h_Q2REC_e");
    TH1D* h_yREC_e = (TH1D*) file->Get("h_yREC_e");
    TH1D* h_energy_REC = (TH1D*) file->Get("h_energy_REC");
    //Epz
    TH1D* h_Epz_REC = (TH1D*) file->Get("h_Epz_REC");
    TH1D* h_trk_Epz_REC = (TH1D*) file->Get("h_trk_Epz_REC");
 	//cluster
    TH1D* h_EoverP_REC = (TH1D*) file->Get("h_EoverP_REC");  
    TH2D* h_emHits_position_REC = (TH2D*) file->Get("h_emHits_position_REC");//default cluster positio
    TH1D* h_energy_calibration_REC = (TH1D*) file->Get("h_energy_calibration_REC");
    TH1D* h_ClusOverHit_REC = (TH1D*) file->Get("h_ClusOverHit_REC");


    TCanvas* c1 = new TCanvas("c1","c1",1,1,1600,800);
    c1->Divide(4,2,0.01,0.01);
    c1->cd(1);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Q2_e->GetXaxis()->SetTitleSize(0.8*h_Q2_e->GetXaxis()->GetTitleSize());
	h_Q2_e->GetXaxis()->SetLabelSize(0.8*h_Q2_e->GetXaxis()->GetLabelSize());
	h_Q2_e->GetYaxis()->SetTitleSize(0.8*h_Q2_e->GetYaxis()->GetTitleSize());
	h_Q2_e->GetYaxis()->SetLabelSize(0.8*h_Q2_e->GetYaxis()->GetLabelSize());
	h_Q2_e->GetXaxis()->SetTitleOffset(1.6*h_Q2_e->GetXaxis()->GetTitleOffset());
	h_Q2_e->GetYaxis()->SetTitleOffset(2.0*h_Q2_e->GetYaxis()->GetTitleOffset());
	h_Q2_e->GetYaxis()->SetTitle("counts");
    h_Q2_e->Draw();
    h_Q2REC_e->SetMarkerStyle(24);
    h_Q2REC_e->Draw("PEsame");
    TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_energy_MC, "MC ", "L");
	w7->AddEntry(h_energy_REC, "RECO", "P");
	w7->Draw("same");

    c1->cd(2);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_y_e->GetXaxis()->SetTitleSize(0.8*h_y_e->GetXaxis()->GetTitleSize());
	h_y_e->GetXaxis()->SetLabelSize(0.8*h_y_e->GetXaxis()->GetLabelSize());
	h_y_e->GetYaxis()->SetTitleSize(0.8*h_y_e->GetYaxis()->GetTitleSize());
	h_y_e->GetYaxis()->SetLabelSize(0.8*h_y_e->GetYaxis()->GetLabelSize());
	h_y_e->GetXaxis()->SetTitleOffset(1.6*h_y_e->GetXaxis()->GetTitleOffset());
	h_y_e->GetYaxis()->SetTitleOffset(2.0*h_y_e->GetYaxis()->GetTitleOffset());
	h_y_e->GetYaxis()->SetTitle("counts");
    h_y_e->Draw();
    h_yREC_e->SetMarkerStyle(24);
    h_yREC_e->Draw("PEsame");
    w7->Draw("same");

    c1->cd(3);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_energy_MC->GetXaxis()->SetTitleSize(0.8*h_energy_MC->GetXaxis()->GetTitleSize());
	h_energy_MC->GetXaxis()->SetLabelSize(0.8*h_energy_MC->GetXaxis()->GetLabelSize());
	h_energy_MC->GetYaxis()->SetTitleSize(0.8*h_energy_MC->GetYaxis()->GetTitleSize());
	h_energy_MC->GetYaxis()->SetLabelSize(0.8*h_energy_MC->GetYaxis()->GetLabelSize());
	h_energy_MC->GetXaxis()->SetTitleOffset(1.6*h_energy_MC->GetXaxis()->GetTitleOffset());
	h_energy_MC->GetYaxis()->SetTitleOffset(2.5*h_energy_MC->GetYaxis()->GetTitleOffset());
	h_energy_MC->GetYaxis()->SetTitle("counts");
    h_energy_MC->Draw();
    h_energy_REC->SetMarkerStyle(24);
    h_energy_REC->Draw("PEsame");
    w7->Draw("same");
    // TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
	// w7->SetLineColor(kWhite);
	// w7->SetFillColor(0);
	// w7->SetTextSize(17);
	// w7->SetTextFont(45);
	// w7->AddEntry(h_energy_MC, "MC ", "L");
	// w7->AddEntry(h_energy_REC, "RECO", "P");
	// w7->Draw("same");

    c1->cd(4);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_Epz_REC->GetXaxis()->SetTitleSize(0.8*h_Epz_REC->GetXaxis()->GetTitleSize());
	h_Epz_REC->GetXaxis()->SetLabelSize(0.8*h_Epz_REC->GetXaxis()->GetLabelSize());
	h_Epz_REC->GetYaxis()->SetTitleSize(0.8*h_Epz_REC->GetYaxis()->GetTitleSize());
	h_Epz_REC->GetYaxis()->SetLabelSize(0.8*h_Epz_REC->GetYaxis()->GetLabelSize());
	h_Epz_REC->GetXaxis()->SetTitleOffset(1.6*h_Epz_REC->GetXaxis()->GetTitleOffset());
	h_Epz_REC->GetYaxis()->SetTitleOffset(2.5*h_Epz_REC->GetYaxis()->GetTitleOffset());
	h_Epz_REC->GetYaxis()->SetTitle("counts");
    h_Epz_REC->SetMarkerStyle(24);
    h_Epz_REC->Draw("PEsame");

    c1->cd(5);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_energy_calibration_REC->GetXaxis()->SetTitleSize(0.8*h_energy_calibration_REC->GetXaxis()->GetTitleSize());
	h_energy_calibration_REC->GetXaxis()->SetLabelSize(0.8*h_energy_calibration_REC->GetXaxis()->GetLabelSize());
	h_energy_calibration_REC->GetYaxis()->SetTitleSize(0.8*h_energy_calibration_REC->GetYaxis()->GetTitleSize());
	h_energy_calibration_REC->GetYaxis()->SetLabelSize(0.8*h_energy_calibration_REC->GetYaxis()->GetLabelSize());
	h_energy_calibration_REC->GetXaxis()->SetTitleOffset(1.6*h_energy_calibration_REC->GetXaxis()->GetTitleOffset());
	h_energy_calibration_REC->GetYaxis()->SetTitleOffset(2.5*h_energy_calibration_REC->GetYaxis()->GetTitleOffset());
	h_energy_calibration_REC->GetYaxis()->SetTitle("counts");
	h_energy_calibration_REC->GetXaxis()->SetTitle("E_{reco} / E_{mc}");
    h_energy_calibration_REC->SetMarkerStyle(24);
    h_energy_calibration_REC->Draw("PEsame");

    c1->cd(6);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_EoverP_REC->GetXaxis()->SetTitleSize(0.8*h_EoverP_REC->GetXaxis()->GetTitleSize());
	h_EoverP_REC->GetXaxis()->SetLabelSize(0.8*h_EoverP_REC->GetXaxis()->GetLabelSize());
	h_EoverP_REC->GetYaxis()->SetTitleSize(0.8*h_EoverP_REC->GetYaxis()->GetTitleSize());
	h_EoverP_REC->GetYaxis()->SetLabelSize(0.8*h_EoverP_REC->GetYaxis()->GetLabelSize());
	h_EoverP_REC->GetXaxis()->SetTitleOffset(1.6*h_EoverP_REC->GetXaxis()->GetTitleOffset());
	h_EoverP_REC->GetYaxis()->SetTitleOffset(2.5*h_EoverP_REC->GetYaxis()->GetTitleOffset());
	h_EoverP_REC->GetYaxis()->SetTitle("counts");
    h_EoverP_REC->SetMarkerStyle(24);
    h_EoverP_REC->Draw("PEsame");

    c1->cd(7);
    gPad->SetLogy(doLog_);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_ClusOverHit_REC->GetXaxis()->SetTitleSize(0.8*h_ClusOverHit_REC->GetXaxis()->GetTitleSize());
	h_ClusOverHit_REC->GetXaxis()->SetLabelSize(0.8*h_ClusOverHit_REC->GetXaxis()->GetLabelSize());
	h_ClusOverHit_REC->GetYaxis()->SetTitleSize(0.8*h_ClusOverHit_REC->GetYaxis()->GetTitleSize());
	h_ClusOverHit_REC->GetYaxis()->SetLabelSize(0.8*h_ClusOverHit_REC->GetYaxis()->GetLabelSize());
	h_ClusOverHit_REC->GetXaxis()->SetTitleOffset(1.6*h_ClusOverHit_REC->GetXaxis()->GetTitleOffset());
	h_ClusOverHit_REC->GetYaxis()->SetTitleOffset(2.5*h_ClusOverHit_REC->GetYaxis()->GetTitleOffset());
	h_ClusOverHit_REC->GetYaxis()->SetTitle("counts");
    h_ClusOverHit_REC->SetMarkerStyle(24);
    h_ClusOverHit_REC->Draw("PEsame");

    c1->cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.15);
	h_emHits_position_REC->GetXaxis()->SetTitleSize(0.8*h_emHits_position_REC->GetXaxis()->GetTitleSize());
	h_emHits_position_REC->GetXaxis()->SetLabelSize(0.8*h_emHits_position_REC->GetXaxis()->GetLabelSize());
	h_emHits_position_REC->GetYaxis()->SetTitleSize(0.8*h_emHits_position_REC->GetYaxis()->GetTitleSize());
	h_emHits_position_REC->GetYaxis()->SetLabelSize(0.8*h_emHits_position_REC->GetYaxis()->GetLabelSize());
	h_emHits_position_REC->GetXaxis()->SetTitleOffset(1.6*h_emHits_position_REC->GetXaxis()->GetTitleOffset());
	h_emHits_position_REC->GetYaxis()->SetTitleOffset(2.*h_emHits_position_REC->GetYaxis()->GetTitleOffset());
	h_emHits_position_REC->GetYaxis()->CenterTitle();
	h_emHits_position_REC->GetYaxis()->SetTitle("y(mm)");
    h_emHits_position_REC->Draw("colzsame");

	c1->Print("../figures/benchmark-phi-DIS-kinematics.pdf");


}