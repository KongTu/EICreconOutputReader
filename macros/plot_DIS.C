#include "RiceStyle.h"
using namespace std;
void plot_DIS(TString energy_config="high"){
	
	TFile* file;
	if(energy_config=="high") file=new TFile("../output/eicrecon-DIS_18x275_Q2-1-output.root");
	if(energy_config=="midd") file=new TFile("../output/eicrecon-DIS_10x100_Q2-1-output.root");
	TH1D* h_Eta_Elect_MC = (TH1D*) file->Get("h_Eta_Elect_MC");
	TH1D* h_Eta_Elect_REC = (TH1D*) file->Get("h_Eta_Elect_REC");
	TH1D* h_Eta_Elect_REC_bkg = (TH1D*) file->Get("h_Eta_Elect_REC_bkg");
	TH1D* h_Eta_Elect_REC_bkg_pfRICH = (TH1D*) file->Get("h_Eta_Elect_REC_bkg_pfRICH");

	//Q2
	TH1D* h_Q2elec_MC = (TH1D*) file->Get("h_Q2elec_MC");
	TH1D* h_Q2elec_REC = (TH1D*) file->Get("h_Q2elec_REC");

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.01);
	gPad->SetLogy(1);
	double ymin=1e0;
	double ymax=1e8;
	TH1D* base1 = makeHist("base1", "", "#eta", "counts", 100,-4.3,1,kBlack);
	base1->GetYaxis()->SetRangeUser(ymin,ymax);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.2,1.2);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(3,6,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	h_Eta_Elect_MC->Draw("HIST SAME");
	h_Eta_Elect_REC->SetMarkerStyle(24);
	h_Eta_Elect_REC->Draw("PESAME");
	h_Eta_Elect_REC_bkg->SetMarkerStyle(25);
	h_Eta_Elect_REC_bkg->Draw("PESAME");
	h_Eta_Elect_REC_bkg_pfRICH->SetMarkerStyle(20);
	h_Eta_Elect_REC_bkg_pfRICH->SetMarkerColor(kRed);
	h_Eta_Elect_REC_bkg_pfRICH->Draw("PESAME");

	TLatex* r42 = new TLatex(0.16,0.91, "EPIC Simulation - Brycecanyon 22.11.2");
	r42->SetNDC();
	r42->SetTextSize(0.04);
	r42->Draw("same");

	TLatex* r43;
	if(energy_config=="high")r43 = new TLatex(0.18, 0.84, "ep DIS 18x275 GeV, 1.0< Q^{2} < 10 GeV^{2}, 0.01 < y < 0.95");
	if(energy_config=="midd")r43 = new TLatex(0.18, 0.84, "ep DIS 10x100 GeV, 1.0< Q^{2} < 10 GeV^{2}, 0.01 < y < 0.95");
	r43->SetNDC();
	r43->SetTextSize(20);
	r43->SetTextFont(43);
	r43->SetTextColor(kBlack);
	r43->Draw("same");

	TLegend *w8 = new TLegend(0.17,0.67,0.46,0.8);
	w8->SetLineColor(kWhite);
	w8->SetFillColor(0);
	w8->SetTextSize(19);
	w8->SetTextFont(45);
	w8->AddEntry(h_Eta_Elect_MC, " Scat' e MC", "L");
	w8->AddEntry(h_Eta_Elect_REC, " Scat' e REC", "P");
	w8->AddEntry(h_Eta_Elect_REC_bkg, " Scat' e background", "P");
	w8->AddEntry(h_Eta_Elect_REC_bkg_pfRICH, " Scat' e background after pfRICH vetos", "P");
	w8->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.01);
	gPad->SetLogy(1);
	TH1D* base2 = makeHist("base2", "", "Q^{2}_{e} (GeV^{2})", "counts", 100,0,12,kBlack);
	base2->GetYaxis()->SetRangeUser(1e2,5e7);
	base2->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base2,1.2,1.2);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
	base2->GetXaxis()->SetNdivisions(3,6,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->Draw();

	h_Q2elec_MC->Draw("SAME");
	h_Q2elec_REC->SetMarkerStyle(24);
	h_Q2elec_REC->Draw("PESAME");
	TLegend *w9 = new TLegend(0.17,0.67,0.46,0.8);
	w9->SetLineColor(kWhite);
	w9->SetFillColor(0);
	w9->SetTextSize(19);
	w9->SetTextFont(45);
	w9->AddEntry(h_Q2elec_MC, "  MC", "L");
	w9->AddEntry(h_Q2elec_REC, "  REC", "P");
	w9->Draw("same");

	r42->Draw("same");
	r43->Draw("same");

}