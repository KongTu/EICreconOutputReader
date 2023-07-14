#include "RiceStyle.h"
using namespace std;
void plot_diffractive_vm_physics_benchmark(TString filename="./output/eicrecon-sartre_coherent_phi_output.root"){
	
	TFile* file = new TFile(filename);
	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	if(filename=="jpsi") {vm_label="J/#psi";daug_label="e^{+}e^{-}";}
	//t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_trk_REC = (TH1D*) file->Get("h_t_trk_REC");
	TH1D* h_t_combo_REC = (TH1D*) file->Get("h_t_combo_REC");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(8e-2, 8e5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	h_t_MC->Draw("same");

	h_t_REC->SetMarkerStyle(20);
	h_t_REC->Draw("PEsame");

	h_t_trk_REC->SetFillColorAlpha(kBlue,0.4);
    h_t_trk_REC->SetFillStyle(1001);
	h_t_trk_REC->SetMarkerStyle(24);
	h_t_trk_REC->SetMarkerColor(kBlue);
	// h_t_trk_REC->Draw("PE3same");

	h_t_combo_REC->SetFillColorAlpha(kRed,0.4);
    h_t_combo_REC->SetFillStyle(1001);
	h_t_combo_REC->SetMarkerStyle(24);
	h_t_combo_REC->SetMarkerColor(kRed);
	// h_t_combo_REC->Draw("PE3same");

	TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.9,0.91, "EPIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.95");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_0 = new TLatex(0.18, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");

	TLatex* r44_2 = new TLatex(0.18, 0.18, ""+vm_label+" #rightarrow "+daug_label);
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TLegend *w7 = new TLegend(0.48,0.68,0.93,0.76);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	if(filename=="MCclusterEnergy")
	{
		w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
		w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true EEMC E ", "P");
		w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only ", "P");

	}
	else if(filename=="MCvmAndelectron")
	{
		w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
		w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true e' ", "P");
		w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" MC w. RECO e' ", "P");
	}
	else{
		w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
		w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. EEMC ", "P");
		// w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only", "P");
	    // w7->AddEntry(h_t_combo_REC, "Sartre "+vm_label+" RECO best", "P");

	}
	

	w7->Draw("same");
	TString outputfilename = filename;
	outputfilename.ReplaceAll("output/","");
	c1->Print("./figures/benchmark-phi-dsigmadt-from-"+outputfilename+".pdf");

}