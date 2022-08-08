#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"

double bethe_bloch_function(double* p, double*par) {
    double m_part_c = 1;
    double l = par[0];
    double n = par[1];
    double r = par[2];
    
    double beta_gamma = p[0]/m_part_c;
    double beta_squared = TMath::Power(beta_gamma,2)/(1+TMath::Power(beta_gamma,2));
    
    double dEdx = 0;
    dEdx = -(l/beta_squared)*(0.5*TMath::Log(n*TMath::Power(beta_gamma,2))-beta_squared)+r;
    return dEdx;
}

void bethe_bloch() {
    //arrays to use -----------------------------------------------------------------------------
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"positive", "negative"};
    float xshift[pbins_num] = {-0.102858, -0.11343, -0.1139, -0.100269, -0.0970618, -0.099067, -0.101021, -0.0971758};
    
    //load data and initialize histograms/functions
    TFile *file = new TFile("root_files/pPb_data.root","read");
    TH1D* distributions[2*pbins_num];
    TF1* pionFit[pbins_num];
    TH1D* positivePionPeaks = new TH1D("positivePionPeaks", "", pbins_num, pbins);
    TH1D* negativePionPeaks = new TH1D("negativePionPeaks", "", pbins_num, pbins);
    TH1D* momentumCenter = new TH1D("momentum_center","", pbins_num, pbins);
    TF1* positiveBBFit;
    TF1* negativeBBFit;
    
    for (int m = 0; m < 2; m++) {
        //define limits for gaussian fit based on data -------------------------------------------
        double low_fit[pbins_num] = {-0.48,-0.45,-0.5,-0.55,-0.55,-0.55,-0.5,-0.5};
        double high_fit[pbins_num] = {0.22,0.18,0.15,0.2,0.25,0.2,0.15,0.15};
        
        for (int p = 0; p < pbins_num; p++) {
            //get data from root file and actually do the fit
            distributions[p] = (TH1D*) file->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            pionFit[p] = new TF1(Form("%s_p%d", p_and_pbar[m].c_str(), p), "gaus",low_fit[p]-xshift[p],high_fit[p]-xshift[p]);
            distributions[p]->Fit(pionFit[p],"RNL");
            cout  << "Chi^2: " << pionFit[p]->GetChisquare() << "\n" << "NDF: " << pionFit[p]->GetNDF() << "\n" << "Chi^2/NDF: " << pionFit[p]->GetChisquare() / pionFit[p]->GetNDF() << "\n" << "Mean: " << pionFit[p]->GetParameter(1) << "\n" << endl;
            
            //plot root data and gaussian fit
            TCanvas *tc = new TCanvas();
            gStyle->SetOptStat(0);
            gPad->SetTicks();
            gPad->SetLogy(1);
            TLatex *latex_data = new TLatex();
            latex_data->SetNDC(kTRUE);
            
            distributions[p]->SetLineWidth(2);
            distributions[p]->SetLineColor(kBlack);
            distributions[p]->Draw();
            pionFit[p]->SetLineWidth(2);
            pionFit[p]->SetLineColor(kRed);
            pionFit[p]->Draw("same");
            tc->SaveAs(Form("dEdx_histograms/quality1/fits/m%d_p%d.pdf",m,p));
            
            // fill histograms with pion peaks
            if (m == 0) {
                double mu = pionFit[p]->GetParameter(1);
                double delta = pionFit[p]->GetParError(1);
                positivePionPeaks->SetBinContent(p+1,exp(mu));
                positivePionPeaks->SetBinError(p+1,TMath::Abs(exp(mu+delta)-exp(mu)));
            }
            else {
                double mu = pionFit[p]->GetParameter(1);
                double delta = pionFit[p]->GetParError(1);
                negativePionPeaks->SetBinContent(p+1,exp(mu));
                negativePionPeaks->SetBinError(p+1,TMath::Abs(exp(mu+delta)-exp(mu)));
            }
        }
    }
    //plot pion peaks as a function of momentum -------------------------------------------------------------
    TCanvas *tc_pm = new TCanvas();
    tc_pm->cd();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(0);
    gStyle->SetErrorX(0);
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    positivePionPeaks->SetMarkerStyle(3);
    positivePionPeaks->SetMarkerColor(kRed);
    positivePionPeaks->SetLineColor(kRed);
    positivePionPeaks->SetMinimum(0.87);
    positivePionPeaks->Draw("P");
    negativePionPeaks->SetMarkerStyle(3);
    negativePionPeaks->SetMarkerColor(kBlue);
    negativePionPeaks->SetLineColor(kBlue);
    negativePionPeaks->Draw("same P");
    positivePionPeaks->SetTitle("dE/dx pion maximum for various momentums");
    positivePionPeaks->SetXTitle("Momentum [GeV/c]");
    positivePionPeaks->SetYTitle("dE/dx [MeV cm^{2} g^{-1}])");
    positivePionPeaks->SetMinimum(0.85);
    positivePionPeaks->SetMaximum(1.11);
    latex->DrawLatex(0.68,0.66,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
    latex->DrawLatex(0.68,0.7,"#scale[0.8]{ATLAS #bf{Internal}}");
    latex->DrawLatex(0.68,0.62,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
    gStyle->SetLegendBorderSize(1);
    
    TLegend *legend = new TLegend(0.68,0.75,0.85,0.85);
    legend->AddEntry(positivePionPeaks, "Positive pions","P");
    legend->AddEntry(negativePionPeaks, "Negative pions", "P");
    legend->Draw("same");
    tc_pm->SaveAs("dEdx_histograms/momentum_vs_maxPion_dEdx.pdf");
    
    const float l = 0.0464;
    const float n = 11;
    const float r = 0.93;
    
    double *p = 0;
    // fit the positive data  ----------------------------------------------
    positiveBBFit = new TF1("positiveBBFit", bethe_bloch_function, 0.3, 1.2, 3.);
    positiveBBFit->SetParameters(l,n,r);
    positiveBBFit->SetParLimits(0,0.04,0.05);
    positiveBBFit->SetParLimits(1,9,13);
    positiveBBFit->SetParLimits(2,0.9,1);
    positivePionPeaks->Fit(positiveBBFit, "RNQ");
    positivePionPeaks->Fit(positiveBBFit, "RNQ");
    positivePionPeaks->Fit(positiveBBFit, "RNQ");
    positivePionPeaks->Fit(positiveBBFit, "RN");

    cout << "Chi^2: " << positiveBBFit->GetChisquare() << "\n NDF: " << positiveBBFit->GetNDF() << "\n Chi^2/NDF: " << positiveBBFit->GetChisquare()/positiveBBFit->GetNDF() << endl;

    //Plot the positive fit to the BB function
    TCanvas *tc_bb_p = new TCanvas();
    tc_bb_p->cd();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(0);
    TLatex *bbLatex = new TLatex();
    bbLatex->SetNDC(kTRUE);
    
    positiveBBFit->SetLineWidth(2);
    positiveBBFit->SetLineColor(kRed);
    positiveBBFit->Draw("");
    positivePionPeaks->SetMarkerStyle(3);
    positivePionPeaks->SetMarkerColor(kBlack);
    positivePionPeaks->Draw("same P");
    
    positiveBBFit->SetTitle("dE/dx pion maximum vs momentums;Momentum [GeV/c] / 0.1395;dE/dx [MeV g^{-1} cm^{2}])");
    bbLatex->DrawLatex(0.68,0.66,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
    bbLatex->DrawLatex(0.68,0.7,"#scale[0.8]{ATLAS #bf{Internal}}");
    bbLatex->DrawLatex(0.68,0.62,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
    TLegend *legend_p = new TLegend(0.68,0.8,0.9,0.9);
    legend_p->AddEntry(positivePionPeaks, "Positive Pions", "P");
    legend_p->AddEntry(positiveBBFit, "Fit to BB function", "l");
    gStyle->SetLegendBorderSize(1);
    legend_p->Draw("same");
    
    tc_bb_p->SaveAs("dEdx_histograms/bbFit_positive_pions.pdf");
    
    // fit the negative data  ----------------------------------------------
    negativeBBFit = new TF1("negativeBBFit", bethe_bloch_function, 0.3, 1.2, 3.);
    negativeBBFit->SetParameters(l,n,r);
    negativeBBFit->SetParLimits(0,0.04,0.07);
    negativeBBFit->SetParLimits(1,9,13);
    negativeBBFit->SetParLimits(2,0.9,1);
    negativePionPeaks->Fit(negativeBBFit, "RNQ");
    negativePionPeaks->Fit(negativeBBFit, "RNQ");
    negativePionPeaks->Fit(negativeBBFit, "RNQ");
    negativePionPeaks->Fit(negativeBBFit, "RN");

    cout << "Chi^2: " << negativeBBFit->GetChisquare() << "\n NDF: " << negativeBBFit->GetNDF() << "\n Chi^2/NDF: " << negativeBBFit->GetChisquare()/negativeBBFit->GetNDF() << endl;
    
    //Plot the negative fit to the BB function
    TCanvas *tc_bb_n = new TCanvas();
    tc_bb_n->cd();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(0);
    negativeBBFit->SetLineWidth(2);
    negativeBBFit->SetLineColor(kBlue);
    negativeBBFit->Draw("");
    negativePionPeaks->SetMarkerStyle(3);
    negativePionPeaks->SetMarkerColor(kBlack);
    negativePionPeaks->Draw("same P");
    
    negativeBBFit->SetTitle("dE/dx pion maximum for vs;Momentum [GeV/c] / 0.1395;dE/dx [MeV g^{-1} cm^{2}])");
    bbLatex->DrawLatex(0.68,0.66,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
    bbLatex->DrawLatex(0.68,0.7,"#scale[0.8]{ATLAS #bf{Internal}}");
    bbLatex->DrawLatex(0.68,0.62,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
    TLegend *legend_n = new TLegend(0.68,0.8,0.9,0.9);
    legend_n->AddEntry(negativePionPeaks, "Negative Pions", "P");
    legend_n->AddEntry(negativeBBFit, "Fit to BB function", "l");
    gStyle->SetLegendBorderSize(1);
    legend_n->Draw("same");
    
    tc_bb_n->SaveAs("dEdx_histograms/bbFit_negative_pions.pdf");
}
