#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text
#include "TPad.h"

double bethe_bloch_function(double* p, double*par) {
    double m_part_c = par[0];
    double na = par[1]; // 1/mol
    double re =par[2]; // cm
    double me_csquared = par[3]; // MeV
    double z = par[4];
    double big_Z = par[5];
    double a = par[6];
    double tMax = par[7];
    double i = par[8];
    
    
    double k = 4*TMath::Pi()*na*TMath::Power(re,2)*me_csquared;
    double l = -k*TMath::Power(z,2)*big_Z/a;
    double n = 2*me_csquared*tMax*TMath::Power(i,-2);
    double dEdx = 0;
    double beta_gamma = p[0]/m_part_c;
    double beta_squared = TMath::Power(beta_gamma,2)/(1+TMath::Power(beta_gamma,2));
    //double TMax_guess = (2*me_csquared*TMath::Power(beta_gamma, 2))/(1 + 2*TMath::Power(TMath::Power(beta_gamma, 2)/beta_squared, 0.5));
    //cout << TMax_guess << endl;
    
    dEdx = l/beta_squared*(0.5*TMath::Log(n*TMath::Power(beta_gamma,2))-beta_squared);//-0.5*delta);
    return dEdx;
}

void bethe_bloch() {
    //arrays to use -----------------------------------------------------------------------------
    const int pbins_num = 9;
    double pbins[pbins_num+1] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"proton", "antiproton"};
    
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
        double low_fit[pbins_num] = {-0.45,-0.48,-0.45,-0.5,-0.55,-0.55,-0.55,-0.5,-0.5};
        double high_fit[pbins_num] = {0.3,0.22,0.18,0.15,0.2,0.25,0.2,0.15,0.15};
        
        for (int p = 0; p < pbins_num; p++) {
            //get data from root file and actually do the fit
            distributions[p] = (TH1D*) file->Get(Form("%ss_p%d", p_and_pbar[m].c_str(), p));
            pionFit[p] = new TF1(Form("%ss_p%d", p_and_pbar[m].c_str(), p), "gaus",low_fit[p],high_fit[p]);
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
            //tc->SaveAs(Form("dEdx_histograms/quality1/fits/m%d_p%d.pdf",m,p));
            
            // fill histograms with pion peaks
            if (m == 0) {
                positivePionPeaks->SetBinContent(p+1,TMath::Power(TMath::E(), pionFit[p]->GetParameter(1)));
            }
            else {
                negativePionPeaks->SetBinContent(p+1,TMath::Power(TMath::E(),pionFit[p]->GetParameter(1)));
            }
        }
    }
    //plot pion peaks as a function of momentum -------------------------------------------------------------
    TCanvas *tc_pm = new TCanvas();
    tc_pm->cd();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(0);
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    positivePionPeaks->SetMarkerStyle(3);
    positivePionPeaks->SetMarkerColor(kRed);
    positivePionPeaks->SetMinimum(0.87);
    positivePionPeaks->Draw("P");
    negativePionPeaks->SetMarkerStyle(3);
    negativePionPeaks->SetMarkerColor(kBlue);
    negativePionPeaks->Draw("same P");
    positivePionPeaks->SetTitle("dE/dx pion maximum for various momentums");
    positivePionPeaks->SetXTitle("Momentum [GeV/c]");
    positivePionPeaks->SetYTitle("dE/dx [MeV cm^{2} g^{-1}])");
    latex->DrawLatex(0.68,0.66,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
    latex->DrawLatex(0.68,0.7,"#scale[0.8]{ATLAS #bf{Internal}}");
    latex->DrawLatex(0.68,0.62,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
    gStyle->SetLegendBorderSize(1);
    
    TLegend *legend = new TLegend(0.6,0.8,0.9,0.9);
    legend->AddEntry(positivePionPeaks, "Positive pions","P");
    legend->AddEntry(negativePionPeaks, "Negative pions", "P");
    legend->Draw("same");
    tc_pm->SaveAs("dEdx_histograms/momentum_vs_maxPion_dEdx.pdf");
    
    //fit the peaks to the BB function - for both the positive and negative functions ------------------------
    // guesses for the fit parameters:
    const float m_part_c = 139.5; //MeV/c
    const float na = TMath::Na(); //1/mol
    const float re = 2.817940326e-13; // cm
    const float me_csquared = 0.510998918; // MeV
    const float z = 1; // no units
    const float big_Z = 14; // no units
    const float a = 28.0855; //g/mol
    const float tMax = 1e-5; //MeV
    const float i = 1; //MeV
    
    positiveBBFit = new TF1("positiveBBFit", bethe_bloch_function, 0.3, 1.2,9.);
    positiveBBFit->SetParameters(m_part_c, na, re, me_csquared, z, big_Z, a, tMax, i);
    positiveBBFit->SetParLimits(0,m_part_c-0.1,m_part_c + 0.1);
    positiveBBFit->FixParameter(1,na);
//    positiveBBFit->SetParLimits(7,1e-6,1e-3);
//    positiveBBFit->SetParLimits(8,1e-6,1e2);
    
    positivePionPeaks->Fit(positiveBBFit, "RNL");
    
    cout << "Chi^2: " << positiveBBFit->GetChisquare() << "\n NDF: " << positiveBBFit->GetNDF() << "\n Chi^2/NDF: " << positiveBBFit->GetChisquare()/positiveBBFit->GetNDF() << endl;
}
