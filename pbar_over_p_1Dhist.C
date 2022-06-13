#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TH1F.h"

void pbar_over_p_1Dhist() {
    //create eta binning
    const int etabin_num = 1;
    double etabins[etabin_num+1] = {-0.8,0.8};
    const int pbin_num = 5;
    double pbins[pbin_num+1] = {0.3,0.32,0.34,0.36,0.38,0.4};
    string p_and_pbar[2] = {"proton","antiproton"};
    
    TLatex *latex = new TLatex;
    latex->SetNDC(kTRUE);
    
    TFile *protonFile_p3_32 = new TFile("root_files/proton_fit_p3_32.root","read");
    TH1D* proton_number_p3_32;
    proton_number_p3_32 = (TH1D*) protonFile_p3_32->Get("totalProtons");
    
    TFile *protonFile_p32_34 = new TFile("root_files/proton_fit_p32_34.root","read");
    TH1D* proton_number_p32_34;
    proton_number_p32_34 = (TH1D*) protonFile_p32_34->Get("totalProtons");
    
    TFile *protonFile_p34_36 = new TFile("root_files/proton_fit_p34_36.root","read");
    TH1D* proton_number_p34_36;
    proton_number_p34_36 = (TH1D*) protonFile_p34_36->Get("totalProtons");
    
    TFile *protonFile_p36_38 = new TFile("root_files/proton_fit_p36_38.root","read");
    TH1D* proton_number_p36_38;
    proton_number_p36_38 = (TH1D*) protonFile_p36_38->Get("totalProtons");
    
    TFile *protonFile_p38_4 = new TFile("root_files/proton_fit_p38_4.root","read");
    TH1D* proton_number_p38_4;
    proton_number_p38_4 = (TH1D*) protonFile_p38_4->Get("totalProtons");
    
    
    TFile *antiprotonFile_p3_32 = new TFile("root_files/antiproton_fit_p3_32.root","read");
    TH1D* antiproton_number_p3_32;
    antiproton_number_p3_32 = (TH1D*) antiprotonFile_p3_32->Get("totalantiprotons");
    
    TFile *antiprotonFile_p32_34 = new TFile("root_files/antiproton_fit_p32_34.root","read");
    TH1D* antiproton_number_p32_34;
    antiproton_number_p32_34 = (TH1D*) antiprotonFile_p32_34->Get("totalantiprotons");
    
    TFile *antiprotonFile_p34_36 = new TFile("root_files/antiproton_fit_p34_36.root","read");
    TH1D* antiproton_number_p34_36;
    antiproton_number_p34_36 = (TH1D*) antiprotonFile_p34_36->Get("totalantiprotons");
    
    TFile *antiprotonFile_p36_38 = new TFile("root_files/antiproton_fit_p36_38.root","read");
    TH1D* antiproton_number_p36_38;
    antiproton_number_p36_38 = (TH1D*) antiprotonFile_p36_38->Get("totalantiprotons");
    
    TFile *antiprotonFile_p38_4 = new TFile("root_files/antiproton_fit_p38_4.root","read");
    TH1D* antiproton_number_p38_4;
    antiproton_number_p38_4 = (TH1D*) antiprotonFile_p38_4->Get("totalantiprotons");
    
    
    TH1D *divide_p3_32 = (TH1D*) antiproton_number_p3_32->Clone("divide_p3_32");
    divide_p3_32->Divide(proton_number_p3_32);
    
    TH1D *divide_p32_34 = (TH1D*) antiproton_number_p32_34->Clone("divide_p32_34");
    divide_p32_34->Divide(proton_number_p32_34);
    
    TH1D *divide_p34_36 = (TH1D*) antiproton_number_p34_36->Clone("divide_p34_36");
    divide_p34_36->Divide(proton_number_p34_36);
    
    TH1D *divide_p36_38 = (TH1D*) antiproton_number_p36_38->Clone("divide_p36_38");
    divide_p36_38->Divide(proton_number_p36_38);
    
    TH1D *divide_p38_4 = (TH1D*) antiproton_number_p38_4->Clone("divide_p38_4");
    divide_p38_4->Divide(proton_number_p38_4);
    double x[5], y[5];
    for (int n = 0; n < 5; n++){
        y[0]=divide_p3_32->GetBinContent(1);
        y[1]=divide_p32_34->GetBinContent(1);
        y[2]=divide_p34_36->GetBinContent(1);
        y[3]=divide_p36_38->GetBinContent(1);
        y[4]=divide_p38_4->GetBinContent(1);
        x[n] = pbins[n];
    }
    auto graph = new TGraph(5, x, y);
    graph->SetTitle("Ratio of antiprotons to protons vs momentum;Momentum [GeV];pbar/p");
    latex->DrawLatex(0.7,0.65,"#scale[0.8]{ATLAS #bf{Internal}}");
    latex->DrawLatex(0.7,0.60,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");
    latex->DrawLatex(0.7,0.55,"#scale[0.6]{#bf{0.3 < p < 0.4}}");
    graph->Draw("AP*");
    graph->SaveAs("pbar_p.pdf");
}
