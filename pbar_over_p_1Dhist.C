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
    const int pbin_num = 6;
    double pbins[pbin_num+1] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5};
    string p_and_pbar[2] = {"proton","antiproton"};
    TCanvas *tc = new TCanvas();
    
    TLatex *latex = new TLatex;
    latex->SetNDC(kTRUE);
    
    TFile *protonFile_p0 = new TFile("root_files/proton_p0_fit.root","read");
    TH1D* proton_number_p0;
    proton_number_p0 = (TH1D*) protonFile_p0->Get("totalprotons_p0");
    
    TFile *protonFile_p1 = new TFile("root_files/proton_p1_fit.root","read");
    TH1D* proton_number_p1;
    proton_number_p1 = (TH1D*) protonFile_p1->Get("totalprotons_p1");
    
    TFile *protonFile_p2 = new TFile("root_files/proton_p2_fit.root","read");
    TH1D* proton_number_p2;
    proton_number_p2 = (TH1D*) protonFile_p2->Get("totalprotons_p2");
    
    TFile *protonFile_p3 = new TFile("root_files/proton_p3_fit.root","read");
    TH1D* proton_number_p3;
    proton_number_p3 = (TH1D*) protonFile_p3->Get("totalprotons_p3");
    /*
    TFile *protonFile_p4 = new TFile("root_files/proton_p4_fit.root","read");
    TH1D* proton_number_p4;
    proton_number_p4 = (TH1D*) protonFile_p4->Get("totalprotons_p4");*/
    
    
    TFile *antiprotonFile_p0 = new TFile("root_files/antiproton_p0_fit.root","read");
    TH1D* antiproton_number_p0;
    antiproton_number_p0 = (TH1D*) antiprotonFile_p0->Get("totalantiprotons_p0");
    
    TFile *antiprotonFile_p1 = new TFile("root_files/antiproton_p1_fit.root","read");
    TH1D* antiproton_number_p1;
    antiproton_number_p1 = (TH1D*) antiprotonFile_p1->Get("totalantiprotons_p1");
    
    TFile *antiprotonFile_p2 = new TFile("root_files/antiproton_p2_fit.root","read");
    TH1D* antiproton_number_p2;
    antiproton_number_p2 = (TH1D*) antiprotonFile_p2->Get("totalantiprotons_p2");
    
    TFile *antiprotonFile_p3 = new TFile("root_files/antiproton_p3_fit.root","read");
    TH1D* antiproton_number_p3;
    antiproton_number_p3 = (TH1D*) antiprotonFile_p3->Get("totalantiprotons_p3");
    /*
    TFile *antiprotonFile_p4 = new TFile("root_files/antiproton_p4_fit.root","read");
    TH1D* antiproton_number_p4;
    antiproton_number_p4 = (TH1D*) antiprotonFile_p4->Get("totalantiprotons_p4");*/
    
    TH1D *divide_p0 = (TH1D*) antiproton_number_p0->Clone("divide_p0");
    divide_p0->Divide(proton_number_p0);
    
    TH1D *divide_p1 = (TH1D*) antiproton_number_p1->Clone("divide_p1");
    divide_p1->Divide(proton_number_p1);
    
    TH1D *divide_p2 = (TH1D*) antiproton_number_p2->Clone("divide_p2");
    divide_p2->Divide(proton_number_p2);
    
    TH1D *divide_p3 = (TH1D*) antiproton_number_p3->Clone("divide_p3");
    divide_p3->Divide(proton_number_p3);
    /*
    TH1D *divide_p4 = (TH1D*) antiproton_number_p4->Clone("divide_p4");
    divide_p4->Divide(proton_number_p4);*/
    double x[5], y[5];
    for (int n = 0; n < 5; n++){
        y[0]=divide_p0->GetBinContent(1);
        y[1]=divide_p1->GetBinContent(1);
        y[2]=divide_p2->GetBinContent(1);
        y[3]=divide_p3->GetBinContent(1);
        //y[4]=divide_p4->GetBinContent(1);
        x[n] = pbins[n];
    }
    auto graph = new TGraph(4, x, y);
    graph->SetTitle("Ratio of antiprotons to protons vs momentum;Momentum [GeV];pbar/p");
    latex->DrawLatex(0.7,0.65,"#scale[0.8]{ATLAS #bf{Internal}}");
    latex->DrawLatex(0.7,0.60,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");
    latex->DrawLatex(0.7,0.55,"#scale[0.6]{#bf{-0.8 < #eta < 0.8}}");
    graph->Draw("AP*");
    tc->SaveAs("pbar_p.pdf");
}
