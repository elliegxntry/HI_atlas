#include "TFile.h" // reads and writes root files
#include "TTree.h" // reads and writes TTrees
#include "TH1D.h" //Makes 1D histograms (double type)
#include "TF1.h" // creates and draws functions
#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text


void pT_dist_1Dhist() {
    TH1D* pT = new TH1D("pT", "Transverse momentum distribution;pT;N_{ch}",100,0,2);
    pT->Sumw2();
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","read");
    
    TTree *tree = (TTree*) file->Get("tree");
    
    int trk_n;
    float trk_pt[100000];
    
    tree->SetBranchAddress("trk_n", &trk_n); //the "&" says it's a single value
    tree->SetBranchAddress("trk_pt", trk_pt); //no "&" means it's an array of variables
    
    for (int e = 0; e < tree->GetEntries(); e++){
        tree->GetEntry(e);
        for (int i = 0; i < trk_n; i++) {
            if (pT[i] > 0.4) {continue;}
            pT->Fill(trk_pt[i]);
        }
    }
    TCanvas *tc = new TCanvas(); //canvas is called tc
    gStyle->SetOptStat(0);         //if 1: Shows information about histogram in box at top of plot, if 0: information not shown
    gPad->SetTicks();
    pT->SetLineWidth(2);
    pT->SetLineColor(kBlue);
    pT->SetMarkerStyle(4);
    pT->SetMarkerColor(kBlue);
    pT->Draw();
    
    tc->SaveAs("pT_Nch.pdf");
}

