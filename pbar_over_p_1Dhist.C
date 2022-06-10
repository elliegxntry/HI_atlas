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

void pbar_over_p_1Dhist() {
    //create eta binning
    const int etabin_num = 1;
    double etabins[etabin_num+1] = {-0.8,0.8};
    const int pbin_num = 1;
    double pbins[pbin_num+1] = {0.3,0.4};
    
    //create canvas and customize
    auto graph = new TGraph();
    gStyle->SetOptStat(0);
    //gPad->SetTicks();
    //gPad->SetLogy(1);
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);

    TFile *proton_file = new TFile("root_files/proton_fit.root","read");
    TFile *antiproton_file = new TFile("root_files/antiproton_fit.root","read");
    
    TH1I* protons[etabin_num*pbin_num];
    TH1I* antiprotons[etabin_num*pbin_num];
    
    //int plotCount = 0;
    
    for (int p = 0; p < 1; p++) {
        /*
        protons[p] = (TH1I*) proton_file->Get("totalProtons");
        int proton_num;
        proton_num = protons[p];
        antiprotons[p] = (TH1I*) antiproton_file->Get("totalAntiProtons");
        int antiproton_num;
        antiproton_num = antiprotons[p];
        int pbar_p;
        pbar_p = antiproton_num/proton_num;*/
        graph->SetPoint(0.,1.,1.);
        graph->SetMarkerColor(kBlue);
        graph->SetMarkerStyle(kFullCircle);
        //tc->Add(graph);
        graph->Draw();
        //tc->SaveAs("pbar_p_test_plot.pdf");
    }
}
