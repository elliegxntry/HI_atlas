#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH2D.h"

void pbar_over_p_1Dhist() {
    //load in the data
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
    TTree *tree = (TTree*) file->Get("tree");
    
    //set addresses to tree variables and initialize arrays to store them
    int trk_n;
    tree->SetBranchAddress("trk_n",&trk_n);
    float trk_pt[10000];
    tree->SetBranchAddress("trk_pt",trk_pt);
    float trk_eta[10000];
    tree->SetBranchAddress("trk_eta",trk_eta);
    float trk_q[10000];
    tree->SetBranchAddress("trk_q",trk_q);
    
    //create arrays for the proton and anti proton data
    float proton[10000];
    float proton_bar[10000];
    
    //create eta binning
    const int eta_bin_num = 1;
    double eta_bins[eta_bin_num+1] = {-0.8,0.8};
    TH1D* histogram[eta_bin_num];
    
    //create a histogram for each bin
    for (int eta = 0; eta < eta_bin_num; eta++){
        ostringstream histName;
        histName << "pbar_over_p_vs_pt_eta" << eta;
        histogram[eta] = new TH1D(histName.str().c_str(),"#frac{#bar{p}}{p} vs pT;pT [GeV/c];#frac{#bar{p}}{p}",100,-1,3);
        histogram[eta]->Sumw2();
    }
    
    std::cout << "Reading the file!" << std::endl;
    
    //get each entry
    for (int entry = 0; entry < tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
    }
    file->Close();
    
}
