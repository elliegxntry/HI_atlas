#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

void pbar_over_p_1Dhist() {
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
    TTree *tree = (TTree*) file->Get("tree");
    
    float proton[10000];
    float proton_bar[10000];
    
    int trk_n;
    tree->SetBranchAddress("trk_n",&trk_n);
    float trk_pt[10000];
    tree->SetBranchAddress("trk_pt",trk_pt);
    float trk_eta[10000];
    tree->SetBranchAddress("trk_eta",trk_eta);
    float trk_q[10000];
    tree->SetBranchAddress("trk_q",trk_q);
    
    for (int entry = 0; entry < Tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
    }
    if (trk_q > 0){
        //XX figure out how to assign these to the array
    }
    if (trk_q < 0){
        XX figure out how to assign these to array
    }
    file->Close();
    
}
