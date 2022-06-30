#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TMath.h"

void dEdx_updated_cuts() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int pbins_num = 9;
    double pbins[pbins_num+1] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    
    string p_and_pbar[2] = {"proton", "antiproton"};
    TFile *outfile = new TFile("root_files/pPb_data.root","recreate");
    TH1D* protons[2*pbins_num];
    TH1D* antiprotons[2*pbins_num];
    //std::cout << "Loading " << p_and_pbar[m] << "s!" << std::endl;
            
    for (int p = 0; p < pbins_num; p++) {
        protons[p] = new TH1D(Form("protons_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        antiprotons[p] = new TH1D(Form("antiprotons_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);

    }
    std::cout << "Reading file!" << std::endl;
    
    TFile *file = new TFile("root_files/user.srdas.29429835.OUTPUT.0.root","READ");
    TTree *tree = (TTree*) file->Get("tree");

    int trk_n;
    int trk_400;
    vector<float> *trk_pt = nullptr;
    vector<float> *trk_eta = nullptr;
    vector<float> *trk_q = nullptr;
    vector<float> *dEdx = nullptr;
    vector<float> *trk_nhits_dEdx = nullptr;
    vector<float> *trk_quality1 = nullptr;

    tree->SetBranchAddress("trk_n", &trk_n);
    tree->SetBranchAddress("trk_400", &trk_400);
    tree->SetBranchAddress("trk_pt", &trk_pt);
    tree->SetBranchAddress("trk_eta", &trk_eta);
    tree->SetBranchAddress("trk_q", &trk_q);
    tree->SetBranchAddress("trk_dEdx", &dEdx);
    tree->SetBranchAddress("trk_nhits_dEdx",&trk_nhits_dEdx);
    tree->SetBranchAddress("trk_quality1",&trk_quality1);

    int nevents = tree->GetEntries();
    for (int e = 0; e < nevents; e++) {
        if (e % 100000 == 0) std::cout << "Loading event " << e << std::endl;
        tree->GetEntry(e);
        for (int i=0; i < trk_n; i++){
            float eta = trk_eta->at(i);
            if (TMath::Abs(eta) > 0.3) {continue;} //restrict to -0.3 < eta < 0.3
            double momentum = TMath::Abs(trk_pt->at(i)*TMath::CosH(trk_eta->at(i)));
            if (trk_nhits_dEdx->at(i) < 3 && trk_nhits_dEdx->at(i) > 3) {continue;} // keep only the nhits = 3
            if (trk_quality1->at(i) == 0) {continue;} // limit to trk_quality1c == true
            for (int p = 0; p < pbins_num; p++) {
                if (pbins[p] <= momentum && momentum < pbins[p+1]) {
                    if (trk_q->at(i) > 0) {
                        protons[p]->Fill(TMath::Log(dEdx->at(i)));
                    }
                    else {
                        antiprotons[p]->Fill(TMath::Log(dEdx->at(i)));
                    }
                }
            }
        }
    }
    file->Close();
    outfile->cd();
    for (int p = 0; p < pbins_num; p++) {
        protons[p]->Scale(1,"width");
        protons[p]->Write();
        protons[p] = nullptr;
        antiprotons[p]->Scale(1,"width");
        antiprotons[p]->Write();
        antiprotons[p] = nullptr;
    }
    outfile->Close();
}
