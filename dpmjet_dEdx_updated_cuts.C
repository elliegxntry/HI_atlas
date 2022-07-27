#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TMath.h"

void dpmjet_dEdx_updated_cuts() {
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {.4,.5,.6,.7,.8,.9,1.0,1.1,1.2};

    string files[2] = {"root_files/user.srdas.061522.DPMJet_gammaA_photNegEta.v1_OUTPUT.root", "root_files/user.srdas.061522.DPMJet_gammaA_photPosEta.v1_OUTPUT.root"};
    
    string p_and_pbar[2] = {"positive", "negative"};
    TFile *outfile = new TFile("root_files/gammaPb_dpmjet.root","recreate");
    TH1D* total[2*pbins_num];
    TH1D* positive[2*pbins_num];
    TH1D* negative[2*pbins_num];
    TH1D* pions[pbins_num*2];
    TH1D* kaons[pbins_num*2];
    TH1D* protons[pbins_num*2];
    float xshift[pbins_num+1] = {0.107718, 0.0918337, 0.0816148, 0.100812, 0.106426, 0.11712, 0.106271, 0.10289};
    
    for (int p = 0; p < pbins_num; p++) {
        total[p] = new TH1D(Form("total_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        positive[p] = new TH1D(Form("positive_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        negative[p] = new TH1D(Form("negative_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        pions[p] = new TH1D(Form("pions_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        kaons[p] = new TH1D(Form("kaons_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        protons[p] = new TH1D(Form("protons_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
    }
    
    for (int a = 0; a < 2; a++) {
        cout << "Reading file " << a << endl;
        TFile *file = new TFile(files[a].c_str(),"READ");
//        TFile *file = new TFile("root_files/user.srdas.061522.DPMJet_gammaA_photNegEta.v1_OUTPUT.root","READ");
        TTree *tree = (TTree*) file->Get("tree");

        int trk_n;
        vector<float> *trk_pt = nullptr;
        vector<float> *trk_eta = nullptr;
        vector<float> *trk_q = nullptr;
        vector<float> *dEdx = nullptr;
        vector<int> *trk_nhits_dEdx = nullptr;
        vector<float> *trk_quality1 = nullptr;
        vector<int> *trk_truthpdgID = nullptr;

        tree->SetBranchAddress("trk_n", &trk_n);
        tree->SetBranchAddress("trk_pt", &trk_pt);
        tree->SetBranchAddress("trk_q", &trk_q);
        tree->SetBranchAddress("trk_eta", &trk_eta);
        tree->SetBranchAddress("trk_dEdx", &dEdx);
        tree->SetBranchAddress("trk_nhits_dEdx", &trk_nhits_dEdx);
        tree->SetBranchAddress("trk_quality1", &trk_quality1);
        tree->SetBranchAddress("trk_truthpdgID", &trk_truthpdgID);

        int nevents = tree->GetEntries();
        for (int e = 0; e < nevents; e++) {
            if (e % 100000 == 0) cout << "Loading event " << e << endl;
            tree->GetEntry(e);
            for (int i=0; i < trk_n; i++){
                float eta = trk_eta->at(i);
                if (TMath::Abs(eta) > 0.3) {continue;} //restrict to -0.3 < eta < 0.3
                double momentum = TMath::Abs(trk_pt->at(i)*TMath::CosH(trk_eta->at(i)));
                if (trk_quality1->at(i) == FALSE) {continue;} // limit to trk_quality1c == true
                if (trk_nhits_dEdx->at(i)!=3) {continue;}
                if (a == 0) {trk_eta->at(i) = -1*trk_eta->at(i);}
                for (int p = 0; p < pbins_num; p++) {
                    if ((pbins[p]) <= momentum && momentum < (pbins[p+1])) {
                        total[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));
                        if (trk_q->at(i) > 0) {
                            positive[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));
                        }
                        else {
                            negative[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));                        }
                        if ((TMath::Abs(trk_truthpdgID->at(i))) == 211) {
                            pions[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));
                        }
                        if ((TMath::Abs(trk_truthpdgID->at(i))) == 321) {
                            kaons[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));
                        }
                        if ((TMath::Abs(trk_truthpdgID->at(i))) == 2212) {
                            protons[p]->Fill(TMath::Log(dEdx->at(i)-xshift[p]));
                        }
                    }
                }
            }
        }
        file->Close();
    }
    outfile->cd();
    for (int p = 0; p < pbins_num; p++) {
        total[p]->Scale(1,"width");
        total[p]->Write();
        total[p] = nullptr;
        positive[p]->Scale(1,"width");
        positive[p]->Write();
        positive[p] = nullptr;
        negative[p]->Scale(1,"width");
        negative[p]->Write();
        negative[p] = nullptr;
        pions[p]->Scale(1,"width");
        pions[p]->Write();
        pions[p] = nullptr;
        kaons[p]->Scale(1,"width");
        kaons[p]->Write();
        kaons[p] = nullptr;
        protons[p]->Scale(1,"width");
        protons[p]->Write();
        protons[p] = nullptr;
    }
    outfile->Close();
}
