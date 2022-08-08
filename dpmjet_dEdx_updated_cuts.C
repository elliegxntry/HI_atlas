#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TMath.h"

void dpmjet_dEdx_updated_cuts() {
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {.4,.5,.6,.7,.8,.9,1.0,1.1,1.2};

    string files[4] = {"root_files/user.srdas.061522.DPMJet_gammaA_photNegEta.v1_OUTPUT.root", "root_files/user.srdas.061522.DPMJet_gammaA_trk2_photNegEta.v1_OUTPUT.root", "root_files/user.srdas.061522.DPMJet_gammaA_photPosEta.v1_OUTPUT.root", "root_files/user.srdas.061522.DPMJet_gammaA_trk2_photPosEta.v1_OUTPUT.root"};
    
    string p_and_pbar[2] = {"positive", "negative"};
    TFile *outfile = new TFile("root_files/gammaPb_dpmjet.root","recreate");
    TH1D* total[2*pbins_num];
    TH1D* positive[2*pbins_num];
    TH1D* negative[2*pbins_num];
    TH1D* pos_pions[pbins_num*2];
    TH1D* pos_kaons[pbins_num*2];
    TH1D* pos_protons[pbins_num*2];
    TH1D* neg_pions[pbins_num*2];
    TH1D* neg_kaons[pbins_num*2];
    TH1D* neg_protons[pbins_num*2];
    
    float xshift[pbins_num+1] = {0.108957, 0.0909291, 0.0816315, 0.0975523, 0.104969, 0.107637, 0.110375, 0.101998};
    
    for (int p = 0; p < pbins_num; p++) {
        total[p] = new TH1D(Form("total_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        positive[p] = new TH1D(Form("positive_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        negative[p] = new TH1D(Form("negative_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))", 100, -2, 3);
        pos_pions[p] = new TH1D(Form("pions_m0_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        pos_kaons[p] = new TH1D(Form("kaons_m0_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        pos_protons[p] = new TH1D(Form("protons_m0_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        neg_pions[p] = new TH1D(Form("pions_m1_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        neg_kaons[p] = new TH1D(Form("kaons_m1_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        neg_protons[p] = new TH1D(Form("protons_m1_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
    }
            
    for (int a = 0; a < 4; a++) {
        cout << "Reading file " << a << endl;
        TFile *file = new TFile(files[a].c_str(),"READ");
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
                if (a == 0 || a == 1) {trk_eta->at(i) = -1*trk_eta->at(i);}
                for (int p = 0; p < pbins_num; p++) {
                    if ((pbins[p]) <= momentum && momentum < (pbins[p+1])) {
                        total[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                        if (trk_q->at(i) > 0) {
                            positive[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 211) {
                                pos_pions[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 321) {
                                pos_kaons[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 2212) {
                                pos_protons[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
                        }
                        else {
                            negative[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 211) {
                                neg_pions[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 321) {
                                neg_kaons[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
                            if ((TMath::Abs(trk_truthpdgID->at(i))) == 2212) {
                                neg_protons[p]->Fill(TMath::Log(dEdx->at(i))-xshift[p]);
                            }
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
        pos_pions[p]->Scale(1,"width");
        pos_pions[p]->Write();
        pos_pions[p] = nullptr;
        pos_kaons[p]->Scale(1,"width");
        pos_kaons[p]->Write();
        pos_kaons[p] = nullptr;
        pos_protons[p]->Scale(1,"width");
        pos_protons[p]->Write();
        pos_protons[p] = nullptr;
        neg_pions[p]->Scale(1,"width");
        neg_pions[p]->Write();
        neg_pions[p] = nullptr;
        neg_kaons[p]->Scale(1,"width");
        neg_kaons[p]->Write();
        neg_kaons[p] = nullptr;
        neg_protons[p]->Scale(1,"width");
        neg_protons[p]->Write();
        neg_protons[p] = nullptr;
    }
    outfile->Close();
}
