#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TMath.h"

void dpmjet_dEdx_updated_cuts() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int pbins_num = 9;
    double pbins[pbins_num+1] = {.3,.4,.5,.6,.7,.8,.9,1.0,1.1,1.2};
    string files[19] = {"35._000001", "35._000002", "35._000003", "35._000004",
        "35._000005", "35._000006", "35._000007", "35._000008", "35._000009", "35._000010", "35._000011", "46._000001", "46._000002", "46._000003", "46._000004", "46._000005", "46._000006", "46._000007", "46._000008",
    };
    
    string p_and_pbar[2] = {"proton", "antiproton"};
    TFile *outfile = new TFile("root_files/gammaPb_dpmjet_data.root","recreate");
    TH1D* positive[2*pbins_num];
    TH1D* negative[2*pbins_num];
    TH1D* pions[2*pbins_num];
            
    for (int p = 0; p < pbins_num; p++) {
        positive[p] = new TH1D(Form("positive_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);
        negative[p] = new TH1D(Form("negative_p%d", p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,3);

    }
    for (int a = 0; a < 19; a++) {
        cout << "Reading file " << a << endl;
        ostringstream filenames;
        filenames << "root_files/gammaA_DPMJET_MC_forEllie/user.bseidlit.mbUPC_43.860046.Starlight_DPMJet_gammaA_trk25_photNegEta.e8077_s3469_r11509_myOutput.root/user.bseidlit.258018" << files[a] << ".myOutput.root";
        TFile *file = new TFile(filenames.str().c_str(),"READ");
        if (file->IsZombie()) {
          cout << "Error opening file number " << a << " - Skipping file!" << endl;
          continue;
        }
        TTree *tree = (TTree*) file->Get("tree");


        int trk_n;
        //int trk_400;
        float trk_pt[100000];
        float trk_eta[100000];
        //float trk_q[100000];
        float dEdx[100000];
        //int trk_nhits_dEdx = nullptr;
        bool trk_quality1[100000];
        int truth_pid[100000];

        tree->SetBranchAddress("trk_n", &trk_n);
        //tree->SetBranchAddress("trk_400", &trk_400);
        tree->SetBranchAddress("trk_pt", trk_pt);
        tree->SetBranchAddress("trk_eta", trk_eta);
        //tree->SetBranchAddress("trk_q", trk_q);
        tree->SetBranchAddress("trk_dEdx", dEdx);
        //tree->SetBranchAddress("trk_nhits_dEdx",&trk_nhits_dEdx);
        tree->SetBranchAddress("trk_quality1", trk_quality1);
        tree->SetBranchAddress("truth_pid", truth_pid);

        int nevents = tree->GetEntries();
        for (int e = 0; e < nevents; e++) {
            if (e % 10000 == 0) std::cout << "Loading event " << e << std::endl;
            tree->GetEntry(e);
            for (int i=0; i < trk_n; i++){
                float eta = trk_eta[i];
                if (TMath::Abs(eta) > 0.3) {continue;} //restrict to -0.3 < eta < 0.3
                double momentum = TMath::Abs(trk_pt[i]*TMath::CosH(trk_eta[i]));
                //if (trk_nhits_dEdx[i]!=3) {continue;} // keep only the nhits = 3
                if (trk_quality1[i] == FALSE) {continue;} // limit to trk_quality1c == true
                for (int p = 0; p < pbins_num; p++) {
                    if (pbins[p] <= momentum && momentum < pbins[p+1]) {
                        if (truth_pid[i] > 0) {
                            positive[p]->Fill(TMath::Log(dEdx[i]));
                        }
                        else {
                            negative[p]->Fill(TMath::Log(dEdx[i]));
                        }
                        if (TMath::Abs(truth_pid) == 211) {
                            pions[p]->Fill(TMath::Log(dEdx[i]));
                        }
                    }
                }
            }
        }
        file->Close();
    }
    outfile->cd();
    for (int p = 0; p < pbins_num; p++) {
        positive[p]->Scale(1,"width");
        positive[p]->Write();
        positive[p] = nullptr;
        negative[p]->Scale(1,"width");
        negative[p]->Write();
        negative[p] = nullptr;
    }
    outfile->Close();
}
