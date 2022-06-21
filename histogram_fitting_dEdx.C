#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"

void histogram_fitting_dEdx() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int nbins_num = 1;
    const int pbins_num = 15;
    double nbins[nbins_num+1] = {25,60};
    double pbins[pbins_num+1] = {0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6};
    
    string p_and_pbar[2] = {"proton", "antiproton"};
    TFile *outfile = new TFile("root_files/data.root","recreate");
    for (int m = 0; m < 2; m++) {
        TH1D* distributions[nbins_num*pbins_num];
        std::cout << "Loading " << p_and_pbar[m] << "s!" << std::endl;
                
        for (int p = 0; p < pbins_num; p++) {
            ostringstream histName;
            histName << p_and_pbar[m] << "_n25_60_p" << p;
            distributions[p] = new TH1D(Form("%s_n25_60_p%d", p_and_pbar[m].c_str(), p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,5);
        }
        std::cout << "Reading file!" << std::endl;
        
        TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
        TTree *tree = (TTree*) file->Get("tree");

        int trk_n;
        int trk_400;
        float eta_sign;
        float trk_pt[10000];
        float trk_eta[10000];
        float trk_q[10000];
        float dEdx[10000];
        float gaps[4];
        
        tree->SetBranchAddress("trk_n", &trk_n);
        tree->SetBranchAddress("trk_400", &trk_400);
        tree->SetBranchAddress("eta_sign", &eta_sign);
        tree->SetBranchAddress("trk_pt", trk_pt);
        tree->SetBranchAddress("trk_eta", trk_eta);
        tree->SetBranchAddress("trk_q", trk_q);
        tree->SetBranchAddress("trk_dEdx", dEdx);
        
        int nevents = tree->GetEntries();
        //int nevents = 10000;
        for (int e = 0; e < nevents; e++) {
            if (e % 100000 == 0) std::cout << "Loading event " << e << std::endl;
            tree->GetEntry(e);
            for (int i=0; i < trk_n; i++){
                float eta = eta_sign*trk_eta[i];
                if (abs(eta) > 0.8) {continue;}
                double momentum = abs(trk_pt[i]*cosh(eta));
                if (m == 0 && (trk_q[i] < 0)) {continue;} //keep only protons
                if (m == 1 && (trk_q[i] > 0)) {continue;} //keep only antiprotons
                for (int p = 0; p < pbins_num; p++) {
                    if (nbins[0] < trk_400 && trk_400 <= nbins[1] && pbins[p] <= momentum && momentum < pbins[p+1]) {
                        distributions[p]->Fill(log(dEdx[i]));
                    }
                }
            }
        }
        file->Close();
        outfile->cd();
        for (int p = 0; p < pbins_num; p++) {
            distributions[p]->Scale(1,"width");
            distributions[p]->Write();
            distributions[p] = nullptr;
        }
    }//m loop
    outfile->Close();
}
