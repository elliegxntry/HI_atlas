#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"

void histogram_fitting_dEdx() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int nbins_num = 1;
    const int pbins_num = 1;
    double nbins[nbins_num+1] = {25,60};
    double pbins[pbins_num+1] = {0.3,0.4};//,0.5,0.6,0.7};//,0.8,0.9,1.0};
    
    string p_and_pbar[2] = {"proton", "antiproton"};
    for (int m = 0; m < 2; m++) {
        TH1D* distributions[nbins_num*pbins_num];
        std::cout << "Loading " << p_and_pbar[m] << "s!" << std::endl;
        
        for (int p = 0; p < pbins_num; p++) {
            distributions[p] = new TH1D(Form("n25_60_p%d",p),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-2,5);
        }
        std::cout << "Reading file!" << std::endl;
        
        TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
        TTree *tree = (TTree*) file->Get("tree");

        int trk_n;
        int trk_400;
        float trk_pt[10000];
        float trk_eta[10000];
        float trk_q[10000];
        float dEdx[10000];
        float gaps[4];
        
        tree->SetBranchAddress("trk_n", &trk_n);
        tree->SetBranchAddress("trk_400", &trk_400);
        tree->SetBranchAddress("trk_pt", trk_pt);
        tree->SetBranchAddress("trk_eta", trk_eta);
        tree->SetBranchAddress("trk_q", trk_q);
        tree->SetBranchAddress("trk_dEdx", dEdx);
        
        int nevents = tree->GetEntries();
        for (int e = 0; e < nevents; e++) {
            if (e % 50000 == 0) std::cout << "Loading event " << e << std::endl;
            tree->GetEntry(e);
            for (int i=0; i < trk_n; i++){
                if (abs(trk_eta[i]) > 0.8) {continue;}
                double momentum = abs(trk_pt[i]*cosh(trk_eta[i]));
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
        for (int p = 0; p < pbins_num; p++) {
            distributions[p]->Scale(1,"width");
        }
        ostringstream filename;
        filename << "root_files/" << p_and_pbar[m] << "_data.root";
        TFile *outfile = new TFile(filename.str().c_str(),"recreate");
        for (int p = 0; p < pbins_num; p++) {
            distributions[p]->Write();
            outfile->Close();
            distributions[p] = nullptr;
        }
    }//m loop
}
