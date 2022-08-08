#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TMath.h"

void scaling() {
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"positive", "negative"};
    
    TFile *data_file = new TFile("root_files/pPb_data.root","read");
    TFile *mc_file = new TFile("root_files/gammaPb_dpmjet.root","UPDATE");
//    TFile *mc_file = new TFile("root_files/gammaPb_dpmjet.root","read");
    TFile *outfile = new TFile("root_files/scaling_factors.root", "recreate");
    
    TH1D* data_distribution[pbins_num];
    TH1D* mc_distribution[pbins_num];
    TH1D* pos_pions[pbins_num];
    TH1D* pos_kaons[pbins_num];
    TH1D* pos_protons[pbins_num];
    TH1D* neg_pions[pbins_num];
    TH1D* neg_kaons[pbins_num];
    TH1D* neg_protons[pbins_num];
    TH1D* scale[pbins_num*2];

    for (int m = 0; m < 2; m++) {
        cout << "Charge bin: " << m << "\n";
        for (int p = 0; p < pbins_num; p++) {
            cout << "Momentum bin: " << p << "\n";
            data_distribution[p] = (TH1D*) data_file->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            mc_distribution[p] = (TH1D*) mc_file->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            pos_pions[p] = (TH1D*) mc_file->Get(Form("pions_m0_p%d", p));
            pos_kaons[p] = (TH1D*) mc_file->Get(Form("kaons_m0_p%d", p));
            pos_protons[p] = (TH1D*) mc_file->Get(Form("protons_m0_p%d", p));
            neg_pions[p] = (TH1D*) mc_file->Get(Form("pions_m1_p%d", p));
            neg_kaons[p] = (TH1D*) mc_file->Get(Form("kaons_m1_p%d", p));
            neg_protons[p] = (TH1D*) mc_file->Get(Form("protons_m1_p%d", p));
            
            //----------------------
            int pi_binnum_data;
            pi_binnum_data = data_distribution[p]->FindBin(0);
            float pi_binval_data;
            pi_binval_data = data_distribution[p]->GetBinContent(pi_binnum_data);
//            cout << pi_binval_data << "\n";

            int pi_binnum_mc;
            pi_binnum_mc = mc_distribution[p]->FindBin(0);
            float pi_binval_mc;
            pi_binval_mc = mc_distribution[p]->GetBinContent(pi_binnum_mc);
//            cout << pi_binval_mc << "\n";

            cout << pi_binval_data/pi_binval_mc << ", ";

            pos_pions[p] = (TH1D*) mc_distribution[p]->Clone(Form("pions_m%d_p%d", m, p));
            pos_pions[p]->GetXaxis()->SetLimits(-2, 0.5);
            pos_pions[p]->Scale(pi_binval_data/pi_binval_mc);
            pos_pions[p]->Write();

            //----------------------

            int ka_binnum_data;
            ka_binnum_data = data_distribution[p]->FindBin(0);
            float ka_binval_data;
            ka_binval_data = data_distribution[p]->GetBinContent(ka_binnum_data);
//            cout << ka_binval_data << "\n";

            int ka_binnum_mc;
            ka_binnum_mc = mc_distribution[p]->FindBin(0);
            float ka_binval_mc;
            ka_binval_mc = mc_distribution[p]->GetBinContent(ka_binnum_mc);
//            cout << ka_binval_mc << "\n";

            cout << ka_binval_data/ka_binval_mc << ", ";

            pos_kaons[p] = (TH1D*) mc_distribution[p]->Clone(Form("kaons_m%d_p%d", m, p));
            pos_kaons[p]->GetXaxis()->SetLimits(-2, 0.5);
            pos_kaons[p]->Scale(ka_binval_data/ka_binval_mc);
            pos_kaons[p]->Write();

            //----------------------

            int pr_binnum_data;
            pr_binnum_data = data_distribution[p]->FindBin(0);
            float pr_binval_data;
            pr_binval_data = data_distribution[p]->GetBinContent(pr_binnum_data);
//            cout << pr_binval_data << "\n";

            int pr_binnum_mc;
            pr_binnum_mc = mc_distribution[p]->FindBin(0);
            float pr_binval_mc;
            pr_binval_mc = mc_distribution[p]->GetBinContent(pr_binnum_mc);
//            cout << pr_binval_mc << "\n";

            cout << pr_binval_data/pr_binval_mc << ", ";

            pos_protons[p] = (TH1D*) mc_distribution[p]->Clone(Form("protons_m%d_p%d", m, p));
            pos_protons[p]->GetXaxis()->SetLimits(-2, 0.5);
            pos_protons[p]->Scale(pr_binval_data/pr_binval_mc);
            pos_protons[p]->Write();

            mc_distribution[p]->Scale(binval_data/binval_mc);
            mc_file->cd();
            mc_distribution[p]->Write("",TObject::kOverwrite);
//
//            outfile->cd();
//            scale[p] = new TH1D(Form("scale_m%d_p%d", m, p),"scale", 1, 0, 2);
//            scale[p]->SetBinContent(1,binval_data/binval_mc);
//            scale[p]->Write();
        }
    }
}
