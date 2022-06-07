#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

void pT_vs_p_2Dhist() {
    //load the data
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","read");
    TTree *tree = (TTree*) file->Get("tree");
    const int eta_bin_num = 10;
    //double eta_bins[eta_bin_num+1] = {-2.5,2.5};
    double eta_bins[eta_bin_num+1] = {-2.5,-2,-1.5,-1,-0.5,0.,0.5,1.0,1.5,2.0,2.5};
    
    //get variables from tree
    int trk_n;
    tree->SetBranchAddress("trk_n", &trk_n);
    
    float trk_pt[10000];
    tree->SetBranchAddress("trk_pt", trk_pt);
    
    float trk_eta[10000];
    tree->SetBranchAddress("trk_eta", trk_eta);
    
    // Set up plots
    TCanvas *tc = new TCanvas();
    gStyle->SetOptStat(0);
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);
    gPad->SetTicks();
    gPad->SetLogz();
    gStyle->SetPalette(62);
    TH2D* pTp_hist[eta_bin_num];
    for (int i=0; i < eta_bin_num;i++) {
        ostringstream histName;
        histName << "pTp_eta" << i;
        pTp_hist[i] = new TH2D(histName.str().c_str(), "Transverse momentum versus total momentum;pT [GeV/c];p [GeV/c]",100,0,50,100,0,50);
    }
    
    for (int entries = 0; entries < tree->GetEntries(); entries++){
        tree->GetEntry(entries);
        for (int tracks = 0; tracks < trk_n; tracks++){
            double momentum = abs(trk_pt[tracks]*cosh(trk_eta[tracks]));
            for (int eta = 0; eta < eta_bin_num; eta++) {
                if (trk_eta[tracks]>eta_bins[eta+1]) continue;
                if (trk_eta[tracks]<eta_bins[eta]) continue;
                pTp_hist[eta]->Fill(trk_pt[tracks],momentum);
            }
        }
    }
    for (int eta = 0; eta < eta_bin_num; eta++) {
        TLegend *legend = new TLegend(0.68,0.3,0.85,0.35);
        pTp_hist[eta]->Draw("colz");
        latex->DrawLatex(0.68,0.50,"#scale[0.8]{ATLAS #bf{Internal}}");   //Always include "ATLAS INTERNAL" on your plots
        latex->DrawLatex(0.68,0.45,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");
        ostringstream eta_label;
        eta_label << "#scale[0.6]{#bf{" << eta_bins[eta] << " < #eta < " << eta_bins[eta+1] << "}}";
        latex->DrawLatex(0.68, 0.40,eta_label.str().c_str());
        ostringstream histName;
        histName << "pT_vs_p_plots/pTp_hist" << "_eta" << eta << ".pdf";
        TF1 *yx;
        yx = new TF1("yx","x",0,50);
        yx->SetLineWidth(2);
        yx->SetLineColor(kBlue);
        yx->SetLineStyle(3);
        legend->AddEntry(yx, "pT = p","l");
        yx->Draw("same");
        legend->Draw("same");
        tc->SaveAs(histName.str().c_str());
    }
}
