#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text
#include "TFile.h"    //allows you to load root files
#include "TH1D.h"     //makes 1D double histograms
#include "TF1.h"    //makes functions
#include "TPad.h"
#include "TH1F.h"

void dEdx_updated_plots_individual() {
    const int pbins_num = 9;
    double pbins[pbins_num+1] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"proton", "antiproton"};
    
    TFile *file = new TFile("root_files/pPb_data.root","read");
    TH1D* distributions[2*pbins_num];
    TH1D* protons[pbins_num];
    TH1D* antiprotons[pbins_num];
    TH1D* ratio[pbins_num];

    for (int m = 0; m < 2; m++) {
        TCanvas *tc_data = new TCanvas();
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        gPad->SetLogy(1);
        TLatex *latex_data = new TLatex();
        latex_data->SetNDC(kTRUE);
        
        TFile *file = new TFile("root_files/pPb_data.root","read");
        TH1D* distributions[2*pbins_num];

        for (int p = 0; p < pbins_num; p++){
            distributions[p] = (TH1D*) file->Get(Form("%ss_p%d", p_and_pbar[m].c_str(), p));
            distributions[p]->SetLineWidth(2);
            distributions[p]->SetLineColor(kBlack);
            distributions[p]->Draw();
            latex_data->DrawLatex(0.68,0.7,"#scale[0.8]{ATLAS #bf{Internal}}");
            latex_data->DrawLatex(0.68,0.66,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
            ostringstream p_label;
            p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
            latex_data->DrawLatex(0.68,0.6, p_label.str().c_str());
            ostringstream particle_type;
            particle_type << "#scale[0.6]{#bf{" <<  p_and_pbar[m] << "s}}";
            latex_data->DrawLatex(0.68,0.55,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
            latex_data->DrawLatex(0.68,0.5, particle_type.str().c_str());
            //tc_data->SaveAs(Form("dEdx_histograms/just_data/m%d_p%d.pdf",m,p));
        }
    }
    
    TCanvas *tc = new TCanvas("tc","New Canvas",50,50,700,600);
    double padSplit = 0.32;
    double leftMargin = 0.08;
    double bottomMargin = 0.15;
    double rightMargin = 0.02;
    double topMargin = 0.02;
    double height = 900;
    double width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);

    TPad* pads[2];
    pads[0] = new TPad("pad0","",0.0,padSplit,1.0,1.0);
    pads[0]->SetTopMargin(topMargin);
    pads[0]->SetRightMargin(rightMargin);
    pads[0]->SetLeftMargin(leftMargin);
    pads[0]->SetBottomMargin(0.001);
    
    tc->cd();
    pads[0]->Draw("same");
    pads[0]->cd();
    tc->cd();

    pads[1] = new TPad("pad1","",0.0,0.0,1.0,padSplit);
    pads[1]->SetTopMargin(0.001);
    pads[1]->SetRightMargin(rightMargin);
    pads[1]->SetLeftMargin(leftMargin);
    pads[1]->SetBottomMargin(bottomMargin);

    //tc->cd();
    pads[1]->Draw("same");
    pads[1]->cd();
    tc->cd();
    
    gStyle->SetOptStat(0);
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);
    for (int p = 0; p < pbins_num; p++) {
        pads[0]->cd();
        gPad->SetTicks(1);
        gPad->SetLogy();
        
        TFile *file = new TFile("root_files/pPb_data.root","read");
        protons[p] = (TH1D*) file->Get(Form("protons_p%d",p));
        antiprotons[p] = (TH1D*) file->Get(Form("antiprotons_p%d",p));
        
        //draw the distributions overlayed
        protons[p]->SetLineWidth(2);
        protons[p]->SetLineColor(kRed);
        protons[p]->Draw();
        antiprotons[p]->SetLineWidth(2);
        antiprotons[p]->SetLineColor(kBlue);
        antiprotons[p]->Draw("same");
        latex->DrawLatex(0.75,0.73,"#scale[0.8]{ATLAS #bf{Internal}}");
        latex->DrawLatex(0.75,0.68,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
        ostringstream p_label;
        p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
        latex->DrawLatex(0.75,0.63, p_label.str().c_str());
        latex->DrawLatex(0.75,0.58,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
        TLegend *legend = new TLegend(0.7,0.8,0.93,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(protons[p], "Positively Charged Particles","l");
        legend->AddEntry(antiprotons[p], "Negatively Charged Particles", "l");
        legend->Draw("same");
        
        tc->cd();
        pads[1] = new TPad("pad1","",0.0,0.0,1.0,padSplit);
        pads[1]->SetTopMargin(0.001);
        pads[1]->SetRightMargin(rightMargin);
        pads[1]->SetLeftMargin(leftMargin);
        pads[1]->SetBottomMargin(bottomMargin);

        //tc->cd();
        pads[1]->Draw("same");
        pads[1]->cd();
        tc->cd();
        
        pads[1]->cd();
        //pads[1]->Clear();
        gPad->SetTicks(1);
        gPad->SetLogy(0);
        ratio[p] = (TH1D*) antiprotons[p]->Clone(Form("antiprotons_%d",p));
        ratio[p]->Divide(protons[p]);
        ratio[p]->SetLineWidth(2);
        ratio[p]->SetLineColor(kBlack);
        ratio[p]->SetMaximum(1.5);
        ratio[p]->SetMinimum(0.8);
        ratio[p]->Draw("same");
        TLine *line = new TLine(-2, 1, 5, 1);
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->Draw("same");
        ratio[p]->GetYaxis()->SetTitleSize(0.06);
        ratio[p]->GetYaxis()->SetTitle("-/+");
        
        ostringstream histName;
        histName << "dEdx_histograms/quality1/p" << p << "_pPb_quality1cut.pdf";
        tc->SaveAs(histName.str().c_str());
        pads[1]->Clear();
    }
}
