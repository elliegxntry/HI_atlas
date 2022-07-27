#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text
#include "TFile.h"    //allows you to load root files
#include "TH1D.h"     //makes 1D double histograms
#include "TF1.h"    //makes functions
#include "TPad.h"
#include "TH1F.h"

void dpmjet_updated_plots() {
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {.4,.5,.6,.7,.8,.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"positive", "negative"};
    
    TFile *file = new TFile("root_files/gammaPb_dpmjet.root","read");
    TH1D* distributions[2*pbins_num];
    TH1D* total[2*pbins_num];
    TH1D* positive[pbins_num];
    TH1D* negative[pbins_num];
    TH1D* ratio[pbins_num];
    TH1D* pions[2*pbins_num];
    TH1D* kaons[pbins_num*2];
    TH1D* protons[pbins_num*2];

    for (int m = 0; m < 2; m++) {
        TCanvas *tc_data = new TCanvas();
        gStyle->SetOptStat(0);
        gPad->SetTicks();
        gPad->SetLogy(1);
        TLatex *latex_data = new TLatex();
        latex_data->SetNDC(kTRUE);
        
        TFile *file = new TFile("root_files/gammaPb_dpmjet.root","read");
        TH1D* distributions[2*pbins_num];

        for (int p = 0; p < pbins_num; p++){
            distributions[p] = (TH1D*) file->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            distributions[p]->SetLineWidth(2);
            distributions[p]->SetLineColor(kBlack);
            distributions[p]->Draw();
            latex_data->DrawLatex(0.68,0.81,"#scale[0.8]{ATLAS #bf{Internal}}");
            latex_data->DrawLatex(0.68,0.77,"#scale[0.6]{#bf{DPMJET gamma+Pb}}");
            ostringstream p_label;
            p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
            latex_data->DrawLatex(0.68,0.73, p_label.str().c_str());
            ostringstream particle_type;
            particle_type << "#scale[0.6]{#bf{" <<  p_and_pbar[m] << "s}}";
            latex_data->DrawLatex(0.68,0.69,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
            latex_data->DrawLatex(0.68,0.65, particle_type.str().c_str());
            tc_data->SaveAs(Form("dEdx_histograms/just_data/dpmjet_m%d_p%d.pdf",m,p));
        }
    }
    
    TCanvas *tc_total = new TCanvas();
    tc_total->cd();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(1);
    TLatex *latex_total = new TLatex();
    latex_total->SetNDC(kTRUE);
    
    TFile *total_file = new TFile("root_files/gammaPb_dpmjet.root","read");
    for (int p = 0; p < pbins_num; p++) {
        total[p] = (TH1D*) total_file->Get(Form("total_p%d", p));
        total[p]->SetLineColor(kBlack);
        total[p]->SetLineWidth(2);
        total[p]->Draw();
        pions[p] = (TH1D*) file->Get(Form("pions_p%d", p));
        pions[p]->SetLineColor(kBlue);
        pions[p]->SetLineWidth(2);
        pions[p]->Draw();
        kaons[p] = (TH1D*) file->Get(Form("kaons_p%d", p));
        kaons[p]->SetLineColor(kRed);
        kaons[p]->SetLineWidth(2);
        kaons[p]->Draw("same");
        protons[p] = (TH1D*) file->Get(Form("protons_p%d", p));
        protons[p]->SetLineColor(kGreen);
        protons[p]->SetLineWidth(2);
        protons[p]->Draw("same");
        latex_total->DrawLatex(0.68,0.66,"#scale[0.8]{ATLAS #bf{Internal}}");
        latex_total->DrawLatex(0.68,0.62,"#scale[0.6]{#bf{DPMJET gamma+Pb}}");
        ostringstream p_label;
        p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
        latex_total->DrawLatex(0.68,0.58, p_label.str().c_str());
        TLegend *legend = new TLegend(0.68,0.7,0.9,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(total[p], "All particles", "l");
        legend->AddEntry(pions[p], "Pions", "l");
        legend->AddEntry(kaons[p], "Kaons", "l");
        legend->AddEntry(protons[p], "Protons", "l");
        legend->Draw("Same");
        tc_total->SaveAs(Form("dEdx_histograms/dpmjet_truth_p%d.pdf",p));
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
        
        TFile *file = new TFile("root_files/gammaPb_dpmjet.root","read");
        positive[p] = (TH1D*) file->Get(Form("positive_p%d",p));
        negative[p] = (TH1D*) file->Get(Form("negative_p%d",p));
        
        //draw the distributions overlayed
        positive[p]->SetLineWidth(2);
        positive[p]->SetLineColor(kRed);
        positive[p]->Draw();
        negative[p]->SetLineWidth(2);
        negative[p]->SetLineColor(kBlue);
        negative[p]->Draw("same");
        latex->DrawLatex(0.75,0.73,"#scale[0.8]{ATLAS #bf{Internal}}");
        latex->DrawLatex(0.75,0.68,"#scale[0.6]{#bf{DPMJET gamma+Pb}}");
        ostringstream p_label;
        p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
        latex->DrawLatex(0.75,0.63, p_label.str().c_str());
        latex->DrawLatex(0.75,0.58,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
        TLegend *legend = new TLegend(0.7,0.8,0.93,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(positive[p], "Positively Charged Particles","l");
        legend->AddEntry(negative[p], "Negatively Charged Particles", "l");
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
        ratio[p] = (TH1D*) negative[p]->Clone(Form("negative_%d",p));
        ratio[p]->Divide(positive[p]);
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
        histName << "dEdx_histograms/quality1/dpmjet_p" << p << "_pPb_quality1cut.pdf";
        tc->SaveAs(histName.str().c_str());
        pads[1]->Clear();
    }
}
