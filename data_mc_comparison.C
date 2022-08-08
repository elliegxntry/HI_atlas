#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"

double bethe_bloch_function(double* p, double*par) {
    double m_part_c = 1;
    double l = par[0];
    double n = par[1];
    double r = par[2];
    
    double beta_gamma = p[0]/m_part_c;
    double beta_squared = TMath::Power(beta_gamma,2)/(1+TMath::Power(beta_gamma,2));
    
    double dEdx = 0;
    dEdx = -(l/beta_squared)*(0.5*TMath::Log(n*TMath::Power(beta_gamma,2))-beta_squared)+r;
    return dEdx;
}

void data_mc_comparison() {
    
    const int pbins_num = 8;
    double pbins[pbins_num+1] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    string p_and_pbar[2] = {"positive", "negative"};
    
    TFile *data_file = new TFile("root_files/pPb_data.root","read");
    TFile *mc_file = new TFile("root_files/gammaPb_dpmjet.root","read");
    TFile *scaled = new TFile("root_files/scaling_factors.root", "recreate");
    
    TH1D* data_distribution[pbins_num];
    TH1D* mc_distribution[pbins_num];
    TH1D* pions[pbins_num];
    TH1D* kaons[pbins_num];
    TH1D* protons[pbins_num];
    TF1* data_pionFit[pbins_num];
    TF1* mc_pionFit[pbins_num];
    TH1D* mc_distribution_new[pbins_num];
    TH1D* ratio[pbins_num];
    
    double difference[pbins_num] = {};
    double y_scale_m0[pbins_num+1] = {6.98579, 7.55825, 8.25687, 9.13616, 11.8121, 12.5909, 13.6024, 14.2701};
    double y_scale_m1[pbins_num+1] = {6.96112, 7.26236, 8.23275, 9.24009, 10.4656, 12.349, 13.6092, 14.8219};

    for (int m = 0; m < 2; m++) {
        double data_low_fit[pbins_num] = {-0.48,-0.45,-0.5,-0.55,-0.55,-0.55,-0.5,-0.5};
        double data_high_fit[pbins_num] = {0.22,0.18,0.15,0.2,0.25,0.2,0.15,0.15};
        double mc_low_fit[pbins_num] = {-0.2,-0.2,-0.2,-0.2,-0.25,-0.25,-0.2,-0.2};
        double mc_high_fit[pbins_num] = {0.35,0.28,0.25,0.3,0.35,0.35,0.35,0.3};
        float data_shift[pbins_num] = {-0.102858, -0.11343, -0.1139, -0.100269, -0.0970618, -0.099067, -0.101021, -0.0971758};
        float mc_shift[pbins_num+1] = {0.108957, 0.0909291, 0.0816315, 0.0975523, 0.104969, 0.107637, 0.110375, 0.101998};
        
        for (int p = 0; p < pbins_num; p++) {
            
            data_distribution[p] = (TH1D*) data_file->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            pions[p] = (TH1D*) scaled->Get(Form("pions_m%d_p%d", m, p));
            kaons[p] = (TH1D*) scaled->Get(Form("kaons_m%d_p%d", m, p));
            protons[p] = (TH1D*) scaled->Get(Form("protons_m%d_p%d", m, p));
            mc_distribution[p] = (TH1D*) scaled->Get(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            
            data_pionFit[p] = new TF1(Form("%s_p%d", p_and_pbar[m].c_str(), p), "gaus", data_low_fit[p] - data_shift[p], data_high_fit[p] - data_shift[p]);
            mc_pionFit[p] = new TF1(Form("%s_p%d", p_and_pbar[m].c_str(), p), "gaus", mc_low_fit[p] - mc_shift[p], mc_high_fit[p] - mc_shift[p]);
            
            data_distribution[p]->Fit(data_pionFit[p],"RNLQ");
//            pions[p]->Fit(mc_pionFit[p], "RNLQ");
            mc_distribution[p]->Fit(mc_pionFit[p],"RNLQ");
            
            difference[p] = data_pionFit[p]->GetParameter(1) - mc_pionFit[p]->GetParameter(1);
            cout << mc_pionFit[p]->GetParameter(1) << ", ";
            
            TCanvas *tc = new TCanvas("tc","New Canvas",50,50,700,600);
            double padSplit = 0.32;
            double leftMargin = 0.07;
            double bottomMargin = 0.15;
            double rightMargin = 0.02;
            double topMargin = 0.08;
            double height = 900;
            double width = height*(1.0 - topMargin*(1.0 - padSplit) - bottomMargin*padSplit)/(1.0 - leftMargin - rightMargin);
            
            TPad* pads[2];
            pads[0] = new TPad("pad0","",0.0,padSplit,1.0,0.98);
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
            
            pads[0]->cd();

            gStyle->SetOptStat(0);
            gPad->SetTicks();
            gPad->SetLogy(1);
            TLatex *latex = new TLatex();
            latex->SetNDC(kTRUE);
            TLegend *legend = new TLegend(0.78,0.92,0.95,0.8);
            gStyle->SetLegendBorderSize(1);

            data_distribution[p]->SetMarkerStyle(3);
            data_distribution[p]->SetMarkerColor(kRed);
            data_distribution[p]->SetLineColor(kRed);
            legend->AddEntry(data_distribution[p], "ATLAS Data", "P");
            mc_distribution[p]->SetMaximum(2e6);
            mc_distribution[p]->SetTitle(";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))");
            mc_distribution[p]->SetMarkerStyle(3);
            mc_distribution[p]->SetMarkerColor(kBlue);
            mc_distribution[p]->SetLineColor(kBlue);
//            pions[p]->SetMaximum(2e6);
//            pions[p]->SetMarkerStyle(3);
//            pions[p]->SetMarkerColor(kBlue);
//            pions[p]->SetLineColor(kBlue);
            legend->AddEntry(mc_distribution[p], "DPMJET Data", "P");
            mc_distribution[p]->Draw();
//            pions[p]->Draw();
            data_distribution[p]->Draw("same");
            
            latex->DrawLatex(0.78,0.76,"#scale[0.8]{ATLAS #bf{Internal}}");
            latex->DrawLatex(0.78,0.72,"#scale[0.6]{#bf{5.02 TeV p+Pb}}");
            latex->DrawLatex(0.78,0.68,"#scale[0.6]{#bf{DPMJET}}");
            latex->DrawLatex(0.78,0.64,"#scale[0.6]{#bf{-0.3 < #eta < 0.3}}");
            ostringstream p_label;
            p_label << "#scale[0.6]{#bf{" << pbins[p] << " < p < " << pbins[p+1] << "}}";
            latex->DrawLatex(0.78,0.6, p_label.str().c_str());
            ostringstream particle;
            particle <<"#scale[0.6]{#bf{" << p_and_pbar[m] << "}}";
            latex->DrawLatex(0.78,0.56, particle.str().c_str());
            
            ostringstream equation;
            if (m == 0) {
                equation << "#scale[0.6]{#bf{DPMJET = (original - " << mc_shift[p] <<  ")#times" << y_scale_m0[p] << "}}";}
            else {
                equation << "#scale[0.6]{#bf{DPMJET = (original - " << mc_shift[p] <<  ") * " << y_scale_m1[p] << "}}";}
            latex->DrawLatex(0.65, 0.52, equation.str().c_str());
            legend->Draw("same");
            tc->cd();
            pads[1] = new TPad("pad1","",0.0,0.0,1.0,padSplit);
            pads[1]->SetTopMargin(0.001);
            pads[1]->SetRightMargin(rightMargin);
            pads[1]->SetLeftMargin(leftMargin);
            pads[1]->SetBottomMargin(bottomMargin);

            pads[1]->Draw("same");
            pads[1]->cd();
            tc->cd();
            
            pads[1]->cd();
            gPad->SetTicks(1);
            gPad->SetLogy(0);
            ratio[p] = (TH1D*) data_distribution[p]->Clone(Form("%s_p%d", p_and_pbar[m].c_str(), p));
            ratio[p]->Divide(mc_distribution[p]);
            ratio[p]->SetLineWidth(2);
            ratio[p]->SetLineColor(kBlack);
            ratio[p]->SetMaximum(3.2);
            ratio[p]->SetMinimum(-0.3);
            ratio[p]->Draw("same");
            TLine *line = new TLine(-2, 1, 5, 1);
            line->SetLineColor(kBlack);
            line->SetLineStyle(kDashed);
            line->Draw("same");
            ratio[p]->GetYaxis()->SetTitleSize(0.06);
            ratio[p]->GetYaxis()->SetTitle("data/mc");
            tc->SaveAs(Form("dEdx_histograms/ratio_comparison_m%d_p%d.pdf",m,p));
            pads[1]->Clear();
            
        }
    }
}
