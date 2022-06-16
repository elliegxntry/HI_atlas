#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text
#include "TFile.h"    //allows you to load root files
#include "TH1D.h"     //makes 1D double histograms
#include "TF1.h"    //makes functions

const int nbins_number = 1;
const int pbins_number = 1;

void dEdx_protonfit_plot() {

    //create canvas and customize
    TCanvas *tc = new TCanvas();
    gStyle->SetOptStat(0);
    gPad->SetTicks();
    gPad->SetLogy(1);
    
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    TFile *data_file = new TFile("root_files/antiproton_data.root","read");
    TFile *file = new TFile("root_files/antiproton_fit.root","read");
    TH1D* myhist[nbins_number*pbins_number];
    TH1D* kaonData[nbins_number*pbins_number];
    TH1D* protonData[nbins_number*pbins_number];
    TF1* pionFit[nbins_number*pbins_number];
    TF1* pionFit_fullRange[nbins_number*pbins_number];
    TF1* kaonFit[nbins_number*pbins_number];
    TF1* kaonFit_fullRange[nbins_number*pbins_number];
    TF1* protonFit[nbins_number*pbins_number];
    TF1* protonFit_fullRange[nbins_number*pbins_number];
    for (int p = 0; p < pbins_number; p++) {
        //load total histogram data
        ostringstream histName;
        histName << "n25_60_p" << p;
//        myhist[p] = (TH1D*) file->Get(histName.str().c_str());
        myhist[p] = (TH1D*) data_file->Get("n25_60_p0");

        //load pion subtracted data
        ostringstream kaondatastream;
        kaondatastream << "kaonData_n25_60_p" << p;
        kaonData[p] = (TH1D*) file->Get(kaondatastream.str().c_str());
        
        //load pion and kaon subtracted data
        ostringstream protondatastream;
        protondatastream << "protonData_n25_60_p" << p;
        protonData[p] = (TH1D*) file->Get(protondatastream.str().c_str());
        
        
        //load the pion fit
        ostringstream pionFitStream;
        pionFitStream << "pion_n25_60_p" << p;
        pionFit[p] = (TF1*) file->Get(pionFitStream.str().c_str());
        ostringstream pionFit_fullrangeName;
        pionFit_fullrangeName << "pion_n25_60_p" << p << "_fullrange";
        pionFit_fullRange[p] = (TF1*) file->Get(pionFit_fullrangeName.str().c_str());
        
        //load the kaon fit
        ostringstream kaonFitStream;
        kaonFitStream << "kaon_n25_60_p" << p;
        kaonFit[p] = (TF1*) file->Get(kaonFitStream.str().c_str());
        ostringstream kaonFit_fullrangeName;
        kaonFit_fullrangeName << "kaon_n25_60_p" << p << "_fullrange";
        kaonFit_fullRange[p] = (TF1*) file->Get(kaonFit_fullrangeName.str().c_str());
        
        //load the proton fit
        ostringstream protonFitStream;
        protonFitStream << "proton_n25_60_p" << p;
        protonFit[p] = (TF1*) file->Get(protonFitStream.str().c_str());
        ostringstream protonFit_fullrangeName;
        protonFit_fullrangeName << "proton_n25_60_p" << p << "_fullrange";
        protonFit_fullRange[p] = (TF1*) file->Get(protonFit_fullrangeName.str().c_str());
    }
    
    for (int p = 0; p < pbins_number; p++) {
        //plot the total histogram
        myhist[p]->SetLineWidth(2);
        myhist[p]->SetLineColor(kBlack);
        myhist[p]->Draw();
        

        // plot the histogram with the pion fit subtracted
        kaonData[p]->SetLineWidth(2);
        kaonData[p]->SetLineColor(kBlue);
        kaonData[p]->Draw("same");
        

        //plot the histogram with the pion and kaon fits subtracted
        protonData[p]->SetLineWidth(2);
        protonData[p]->SetLineColor(kMagenta+1);
        protonData[p]->Draw("same");
        
        //plot the fits:
         
        //plot the pion fit
        pionFit_fullRange[p]->SetLineWidth(2);
        pionFit_fullRange[p]->SetLineColor(kPink-1);
        pionFit_fullRange[p]->SetLineStyle(3);
        pionFit_fullRange[p]->Draw("same");
        
        pionFit[p]->SetLineWidth(2);
        pionFit[p]->SetLineColor(kRed);
        pionFit[p]->SetLineStyle(kDashed);
        pionFit[p]->Draw("same");
        
        

        //plot the kaon fit
        kaonFit_fullRange[p]->SetLineWidth(2);
        kaonFit_fullRange[p]->SetLineColor(kOrange-4);
        kaonFit_fullRange[p]->SetLineStyle(3);
        kaonFit_fullRange[p]->Draw("same");
        kaonFit[p]->SetLineWidth(2);
        kaonFit[p]->SetLineColor(kOrange-3);
        kaonFit[p]->SetLineStyle(kDashed);
        kaonFit[p]->Draw("same");
        

        //plot the proton fit
        protonFit_fullRange[p]->SetLineWidth(2);
        protonFit_fullRange[p]->SetLineStyle(3);
        protonFit_fullRange[p]->SetLineColor(kGreen-7);
        protonFit_fullRange[p]->Draw("same");
        protonFit[p]->SetLineWidth(2);
        protonFit[p]->SetLineColor(kGreen-3);
        protonFit[p]->SetLineStyle(kDashed);
        protonFit[p]->Draw("same");
        
        // build the legend
        TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(myhist[p], "Data","l");
        legend->AddEntry(kaonData[p], "Pion subtracted data", "l");
        legend->AddEntry(protonData[p], "Pion and Kaon subtracted data", "l");
        
        legend->AddEntry(pionFit[p], "Pion Fit","l");
        legend->AddEntry(pionFit_fullRange[p], "Extended pion Fit", "l");
        legend->AddEntry(kaonFit[p], "Kaon Fit","l");
        legend->AddEntry(kaonFit_fullRange[p],"Extended Kaon Fit","l");
        legend->AddEntry(protonFit[p],"Proton Fit","l");
        legend->AddEntry(protonFit_fullRange[p], "Extended Proton Fit","l");
        legend->Draw("same");
        
        latex->DrawLatex(0.68,0.65,"#scale[0.8]{ATLAS #bf{Internal}}");
        latex->DrawLatex(0.68,0.60,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");
        latex->DrawLatex(0.68,0.55,"#scale[0.6]{#bf{0.3 < p < 0.4}}");
        latex->DrawLatex(0.68,0.5,"#scale[0.6]{#bf{-0.8 < #eta < 0.8}}");
        latex->DrawLatex(0.68,0.45,"#scale[0.6]{#bf{25 < N_{ch} < 60}}");
        tc->SaveAs(Form("dEdx_histograms/antiproton_n25_60_p%d.pdf",p));
    }
}
