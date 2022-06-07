
//ROOT header files

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
    
    //create Latex object
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);

    //get histograms from the root file you output in analysis macro
    TFile *file = new TFile("root_files/fit.root","read");
    TH1D* myhist[nbins_number*pbins_number];
    TH1D* kaonData[nbins_number*pbins_number];
    TH1D* protonData[nbins_number*pbins_number];
    TF1* pionFit[nbins_number*pbins_number];
    TF1* pionFit_fullRange[nbins_number*pbins_number];
    TF1* kaonFit[nbins_number*pbins_number];
    TF1* kaonFit_fullRange[nbins_number*pbins_number];
    TF1* protonFit[nbins_number*pbins_number];
    TF1* protonFit_fullRange[nbins_number*pbins_number];
    for (int n = 0; n < nbins_number*pbins_number; n++){
        //load total histogram data
        ostringstream histName;
        histName << "n" << n << "p0";// << p;
        myhist[n] = (TH1D*) file->Get(histName.str().c_str());

        //load pion subtracted data
        ostringstream kaondatastream;
        kaondatastream << "kaonData_n" << n << "p0";// << p;
        kaonData[n] = (TH1D*) file->Get(kaondatastream.str().c_str());
        
        //load pion and kaon subtracted data
        ostringstream protondatastream;
        protondatastream << "protonData_n" << n << "p0";// << p;
        protonData[n] = (TH1D*) file->Get(protondatastream.str().c_str());
        
        
        //load the pion fit
        ostringstream pionFitStream;
        pionFitStream << "pion_n" << n << "p0";// << p;
        pionFit[n] = (TF1*) file->Get(pionFitStream.str().c_str());
        pionFit_fullRange[n] = (TF1*) file->Get("pion_n0p0_fullrange");
        
        //load the kaon fit
        ostringstream kaonFitStream;
        kaonFitStream << "kaon_n" << n << "p0";// << p;
        kaonFit[n] = (TF1*) file->Get(kaonFitStream.str().c_str());
        kaonFit_fullRange[n] = (TF1*) file->Get("kaon_n0p0_fullrange");
        
        //load the proton fit
        ostringstream protonFitStream;
        protonFitStream << "proton_n" << n << "p0";// << p;
        protonFit[n] = (TF1*) file->Get(protonFitStream.str().c_str());
        protonFit_fullRange[n] = (TF1*) file->Get("proton_n0p0_fullrange");
    }
    
    for (int n=0; n < nbins_number; n++) {
        //plot the total histogram
        myhist[n]->SetLineWidth(2);
        myhist[n]->SetLineColor(kBlack);
        myhist[n]->Draw();
        

        // plot the histogram with the pion fit subtracted
        kaonData[n]->SetLineWidth(2);
        kaonData[n]->SetLineColor(kBlue);
        kaonData[n]->Draw("same");
        

        //plot the histogram with the pion and kaon fits subtracted
        protonData[n]->SetLineWidth(2);
        protonData[n]->SetLineColor(kMagenta+1);
        protonData[n]->Draw("same");
        
        //plot the fits:
         
        //plot the pion fit
        pionFit_fullRange[n]->SetLineWidth(2);
        pionFit_fullRange[n]->SetLineColor(kPink-1);
        pionFit_fullRange[n]->SetLineStyle(3);
        pionFit_fullRange[n]->Draw("same");
        
        pionFit[n]->SetLineWidth(2);
        pionFit[n]->SetLineColor(kRed);
        pionFit[n]->SetLineStyle(kDashed);
        pionFit[n]->Draw("same");
        
        

        //plot the kaon fit
        kaonFit_fullRange[n]->SetLineWidth(2);
        kaonFit_fullRange[n]->SetLineColor(kOrange-4);
        kaonFit_fullRange[n]->SetLineStyle(3);
        kaonFit_fullRange[n]->Draw("same");
        kaonFit[n]->SetLineWidth(2);
        kaonFit[n]->SetLineColor(kOrange-3);
        kaonFit[n]->SetLineStyle(kDashed);
        kaonFit[n]->Draw("same");
        

        //plot the proton fit
        protonFit_fullRange[n]->SetLineWidth(2);
        protonFit_fullRange[n]->SetLineStyle(3);
        protonFit_fullRange[n]->SetLineColor(kGreen-7);
        protonFit_fullRange[n]->Draw("same");
        protonFit[n]->SetLineWidth(2);
        protonFit[n]->SetLineColor(kGreen-3);
        protonFit[n]->SetLineStyle(kDashed);
        protonFit[n]->Draw("same");
        
        // build the legend
        TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(myhist[n], "Data","l");
        legend->AddEntry(kaonData[n], "Pion subtracted data", "l");
        legend->AddEntry(protonData[n], "Pion and Kaon subtracted data", "l");
        
        legend->AddEntry(pionFit[n], "Pion Fit","l");
        legend->AddEntry(pionFit_fullRange[n], "Extended pion Fit", "l");
        legend->AddEntry(kaonFit[n], "Kaon Fit","l");
        legend->AddEntry(kaonFit_fullRange[n],"Extended Kaon Fit","l");
        legend->AddEntry(protonFit[n],"Proton Fit","l");
        legend->AddEntry(protonFit_fullRange[n], "Extended Proton Fit","l");
        legend->Draw("same");
        
        latex->DrawLatex(0.7,0.65,"#scale[0.8]{ATLAS #bf{Internal}}");
        latex->DrawLatex(0.7,0.60,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");
        latex->DrawLatex(0.7,0.55,"#scale[0.6]{#bf{0.3 < p < 0.4}}");
        tc->SaveAs(Form("dEdx_histograms/n%d_p0.pdf",n));
    }
}
