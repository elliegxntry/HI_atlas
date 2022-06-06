
//ROOT header files

#include "TCanvas.h"  //allows you to create "canvases" to "draw" histograms on
#include "TStyle.h"   //allows you to customize style settings
#include "TLegend.h"  //allows you to make legends
#include "TLatex.h"   //allows you to make Latex-style text
#include "TFile.h"    //allows you to load root files
#include "TH1D.h"     //makes 1D double histograms
#include "TF1.h"    //makes functions

const int nbins_number =1;
const int pbins_number = 1;

void dEdx_protonfit_plot() {

    //create canvas and customize
    TCanvas *tc = new TCanvas(); //canvas is called tc
    gStyle->SetOptStat(0);         //if 1: Shows information about histogram in box at top of plot, if 0: information not shown
    gPad->SetTicks();              //Adds tick marks to the top and right edges of the plot
    gPad->SetLogy(1);              //Set the y-axis to a log-scale. Better readability for dE/dx data. To turn off, comment the line out or replace 1 with 0

    //create Latex object
    TLatex *latex = new TLatex();
    latex->SetNDC(kTRUE);          //allows you to place text according to relative position on canvas. To change to placing text according to plot coordinates, comment out line or replace TRUE with FALSE


    //get histograms from the root file you output in analysis macro
    TFile *file = new TFile("root_files/fit.root","read");
    TH1D* myhist[nbins_number*pbins_number];
    TH1D* kaonData[nbins_number*pbins_number];
    TH1D* protonData[nbins_number*pbins_number];
    TF1* pionFit[nbins_number*pbins_number];
    TF1* kaonFit[nbins_number*pbins_number];
    TF1* protonFit[nbins_number*pbins_number];
    for (int n = 0; n < nbins_number*pbins_number; n++){
        //load total histogram data
        ostringstream histName;
        histName << "n" << n << "p0";// << p;
        //myhist[histCount] = (TH1D*) file->Get(Form("n%ip%i",n,p));
        myhist[n] = (TH1D*) file->Get(Form("n%ip0",n));

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
        
        //load the kaon fit
        ostringstream kaonFitStream;
        kaonFitStream << "kaon_n" << n << "p0";// << p;
        kaonFit[n] = (TF1*) file->Get(kaonFitStream.str().c_str());
        
        //load the proton fit
        ostringstream protonFitStream;
        protonFitStream << "proton_n" << n << "p0";// << p;
        protonFit[n] = (TF1*) file->Get(protonFitStream.str().c_str());
    }

    int histCount = 0;
    
    for (int n=0; n < nbins_number; n++) {
        //for (int p=0; p < pbins_number; p++){
            
            //plot the data and the isolated peaks:
            
        //plot the total histogram
        myhist[histCount]->SetLineWidth(2);
        myhist[histCount]->SetLineColor(kBlack);
        myhist[histCount]->Draw();
        

        // plot the histogram with the pion fit subtracted
        kaonData[histCount]->SetLineWidth(2);
        kaonData[histCount]->SetLineColor(kBlue);
        kaonData[histCount]->Draw("same");
        

        //plot the histogram with the pion and kaon fits subtracted
        protonData[histCount]->SetLineWidth(2);
        protonData[histCount]->SetLineColor(kMagenta+1);
        protonData[histCount]->Draw("same");
        
        //plot the fits:
        

        //plot the pion fit
        pionFit[histCount]->SetLineWidth(2);
        pionFit[histCount]->SetLineColor(kRed);
        pionFit[histCount]->SetLineStyle(kDashed);
        pionFit[histCount]->Draw("same");
        

        //plot the kaon fit
        kaonFit[histCount]->SetLineWidth(2);
        kaonFit[histCount]->SetLineColor(kOrange-3);
        kaonFit[histCount]->SetLineStyle(kDashed);
        kaonFit[histCount]->Draw("same");
        

        //plot the proton fit
        protonFit[histCount]->SetLineWidth(2);
        protonFit[histCount]->SetLineColor(kGreen-3);
        protonFit[histCount]->SetLineStyle(kDashed);
        protonFit[histCount]->Draw("same");
        
        // build the legend
        TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
        gStyle->SetLegendBorderSize(1);
        legend->AddEntry(myhist[histCount], "Data","l");
        legend->AddEntry(kaonData[histCount], "Pion subtracted data", "l");
        legend->AddEntry(protonData[histCount], "Pion and Kaon subtracted data", "l");
        
        legend->AddEntry(pionFit[histCount], "Pion Fit","l");
        legend->AddEntry(kaonFit[histCount], "Kaon Fit","l");
        legend->AddEntry(protonFit[histCount],"Proton Fit","l");
        legend->Draw("same");
        
        latex->DrawLatex(0.7,0.65,"#scale[0.8]{ATLAS #bf{Internal}}");   //Always include "ATLAS INTERNAL" on your plots
        latex->DrawLatex(0.7,0.60,"#scale[0.6]{#bf{0nXn 5.02 TeV Pb+Pb}}");  //information about your data set
        //latex->DrawLatex(0.7,0.55,"#scale[0.6]{#bf{0.2 < p < 0.3 GeV}}");
                    
        //tc->SaveAs(Form("dEdx_histograms/n%d_p_%d.pdf",n,p));
        tc->SaveAs(Form("dEdx_histograms/n%d_p0_test.pdf",n));

        histCount++;
        //}
    }
}
