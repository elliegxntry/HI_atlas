#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1I.h"

double dscb(double* x, double*par) {
    //this names the parameters of the fit funtion
    double N = par[0];
    double mean = par[1];
    double sigma = par[2];

    double aL = par[3];
    double nL = par[4];
    double aH = par[5];
    double nH = par[6];
    
    //This is using the base parameters to create new parameters
    double u = (x[0] - mean) / sigma;
    double cL = TMath::Exp(-0.5*TMath::Power(aL,2));
    double cH = TMath::Exp(-0.5*TMath::Power(aH,2));

    double result = 0;
    
    //This defines when to use the Gaussian fit and when to use the exponential fit
    if (-aL <= u && u <= aH) {
        result = TMath::Exp(-0.5*TMath::Power(u,2));}
    if (u < -aL) {
        result = cL * TMath::Power(((aL/nL)*((nL/aL)-aL-u)),-nL);}
    if (u > aH) {
        result = cH * TMath::Power(((aH/nH)*((nH/aH)-aH+u)),-nH);}

    return N*result;
}

void calculations_fitting_dEdx() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int nbins_num = 1;
    const int pbins_num = 15;
    double nbins[nbins_num+1] = {25,60};
    double pbins[pbins_num+1] = {0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6};
    
    string p_and_pbar[2] = {"proton", "antiproton"};
    
    TFile *file = new TFile("root_files/data.root","read");
    TFile *outfile = new TFile("root_files/fit.root","recreate");
    
    for (int m = 0; m < 2; m++) {
        //initialize the histograms
        TH1D *kData[nbins_num*pbins_num];
        TH1D *protonData[nbins_num*pbins_num];
        TH1D* distributions[nbins_num*pbins_num];
        TH1D* proton_histogram[nbins_num*pbins_num];
            
            // ///////////////////////////////////Fitting ///////////////////////////////////
        const float pionN[1] = {770000};
        const float pionMean = -0.19;
        const float pionSigma = 0.18;
        const float pionAlphaL = 1.2;
        const float pionnL = 1000000;
        const float pionAlphaH = 1.1;
        const float pionnH = 100000;
        const float pionLowLim = -0.7;
        const float pionUpLim = 0.35;
         
        const float kaonN = 18000;
        const float kaonMean = 0.688;
        const float kaonSigma = 0.16;
        const float kaonAlphaL = 1.5;
        const float kaonnL = 10;
        const float kaonAlphaH = 1.1;
        const float kaonnH = 1000000;
        const float kaonLowLim = pionUpLim;
        const float kaonUpLim = 1.15;
         
        const float protonN = 10000;
        const float protonMean = 1.5;
        const float protonSigma = .16;
        const float protonAlphaL = 1.2;
        const float protonnL = 11;
        const float protonAlphaH = 1.2;
        const float protonnH = 1e6;
        const float protonLowLim = kaonUpLim;
        const float protonUpLim = 2.1;
        
        //initialize the fit functions
        TF1* pionFit[nbins_num*pbins_num];
        TF1* pionFit_fullrange[nbins_num*pbins_num];
        TF1* kaonFit[nbins_num*pbins_num];
        TF1* kaonFit_fullrange[nbins_num*pbins_num];
        TF1* protonFit[nbins_num*pbins_num];
        TF1* protonFit_fullrange[nbins_num*pbins_num];

        for (int p = 0; p < pbins_num; p++) {
            ostringstream output_file_name;
            output_file_name << p_and_pbar[m] << "_p" << p << "_UPCFitStatus.txt";
            outfile->cd();
            
            ostringstream histName;
            histName << p_and_pbar[m] << "_n25_60_p" << p;
            distributions[p] = (TH1D*) file->Get(histName.str().c_str());
            distributions[p]->Sumw2();
            proton_histogram[p] = new TH1D(Form("total%ss_p%d",p_and_pbar[m].c_str(),p),"Proton count;N_{ch};p{+/-}",1,-2,5);
            std::cout << "starting analysis on " << pbins[p] << std::endl;
            pionFit[p] = new TF1(Form("%s_pion_n25_60_p%d",p_and_pbar[m].c_str(),p), dscb, pionLowLim, pionUpLim,7.);
            pionFit_fullrange[p] = new TF1(Form("%s_pion_n25_60_p%d_fullrange",p_and_pbar[m].c_str(),p), dscb, -2,5,7.);
            pionFit[p]->SetParameters(pionN[0], pionMean, pionSigma, pionAlphaL, pionnL, pionAlphaH, pionnH);
            
            distributions[p]->Fit(pionFit[p],"RNLQ");
            distributions[p]->Fit(pionFit[p],"RNLQ");
            distributions[p]->Fit(pionFit[p],"RNLQ");
            distributions[p]->Fit(pionFit[p],"RNL");
            
            pionFit_fullrange[p]->SetParameters(pionFit[p]->GetParameter(0),pionFit[p]->GetParameter(1),pionFit[p]->GetParameter(2),pionFit[p]->GetParameter(3),pionFit[p]->GetParameter(4),pionFit[p]->GetParameter(5),pionFit[p]->GetParameter(6));

            std::cout  << "Chi^2: " << pionFit[p]->GetChisquare() << "\n" << "NDF: " << pionFit[p]->GetNDF() << "\n" << "Chi^2/NDF: " << pionFit[p]->GetChisquare() / pionFit[p]->GetNDF() << "\n" << std::endl;

            // kaons
            kaonFit[p] = new TF1(Form("%s_kaon_n25_60_p%d",p_and_pbar[m].c_str(),p), dscb, kaonLowLim,kaonUpLim,7.);
            kaonFit_fullrange[p] = new TF1(Form("%s_kaon_n25_60_p%d_fullRange",p_and_pbar[m].c_str(),p), dscb, -2,5,7.);
            kaonFit[p]->SetParameters(kaonN, kaonMean, kaonSigma, kaonAlphaL, kaonnL, kaonAlphaH, kaonnH);
            kaonFit[p]->SetParLimits(4,1.1,20.1);
            //kaonFit[p]->SetParLimits(6,0.1,1000.1);
            
            //subtract pion fit from the histogram:
            kData[p] = new TH1D(Form("%s_kaonData_n25_60_p%d",p_and_pbar[m].c_str(),p),"",100,-2,5);
            for (int i = 0; i <= 100; i++) {
                double hist = distributions[p]->GetBinContent(i);
                double bin = distributions[p]->GetBinCenter(i);
                double fit = pionFit[p]->Eval(bin);
                double difference = hist - fit;
                kData[p]->SetBinContent(i,difference);
            }
            //kData[p]->Fit(kaonFit[p],"RNLQ");
            //kData[p]->Fit(kaonFit[p],"RNLQ");
            //kData[p]->Fit(kaonFit[p],"RNLQ");
            kData[p]->Fit(kaonFit[p],"RNL");
            
            kaonFit_fullrange[p]->SetParameters(kaonFit[p]->GetParameter(0),kaonFit[p]->GetParameter(01),kaonFit[p]->GetParameter(02),kaonFit[p]->GetParameter(03),kaonFit[p]->GetParameter(04),kaonFit[p]->GetParameter(05),kaonFit[p]->GetParameter(06));
            int integral_kaons = abs(kaonFit_fullrange[p]->Integral(-0.5,0)) + kaonFit_fullrange[p]->Integral(0,4);
            std::cout << "Number of kaons: " << integral_kaons << std::endl;

            std::cout  << "Chi^2: " << kaonFit[p]->GetChisquare() << "\n" << "NDF: " << kaonFit[p]->GetNDF() << "\n" << "Chi^2/NDF: " << kaonFit[p]->GetChisquare() / kaonFit[p]->GetNDF()  << "\n" << std::endl;
            
            //protons:
            protonFit[p] = new TF1(Form("%s_proton_n25_60_p%d",p_and_pbar[m].c_str(),p), dscb, protonLowLim, protonUpLim, 7);
            protonFit_fullrange[p] = new TF1(Form("%s_proton_n25_60_p%d",p_and_pbar[m].c_str(),p), dscb, -2, 5, 7);

            protonFit[p]->SetParameters(protonN, protonMean, protonSigma, protonAlphaL, protonnL, protonAlphaH, protonnH);
            protonFit[p]->SetParLimits(4,1.1,20.1);
            protonFit[p]->SetParLimits(6,1e5,1e6);
            
            //subtract kaon and pion fit from histogram:
            protonData[p] = new TH1D(Form("%s_protonData_n25_60_p%d",p_and_pbar[m].c_str(),p),"",100,-2,5);
            for (int i = 0; i <= 100; i++) {
                double hist = kData[p] ->GetBinContent(i);
                double bin = kData[p]->GetBinCenter(i);
                double fit = kaonFit[p]->Eval(bin);
                double difference = hist - fit;
                
                protonData[p]->SetBinContent(i,difference);
            }
            //protonData[p]->Fit(protonFit[p],"RNLQ");
            //protonData[p]->Fit(protonFit[p],"RNLQ");
            //protonData[p]->Fit(protonFit[p],"RNLQ");
            protonData[p]->Fit(protonFit[p],"RNL");
            
            protonFit_fullrange[p]->SetParameters(protonFit[p]->GetParameter(0),protonFit[p]->GetParameter(01),protonFit[p]->GetParameter(02),protonFit[p]->GetParameter(03),protonFit[p]->GetParameter(04),protonFit[p]->GetParameter(05),protonFit[p]->GetParameter(06));
            int integral = abs(protonFit_fullrange[p]->Integral(-0.5,0)) + protonFit_fullrange[p]->Integral(0,4);
            std::cout << "Number of " << p_and_pbar[m] << "s: " << integral << std::endl;
            proton_histogram[p]->AddBinContent(p+1,integral);

            std::cout  << "Chi^2: " << protonFit[p]->GetChisquare() << "\n" << "NDF: " << protonFit[p]->GetNDF() << "\n" << "Chi^2/NDF: " << protonFit[p]->GetChisquare() / protonFit[p]->GetNDF()  << "\n" << std::endl;
            
            distributions[p]->Write();
            pionFit[p]->Write();
            pionFit_fullrange[p]->Write();
            kData[p]->Write();
            kaonFit[p]->Write();
            kaonFit_fullrange[p]->Write();
            protonData[p]->Write();
            protonFit[p]->Write();
            protonFit_fullrange[p]->Write();
            proton_histogram[p]->Write();
            distributions[p] = nullptr;
            proton_histogram[p] = nullptr;
        }//p loop
    }//m loop
    outfile->Close();
}

