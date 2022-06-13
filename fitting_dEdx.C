d#include "TFile.h"
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

void fitting_dEdx() {
    // ///////////////////////////////////Load in the data!///////////////////////////////////
    const int nbins_num = 1;
    const int pbins_num = 1;
    //double nbins[nbins_num+1] = {25,30,35,40,45,50,55,60};
    double nbins[nbins_num+1] = {25,60};
    double pbins[pbins_num+1] = {0.38,0.4};
    TH1D* distributions[nbins_num*pbins_num];
    
    for (int n = 0; n < nbins_num; n++) {
        ostringstream histName;
        histName << "n" << n << "p0";
        distributions[n] = new TH1D(histName.str().c_str(),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-1,4);
        distributions[n]->Sumw2();
    }
    std::cout << "Reading file!" << std::endl;
    
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
    if (file->IsZombie()){
        std::cout << "Error opening file - skipping file!" << std::endl;
    }
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
        
    for (int e = 0; e < tree->GetEntries(); e++) {
        if (e % 100000 == 0) std::cout << "Loading event " << e << std::endl;
        tree->GetEntry(e);
        for (int i=0; i < trk_n; i++){
            if (abs(trk_eta[i]) > 0.8) {continue;}
            double momentum = abs(trk_pt[i]*cosh(trk_eta[i]));
            if (trk_q[i] > 0) {continue;}
            for (int n = 0; n < nbins_num; n++) {
                if (nbins[n] < trk_400 && trk_400 <= nbins[n+1] && pbins[0] <= momentum && momentum < pbins[1]) {
                    distributions[n]->Fill(log(dEdx[i]));
                }
            }
        }
    }
    file->Close();
    for (int n = 0; n < nbins_num; n++) {
        distributions[n]->Scale(1,"width");
    }
    TFile *outfile = new TFile("root_files/antiproton_fit_p38_4.root","recreate");
    gSystem->RedirectOutput("antiproton_p38_4_UPCFitStatus.txt","w");
        
        // ///////////////////////////////////Fitting ///////////////////////////////////
        // Fit parameters: (N, mean, Sigma, AlphaL, nL, AlphaH, nH) - defined in dscb function
    // set the parameters of the fit - aka first guesses!
    const float pionN = 770000;
    const float pionMean = -0.19;
    const float pionSigma = 0.18;
    const float pionAlphaL = 1.2;
    const float pionnL = 1000000;
    const float pionAlphaH = 1.1;
    const float pionnH = 100000;
    const float pionLowLim = -0.7;
    const float pionUpLim = 0.35;
    
    //Set the parameters of the fit - initial guesses
    const float kaonN = 18000;
    const float kaonMean = 0.688;
    const float kaonSigma = 0.16;
    const float kaonAlphaL = 1.5;
    const float kaonnL = 100000;
    const float kaonAlphaH = 1.1;
    const float kaonnH = 1000000;
    const float kaonLowLim = pionUpLim;
    const float kaonUpLim = 1.15;
    
    //Set the parameters of the fit - initial guesses:
    const float protonN = 14000;
    const float protonMean = 1.5;
    const float protonSigma = .17;
    const float protonAlphaL = 1.2;
    const float protonnL = 1.4;
    const float protonAlphaH = 1.2;
    const float protonnH = 1000000;
    const float protonLowLim = kaonUpLim;
    const float protonUpLim = 2.1;
    
    //initialize the fit functions
    TF1* pionFit[nbins_num*pbins_num];
    TF1* pionFit_fullrange[nbins_num*pbins_num];
    TF1* kaonFit[nbins_num*pbins_num];
    TF1* kaonFit_fullrange[nbins_num*pbins_num];
    TF1* protonFit[nbins_num*pbins_num];
    TF1* protonFit_fullrange[nbins_num*pbins_num];
    
    //initialize the histograms
    TH1D *kData[nbins_num*pbins_num];
    TH1D *protonData[nbins_num*pbins_num];
    TH1D *totalantiprotons = new TH1D("totalantiprotons","Proton count;N_{ch};p{+/-}",1,-1,4);

    for (int n = 0; n < nbins_num; n++) {
        
        std::cout << "starting analysis" << std::endl;
        // pions - note: the "1000000"s are the power law, so it's just arbitrarily big
        ostringstream pionFitName;
        pionFitName << "pion_n" << n << "p0";
        pionFit[n] = new TF1(pionFitName.str().c_str(), dscb, pionLowLim, pionUpLim,7.);
        pionFit_fullrange[n] = new TF1("pion_n0p0_fullrange", dscb, -1,4,7.);
        pionFit[n]->SetParameters(pionN, pionMean, pionSigma, pionAlphaL, pionnL, pionAlphaH, pionnH);
        //pionFit[n]->SetParLimits(3,0,5);
        //pionFit[n]->SetParLimits(4,10000,10000000);
        //pionFit[n]->SetParLimits(6,10000,10000000);

        distributions[n]->Fit(pionFit[n],"RNLQ");
        distributions[n]->Fit(pionFit[n],"RNLQ");
        distributions[n]->Fit(pionFit[n],"RNLQ");
        distributions[n]->Fit(pionFit[n],"RNL");
        
        pionFit_fullrange[n]->SetParameters(pionFit[n]->GetParameter(0),pionFit[n]->GetParameter(1),pionFit[n]->GetParameter(2),pionFit[n]->GetParameter(3),pionFit[n]->GetParameter(4),pionFit[n]->GetParameter(5),pionFit[n]->GetParameter(6));

        std::cout  << "Chi^2: " << pionFit[n]->GetChisquare() << "\n" << "NDF: " << pionFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << pionFit[n]->GetChisquare() / pionFit[n]->GetNDF() << "\n" << std::endl;

        // kaons
        ostringstream kaonFitName;
        kaonFitName << "kaon_n" << n << "p0";
        kaonFit[n] = new TF1(kaonFitName.str().c_str(), dscb, kaonLowLim,kaonUpLim,7.);
        kaonFit_fullrange[n] = new TF1("kaon_n0p0_fullrange", dscb, -1,4,7.);
        kaonFit[n]->SetParameters(kaonN, kaonMean, kaonSigma, kaonAlphaL, kaonnL, kaonAlphaH, kaonnH);
        //kaonFit[n]->SetParLimits(4,0,100000000);
        //kaonFit[n]->SetParLimits(5,0,5);
        //kaonFit[n]->SetParLimits(6,100000,10000000);
        
        //subtract pion fit from the histogram:
        ostringstream kDataName;
        kDataName << "kaonData_n" << n << "p0";
        kData[n] = new TH1D(kDataName.str().c_str(),"",100,-1,4);
        for (int i = 0; i <= 100; i++) {
            double hist = distributions[n]->GetBinContent(i);
            double bin = distributions[n]->GetBinCenter(i);
            double fit = pionFit[n]->Eval(bin);
            double difference = hist - fit;
            kData[n]->SetBinContent(i,difference);
        }
        //kData[n]->Fit(kaonFit[n],"RNLQ");
        kData[n]->Fit(kaonFit[n],"RNLQ");
        kData[n]->Fit(kaonFit[n],"RNLQ");
        kData[n]->Fit(kaonFit[n],"RNL");
        
        kaonFit_fullrange[n]->SetParameters(kaonFit[n]->GetParameter(0),kaonFit[n]->GetParameter(01),kaonFit[n]->GetParameter(02),kaonFit[n]->GetParameter(03),kaonFit[n]->GetParameter(04),kaonFit[n]->GetParameter(05),kaonFit[n]->GetParameter(06));

        std::cout  << "Chi^2: " << kaonFit[n]->GetChisquare() << "\n" << "NDF: " << kaonFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << kaonFit[n]->GetChisquare() / kaonFit[n]->GetNDF()  << "\n" << std::endl;
        
        //protons:
        ostringstream protonFitName;
        protonFitName << "proton_n" << n << "p0";
        protonFit_fullrange[n] = new TF1("proton_n0p0_fullrange", dscb, -1,4,7.);
        protonFit[n] = new TF1(protonFitName.str().c_str(), dscb, protonLowLim, protonUpLim, 7);

        protonFit[n]->SetParameters(protonN, protonMean, protonSigma, protonAlphaL, protonnL, protonAlphaH, protonnH);
        //protonFit[n]->SetParLimits(3,0,100);
        //protonFit[n]->SetParLimits(4,0,100);
        //protonFit[n]->SetParLimits(6,100000,10000000);
        
        //subtract kaon and pion fit from histogram:
        ostringstream pDataName;
        pDataName << "protonData_n" << n << "p0";
        protonData[n] = new TH1D(pDataName.str().c_str(),"",100,-1, 4);
        for (int i = 0; i <= 100; i++) {
            double hist = kData[n] ->GetBinContent(i);
            double bin = kData[n]->GetBinCenter(i);
            double fit = kaonFit[n]->Eval(bin);
            double difference = hist - fit;
            
            protonData[n]->SetBinContent(i,difference);
        }
        //protonData[n]->Fit(protonFit[n],"RNLQ");
        protonData[n]->Fit(protonFit[n],"RNLQ");
        protonData[n]->Fit(protonFit[n],"RNLQ");
        protonData[n]->Fit(protonFit[n],"RNL");
        
        protonFit_fullrange[n]->SetParameters(protonFit[n]->GetParameter(0),protonFit[n]->GetParameter(01),protonFit[n]->GetParameter(02),protonFit[n]->GetParameter(03),protonFit[n]->GetParameter(04),protonFit[n]->GetParameter(05),protonFit[n]->GetParameter(06));
        int integral = abs(protonFit_fullrange[n]->Integral(-0.5,0)) + protonFit_fullrange[n]->Integral(0,4);
        std::cout << "Number of protons: " << integral << std::endl;
        totalantiprotons->AddBinContent(n+1,integral);

        std::cout  << "Chi^2: " << protonFit[n]->GetChisquare() << "\n" << "NDF: " << protonFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << protonFit[n]->GetChisquare() / protonFit[n]->GetNDF()  << "\n" << std::endl;
        
        distributions[n]->Write();
        pionFit[n]->Write();
        pionFit_fullrange[n]->Write();
        kData[n]->Write();
        kaonFit[n]->Write();
        kaonFit_fullrange[n]->Write();
        protonData[n]->Write();
        protonFit[n]->Write();
        protonFit_fullrange[n]->Write();
        totalantiprotons->Write();
        outfile->Close();
    }
}
