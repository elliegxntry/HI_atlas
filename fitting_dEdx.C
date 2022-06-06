#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

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
    double pbins[pbins_num+1] = {0.3,0.4};
    TH1D* distributions[nbins_num*pbins_num];
    
    for (int n = 0; n < nbins_num; n++) {
        ostringstream histName;
        histName << "n" << n << "p0";
        distributions[n] = new TH1D(histName.str().c_str(),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-1,3);
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
    float dEdx[10000];
    float gaps[4];
    
    tree->SetBranchAddress("trk_n", &trk_n);
    tree->SetBranchAddress("trk_400", &trk_400);
    tree->SetBranchAddress("trk_pt", trk_pt);
    tree->SetBranchAddress("trk_eta", trk_eta);
    tree->SetBranchAddress("trk_dEdx", dEdx);
        
    for (int e = 0; e < tree->GetEntries(); e++) {
        tree->GetEntry(e);
        for (int i=0; i < trk_n; i++){;
            double momentum = abs(trk_pt[i]*cosh(trk_eta[i]));
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
    gSystem->RedirectOutput("UPCFitStatus.txt","w");
        
        // ///////////////////////////////////Fitting ///////////////////////////////////
        // Fit parameters: (N, mean, Sigma, AlphaL, nL, AlphaH, nH) - defined in dscb function
    // set the parameters of the fit - aka first guesses!
    const float pionN = 2000000;
    const float pionMean = -0.15;
    const float pionSigma = 0.2;
    const float pionAlphaL = 1.2;
    const float pionnL = 1000000;
    const float pionAlphaH = 1.2;
    const float pionnH = 1000000;
    
    //Set the parameters of the fit - initial guesses
    const float kaonN = 40000;
    const float kaonMean = 0.75;
    const float kaonSigma = 0.19;
    const float kaonAlphaL = 1.2;
    const float kaonnL = 3.;
    const float kaonAlphaH = 1.1;
    const float kaonnH = 1000000;
    
    //Set the parameters of the fit - initial guesses:
    const float protonN = 12500;
    const float protonMean = 1.58;
    const float protonSigma = .19;
    const float protonAlphaL = 1.2;
    const float protonnL = 2.5;
    const float protonAlphaH = 1.2;
    const float protonnH = 1000000;
    
    //initialize the fit functions
    TF1* pionFit[nbins_num*pbins_num];
    TF1* kaonFit[nbins_num*pbins_num];
    TF1* protonFit[nbins_num*pbins_num];
    
    //initialize the histograms
    TH1D *kData[nbins_num*pbins_num];
    TH1D *protonData[nbins_num*pbins_num];


    for (int n = 0; n < nbins_num; n++) {
        // pions - note: the "1000000"s are the power law, so it's just arbitrarily big
        ostringstream pionFitName;
        pionFitName << "pion_n" << n << "p0";
        pionFit[n] = new TF1(pionFitName.str().c_str(), dscb, -0.7, 0.5,7.);
            
        pionFit[n]->SetParameters(pionN, pionMean, pionSigma, pionAlphaL, pionnL, pionAlphaH, pionnH);
        
        distributions[n]->Fit(pionFit[n],"RNL");
        
        std::cout  << "Chi^2: " << pionFit[n]->GetChisquare() << "\n" << "NDF: " << pionFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << pionFit[n]->GetChisquare() / pionFit[n]->GetNDF() << "\n" << std::endl;

        // kaons
        ostringstream kaonFitName;
        kaonFitName << "kaon_n" << n << "p0";
        kaonFit[n] = new TF1(kaonFitName.str().c_str(), dscb, 0.5,1.2,7.);

        kaonFit[n]->SetParameters(kaonN, kaonMean, kaonSigma, kaonAlphaL, kaonnL, kaonAlphaH, kaonnH);
        
        //subtract pion fit from the histogram:
        ostringstream kDataName;
        kDataName << "kaonData_n" << n << "p0";
        kData[n] = new TH1D(kDataName.str().c_str(),"",100,-1,3);
        for (int i = 0; i <= 100; i++) {
            double hist = distributions[n]->GetBinContent(i);
            double bin = distributions[n]->GetBinCenter(i);
            double fit = pionFit[n]->Eval(bin);
            double difference = hist - fit;
            kData[n]->SetBinContent(i,difference);
        }
        
        kData[n]->Fit(kaonFit[n],"RNL");
        
        std::cout  << "Chi^2: " << kaonFit[n]->GetChisquare() << "\n" << "NDF: " << kaonFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << kaonFit[n]->GetChisquare() / kaonFit[n]->GetNDF()  << "\n" << std::endl;
        
        //protons:
        ostringstream protonFitName;
        protonFitName << "proton_n" << n << "p0";
        protonFit[n] = new TF1(protonFitName.str().c_str(), dscb, 1.2, 2.2, 7);

        protonFit[n]->SetParameters(protonN, protonMean, protonSigma, protonAlphaL, protonnL, protonAlphaH, protonnH);
        
        //subtract kaon and pion fit from histogram:
        ostringstream pDataName;
        pDataName << "protonData_n" << n << "p0";
        protonData[n] = new TH1D(pDataName.str().c_str(),"",100,-1, 3);
        for (int i = 0; i <= 100; i++) {
            double hist = kData[n] ->GetBinContent(i);
            double bin = kData[n]->GetBinCenter(i);
            double fit = kaonFit[n]->Eval(bin);
            double difference = hist - fit;
            
            protonData[n]->SetBinContent(i,difference);
        }
        
        protonData[n]->Fit(protonFit[n],"RNL");
        
        std::cout  << "Chi^2: " << protonFit[n]->GetChisquare() << "\n" << "NDF: " << protonFit[n]->GetNDF() << "\n" << "Chi^2/NDF: " << protonFit[n]->GetChisquare() / protonFit[n]->GetNDF()  << "\n" << std::endl;
        
        TFile *outfile = new TFile("root_files/fit.root","recreate");
        distributions[n]->Write();
        pionFit[n]->Write();
        kData[n]->Write();
        kaonFit[n]->Write();
        protonData[n]->Write();
        protonFit[n]->Write();
        outfile->Close();
    }
}
