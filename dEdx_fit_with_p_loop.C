// ROOT header files - calls implicit root functions
#include "TFile.h" // reads and writes root files
#include "TTree.h" // reads and writes TTrees
#include "TH1D.h" //Makes 1D histograms (double type)
#include "TF1.h" // creates and draws functions

//TF1 Fit Function: Double-Sided Crystal Ball - morgan created this!
// The DSCB function creates a gaussian function with exponential decay on both sides

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

void dEdx_fit_with_p_loop() {
    //here is where we load the data
    //bins!! This defines all of the cutoffs for the bins - I'm pretty sure this is the upper bound but check on this XX
    const int nbins_number = 1;
    const int pbins_number = 1;
    double nbins[nbins_number+1] = {25,60};
    double pbins[pbins_number+1] = {0.3,0.4};
    TH1D* distributions[nbins_number*pbins_number];
    
    // create the histograms! The x axis should be pT later but for now I'll keep dEdx XX
    int histCount = 0;
    for (int n = 0; n < nbins_number; n++) { // iterate over every n bin
        for (int p = 0; p < pbins_number; p++){ // iterate over every momentum bin
            ostringstream histName;
            histName << "n" << n << "p" << p;
            distributions[histCount] = new TH1D(histName.str().c_str(),";ln(dE/dx [MeV g^{-1} cm^{-1}]);dN^{trk} / d(ln(dE/dx))",100,-1,3);
            distributions[histCount] -> SetLineWidth(2);
            distributions[histCount] -> Sumw2();
            histCount++;
        }
    }
    //  for (int a = 0; a < numTrees; a++) {
    std::cout << "Reading file" << std::endl; //displays which file is being loaded

    // If the file is ~broken~ (ie empty, won't load, etc), it will just skip the file instead of shutting down the whole thing
    //TFile *file = new TFile(filenames[a].c_str(),"READ");
    TFile *file = new TFile("root_files/slimmed_user.srdas.29054111.OUTPUT._0.root","READ");
    if (file->IsZombie()) {
      std::cout << "Error opening file - Skipping file!" << std::endl;
    }
    //retreives TTrees from the file
    TTree *tree = (TTree*) file->Get("tree");

    //Declaring TTree variables:
    int trk_n; //number of particles in the event
    int trk_400; //a different ways to count particcles in the event
    float trk_pt[10000]; //transverse momentum
    float trk_eta[10000]; //pseudorapidity - angle of where most of the explosion goes, in a log of tanh scale
    float dEdx[10000]; //ionization energy across each cluster (set of pixels)
    float gaps[4]; //chooses if it's a upc or peripheral event - always 4

    //associated TTree branches - this connects the declared variables with the titles in the tree file
    tree->SetBranchAddress("trk_n", &trk_n); //the "&" says it's a single value
    tree->SetBranchAddress("trk_400", &trk_400);
    tree->SetBranchAddress("trk_pt", trk_pt); //no "&" means it's an array of variables
    tree->SetBranchAddress("trk_eta", trk_eta);
    tree->SetBranchAddress("trk_dEdx", dEdx);
    //tree->SetBranchAddress("Gaps_PsAsPeAe", gaps); //"<name of tree variable>,<local name>"
    
    for (int e = 0; e < tree->GetEntries(); e++) {
        tree->GetEntry(e);
        //iterate over all particles
        for (int i = 0; i < trk_n; i++) {
            //double momentum = trk_pt[i];
            double momentum = abs(trk_pt[i]*cosh(trk_eta[i]));
            int histNumber = 0;
            for (int n = 0; n < nbins_number; n++) {
                for (int p = 0; p < pbins_number; p++) {
                    if (nbins[n] < trk_400 && trk_400 <= nbins[n+1] && pbins[p] <= momentum && momentum < pbins[p+1]){distributions[histNumber]->Fill(log(dEdx[i]));}    //test if the particle belongs in that bin
                        histNumber++;
                    }
                }
            }
        }
    file->Close();
    
    //Normalize histograms
    int histNumber = 0;
    for (int n = 0; n < nbins_number; n++) {
        for (int p = 0; p < pbins_number; p++) {
            distributions[histNumber]->Scale(1,"width");
            histNumber++;
      }
    }
    //Redirect output information to a new txt file
    gSystem->RedirectOutput("UPCFitStatus.txt","w");
    
    //------------------------------Fitting ------------------------------------
    
    //create arrays to store TF1s
    // format: TF1 *<NAME> = new TF1("<NAME>",<name of fit function>,<start of fit on x axis>, <number of fit parameters>);
    TF1 *PionFits[nbins_number*pbins_number];
    TF1 *kaonFits[nbins_number*pbins_number];
    TF1 *ProtonFits[nbins_number*pbins_number];
    
    int PionFitCount = 0;
    int kaonFitCount = 0;
    int protonFitCount = 0;
    
    
    //Pion Fitting:
    //  initial guesses for parameters:
    float PionUpLim[8] = {0.75,0.65,0.55,0.5,0.4,0.35,0.25,0.15};         //Note: These values are fixed - they're how far along the x-axis the fit should go, so you have to set this!
    float PionN[8] = {2000000,2500000,3000000,3000000,4000000,4000000,4500000,6000000};
    float PionMean[8] = {0.0641,-0.01699,-0.07105,-0.1085,-0.1374,-0.1595,-0.1735,-0.1826};
    float PionSigma[8] = {0.1816,0.1736,0.1707,0.1701,0.1710,0.1714,0.1712,0.1710};
    float PionAlphaL[8] = {1.179,1.179,1.175,1.183,1.194,1.192,1.185,1.179};
    float PionAlphaH[8] = {1.29,1.24,1.23,1.22,1.22,1.20,1.16,1.07};
    
    TH1D* PionYields[8];
    for (int p = 0; p < pbins_number; p++) {
        ostringstream PionYieldName;
        PionYieldName << "PionYieldRange_p" << p;
        PionYields[p] = new TH1D(PionYieldName.str().c_str(),"Charged Pion Counts;N_{ch};#pi^{+/-} Counts",1,nbins);
    }

    for (int n = 0; n < nbins_number; n++) {
        for (int p = 0; p < pbins_number; p++) {
            ostringstream fitName;
            fitName << "Pion_n" << n << "_p" << p;
            PionFits[PionFitCount] = new TF1(fitName.str().c_str(),dscb,-0.7,PionUpLim[p],7);
            PionFits[PionFitCount]->SetParameters(PionN[p],PionMean[p],PionSigma[p],PionAlphaL[p],1000000,PionAlphaH[p],1000000);
            // The fixed parameters are from morgan's knowledge of what they should be - XX add your own later, leave out for now
            PionFits[PionFitCount]->FixParameter(4,1000000);
            PionFits[PionFitCount]->FixParameter(6,1000000);
            PionFits[PionFitCount]->FixParameter(1,PionMean[p]);
            PionFits[PionFitCount]->FixParameter(2,PionSigma[p]);
            PionFits[PionFitCount]->FixParameter(3,PionAlphaL[p]);
            PionFits[PionFitCount]->FixParameter(5,PionAlphaH[p]);
    
            //fit the dE/dx distribution
            distributions[PionFitCount]->Fit(PionFits[PionFitCount],"RNL");
            
            //values characterizing how good the fit is:
            std::cout << "Pion " << nbins[n] << " < N_{ch} <= " << nbins[n+1] << " ," << pbins[p] << " <= p < " << pbins[p+1] << " Chi ^2: " << PionFits[PionFitCount]->GetChisquare() << " NDF: " << PionFits[PionFitCount]->GetNDF() << " Chi^2/NDF: " << PionFits[PionFitCount]->GetChisquare() / PionFits[PionFitCount]->GetNDF() << std::endl;
            
            //calculate the integral of the fit function and save it in the PionYields histogram:
            int integral = abs(PionFits[PionFitCount]->Integral(-0.5,0)) + PionFits[PionFitCount]->Integral(0,3);
            PionYields[p]->SetBinContent(n+1,integral);
            PionYields[p]->SetBinError(n+1,sqrt(integral));
            std::cout << "Integral: " << integral << std::endl;
            
            PionFitCount++;
        }
    }

    //Kaon fitting:
    TH1D *kdata[nbins_number*pbins_number];
    TH1D *KaonYields[8];
    for (int p = 0; p < pbins_number; p++) {
        ostringstream kaonyields;
        kaonyields << "KaonYieldRange" << p;
        KaonYields[p] = new TH1D(kaonyields.str().c_str(),"Charged Kaon Counts;N_{ch};K^{+/-} Count",1,nbins);
    }
    
    float KaonLowLim[8] = {1,0.85,0.8,0.65,0.5,0.45,0.335,0.3};
    float KaonUpLim[8] = {1.9,1.7,1.55,1.3,1.2,1.,0.9,0.75};
    
    // initial guesses for each parameter of the fit:
    float KaonN[8] = {3000,6000,8000,10000,15000,17000,20000,30000}; //number of kaons
    float KaonMeans[8] = {1.3924,1.2379,1.0875,0.9510,0.8141,0.6828,0.5643,0.4522}; // averages
    float KaonAlphaL[8] = {1.438,1.498,1.349,1.167,1.879,1.534,1.571,1.585}; // Parameter from crystal ball fit
    float KaonNL[8] = {2.87,3.18,3.99,6.67,3.1,4.15,4.76,5.13}; // parameter from crystal ball fit
    
    for (int n = 0; n < nbins_number; n++){
        for (int p = 0; p < pbins_number; p++){
            
            //This part isolates the kaon peak from the pion peak
            ostringstream kaondatastream;
            kaondatastream << "KaonData_n" << n << "_p" << p;
            kdata[kaonFitCount] = new TH1D(kaondatastream.str().c_str(),"",100,-1,3);
            
            // subtracts fit for pion graph from data histogram
            for (int i = 0; i <= 100; i++){
                double hist = distributions[kaonFitCount]->GetBinContent(i);
                double bin = distributions[kaonFitCount]->GetBinCenter(i);
                double fit = PionFits[kaonFitCount]->Eval(bin);
                double difference = hist - fit;
                kdata[kaonFitCount]->SetBinContent(i,difference);
            }
            
            //fit
            ostringstream fitName;
            fitName << "Kaon_n" << n << "_p" << p;
            kaonFits[kaonFitCount] = new TF1(fitName.str().c_str(),dscb,KaonLowLim[p],KaonUpLim[p],7);
            kaonFits[kaonFitCount]->SetParameters(KaonN[p],KaonMeans[p],0.2,KaonAlphaL[p],KaonNL[p],1.1,1000000);
            // fixed parameters - don't use for now:
            kaonFits[kaonFitCount]->FixParameter(3,KaonAlphaL[p]);
            //kaonFits[kaonFitCount]->SetParLimits(3,aL, aH)
            kaonFits[kaonFitCount]->FixParameter(4,KaonNL[p]);
            kaonFits[kaonFitCount]->FixParameter(6,1000000);
            kaonFits[kaonFitCount]->FixParameter(1,KaonMeans[p]);
            kaonFits[kaonFitCount]->FixParameter(2,0.158);
            kaonFits[kaonFitCount]->FixParameter(5,0.772);
            
            kdata[kaonFitCount]->Fit(kaonFits[kaonFitCount],"RNL");
            
            std::cout << "Kaon " << nbins[n] << " < N_{ch} <= " << nbins[n+1] << " ," << pbins[p] << " <= p < " << pbins[p+1] << " Chi ^2: " << kaonFits[kaonFitCount]->GetChisquare() << " NDF: " << kaonFits[kaonFitCount]->GetNDF() << " Chi^2/NDF: " << kaonFits[kaonFitCount]->GetChisquare() / kaonFits[kaonFitCount]->GetNDF() << std::endl;
            
            int integral = abs(kaonFits[kaonFitCount]-> Integral(-0.5,0) + kaonFits[kaonFitCount]->Integral(0,3));
            KaonYields[p]->SetBinContent(n+1,integral);
            KaonYields[p]->SetBinError (n+1, sqrt(integral));
            std::cout << "Integral: " << integral << std::endl;
            
            kaonFitCount++;
        }
    }
    
    //Proton Fitting:
    TH1D *ProtonData[nbins_number*pbins_number];
    TH1D *ProtonYields[8];
    
    for (int p = 0; p < pbins_number; p++) {
        ostringstream ProtonCountName;
        ProtonCountName << "PbPb_YieldsRange" << p;
        ProtonYields[p] =  new TH1D(ProtonCountName.str().c_str(),"Charged Proton Counts;N_{ch}; p + #bar{p} Counts",1,nbins);
    }
    
    // Lower limits on each bin
    float protonLowLim[6] = {1.55,1.5,1.325,1.2,1.05,0.85};
    float protonUpLim[6] = {2.2,2.15,2.15,2.1,1.85,1.8};
    
    // Initial guesses
    float protonMean[6] = {1.854,1.721,1.642,1.521,1.384,1.217};
    float protonN[6] = {500,1000,1700,8000,10000,20000};
    float protonAlphaL[6] = {0.611,0.997,1.076,1.32,1.54,1.65};
    float protonnL[6] = {3.03,2.47,2.66,2.53,2.28,2.20};
    
    //get the isolated proton peak by subtracting the kaon fit from the isolated kaon peak
    for (int n = 0; n < nbins_number; n++) {
      for (int p = 0; p < pbins_number; p++) {
        ostringstream protondatastream;
        protondatastream << "ProtonData_n" << n << "_p" << p;
        ProtonData[protonFitCount] = new TH1D(protondatastream.str().c_str(),"",100,-1,3);
        for (int i = 0; i <= 100; i++) {
          double hist = kdata[protonFitCount]->GetBinContent(i);
          double bin = kdata[protonFitCount]->GetBinCenter(i);
          double fit = kaonFits[protonFitCount]->Eval(bin);
          double difference = hist - fit;

          ProtonData[protonFitCount]->SetBinContent(i,difference);
        }

        //fitting
        //if (p >= 2) {

            std::cout << "Fit Function: n" << n << "_p" << p << std::endl;

            ostringstream fitName;
            fitName << "Proton_n" << n << "_p" << p;
            ProtonFits[protonFitCount] = new TF1(fitName.str().c_str(),dscb,protonLowLim[p],protonUpLim[p],7);
            ProtonFits[protonFitCount]->SetParameter(0,protonN[p]);
            ProtonFits[protonFitCount]->SetParLimits(0,10,1000000); //the height can't be negative
            ProtonFits[protonFitCount]->FixParameter(1,protonMean[p-2]);
            ProtonFits[protonFitCount]->FixParameter(2,0.169);
            ProtonFits[protonFitCount]->FixParameter(3,protonAlphaL[p-2]);
            ProtonFits[protonFitCount]->FixParameter(4,protonnL[p-2]);
            ProtonFits[protonFitCount]->FixParameter(5,1.222);
            ProtonFits[protonFitCount]->FixParameter(6,1000000);
        
            ProtonData[protonFitCount]->Fit(ProtonFits[protonFitCount],"RNL");

            std::cout << "Proton " << nbins[n] << " < Nch <= " << nbins[n+1] << ", " << pbins[p] << " <= p < " << pbins[p+1] << " Chi ^2: " << ProtonFits[protonFitCount]->GetChisquare() << " NDF: " << ProtonFits[protonFitCount]->GetNDF() << " Chi^2/NDF: " << ProtonFits[protonFitCount]->GetChisquare() / ProtonFits[protonFitCount]->GetNDF() << std::endl;

            int integral = abs(ProtonFits[protonFitCount]->Integral(-0.5,0)) + ProtonFits[protonFitCount]->Integral(0,3);
            ProtonYields[p]->SetBinContent(n+1,integral);
            ProtonYields[p]->SetBinError(n+1,sqrt(integral));
            std::cout << "Integral: " << integral << std::endl;
          //}
          protonFitCount++;
      }
    }
    
    
    //Particle Yields as a function of number of charged particles:
    TH1D *TotalPions = new TH1D("TotalPions","Pion Count;N_{ch};#pi^{+/-}",1,nbins);
    TH1D *TotalKaons = new TH1D("TotalKaons","Kaon Count;N_{ch};K^{+/-}",1,nbins);
    TH1D *TotalProtons = new TH1D("TotalProtons","Proton Count;N_{ch};(p + #bar{p})",1,nbins);

    for (int n = 1; n < nbins_number; n++) {
      double PionTotal = 0;
      double KaonTotal = 0;
      double ProtonTotal = 0;
      for (int p = 0; p < pbins_number; p++) {
        PionTotal = PionTotal + PionYields[p]->GetBinContent(n);
        KaonTotal = KaonTotal + KaonYields[p]->GetBinContent(n);
      }
      for (int p = 0; p < pbins_number; p++) {
        ProtonTotal = ProtonTotal + ProtonYields[p]->GetBinContent(n);
      }
      TotalPions->SetBinContent(n,PionTotal);
      TotalPions->SetBinError(n,sqrt(PionTotal));
      TotalKaons->SetBinContent(n,KaonTotal);
      TotalKaons->SetBinError(n,sqrt(KaonTotal));
      TotalProtons->SetBinContent(n,ProtonTotal);
      TotalProtons->SetBinError(n,sqrt(ProtonTotal));
    }

    //Particle ratios as a function of the number of charged particles:
    TH1D *TotalKPiRatio = (TH1D*) TotalKaons->Clone("TotalKPiRatio");
    TotalKPiRatio->Divide(TotalPions);
    TotalKPiRatio->SetTitle(";N_{ch};(K^{+/-}) / #pi^{+/-}");

    TH1D *TotalpPiRatio = (TH1D*) TotalProtons->Clone("TotalpPiRatio");
    TotalpPiRatio->Divide(TotalPions);
    TotalpPiRatio->SetTitle(";N_{ch};(p + #bar{p}) / #pi^{+/-}");


    //Create a file for all of the histograms and fit functions
    TFile *outfile = new TFile("root_files/protonfit.root","recreate");
    int counter = 0;
    for (int n = 0; n < nbins_number; n++) {
      for (int p = 0; p < pbins_number; p++) {
        distributions[counter]->Write();
        PionFits[counter]->Write();
        kdata[counter]->Write();
        kaonFits[counter]->Write();
        ProtonData[counter]->Write();
        ProtonFits[counter]->Write();

        counter++;
      }
    }

    for (int p = 0; p < pbins_number; p++) {
        PionYields[p]->Write();
        KaonYields[p]->Write();
        ProtonYields[p]->Write();
        outfile->Close();
    }
}
