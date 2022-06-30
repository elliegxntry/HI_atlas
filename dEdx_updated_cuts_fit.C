#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1F.h"

double bethe_bloch(double* p, double* par) {
    double na = TMath::Na(); // 1/mol
    double re = 2.817940325e28; // fm
    double me_csquared = 0.510998918e44; // MeV
    double z = 1;
    double big_Z = 14;
    double a = 28.0855;// g/mol
    double tMax = 1; // change this later
    double i = 1; // also change this later
    double delta = 0.14; //also this is probably wrong
    
    double k = 4*TMath::Pi()*na*TMath::Power(re,2)*me_csquared;
    double l = -k*TMath::Power(z,2)*big_Z/a;
    double n = 2*me_csquared*tMax*TMath::Power(i,-2);
    
    double dEdx = 0;
    
    double beta = par[0];
    double gamma = par[1];
    
    dEdx = l*TMath::Power(beta,-2)*(0.5*TMath::Log(n*TMath::Power(beta*gamma,2))-TMath::Power(beta,2)-0.5*delta);
    return dEdx;
}

void dEdx_updated_cuts_fit() {
    TFile *file = new TFile("root_files/updated_data.root","read");
    const int pbins_num = 9;
    double pbins[pbins_num+1] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
    TH1D* distributions[2*pbins_num];
    for (int m = 0; m < 2; m++){
        //distributions[p] = (TH1D*) file->Get(Form("%_p%d", p_and_pbar[m].c_str(), p));
    }
}
