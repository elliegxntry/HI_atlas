#include "TGraph.h"

void plotting_test() {
    double x[2], y[2];
    x[0] = 1;
    x[1] = 5;
    y[0] = 3;
    y[1] = 10;
    auto g = new TGraph(2,x,y);
    g->Draw("AC+");
    
    
}
