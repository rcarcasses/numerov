
#include <Rcpp.h>
#include "numerov.hpp"

using namespace Rcpp;

Numerov n;

// [[Rcpp::export]]
List rcpp_numerov() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

// [[Rcpp::export]]
void setPotential(NumericVector px = NumericVector(), NumericVector py = NumericVector()) {
    if(px.size() != py.size()) {
        cout << "Please pass two columns with the same size for the potential" << endl;
        return;
    }
    // translate the potential and feed it to numerov
	vector<Point> potential(px.size());
    for(int i = 0; i < px.size(); i++) {
        potential[i].x = px(i);
        potential[i].y = py(i);
    }
    
    n.setPotential(potential);
}

// [[Rcpp::export]]
List getPotential() {
    vector<Point> p = n.getPotential();
    int length = p.size();
    NumericVector x;
    NumericVector y;
    for(int i = 0; i < length; i++) {
        x.push_back(p[i].x);
        y.push_back(p[i].y);
    }
    
    return List::create(Named("x") = x, Named("y") = y);
}

// [[Rcpp::export]]
void computeSpectrum(int nEigen, double dE = 0.1, double tol = 1e-9) {
    n.dEmin = dE;
    n.tol = tol;
    n.findSpectrum(nEigen);
}

// [[Rcpp::export]]
NumericVector getEnergies() {
    vector<double> energies = n.getSpectrum().getEnergies();
    NumericVector x(energies.begin(), energies.end());

    return x;
}

// [[Rcpp::export]]
List getWavefunctions() {
    List wfs;
    vector<vector<Point>> WFs = n.getSpectrum().getWavefunctions();
    
    for(int i = 0; i < WFs.size(); i++) {
        int length = WFs[i].size();
        NumericVector x;
        NumericVector y;
        for(int j = 0; j < length; j++) {
            x.push_back(WFs[i][j].x);
            y.push_back(WFs[i][j].y);
        }
        List wf = List::create(Named("x") = x, Named("y") = y);
        wfs.push_back(wf);
    }
    
    return wfs;
}