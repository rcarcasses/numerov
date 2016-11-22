#ifndef _NUMEROV_
#define _NUMEROV_
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <functional>

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Point {
    double x = 0;
    double y = 0;
    Point() {
        x = 0;
        y = 0;
    }
    Point(double xx, double yy) {
        x = xx;
        y = yy;
    }
};

struct Range {
	double eMin = 0;
	double eMax = 1;
	Range(double m0, double m1) {
	    eMin = m0;
	    eMax = m1;
	}
};

struct Mode {
	double energy;
    int index;
	vector<Point> wavefunction;
	Mode() {   // default constructor for non good modes
	    index = -1;   
	}
	Mode(double e, vector<Point> f, int n){
		energy = e;
		wavefunction = f;
		index = n;
	}
};

struct Spectrum {
	vector<Mode>  modes;
	vector<Point> potential;

	void addMode(Mode m) {
	   	modes.push_back(m); 
	}
	
	void clear() {
	   	modes.clear();
		potential.clear();
   	}
       
    vector<double> getEnergies() {
        vector<double> energies;
        for(int i = 0; i < modes.size(); i++) 
            energies.push_back(modes[i].energy);
        
        return energies;
    }
    
    vector<vector<Point>> getWavefunctions() {
        vector<vector<Point>> wfs;
        for(int i = 0; i < modes.size(); i++) 
            wfs.push_back(modes[i].wavefunction);
        
        return wfs;
    }
    
};

class Numerov {
	private:
		vector<Point> potential;
		vector<Range> zeros;
		vector<Point> solLR;
		vector<Point> solRL;
		vector<Point> sol;
  	    double *A = NULL, *B = NULL, *C = NULL, *D = NULL;   // spline variables
		Point minPot;
		Spectrum spectrum;
		std::function<double(double)> potFunc;
		function<double(double)> diffFunc = [this](double E){
			return diff(E);
		};

		double bisection(function<double(double)> diffFunc, double min, double max);
		void buildSol() ;
		int getSolIndex();
		double diff(double E);
		void displaySol(vector<Point> sol);
		void findMinPot();
		//void plotXY(std::vector<double> x, std::vector<double> y);
		void savePotential() ;
		void scanForZeroRanges(int nZeros) ;
		void wait_for_key();
	 	double zbrent(function<double(double)>& func, double x1, double x2, double tol, bool);
	
	public:
		// these values can be overrided later
		double h     = 0.01;
		double xMin  = -15;
		double xMax  = 15;
		double tol   = 1e-9;
		double dEmin = 0.1;
		int nPoints  = 1000;

		void setPotential(vector<Point>);
		vector<Point> getPotential();
		void findSpectrum(int nEigen);
		Spectrum getSpectrum();

};


vector<Point> Numerov::getPotential() {
    return potential;
}

Spectrum Numerov::getSpectrum() {
	return spectrum;
}

// Van Wijngaarden–Dekker–Brent Method for finding root, from NR
#define ITMAX 600
#define EPS 1e-9
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double Numerov::zbrent(function<double(double)>& func, double x1, double x2, double tol, bool silent = false)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	
	if (!silent && ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)))
		cout << "?";
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c  = a;			//Rename a, b, c and adjust bounding interval d.
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.0*EPS*fabs(b)+0.5*tol; //Convergence check.
		xm = 0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;		//Attempt inverse quadratic interpolation.
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;		//Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				//Accept interpolation.
				d = p/q;
			} else {
				d = xm;
				//Interpolation failed, use bisection.
				e = d;
			}
		} else {		//Bounds decreasing too slowly, use bisection.
			d = xm;
			e = d;
		}
		a = b;			//Move last best guess to a.
		fa = fb;
		if (fabs(d) > tol1)	//Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = func(b);
	}
	cout<<"[WARN] Maximum number of iterations exceeded in zbrent, returning biggest value"<<endl;
	return x2;
}

// just returns the overall minima of the potential
void Numerov::findMinPot(){
	minPot = potential[0];

	for (int i = 0; i < potential.size(); i++) 
		if (minPot.y > potential[i].y) 
			minPot = potential[i];
		
    //cout << "minPlot " << minPot.y << " at " << minPot.x <<endl;
}

void Numerov::setPotential(vector<Point> v = vector<Point>(1)){
	// the potential can be read at this point from a file or similar, or can be hard coded.
    if(v.size() > 1) {
        potential = v;
    } else { // show an harmonic oscillator data by default
    	potential.resize(nPoints);
    	for (int i = 0; i < nPoints; i++) {
    		potential[i].x = xMin + i * (xMax - xMin) / (nPoints - 1);
    		potential[i].y = potential[i].x * potential[i].x;  // harmonic oscillator for testing
    	}
        cout << "[DEBUG] Using harmonic oscillator potential..." << endl;
    }

    nPoints = potential.size();
    xMin    = potential.front().x;
    xMax    = potential.back().x;
  	h = (xMax - xMin) / nPoints;

	 	int N = nPoints - 1;
	  // delete previous and create new ones
	  delete[] A; A = new double[N];
	  delete[] B; B = new double[N];
	  delete[] C; C = new double[N];
	  delete[] D; D = new double[N];

	  double w[N];
	  double h[N];
	  double ftt[N+1];

	  for (int i=0; i<N; i++) {
  	    w[i] = (potential[i + 1].x-potential[i].x);
  	    h[i] = (potential[i + 1].y-potential[i].y) / w[i];
	  }

	  ftt[0] = 0;
	  for (int i=0; i<N-1; i++)
	      ftt[i+1] = 3*(h[i+1]-h[i])/(w[i+1]+w[i]);

	  ftt[N] = 0;

	  for (int i=0; i<N; i++) {
  	    A[i] = (ftt[i+1]-ftt[i])/(6*w[i]);
  	    B[i] = ftt[i]/2;
  	    C[i] = h[i]-w[i]*(ftt[i+1]+2*ftt[i])/6;
  	    D[i] = potential[i].y;
	  }

  	// here we use the nice feature of c++ 11, the lambda functions
  	// oh dear, this is soo cool!
  	potFunc = [this] (double x){
  		// here we have to interpolate the potential since all the testing points are not in the grid!
  		// first get the closest left point
  		int index = 0;
  		for (int i = 0; i < potential.size(); i++){
    			if (potential[i].x > x) {
    				index = i - 1;
    				break;
    			}
    			// if we arrived to the end then it because x is the end point.
    			index = potential.size() - 2;
  		}

  		return A[index] * pow(x - potential[index].x, 3) +
  		       B[index] * pow(x - potential[index].x, 2) +
  		       C[index] * (x - potential[index].x) +
  		       D[index];
	};
}


void Numerov::displaySol(vector<Point> sol){
	for (int i = 0; i < sol.size(); ++i) {
		cout << i <<" "<< sol[i].x<< ", " << sol[i].y << endl;
	}
}

double Numerov::diff(double E){
	// first we find the matching point
    function<double(double)> shiftedPot = [this, &E](double x) {
		//cout<<"E=" << E << " shiftedPot(" << x <<")=" << E - potFunc(x) << endl;
		return E - potFunc(x);
	};
	
	// since we want to find preferently the right turning point we
	// look from the minimum position to the right
	double matchPoint = (xMax + xMin) / 2;//zbrent(shiftedPot, minPot.x, xMax, tol, true);
	//cout << "match point at " << matchPoint << " for E = " << E << endl;
	int matchPointIndex = int((matchPoint - xMin) / h);
	matchPointIndex = max(4, matchPointIndex);
	matchPointIndex = min(nPoints - 6, matchPointIndex);
	// for our specific problem it may be convenient just to set the match point close to the end
	matchPointIndex = nPoints - 6;
	// now we have to propagate left and right solutions until the matching point
	solLR.clear();
    solLR.resize(matchPointIndex + 1);
    solLR[0].x = xMin;
    solLR[0].y = 0;
    solLR[1].x = xMin + h;
    solLR[1].y = 0.0001;  // this number should not be important, but it has to be small to avoid overflow

	//cout << "matchPointIndex = " << matchPointIndex << " solLR.size() = " << solLR.size() << endl;
	double h2 = h * h;
    for(int i = 2; i < solLR.size(); i++)
    {
        double x2  = xMin + (i - 2) * h; 
        double x1  = xMin + (i - 1) * h;
        double x   = xMin + i * h;     
        double p1  = 2 - 0.83333333333333 * h2 * shiftedPot(x1);
        double p2  = 1 + 0.08333333333333 * h2 * shiftedPot(x2);
        double p3  = 1 + 0.08333333333333 * h2 * shiftedPot(x);
        solLR[i].y = (p1 * solLR[i-1].y - p2 * solLR[i-2].y) / p3;
        solLR[i].x = x;
    }
 
	solRL.clear(); 
	solRL.resize(potential.size() - solLR.size());
    //fill the array, propagating the solRLution from the Right to the Left
    solRL[solRL.size() - 1].x = xMax;
    solRL[solRL.size() - 1].y = 0;
    solRL[solRL.size() - 2].x = xMax - h;
    solRL[solRL.size() - 2].y = 0.0001;

    for(int i = solRL.size() - 3; i > -1; i--)
    {
        int j = solRL.size() - 1 - i;  // distance to the end
        double x2  = xMax - (j - 2) * h;   
        double x1  = xMax - (j - 1) * h;  
        double x   = xMax - j * h;       
        double p1  = 2 - 0.83333333333333 * h2 * shiftedPot(x1);
        double p2  = 1 + 0.08333333333333 * h2 * shiftedPot(x2);
        double p3  = 1 + 0.08333333333333 * h2 * shiftedPot(x);
        solRL[i].y = (p1 * solRL[i+1].y - p2 * solRL[i+2].y) / p3;
        solRL[i].x = x;
    }

	// now we have to find the log derivative at the matching point
	double v1 = solRL[0].y;
	double d2 = ((25/12)*solLR[matchPointIndex].y-4*solLR[matchPointIndex-1].y+3*solLR[matchPointIndex-2].y-(4/3)* solLR[matchPointIndex-3].y+(1/4)* solLR[matchPointIndex-4].y);
	double v2 = solLR[matchPointIndex].y;
	double d1 = -((25/12)*solRL[0].y-4*solRL[1].y+3*solRL[2].y-(4/3)*solRL[3].y+(1/4) * solRL[4].y);
	
	double logDer = (d1 * v2 - d2 * v1) / h;
	//logDer /= (ld1 * ld1 + ld2 * ld2);
	return logDer;
}



/**
 * Starting from the potential minima with a step of dEmin
 * computes diff(E), if a change of sign appears it means
 * that there is a zero thus store the range in the zeros vector.
 * These ranges are later used to find the root in a more refined phase.
 */
void Numerov::scanForZeroRanges(int nZeros) {
	findMinPot();
	double E = minPot.y;
	zeros.clear();
	double lastDiff = diff(E), newDiff;
	while(zeros.size() < nZeros){
		newDiff = diff(E);
		if(newDiff * lastDiff < 0) {
			Range range(E - dEmin, E);
			zeros.push_back(range);
			//cout << "zero in [" << range.eMin << ", " << range.eMax << "]" << endl;
		}
		lastDiff = newDiff;
		// just in case we hit a zero, which is not so unlikely to happen
		if(abs(lastDiff) < 1e-8) {
		    Range range(E - 1e-8, E + 1e-8);
		    zeros.push_back(range);
		    E += 1e-8;
		    //cout << "zero hit at " << E << endl;
		}
		
		E += dEmin;
	}
}

void Numerov::buildSol() {
	sol.clear();
	double scale = solRL[0].y / solLR[solLR.size() - 1].y;
	for (int i = 0; i < solLR.size(); i++) {
		Point p;
		p.x = solLR[i].x;
		p.y = scale * solLR[i].y;
		sol.push_back(p);
	}
	for (int i = 0; i < solRL.size(); i++) {
		sol.push_back(solRL[i]);
	}

	// normalize it
	double c = 0;
	for (int i = 1; i < sol.size(); i++) 
		c += sol[i].y * sol[i].y * (sol[i].x - sol[i - 1].x);
	
	for (int i = 0; i < sol.size(); i++) 
		sol[i].y /= sqrt(c);
	
	// we want to impose that the first maxima is alway positive
	// to avoid "jumps" while changing parameters in the potential
    auto derivative = [&](int i) {
        //https://en.wikipedia.org/wiki/Five-point_stencil
        return -sol[i + 2].y + 8 * sol[i + 1].y  - 8 * sol[i - 1].y + sol[i - 2].y;
    };
    
    double der = derivative(2);
    if(der < 0) // in this case change the overall sign
        for (int i = 0; i < sol.size(); i++) 
            sol[i].y *= -1;
}

// given a solution this find its number
// of "nodes" which actually labels the solution
int Numerov::getSolIndex() {
    int mp = solLR.size();
    auto derivative = [&](int i) {
        //https://en.wikipedia.org/wiki/Five-point_stencil
        return -sol[i + 2].y + 8 * sol[i + 1].y  - 8 * sol[i - 1].y + sol[i - 2].y;
    };
    
    int flips = -1;
    //cout << "node found at ";
    double lastDer = derivative(2);
    for(int i = 3; i < sol.size() - 6; i++) {
        double newDer = derivative(i);
        
        if(lastDer * newDer < 0) {
            lastDer = newDer;
            // do not count the artifical nodes near the matching point
            if(mp - 50 < i && i < mp + 50)
                continue;
            //cout << sol[i].x << " ";
            flips++;
        }
    }
    //cout << endl;
 
    return flips;   
}

void Numerov::findSpectrum(int nEigen){
    if(potential.size() < 2) {
        cout << "Please use the setPotential() function before using this one." << endl;
        return;
    }
    
	spectrum.clear();
	spectrum.potential = potential;
	scanForZeroRanges(nEigen);
	vector<Mode> modes;
	// check if a given index has been already computed
	auto hasBeenComputed = [&] (int index) {
		for (int i = 0; i < modes.size(); i++) {
		    if(modes[i].index == index)
		        return true;
		}
		  
		return false;
	};
	
	const int MAX_DIVISIONS = 100;
	std::function<Range(double, double, int)> explore = [&](double from, double to, int N = 10) {
	    //cout << "[DEBUG] exploring interval [" << from <<", " << to << ") N = " << N << endl;
	    // this is too much, there should be something wrong (root not in the interval)
	    if(N > MAX_DIVISIONS) {
	        cout << "[WARN] The number of divisions " << N << " has gone beyond the limit for interval [" << from <<", " << to << "], aborting..." << endl;
	        throw std::runtime_error("run time error");
	    }
	    
	    // divide the interval in N and look for sign changes
	    // if there is no success then attempts again with a finer grid
	    double h = (to - from) / N;
	    // find the first sign change
	    Point lastVal(from, diffFunc(from));
	    for(int i = 1; i < N; i++) {
	        double p = from + i * h;
	        Point newVal(p, diffFunc(p));
	        if(lastVal.y * newVal.y < 0) {
	            //cout << "New range for root [" << lastVal.x << ", "  << newVal.x << "]" << endl;
	            Range r(lastVal.x, newVal.x);
	            return r;
	        }
	        lastVal = newVal;
	    }
	        
	    // if we reach this point we need to do a more refined search
	    return explore(from, to, N + 50);
	};
	
	// find a solution in a range
	auto findSol = [&](Range r) {
	    //cout << "finding sol in [" << r.eMin << ", " << r.eMax << "]..." << endl; 
	    double E = zbrent(diffFunc, r.eMin, r.eMax, tol);
	    buildSol();
	    int n = getSolIndex();
	    //cout << "E(" << n << ") = " << E << endl;
	    Mode mode(E, sol, n);
	    return mode;
	};
	
	// using this for we get the right lowest nEigen eigenfunctions
	// most of the times, but sometimes we may skip some, see below.
	int attempts = 0;
	int const MAX_ATTEMPTS = 1;
	for (int i = 0; i < nEigen; i++) {
	    Mode m = findSol(zeros[i]);
	    modes.push_back(m);
	    if(m.index >= nEigen) {
	        //cout << "Possible bad index for eigenvalue " << m.index << endl;
	        //cout << "Bound found, index " << m.index << endl;
	        //break;
	    }
	}
	
	auto findSpecSol = [&](int index) {
	    int attempts = 0;
	    while(!hasBeenComputed(index)) {
	        int j = 0;
	        // get the upper bound
	        for(int i = 0; i < modes.size(); i++)
	            if(modes[i].index > index) {
	                j = i;
	                break;
	            }
	   
	        double lowLim = j == 0 ? minPot.y : modes[j - 1].energy;   // default is the lower limit, works for the ground state
	        double minDiff = 1e-7;//(modes[i].energy - lowLim) / 100;
	        // get the right interval to look for another root
	        Range r = explore(lowLim + minDiff, modes[j].energy - minDiff, 20);
	        Mode m = findSol(r);
	        cout << "new mode found " << m.index << endl;
	        modes.push_back(m);
	        // we need to sort the modes and remove repeated after this insertion
	        auto comp = [&](Mode m1, Mode m2) -> bool {
	           return m1.energy < m2.energy;
	        };
	        auto rm = [&](Mode m1, Mode m2) -> bool {
	            return m1.index == m2.index;
	        };
	        std::sort(modes.begin(), modes.end(), comp);
	        std::unique(modes.begin(), modes.end(), rm);
	       
	        attempts++;
	        if(attempts > 10) {
	            cout << "[WARN] Too many attempts while finding eigenvalue " << index << ", computation is compromised." << endl;
	            break;
	        }
	    }  
	};
	
	// second recovery strategy
	// now we need to check if we get the right eigenfunctions
	// if we were asked to return the first 4 we don't want to return the 1,2,3 and 5!
	for (int i = 0; i < nEigen; i++) 
	    //TODO: finish this nice stuff
	    if(false && !hasBeenComputed(i)) { // we need to look for a missing mode
	        cout << "Looking back for E(" << i << ")..." << endl;
	        try {
	            findSpecSol(i);
	        }catch(std::exception &e) {
	            cout<<"Caught exception: "<<e.what()<<"\n";
	            break;
	        }
	    }   
	
	// finally safelly add all the modes found to the spectrum (already in a nice way)
	for (int i = 0; i < nEigen; i++)
	    if(i < modes.size())
            spectrum.addMode(modes[i]);
	    else{
	        Mode m;
	        spectrum.addMode(m);
	    }
}

void Numerov::savePotential() {
    ofstream f("potential.dat");
	for (double x = xMin; x < xMax; x += 0.01) {
		f << x << " " << potFunc(x) << endl;
	}
	f.close();
}

double Numerov::bisection(function<double(double)> diffFunc, double min, double max){
	double E;
	int i=0;			     /* counter for iterations */
	do{
		i++;  
		E=(max+min)/2.0;		     /* divide energy range */
		if (diffFunc(max)*diffFunc(E)>0) 
			max=E;   /* the bisection algorithm */
		else 
			min=E;
	}while(fabs(diffFunc(E))>EPS);
	return E;
}

#endif
