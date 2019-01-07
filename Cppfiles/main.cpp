#include "random.h"
#include <iostream>
#include <vector>

struct parameters{
    double b1;
    double b2;
    double d1;
    double d2;
    double k1;
    double k2;
};

struct popstate{
    popstate(const int &n1, const int &n2, const parameters &pars):n1(n1),n2(n2){CalcRates(pars);};
    int n1, n2;
    double b1, b2, d1, d2;

    private:
    void CalcRates(const parameters &pars){
            b1 = (n1>0 && n1+n2-1<pars.k1) ? 
                pars.b1*(double)n1*(1.0-((double)(n1+n2))/pars.k2) : 0.0;
            b2 = (n2>0 && n1+n2-1<pars.k2) ? 
                pars.b2*(double)n2*(1.0-((double)(n1+n2))/pars.k1) : 0.0;
            d1 = pars.d1*n1;
            d2 = pars.d2*n2;
        }
};

class Gillespie{
    Gillespie(const int &n1, const int &n2, const parameters &pars):pars(pars){
        popstate init(n1,n2,pars);
        timelaps.push_back(init);
    }

    void MCstep(){
        timelaps.back();
    }

    private:
    const parameters pars;
    
    std::vector<popstate> timelaps;
};

int main(){
    parameters pars;
    pars.b1 = 0.2;
    pars.b2 = 0.6;
    pars.d1 = 0.05;
    pars.d2 = 0.05;
    pars.k1 = 50;
    pars.k2 = 50;

    return 0;
}

