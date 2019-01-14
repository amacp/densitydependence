#include "random.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
//#include "WolframLibrary.h"
#include <iostream>

//EXTERN_C DLLEXPORT mint WolframLibrary_getVersion(){return WolframLibraryVersion;}
//EXTERN_C DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {return 0;}
//EXTERN_C DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {}

//include openmp

enum cases { birth1, birth2, death1, death2};

struct parameters{
    public:
    double b1, b2, d1, d2, k1, k2;
    int tmax;
    int nrep;
};

struct popstate{
    public:
    popstate(const int &n1, const int &n2, const double &t, const parameters &pars):n1(n1),n2(n2),time(t) {CalcRates(pars);};
    double time;
    int n1, n2;

    inline double gettotalrate() {return totalrate;}
    inline double getdist(const int &type) {return dist[type];}
void CalcRates(const parameters &pars){
        dist[birth1] = (n1>0 && n1+n2-1<pars.k2) ? 
            pars.b1*(double)n1*(1.0-((double)(n1+n2))/pars.k2) : 0.0;
        dist[birth2] = (n2>0 && n1+n2-1<pars.k1) ? 
            pars.b2*(double)n2*(1.0-((double)(n1+n2))/pars.k1) : 0.0;
        dist[death1] = pars.d1*(double)n1;
        dist[death2] = pars.d2*(double)n2;
        totalrate = dist[birth1]+dist[birth2]+dist[death1]+dist[death2];
    };

    private:
    double totalrate;
    double dist[4];
};

class Gillespie{
    public:
    Gillespie(const int &n1, const int &n2, const parameters &pars, int** counts_type1, int** counts_type2):pars(pars), counts_type1(counts_type1), counts_type2(counts_type2){
        popstate init(n1,n2,0,pars);
        timelaps.push_back(init);
        tmin = 0;
    }

    bool MCstep(){
        static int counter = 0;
        popstate focal = timelaps.back();
        rnd::discrete_distribution dist(4);
        dist[death1] = focal.getdist(death1);
        dist[death2] = focal.getdist(death2);
        dist[birth1] = focal.getdist(birth1);
        dist[birth2] = focal.getdist(birth2);
        const int n1 = focal.n1;
        const int n2 = focal.n2;
        std::cout << counter << " ";
        counter++;
        if(focal.n1 + focal.n2 == 0){
            writedata(n1,n2,pars.tmax);
            return false;
        }
        double t = focal.time+rnd::exponential(1.0/focal.gettotalrate());
        if(t > pars.tmax){
            t = pars.tmax;
            writedata(n1,n2,t);
            return false;
        } 

        writedata(n1,n2,t);
        // determine event
        popstate next(n1,n2,t,pars);
        const int event = dist.sample();

        switch(event) {
            case birth1 : next.n1+=1; break;
            case birth2 : next.n2+=1; break;
            case death1 : next.n1-=1; break;
            case death2 : next.n2-=1; break;
        } 

        next.CalcRates(pars);
        timelaps.push_back(next);
         
        return true;
    }

    private:
    std::vector<popstate> timelaps;
    const parameters pars;
    int** counts_type1;
    int** counts_type2;
    int tmin;
    void writedata(const int &n1, const int &n2, const double &tmax){
        const int inttmax = (int)tmax;
        if(tmin < inttmax){
            for(;tmin < inttmax; ++tmin){
                ++counts_type1[tmin][n1];
                ++counts_type2[tmin][n2];
            }
            tmin = inttmax;
        }
    }
};

int main(){
    rnd::set_seed();
    const int n1_init = 1;
    const int n2_init = 1;
    parameters pars;
    pars.b1 = 0.2;
    pars.b2 = 0.6;
    pars.d1 = 0.05;
    pars.d2 = 0.05;
    pars.k1 = 10;
    pars.k2 = 2;
    pars.tmax = 100;
    pars.nrep = 10000;
    assert(n1_init <= pars.k1);
    assert(n2_init <= pars.k2);
    
    const int rowCount = pars.tmax;
    const int colCount1 = pars.k1;
    const int colCount2 = pars.k2;
    int** counts_type1 = new int*[rowCount];
    int** counts_type2 = new int*[rowCount];
    for(int i = 0; i < rowCount; ++i){
        counts_type1[i] = new int[colCount1];
        counts_type2[i] = new int[colCount2];
    }

    for(int i = 0; i < pars.k1; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            counts_type1[j][i] =0;
        }
    }
    for(int i = 0; i < pars.k2; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            counts_type2[j][i] =0;
        }
    }

    // Run simulations (parallel)?
    #pragma omp parallel
    #pragma omp for
    for(int i = 0; i < pars.nrep; ++i){
        Gillespie test(n1_init,n2_init,pars, counts_type1, counts_type2);  
        while(test.MCstep());
    }
    
    std::ofstream outputcsv1("data_type1.csv");
    outputcsv1.fill(',');
    std::ofstream outputcsv2("data_type2.csv");
    outputcsv2.fill(',');

    for(int i = 0; i < pars.k1; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            outputcsv1 << (double)counts_type1[j][i]/(double)pars.nrep << outputcsv1.fill();
        }

    outputcsv1 << std::endl;
    }

    for(int i = 0; i < pars.k2; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            outputcsv2 << (double)counts_type2[j][i]/(double)pars.nrep << outputcsv2.fill();
        }
    outputcsv2 << std::endl;
    }

    return 0;
}

/*
EXTERN_C DLLEXPORT int myFunction2(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
    
    int out = MArgument_getInteger(Args[0]) + MArgument_getInteger(Args[1]);
    MArgument_setInteger(Res,out);
    return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int myFunction(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){

    rnd::set_seed();
    const int n1_init = MArgument_getInteger(Args[0]);
    const int n2_init = MArgument_getInteger(Args[0]);
    parameters pars;
    pars.b1 = 0.2;
    pars.b2 = 0.6;
    pars.d1 = 0.05;
    pars.d2 = 0.05;
    pars.k1 = 10;
    pars.k2 = 10.0/3.0;
    pars.tmax = 100;
    pars.nrep = 100;
    assert(n1_init <= pars.k1);
    assert(n2_init <= pars.k2);
    
    const int rowCount = pars.tmax;
    const int colCount1 = pars.k1;
    const int colCount2 = pars.k2;
    int** counts_type1 = new int*[rowCount];
    int** counts_type2 = new int*[rowCount];
    for(int i = 0; i < rowCount; ++i){
        counts_type1[i] = new int[colCount1];
        counts_type2[i] = new int[colCount2];
    }

    for(int i = 0; i < pars.k1; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            counts_type1[j][i] =0;
        }
    }
    for(int i = 0; i < pars.k2; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            counts_type2[j][i] =0;
        }
    }

    // Run simulations (parallel)?
    for(int i = 0; i < pars.nrep; ++i){
        Gillespie test(n1_init,n2_init,pars, counts_type1, counts_type2);  
        while(test.MCstep());
    }


    MArgument_setInteger(Res,100);

 return LIBRARY_NO_ERROR;

}


 Graveyard 


    for(int j = 0; j < pars.tmax; ++j){
        double avg = 0.0; 
        for(int i = 0; i < pars.k1; ++i){
            avg += i*counts_type1[j][i];
            outputcsv1 << (double)counts_type1[j][i]/(double)pars.nrep << outputcsv1.fill();
        }

    outputcsv1 << std::endl;
    }

    std::cout << "\n \n";

    for(int i = 0; i < pars.k2; ++i){
        for(int j = 0; j < pars.tmax; ++j){
            outputcsv2 << (double)counts_type2[j][i]/(double)pars.nrep << outputcsv2.fill();
        }
    outputcsv2 << std::endl;
    }

    return 0;
*/

