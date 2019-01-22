#include "random.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <thread>
#include "WolframLibrary.h"

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion(){return WolframLibraryVersion;}
EXTERN_C DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {return 0;}
EXTERN_C DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {}

//include openmp
enum cases {birth1, birth2, death1, death2};

struct parameters{
    public:
    parameters(MArgument *Args):
        b1(MArgument_getReal(Args[2])),
        b2(MArgument_getReal(Args[3])),
        d1(MArgument_getReal(Args[4])),
        d2(MArgument_getReal(Args[5])),
        k1(MArgument_getInteger(Args[6])),
        k2(MArgument_getInteger(Args[7])),
        tmax(MArgument_getInteger(Args[8])),
        nrep(MArgument_getInteger(Args[9])){}
    const double b1, b2, d1, d2, k1, k2;
    const int tmax;
    const int nrep;
};

struct moments{
    public:
    moments():firstM(2,0), secondM(3,0), thirdM(4,0), firstM_moment(2,0.0), secondM_Cmoment(3,0.0), thirdM_Cmoment(5,0.0)
    {
    };

    void update(const int &n1, const int &n2){
        firstM[0] += n1;
        firstM[1] += n2;
        secondM[0] += n1*n1;
        secondM[1] += n1*n2;
        secondM[2] += n2*n2;
        thirdM[0] += n1*n1*n1;
        thirdM[1] += n1*n1*n2;
        thirdM[2] += n1*n2*n2;
        thirdM[3] += n2*n2*n2;
        ++counter;
    }

    void calcCentralMoments(){
        std::vector<double> secondM_moment(3,0.0);
        std::vector<double> thirdM_moment(4,0.0);

        firstM_moment[0] = (long double)firstM[0]/(long double)counter;
        firstM_moment[1] = (long double)firstM[1]/(long double)counter;
        secondM_moment[0] = (long double)secondM[0]/(long double)counter;
        secondM_moment[1] = (long double)secondM[1]/(long double)counter;
        secondM_moment[2] = (long double)secondM[2]/(long double)counter;
        thirdM_moment[0] = (long double)thirdM[0]/(long double)counter;
        thirdM_moment[1] = (long double)thirdM[1]/(long double)counter;
        thirdM_moment[2] = (long double)thirdM[2]/(long double)counter;
        thirdM_moment[3] = (long double)thirdM[3]/(long double)counter;

        secondM_Cmoment[0] = secondM_moment[0] - firstM_moment[0]*firstM_moment[0];
        secondM_Cmoment[1] = secondM_moment[1] - firstM_moment[0]*firstM_moment[1];
        secondM_Cmoment[2] = secondM_moment[2] - firstM_moment[1]*firstM_moment[1];
        thirdM_Cmoment[0] = thirdM_moment[0] - firstM_moment[0]*firstM_moment[0]*firstM_moment[0];
        thirdM_Cmoment[1] = thirdM_moment[1] - secondM_moment[0]*firstM_moment[1]+2*firstM_moment[0]*firstM_moment[0]*firstM_moment[1]-2*firstM_moment[0]*secondM_moment[1];
        thirdM_Cmoment[2] = thirdM_moment[2] - secondM_moment[2]*firstM_moment[0]+2*firstM_moment[0]*firstM_moment[1]*firstM_moment[1]-2*firstM_moment[1]*secondM_moment[1];
        thirdM_Cmoment[3] = thirdM_moment[3] - firstM_moment[1]*firstM_moment[1]*firstM_moment[1];
    }

    double returnmean(const int &i){
        return firstM_moment[i];
    }
    double retursecondCM(const int &i){
        return secondM_Cmoment[i];
    }
    double returnthirdCM(const int &i){
        return thirdM_Cmoment[i];
    }

    private:
    int counter = 0;
    std::vector<long int> firstM;
    std::vector<long int> secondM;
    std::vector<long int> thirdM;
    std::vector<double> firstM_moment;
    std::vector<double> secondM_Cmoment;
    std::vector<double> thirdM_Cmoment;
};

class Gillespie{
    public:
    Gillespie(const int &n1, const int &n2, const parameters &pars): pars(pars), n1(n1),n2(n2){
    t=0.0; 
    twrite = 0; 
    }

    bool MCstep(std::vector<std::vector<int> > &data1, std::vector<std::vector<int> > &data2){
        //assert(n1+n2>0);
        //assert(pars.k2>0); assert(pars.k1>0);
        rnd::discrete_distribution dist(4);
        dist[birth1] = (n1>0 && n1+n2-1<pars.k2) ? pars.b1*(double)n1*(1.0-((double)(n1+n2))/pars.k2) : 0.0;
        dist[birth2] = (n2>0 && n1+n2-1<pars.k1) ? pars.b2*(double)n2*(1.0-((double)(n1+n2))/pars.k1) : 0.0;
        dist[death1] = pars.d1*(double)n1;
        dist[death2] = pars.d2*(double)n2;

        tprime =t+ rnd::exponential(1.0/(dist[death1]+dist[death2]+dist[birth1]+dist[birth2]));
        
        if(tprime >= pars.tmax){
            writedata(n1,n2,data1, data2, pars.tmax);
            return false;
        } 
        else{
            switch(dist.sample()) {
                case birth1 : n1+=1; break;
                case birth2 : n2+=1; break;
                case death1 : n1-=1; break;
                case death2 : n2-=1; break;
            }

            if(n1+n2==0){               
                writedata(n1,n2,data1, data2, pars.tmax);
                return false;
            }

            if(std::floor(t) < std::floor(tprime)){
                writedata(n1,n2,data1, data2, (int)std::floor(tprime));
            }
        }
        t = tprime;
        
        return true;
    }

    private:
    double t;
    double tprime;
    int n1;
    int n2;
    const parameters pars;
    int twrite;

    void writedata(const int &n1, const int &n2,std::vector<std::vector<int> > &data1, std::vector<std::vector<int> > &data2, const int &tmax){
        for(;twrite < tmax; ++twrite){
            ++data1[twrite][n1];
            ++data2[twrite][n2];
        }
    }

    // int** counts_type1;
    // int** counts_type2;
};

class GillespieMoments{
    public:
    GillespieMoments(const int &n1, const int &n2, const parameters &pars): pars(pars), n1(n1),n2(n2){
    t=0.0; 
    twrite = 0; 
    }

    bool MCstep(std::vector<moments> &M){
        //assert(n1+n2>0);
        //assert(pars.k2>0); assert(pars.k1>0);
        rnd::discrete_distribution dist(4);
        dist[birth1] = (n1>0 && n1+n2-1<pars.k2) ? pars.b1*(double)n1*(1.0-((double)(n1+n2))/pars.k2) : 0.0;
        dist[birth2] = (n2>0 && n1+n2-1<pars.k1) ? pars.b2*(double)n2*(1.0-((double)(n1+n2))/pars.k1) : 0.0;
        dist[death1] = pars.d1*(double)n1;
        dist[death2] = pars.d2*(double)n2;

        tprime =t+ rnd::exponential(1.0/(dist[death1]+dist[death2]+dist[birth1]+dist[birth2]));
        
        if(tprime >= pars.tmax){
            writedata(n1,n2, M, pars.tmax);
            return false;
        } 
        else{
            switch(dist.sample()) {
                case birth1 : n1+=1; break;
                case birth2 : n2+=1; break;
                case death1 : n1-=1; break;
                case death2 : n2-=1; break;
            }

            if(n1+n2==0){               
                writedata(n1,n2, M, pars.tmax);
                return false;
            }

            if(std::floor(t) < std::floor(tprime)){
                writedata(n1,n2, M, (int)std::floor(tprime));
            }
        }
        t = tprime;
        
        return true;
    }
    private:
    double t;
    double tprime;
    int n1;
    int n2;
    const parameters pars;
    int twrite;

    void writedata(const int &n1, const int &n2, std::vector<moments> &M , const int &tmax){
        for(;twrite < tmax; ++twrite){
            M[twrite].update(n1,n2);
        }
    }
};

bool verifypars(const int &n1_init, const int &n2_init, const parameters &pars, WolframLibraryData &libData){
    if(pars.k1 < n1_init){
        libData->Message("k1<n1");
        return true;
    };
    if(pars.k2 < n2_init){
        libData->Message("k2<n2");
        return true;
    };
    if(pars.b1 <= 0.0){
        libData->Message("b1<=0");
        return true;
    };
    if(pars.b2 <= 0.0){
        libData->Message("b2<=0");
        return true;
    };
    if(pars.d1 <= 0.0){
        libData->Message("d1<=0");
        return true;
    };
    if(pars.d2 <= 0.0){
        libData->Message("d2<=0");
        return true;
    };
    if(pars.tmax < 1){
        libData->Message("tmax_lowerbound");
        return true;
    };
    return false; 
}

EXTERN_C DLLEXPORT int myFunction(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
    rnd::set_seed();
    parameters pars(Args);

    const int n1_init = MArgument_getInteger(Args[0]);
    const int n2_init = MArgument_getInteger(Args[1]);

    if(verifypars(n1_init,n2_init, pars, libData)){
        return LIBRARY_FUNCTION_ERROR;
    };

    std::vector<moments> M(pars.tmax);
    //std::thread t[2];
    for(int i = 0; i < pars.nrep; ++i){
        GillespieMoments test(n1_init,n2_init, pars); 
        while(test.MCstep(M));
    }

    //Create output for mathematica
    MTensor T0;
    mint type = MType_Real;
    mint dims[2];
    mint rank = 2;
    int err;
    dims[0] = pars.tmax;
    dims[1] = 10;
    err = libData->MTensor_new( type, rank, dims, &T0);

    for ( mint i = 1; i <= pars.tmax && !err; i++) {
        M[i-1].calcCentralMoments();
        mint lens[2];
        lens[0] = i;

        lens[1] = 1;
        err = libData->MTensor_setReal( T0, lens, (double)i);
        lens[1] = 2;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnmean(0));
        lens[1] = 3;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnmean(1));
        lens[1] = 4;
        err = libData->MTensor_setReal( T0, lens, M[i-1].retursecondCM(0));
        lens[1] = 5;
        err = libData->MTensor_setReal( T0, lens, M[i-1].retursecondCM(1));
        lens[1] = 6;
        err = libData->MTensor_setReal( T0, lens, M[i-1].retursecondCM(2));
        lens[1] = 7;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnthirdCM(0));
        lens[1] = 8;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnthirdCM(1));
        lens[1] = 9;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnthirdCM(2));
        lens[1] = 10;
        err = libData->MTensor_setReal( T0, lens, M[i-1].returnthirdCM(3));
    }

    MArgument_setMTensor(Res, T0);

    return LIBRARY_NO_ERROR; 
}
