#include "random.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <thread>
#include <assert.h>     /* assert */
#include <Eigen/Dense>

using namespace Eigen;

enum cases {birth1, birth2, death1, death2};

struct parameters{
    public:
    parameters(double b1, double b2, double d1, double d2, int k1, int k2, int tmax, int nrep):
        b1(b1),
        b2(b2),
        d1(d1),
        d2(d2),
        k1(k1),
        k2(k2),
        tmax(tmax),
        nrep(nrep){}
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

struct individual1{
    individual1():x1(rnd::uniform()),x2(rnd::uniform()){}
    individual1(const int &x1, const int &x2):x1(x1),x2(x2){};
    const double x1,x2;
};
struct individual2{
    individual2():x1(rnd::uniform()),x2(rnd::uniform()){}
    individual2(const int &x1, const int &x2):x1(x1),x2(x2){};
    const double x1,x2;
};

MatrixXd distanceSquared(MatrixXd n1, MatrixXd n2){
    const int n1cols = n1.cols(), n1rows = n1.rows();
    const int n2cols = n2.cols(), n2rows = n2.rows();
    assert(n1cols == n2cols);
    return  (n1.array()*n1.array()).matrix()*MatrixXd::Ones(n1cols,n2rows)+ 
            MatrixXd::Ones(n1rows,n1cols)*(n2.transpose().array()*n2.transpose().array()).matrix() - 
            2*n1*n2.transpose();
}

class GillespieMoments{
    public:
    GillespieMoments(std::vector<individual1> t1, std::vector<individual2) t2: pars(pars){
        t=0.0; 
        twrite = 0; 
        n1 = n1types;
        n2 = n2types;
        // type1 = ((MatrixXd::Random(n1types,2).array()+1.0)/2.0).matrix(); 
        // type2 = ((MatrixXd::Random(n2types,2).array()+1.0)/2.0).matrix();
    }

    bool MCstep(std::vector<moments> &M){

        distanceSquared(type1, type2);
        distanceSquared(type1, type1);
        distanceSquared(type2, type2);

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
    std::vector<individual1> type1;
    std::vector<individual2> type2;
    int n1, n2;
    const parameters pars;
    int twrite;

    void writedata(const int &n1, const int &n2, std::vector<moments> &M , const int &tmax){
        for(;twrite < tmax; ++twrite){
            M[twrite].update(n1,n2);
        }
    }
};

int main(){
    rnd::set_seed();
    parameters pars(0.1,0.1,0.05,0.05,100,100,100,1);

    const int n1_init = 1;
    const int n2_init = 1;

    std::vector<moments> M(pars.tmax);
    for(int i = 0; i < pars.nrep; ++i){
        GillespieMoments test(n1_init,n2_init, pars); 
        while(test.MCstep(M));
    }
    return 0;
}