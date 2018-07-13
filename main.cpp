#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <vector>
#include "random.h"
#include "utils.h"

using namespace std;
//global variables
char filebase[200];
char filename[200];
//output files
ofstream out_Pars;
ofstream out_Data;
ofstream out_xy;
// Classes
class ind;
class gamete;
class mPop
{
    public:
        int sel,geo,init,nLoci,gamma,n0;
        double *aVec,*dVec,mu,r,rho,delta,sigmaEta,sigmaMate,sigmaMig,zMax,zMin,zBar,*pVec;
        vector<ind*> inds;
        vector<gamete*> gametes;
        //functions
        mPop();
        void input(int* gmaxPtr,int*sampFreqPtr);
        void initialize();
        void calcZlimit();
        void reproduction();
        void survivalDI();
        void survivalDD();
        void census(int gen);
};
class ind
{
    public:
    int**genome;
    double z,x,y,eta,Wrepro,gammai,ri,deltai;
    vector<double> neighbors;
    //vector<double> pMateVec;
    //function
    ind();
    ind(mPop metaPop,double x, double y);
    ind(gamete gam1,gamete gam2,mPop metaPop);
    void findNeighbors(mPop metaPop);
    void calcEta(mPop metaPop);
    void calcZdelta(mPop metaPop);
    void calcWrepro(mPop metaPop);
    //void calcpMate(mPop metaPop);
};
class gamete
{
    public:
    double x,y;
    vector<double> neighbors;
    vector<double> pMateVec;
    int*genome;
    //functions
    gamete(mPop metaPop,int parent);
    int findMate(mPop metaPop);
};
double randomNormal(double mu, double sigma);
int* shuffleList(int length);

//MAIN
int main()
{
   
    //random seed
    srand((unsigned int)time(NULL));//Initializing the random seed.
    //setting up output files
    sprintf(filebase,"/Users/ailenemacpherson/Documents/VisualStudio/density_dependence");
    sprintf(filename,"%s/out.csv",filebase);
    out_Pars.open(filename);
    sprintf(filename,"%s/data.csv",filebase);
    out_Data.open(filename);
    sprintf(filename,"%s/xy.csv",filebase);
    out_xy.open(filename);
    //Variables
    mPop metaPop;
    int gmax,sampFreq;
    //Running
    metaPop.input(&gmax,&sampFreq);
    metaPop.initialize();
    for(int g=0;g<=gmax;g++)
    {
        if((g%sampFreq)==0) metaPop.census(g);
        //cout<<"Reproduction"<<endl;
        metaPop.reproduction();
        //cout<<"Survival"<<endl;
        if(metaPop.sel==1) metaPop.survivalDI();//density-independent
        else if(metaPop.sel==2) metaPop.survivalDD();//density-dependent
    }
    out_Pars.close();
    out_Data.close();
    out_xy.close();
    cout<<"done!"<<endl;
}
mPop::mPop()
{
    aVec=NULL;dVec=NULL;
}
void mPop::input(int* gmaxPtr,int*sampFreqPtr)
{
	ifstream inFile;
	sprintf(filename,"%s/Input.txt",filebase);//Makes the parameter file name.
	inFile.open(filename);
	string myString;
	myString="Nothing";
	while(myString != "(gmax):" && inFile.good())	{inFile>>myString;}
	inFile>>(*gmaxPtr);
	out_Pars<<"gmax:, "<<(*gmaxPtr)<<endl;
	while(myString != "(sampFreq):" && inFile.good())	{inFile>>myString;}
	inFile>>(*sampFreqPtr);
	out_Pars<<"sampFreq:, "<<(*sampFreqPtr)<<endl;
	while(myString != "(sel):" && inFile.good())	{inFile>>myString;}
	inFile>>sel;
	out_Pars<<"sel:, "<<sel<<endl;    
	while(myString != "(geo):" && inFile.good())	{inFile>>myString;}
	inFile>>geo;
	out_Pars<<"geo:, "<<geo<<endl; 
	while(myString != "(init):" && inFile.good())	{inFile>>myString;}
	inFile>>init;
	out_Pars<<"init:, "<<init<<endl; 
	while(myString != "(n0):" && inFile.good())	{inFile>>myString;}
	inFile>>n0;
	out_Pars<<"n0:, "<<n0<<endl;
	while(myString != "(nLoci):" && inFile.good())	{inFile>>myString;}
	inFile>>nLoci;
	out_Pars<<"nLoci:, "<<nLoci<<endl;
	while(myString != "(aVec):" && inFile.good())	{inFile>>myString;}
    aVec=new double[nLoci];
    out_Pars<<"aVec: ,";
    for(int l=0;l<nLoci;l++)
    {
       inFile>>aVec[l]; 
       out_Pars<<aVec[l]<<",";
    }
    out_Pars<<endl;
	while(myString != "(dVec):" && inFile.good())	{inFile>>myString;}
    dVec=new double[nLoci];
    out_Pars<<"dVec: ,";
    for(int l=0;l<nLoci;l++)
    {
       inFile>>dVec[l]; 
       out_Pars<<dVec[l]<<",";
    }
    out_Pars<<endl;
	while(myString != "(mu):" && inFile.good())	{inFile>>myString;}
	inFile>>mu;
	out_Pars<<"mu:, "<<mu<<endl;
	while(myString != "(r):" && inFile.good())	{inFile>>myString;}
	inFile>>r;
	out_Pars<<"r:, "<<r<<endl;
	while(myString != "(rho):" && inFile.good())	{inFile>>myString;}
	inFile>>rho;
	out_Pars<<"rho:, "<<rho<<endl;
	while(myString != "(gamma):" && inFile.good())	{inFile>>myString;}
	inFile>>gamma;
	out_Pars<<"gamma:, "<<gamma<<endl;
	while(myString != "(delta):" && inFile.good())	{inFile>>myString;}
	inFile>>delta;
	out_Pars<<"delta:, "<<delta<<endl;
	while(myString != "(sigmaEta):" && inFile.good())	{inFile>>myString;}
	inFile>>sigmaEta;
	out_Pars<<"sigmaEta:, "<<sigmaEta<<endl;
	while(myString != "(sigmaMate):" && inFile.good())	{inFile>>myString;}
	inFile>>sigmaMate;
	out_Pars<<"sigmaMate:, "<<sigmaMate<<endl;
	while(myString != "(sigmaMig):" && inFile.good())	{inFile>>myString;}
	inFile>>sigmaMig;
	out_Pars<<"sigmaMig:, "<<sigmaMig<<endl;
}
void mPop::initialize()
{
    double x,y;
    if(geo==1)//Metapopulation
    {
        for(int i=0;i<n0;i++)
        {
            x=rand()/((double)RAND_MAX);
            y=rand()/((double)RAND_MAX);
            inds.push_back(new ind((*this),x,y));
        }
    }
    else{cout<<"Not set up yet!"<<endl; getchar();}
    calcZlimit();
    pVec=new double[nLoci];
}
void mPop::calcZlimit()
{
    zMax=0;zMin=0;
    for(int l=0;l<nLoci;l++)
    {
        zMax+=max(max(0.0,aVec[l]),2.0*aVec[l]+dVec[l]);
        zMin+=min(min(0.0,aVec[l]),2.0*aVec[l]+dVec[l]);
    }
    out_Pars<<"zMax: "<<zMax<<", zMin: "<<zMin<<endl;
}
void mPop::census(int gen)
{
    cout<<"gen: "<<gen<<endl;
    //calculating allele frequency
    for(int l=0;l<nLoci;l++) pVec[l]=0;
    for(int i=0;i<inds.size();i++)
    {
        for(int l=0;l<nLoci;l++)
        {
            pVec[l]+=inds[i]->genome[0][l]+inds[i]->genome[1][l];
        }
    }
    for(int l=0;l<nLoci;l++)
    {
        pVec[l]/=(2.0*inds.size());
    }
    for(int l=0;l<nLoci-1;l++)
    {
        out_Data<<pVec[l]<<",";
    }
    out_Data<<pVec[nLoci-1]<<endl;
    //Printing x,y coordinates and z values
    //cout<<"Number of individuals: "<<inds.size()<<endl;getchar();
    for(int i=0;i<inds.size()-1;i++) out_xy<<inds[i]->x<<",";
    out_xy<<inds[inds.size()-1]->x<<endl;
    for(int i=0;i<inds.size()-1;i++) out_xy<<inds[i]->y<<",";
    out_xy<<inds[inds.size()-1]->y<<endl;
    for(int i=0;i<inds.size()-1;i++) out_xy<<(inds[i]->z-zMin)/(zMax-zMin)<<",";
    out_xy<<(inds[inds.size()-1]->z-zMin)/(zMax-zMin)<<endl;
}
void mPop::reproduction()
{
    //Making Gametes
    int nGametes,rand1,rand2;
    for(int i=0;i<inds.size();i++)
    {
        inds[i]->findNeighbors((*this));
        inds[i]->calcWrepro((*this));
        nGametes=rnd::poisson(2.0*(inds[i]->Wrepro));
        for(int g=0;g<nGametes;g++)
        {
            gametes.push_back(new gamete((*this),i));
        }
    }
    //Pairing gametes
    inds.clear();
    while(gametes.size()>2)
    {
        nGametes=gametes.size();
        rand1=rand()%gametes.size();
        rand2=gametes[rand1]->findMate((*this));
       // cout<<nGametes<<" ,"<<rand1<<", "<<rand2<<endl;
        //cout<<"making ind"<<endl;
        inds.push_back(new ind((*gametes[rand1]),(*gametes[rand2]),(*this)));
        //removing gametes
        //cout<<"removing gametes"<<endl;
        gametes.erase(gametes.begin()+rand1);
        if(rand1!=rand2) gametes.erase(gametes.begin()+(rand2-1));
    }
}
void mPop::survivalDI()
{
    for(int i=0;i<inds.size();i++)
    {
        inds[i]->calcZdelta((*this));
        //cout<<"delti: "<<inds[i]->deltai<<endl;getchar();
        if(rand()/((double)RAND_MAX)<inds[i]->deltai)
        {
            inds.erase(inds.begin()+i);
            i--;
        }
    }
}
void mPop::survivalDD()
{
    for(int i=0;i<inds.size();i++)
    {
        //inds[i]->calcZdelta((*this));
        if(rand()/((double)RAND_MAX)<delta)
        {
            inds.erase(inds.begin()+i);
            i--;
        }
    }
}
ind::ind()
{
    genome=NULL;
}
ind::ind(mPop metaPop,double xIn, double yIn)
{
    x=xIn;y=yIn;
    genome=new int*[2];
    for(int c=0;c<2;c++) 
    {
        genome[c]=new int[metaPop.nLoci];
        for(int l=0;l<metaPop.nLoci;l++)
        {
            //genome[c][l]=0;
            genome[c][l]=rand()%2;
        }
    }
    calcZdelta(metaPop);
}
ind::ind(gamete gam1,gamete gam2,mPop metaPop)
{
    genome=new int*[2];
    //x=(gam1.x+gam2.x)/2.0;y=(gam1.y+gam2.y)/2.0;
    x=gam1.x;y=gam1.y;
    x=x+randomNormal(0.0,metaPop.sigmaMig);
    y=y+randomNormal(0.0,metaPop.sigmaMig);
    for(int c=0;c<2;c++) genome[c]=new int[metaPop.nLoci];
    for(int l=0;l<metaPop.nLoci;l++)
    {
        genome[0][l]=gam1.genome[l];
        genome[1][l]=gam2.genome[l];
    }
}
void ind::findNeighbors(mPop metaPop)
{
    double dist;
    neighbors.clear();
    for(int i=0;i<metaPop.inds.size();i++)
    {
        dist=sqrt(pow((x-metaPop.inds[i]->x),2)+pow(y-metaPop.inds[i]->y,2));
        neighbors.push_back(dist);
    }
}
void ind::calcEta(mPop metaPop)
{
    eta=0;
    for(int i=0;i<metaPop.inds.size();i++)
    {
        eta+=exp(-1.0*pow(neighbors[i],2.0)/pow(metaPop.sigmaEta,2.0));
    }
}
void ind::calcZdelta(mPop metaPop)
{
    z=0;
    for(int l=0;l<metaPop.nLoci;l++)
    {
        z+=metaPop.aVec[l]*(genome[0][l]+genome[1][l])+metaPop.dVec[l]*(genome[0][l]*genome[1][l]);
    }
    deltai=metaPop.delta*(z-metaPop.zMin)/(metaPop.zMax-metaPop.zMin);
}
void ind::calcWrepro(mPop metaPop)
{
    calcEta(metaPop);
    calcZdelta(metaPop);
    if(metaPop.sel==1)//density-independent
    {
       Wrepro=1.0+metaPop.rho*(1.0-eta/metaPop.gamma);
    }
    else if(metaPop.sel==2) //density-dependent
    {
        ri=metaPop.delta+metaPop.rho-(1.0-metaPop.rho)*deltai;
        gammai=metaPop.gamma*(metaPop.delta+metaPop.rho-(1.0-metaPop.rho)*deltai)/(metaPop.rho*(1.0-deltai));
        Wrepro=1.0+ri*(1.0-eta/gammai);
    }
    else
    {
        cout<<"error in sel"<<endl;getchar();
    }
}
gamete::gamete(mPop metaPop,int parent)
{
    int wrkPar=0;
    x=metaPop.inds[parent]->x;y=metaPop.inds[parent]->y;
    //cout<<"gamete coordinates: "<<x<<", "<<y<<endl;getchar();
    genome=new int[metaPop.nLoci];
    for(int l=0;l<metaPop.nLoci;l++)
    {
        genome[l]=metaPop.inds[parent]->genome[wrkPar][l];
        if(rand()/((double)RAND_MAX)<metaPop.mu) genome[l]=(genome[l]+1)%2;
        if(rand()/((double)RAND_MAX)<metaPop.r) wrkPar=(wrkPar+1)%2;
    }
}
int gamete::findMate(mPop metaPop)
{
    int out;
    double dist=0,total,random;
    total=0;
    for(int g=0;g<metaPop.gametes.size();g++)
    {
        dist=sqrt(pow(x-metaPop.gametes[g]->x,2.0)+pow(y-metaPop.gametes[g]->y,2.0));
        pMateVec.push_back(exp(-pow(dist,2.0)/pow(metaPop.sigmaMate,2.0)));
        total+=exp(-pow(dist,2.0)/pow(metaPop.sigmaMate,2.0));
    }
    random=rand()/((double)RAND_MAX)*total;
    out=-1;total=0;
    while(total<random)
    {
        out++;
        total+=pMateVec[out];
    }
    return out;
}
int* shuffleList(int length)
{
	int* list=new int[length];
	int r,temp;
	for(int e=0;e<length;e++) list[e]=e;
	for(int e=0;e<length;e++)
	{
		r=e+rand()%(length-e);
		temp=list[r];
		list[r]=list[e];list[e]=temp;
	}
	return list;
}
double randomNormal(double mu, double sigma)
{
    double u,v,x;
    u=rand()/((double)RAND_MAX);
    v=rand()/((double)RAND_MAX);
    x=sqrt(-2.0*log(u))*cos(2.0*3.1415*v);
    x=x*sigma+mu;
    //y=sqrt(-2.0*log(u))*sin(2.0*3.1415*v); //Also distributed as a standard normal
    return x;
}
//Make census
//Make output folders and files