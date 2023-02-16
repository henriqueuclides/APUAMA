#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <iostream>
#include <algorithm>
#include <vector>

//constants
const int xmep=10;
const double psi=0.9705646;
const double h=6.62607004*pow(10,-34);
const double c=2.99793*pow(10,10);
const double mc=1.6605402*pow(10,-27);
const double hJ=4.3597483*pow(10,-18);
const double pi=4.0*atan(1.0);
const double press=1.013*pow(10,5);
const double r=1.9872;
const double kB=1.38064852*pow(10,-23);
const double xfb=mc*pow(10,-20);
const double A=6.022*pow(10,23);
const double hl=h/(2.0*pi);
const double fca1a2=(2.0*pi*hJ)/(627.5095*h*c);

typedef struct{
    char name[50];
    int nA,nFreq,type,symmetry,v,J;
    double *mass,*x,*y,*z,*freq;
    double Mass,energy,energyAux,freqN,zep;
    double xT7_1[10],xT7_2[10],xT9_1[12],xT9_2[12];
    double eao[60],sto[60],cpo[60];
    double di[5];
    double Eref,De,Req;
    double DEn[50][50],En[50][50],Y[11][11];
}TypeSpecie;

typedef struct{
    char react1[50],react2[50],ts[50],prod1[50],prod2[50];
    double Temp[60],temp[60],qt[50][45];
    double beta1,beta2,enthalpy,vg,vag,veff,vref1,vref2,Mass,yc,energyRC,energyPC,Psi;
    double K[30],K_Wigner[30],K_Eckart[30],K_TSC[30],arh[22];
    double vagg[101],vmep[101],s[101];
}TypeReaction;

double setAtom(StringVector atoms, int i);

void extractSpecies(int n, StringVector names, NumericVector natoms, NumericVector nfreqs, NumericVector symmetry, StringVector types, NumericVector energy, StringVector atoms, NumericVector xs, NumericVector ys, NumericVector zs, NumericVector freqs);

void extractReactions(int n, StringVector react1, StringVector react2, StringVector ts, StringVector prod1, StringVector prod2, NumericVector Psi);

void calculateRate();

void partition_functions_enthalpy();

void reaction_rate_Wigner();

void MEP_Eckart();

void tunneling_small_curvature();

void arrhenius();

void rovibrational_levels(int I);

void print_results();

void sys_linear10(double M[10][10], double b[], double x[]);

void sys_linear12(double M[12][12], double b[], double x[]);

