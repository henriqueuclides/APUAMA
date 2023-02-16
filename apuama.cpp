#include "apuama.h"

TypeSpecie *Specie;
TypeReaction *Reaction;
int Nspecies, Nreactions;

// [[Rcpp::export]]
void extractSpecies(int n, StringVector names, NumericVector natoms, NumericVector nfreqs, NumericVector symmetry, StringVector types, NumericVector energy, StringVector atoms, NumericVector xs, NumericVector ys, NumericVector zs, NumericVector freqs){
    Nspecies = n;
    Specie = (TypeSpecie*)malloc(Nspecies*sizeof(TypeSpecie));
    int i, j, cont1=0, cont2=0;

    for(i = 0; i < Nspecies; i++){
        strcpy(Specie[i].name,names[i]);
        Specie[i].nA = natoms[i];
        Specie[i].nFreq = nfreqs[i];
        Specie[i].symmetry = symmetry[i];
        if(strcmp(types[i],"atom") == 0)
            Specie[i].type = 1;
        if(strcmp(types[i],"linear") == 0)
            Specie[i].type = 2;
        if(strcmp(types[i],"nonlinear") == 0)
            Specie[i].type = 3;
        Specie[i].mass=(double*)malloc(Specie[i].nA*sizeof(double));
        Specie[i].x=(double*)malloc(Specie[i].nA*sizeof(double));
        Specie[i].y=(double*)malloc(Specie[i].nA*sizeof(double));
        Specie[i].z=(double*)malloc(Specie[i].nA*sizeof(double));
        Specie[i].Mass=0.0;
        Specie[i].energy = energy[i];
        Specie[i].energyAux=Specie[i].energy;

        for(j=0; j < Specie[i].nA; j++){
            Specie[i].mass[j] = setAtom(atoms,cont1);
            Specie[i].x[j] = xs[cont1];
            Specie[i].y[j] = ys[cont1];
            Specie[i].z[j] = zs[cont1];
            Specie[i].Mass+=Specie[i].mass[j];
            cont1++;
        }
        Specie[i].zep=0.0;
        Specie[i].freq=(double*)malloc(Specie[i].nFreq*sizeof(double));

        if (freqs[cont2] < 0 && Specie[i].nFreq > 0){
            Specie[i].freqN = freqs[cont2];
            Specie[i].nFreq = Specie[i].nFreq - 1;
            cont2++;
        }
        for(j=0; j < Specie[i].nFreq; j++){
            Specie[i].freq[j] = freqs[cont2];
            Specie[i].zep+=(1.0*pow(10,-3))*(0.5*A*pow(10,3)*h*c*Specie[i].freq[j])/4.184;
            cont2++;
        }
        Specie[i].v=Specie[i].J=0;
    }
}

// [[Rcpp::export]]
NumericMatrix setSpecies(int I){
  int i;
  i = Specie[I].nA;
  if(Specie[I].nA < Specie[I].nFreq)
    i = Specie[I].nFreq;
  NumericMatrix auxSpecie(6,i+1);

  auxSpecie(0,0) = Specie[I].energy;
  for(i = 0; i < Specie[I].nA; i++){
    auxSpecie(1,i) = Specie[I].mass[i];
    auxSpecie(2,i) = Specie[I].x[i];
    auxSpecie(3,i) = Specie[I].y[i];
    auxSpecie(4,i) = Specie[I].z[i];
  }
  if(Specie[I].freqN < 0){
    auxSpecie(5,0) = Specie[I].freqN;
    for(i = 0; i < Specie[I].nFreq; i++)
      auxSpecie(5,i+1) = Specie[I].freq[i];
  }
  else
    for(i = 0; i < Specie[I].nFreq; i++)
      auxSpecie(5,i) = Specie[I].freq[i];

  return auxSpecie;
}

// [[Rcpp::export]]
void extractReactions(int n, StringVector react1, StringVector react2, StringVector ts, StringVector prod1, StringVector prod2, NumericVector Psi){
    Nreactions = n;
    Reaction = (TypeReaction*)malloc(Nreactions*sizeof(TypeReaction));

    for(int i=0; i < Nreactions; i++){
        strcpy(Reaction[i].react1,react1[i]);
        strcpy(Reaction[i].react2,react2[i]);
        strcpy(Reaction[i].ts,ts[i]);
        strcpy(Reaction[i].prod1,prod1[i]);
        strcpy(Reaction[i].prod2,prod2[i]);
        if(Psi[i] == 0)
            Reaction[i].Psi = psi;
        else
            Reaction[i].Psi = Psi[i];
    }
}

// [[Rcpp::export]]
void extractRov(int I, NumericVector di, double Eref, double De, double Req){
  int i;
  for(i=0; i<5; i++)
    Specie[I].di[i] = di[i];
  Specie[I].Eref = Eref;
  Specie[I].De = De;
  Specie[I].Req = Req;
  rovibrational_levels(I);
}

// [[Rcpp::export]]
void calculateRate(){
    int i,j;

    for(i=0; i < Nreactions; i++){
        //define temperature rage
        for(j=0; j < 3; j++)
            if(j == 0)
                Reaction[i].temp[j]=50.0;
            else
                Reaction[i].temp[j]=100.0*j;
        Reaction[i].temp[3]=250.0;
        Reaction[i].temp[4]=298.15;
        for(j=5; j < 10; j++)
            Reaction[i].temp[j]=300.0+50.0*(j-5);
        for(j=10; j < 60; j++)
            Reaction[i].temp[j]=600.0+100.0*(j-10);
    }

    //Calculate partition functions and enthalpy
    partition_functions_enthalpy();

    //Calculate reaction rate and Wigner correction
    reaction_rate_Wigner();

    //Calculate MEP and Eckart correction
    MEP_Eckart();

    //Calculate tunneling small curvature
    tunneling_small_curvature();

    //set Ahrrenius
    arrhenius();

    //print some results
    //print_results();
}

// [[Rcpp::export]]
NumericMatrix setRate(int I){
    int j;
    NumericMatrix Rate(24,5);
    for(j=0; j <24; j++){
        Rate(j,0)=10000/Reaction[I].temp[j];
        Rate(j,1)=log10(Reaction[I].K[j]);
        Rate(j,2)=log10(Reaction[I].K_Wigner[j]);
        Rate(j,3)=log10(Reaction[I].K_Eckart[j]);
        Rate(j,4)=log10(Reaction[I].K_TSC[j]);
    }
    return Rate;
}

// [[Rcpp::export]]
NumericMatrix setMEP(int I){
  int j;
  NumericMatrix MEP(101,3);
  for(j=0; j <101; j++){
    MEP(j,0)=Reaction[I].s[j];
    MEP(j,1)=Reaction[I].vmep[j];
    MEP(j,2)=Reaction[I].vagg[j];
  }
  return MEP;
}

// [[Rcpp::export]]
NumericMatrix setThermProp(int I){
  int i,j;
  NumericMatrix ThermProp(60,10);
  for(i=0; i <60; i++){
    ThermProp(i,0)=Reaction[0].Temp[i];
    ThermProp(i,1)=Specie[I].eao[i]; //enthalpy
    ThermProp(i,2)=Specie[I].sto[i]; //entropy
    ThermProp(i,3)=Specie[I].cpo[i]; //heat capacity
  }
  j=0;
  for(i=0; i<5; i++){
    ThermProp(0,i+5)=Specie[I].xT7_1[j];
    j++;
  }
  for(i=0; i<5; i++){
    if(j<7)
      ThermProp(1,i+5)=Specie[I].xT7_1[j];
    else{
      j=0;
      ThermProp(1,i+5)=Specie[I].xT7_2[j];
    }
    j++;
  }
  for(i=0; i<5; i++){
    if(j<7)
      ThermProp(2,i+5)=Specie[I].xT7_2[j];
    else
      ThermProp(2,i+5)=0.0;
    j++;
  }
  j=0;
  for(i=0; i<5; i++){
    ThermProp(3,i+5)=Specie[I].xT9_1[j];
    j++;
  }
  for(i=0; i<5; i++){
    if(j < 7)
      ThermProp(4,i+5)=Specie[I].xT9_1[j];
    if(j == 7)
      ThermProp(4,i+5)=0.0;
    if(j > 7)
      ThermProp(4,i+5)=Specie[I].xT9_1[j-1];
    j++;
  }
  j=0;
  for(i=0; i<5; i++){
    ThermProp(5,i+5)=Specie[I].xT9_2[j];
    j++;
  }
  for(i=0; i<5; i++){
    if(j < 7)
      ThermProp(6,i+5)=Specie[I].xT9_2[j];
    if(j == 7)
      ThermProp(6,i+5)=0.0;
    if(j > 7)
      ThermProp(6,i+5)=Specie[I].xT9_2[j-1];
    j++;
  }

  return ThermProp;
}

// [[Rcpp::export]]
NumericMatrix setArrhenius(int I){
  int i;
  NumericMatrix Arrhenius(4,3);
  for(i=0; i<3; i++){
    Arrhenius(0,i) = Reaction[I].arh[i];
    Arrhenius(1,i) = Reaction[I].arh[i+4];
    Arrhenius(2,i) = Reaction[I].arh[i+8];
    Arrhenius(3,i) = Reaction[I].arh[i+12];
  }
  return Arrhenius;
}

// [[Rcpp::export]]
NumericMatrix setRvl(int I){
  int i,j;
  NumericMatrix Rvl(50,50);
  for(i=0; i<50; i++)
    for(j=0; j<50; j++)
      Rvl(i,j) = Specie[I].En[i][j];

  return Rvl;
}

// [[Rcpp::export]]
NumericMatrix setSpecConst(int I){
  int i,j;
  NumericMatrix Spec(11,11);
  for(i=0; i<11; i++)
    for(j=0; j<11; j++)
      Spec(i,j) = Specie[I].Y[i][j];

  return Spec;
}

// [[Rcpp::export]]
double addRvl(int I, int v, int J){
  Specie[I].energy = Specie[I].energyAux;
  Specie[I].energy += Specie[I].DEn[v][J];
  return Specie[I].energy;
}

void print_results(){
  int I,i,ts;
  FILE *f1,*f2;
  char directory[200],aux[200];
  for(I=0; I < Nreactions; I++){
    for(i=0; i < Nspecies; i++)
      if(strcmp(Reaction[I].ts,Specie[i].name) == 0)
        ts=i;
      strcpy(directory,"");
      strcpy(directory,Specie[ts].name);
      mkdir(directory,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      strcat(directory,"/");
      strcpy(aux,directory);
      f1=fopen(strcat(aux,"kt.dat"),"w");
      fprintf(f1,"k(T)\tkW(T)\tkE(T)\tkSCT(T)\n");
      for(i=0 ;i<24; i++){
        fprintf(f1,"%e\t",(10000.0/Reaction[I].temp[i]));
        fprintf(f1,"%e\t",log10(Reaction[I].K[i]));
        fprintf(f1,"%e\t",log10(Reaction[I].K_Wigner[i]));
        fprintf(f1,"%e\t",log10(Reaction[I].K_Eckart[i]));
        fprintf(f1,"%e\t",log10(Reaction[I].K_TSC[i]));
        fprintf(f1,"\n");
      }
      strcpy(aux,directory);
      f2=fopen(strcat(aux,"arrh.dat"),"w");
      fprintf(f2,"k\t\t A_tst\t\t n_tst\t\t Ea_tst\t\t diff\n");
      fprintf(f2,"k(T)\t\t");
      for(i=0; i <= 15 ;i++){
        fprintf(f2,"%e\t",Reaction[I].arh[i]);
        if(i==3)
          fprintf(f2,"\nkW(T)\t\t");
        if(i==7)
          fprintf(f2,"\nkE(T)\t\t");
        if(i==11)
          fprintf(f2,"\nkSCT(T)\t\t");
      }
      fclose(f1);
      fclose(f2);
  }
}

//set atom mass
double setAtom(StringVector atoms, int i){
    if(strcmp(atoms[i],"H") == 0)
        return 1.007831;
    if(strcmp(atoms[i],"He") == 0)
        return 4.002602;
    if(strcmp(atoms[i],"Li") == 0)
        return 6.941;
    if(strcmp(atoms[i],"Be") == 0)
        return 9.012182;
    if(strcmp(atoms[i],"B") == 0)
        return 10.811;
    if(strcmp(atoms[i],"C") == 0)
        return 12.0;
    if(strcmp(atoms[i],"N") == 0)
        return 14.0067;
    if(strcmp(atoms[i],"O") == 0)
        return 15.99491;
    if(strcmp(atoms[i],"F") == 0)
        return 18.9989518;
    if(strcmp(atoms[i],"Ne") == 0)
        return 20.1797;
    if(strcmp(atoms[i],"Na") == 0)
        return 22.98976928;
    if(strcmp(atoms[i],"Mg") == 0)
        return 24.3050;
    if(strcmp(atoms[i],"Al") == 0)
        return 26.9815386;
    if(strcmp(atoms[i],"Si") == 0)
        return 28.0855;
    if(strcmp(atoms[i],"P") == 0)
        return 30.973762;
    if(strcmp(atoms[i],"S") == 0)
        return 32.065;
    if(strcmp(atoms[i],"Cl") == 0)
        return 35.453;
    if(strcmp(atoms[i],"Ar") == 0)
        return 39.948;
    if(strcmp(atoms[i],"K") == 0)
        return 39.0983;
    if(strcmp(atoms[i],"Ca") == 0)
        return 40.078;
    if(strcmp(atoms[i],"Sc") == 0)
        return 44.955912;
    if(strcmp(atoms[i],"Ti") == 0)
        return 47.867;
    if(strcmp(atoms[i],"V") == 0)
        return 50.9415;
    if(strcmp(atoms[i],"Cr") == 0)
        return 51.9961;
    if(strcmp(atoms[i],"Mn") == 0)
        return 54.938;
    if(strcmp(atoms[i],"Fe") == 0)
        return 55.845;
    if(strcmp(atoms[i],"Co") == 0)
        return 58.9332;
    if(strcmp(atoms[i],"Ni") == 0)
        return 58.6934;
    if(strcmp(atoms[i],"Cu") == 0)
        return 63.546;
    if(strcmp(atoms[i],"Zn") == 0)
        return 65.39;
    if(strcmp(atoms[i],"Ga") == 0)
        return 69.723;
    if(strcmp(atoms[i],"Ge") == 0)
        return 72.64;
    if(strcmp(atoms[i],"As") == 0)
        return 74.9216;
    if(strcmp(atoms[i],"Se") == 0)
        return 78.96;
    if(strcmp(atoms[i],"Br") == 0)
        return 79.904;
    if(strcmp(atoms[i],"Kr") == 0)
        return 83.8;
    if(strcmp(atoms[i],"Rb") == 0)
        return 85.4678;
    if(strcmp(atoms[i],"Sr") == 0)
        return 87.62;
    if(strcmp(atoms[i],"Y") == 0)
        return 88.9059;
    if(strcmp(atoms[i],"Zr") == 0)
        return 91.224;
    if(strcmp(atoms[i],"Nb") == 0)
        return 92.9064;
    if(strcmp(atoms[i],"Mo") == 0)
        return 95.94;
    if(strcmp(atoms[i],"Tc") == 0)
        return 98.0;
    if(strcmp(atoms[i],"Ru") == 0)
        return 101.07;
    if(strcmp(atoms[i],"Rh") == 0)
        return 102.9055;
    if(strcmp(atoms[i],"Pd") == 0)
        return 106.42;
    if(strcmp(atoms[i],"Ag") == 0)
        return 107.8682;
    if(strcmp(atoms[i],"Cd") == 0)
        return 112.411;
    if(strcmp(atoms[i],"In") == 0)
        return 114.818;
    if(strcmp(atoms[i],"Sn") == 0)
        return 118.71;
    if(strcmp(atoms[i],"Sb") == 0)
        return 121.76;
    if(strcmp(atoms[i],"Te") == 0)
        return 127.6;
    if(strcmp(atoms[i],"I") == 0)
        return 126.9045;
    if(strcmp(atoms[i],"Xe") == 0)
        return 131.293;
    return 0;
}

//calculate partition functions, thermodynamic properties and NASA polynomials
void partition_functions_enthalpy(){
    int i,k,K,j,I;
    double sumx,sumy,sumz,sumx2,sumy2,sumz2,sumxy,sumxz,sumyz,sumxyq,sumxzq,sumyzq,ener298;
    double ss,sumfim,xa,xb,xc,xd,xe,xf,xcm,ycm,zcm,qt[60],pia,qr,qv,fac,qtt[Nspecies][60],tk[45];
    double x,ear,str,cpr,eav,stv,cpv,eat[60],stt[60],cpt[60];
    double tT[12][12][2],vT[12][2];
    double tT10_1[10][10],vT10_1[10],tT10_2[10][10],vT10_2[10];
    double tT_1[12][12],vT_1[12],tT_2[12][12],vT_2[12];

     //FILE *f=fopen("St_Cp_Ht.dat","w");
     //FILE *f1=fopen("Qtot.dat","w");

    for(I=0; I < Nspecies; I++){
        ener298=0;
        for(k=0; k<60; k++){
            qt[k]=pow((2.0*pi*mc*Specie[I].Mass*kB*Reaction[0].temp[k]/(h*h)),(1.5));
            stt[k]=37.0+3/2.0*r*log(Specie[I].Mass/40.0)+5/2.0*r*log(Reaction[0].temp[k]/298.0);
            cpt[k]=1.5*r;
            eat[k]=1.5*r*Reaction[0].temp[k];
        }
        ener298+=eat[4];
        //calculation of mass center
        sumx=sumy=sumz=0.0;
        for(j=0; j<Specie[I].nA; j++){
            sumx=sumx+Specie[I].mass[j]*Specie[I].x[j];
            sumy=sumy+Specie[I].mass[j]*Specie[I].y[j];
            sumz=sumz+Specie[I].mass[j]*Specie[I].z[j];
        }
        xcm=sumx/Specie[I].Mass;
        ycm=sumy/Specie[I].Mass;
        zcm=sumz/Specie[I].Mass;

        //calculation of inertial moment
        sumx=sumy=sumz=0.0;
        sumx2=sumy2=sumz2=sumxy=sumxz=sumyz=sumxyq=sumxzq=sumyzq=0.0;
        for(j=0; j<Specie[I].nA; j++){
            sumx=sumx+Specie[I].mass[j]*(Specie[I].x[j]);
            sumy=sumy+Specie[I].mass[j]*(Specie[I].y[j]);
            sumz=sumz+Specie[I].mass[j]*(Specie[I].z[j]);
            sumx2=sumx2+Specie[I].mass[j]*pow((Specie[I].x[j]-xcm),2);
            sumy2=sumy2+Specie[I].mass[j]*pow((Specie[I].y[j]-ycm),2);
            sumz2=sumz2+Specie[I].mass[j]*pow((Specie[I].z[j]-zcm),2);
            if(sumx2 != 0.0)
                sumfim=sumx2;
            if(sumy2 != 0.0)
                sumfim=sumy2;
            if(sumz2 != 0.0)
                sumfim=sumz2;
            sumxy=sumxy+Specie[I].mass[j]*(Specie[I].x[j])*(Specie[I].y[j]);
            sumxz=sumxz+Specie[I].mass[j]*(Specie[I].x[j])*(Specie[I].z[j]);
            sumyz=sumyz+Specie[I].mass[j]*(Specie[I].y[j])*(Specie[I].z[j]);
            sumxyq=sumxyq+Specie[I].mass[j]*(pow((Specie[I].x[j]),2)+pow((Specie[I].y[j]),2));
            sumxzq=sumxzq+Specie[I].mass[j]*(pow((Specie[I].x[j]),2)+pow((Specie[I].z[j]),2));
            sumyzq=sumyzq+Specie[I].mass[j]*(pow((Specie[I].y[j]),2)+pow((Specie[I].z[j]),2));
        }
        xa=sumyzq-(pow(sumy,2)+pow(sumz,2))/Specie[I].Mass;
        xb=sumxzq-(pow(sumx,2)+pow(sumz,2))/Specie[I].Mass;
        xc=sumxyq-(pow(sumx,2)+pow(sumy,2))/Specie[I].Mass;
        xd=sumxy-sumx*sumy/Specie[I].Mass;
        xe=sumxz-sumx*sumz/Specie[I].Mass;
        xf=sumyz-sumy*sumz/Specie[I].Mass;
        if(Specie[I].type == 1)
            pia=0.0;
        if(Specie[I].type == 2)
            pia=sumfim;
        if(Specie[I].type == 3){
            pia=(xa*xb*xc-xa*xf*xf-xc*xd*xd-2.0*xd*xe*xf-xb*xe*xe);
            pia=sumyzq*sumxzq*sumxyq-sumxy*sumyz*sumxz-sumxz*sumxy*sumyz-pow(sumxz,2)*sumxzq-pow(sumyz,2)*sumyzq-pow(sumxy,2)*sumxyq;
        }
        //calculation of ener298
        if(Specie[I].type == 1)
            ear=0.0;
        else{
            if(Specie[I].type == 2)
                ear=r*Reaction[0].temp[4];
            else
                ear=1.5*r*Reaction[0].temp[4];
        }
        eav=0.0;
        fac=h*c/(kB*Reaction[0].temp[4]);
        for(j=0; j<Specie[I].nFreq; j++){
            x=fac*Specie[I].freq[j];
            eav=eav+r*Reaction[0].temp[4]*x/(exp(x)-1.0);
        }
        ener298+=ear+eav;

        //calculation of the rotation partition function
        //fprintf(f,"\n%s entropia, Calor e., entalpia\n",Specie[I].name);
        //fprintf(f1,"\n%s\nIxx: %e, Iyy: %e, Izz: %e, Ixy: %e, Ixz: %e, Iyz: %e\n",Specie[I].name,sumyzq,sumxzq,sumxyq,sumxy,sumxz,sumyz);
        for(k=0; k<60; k++){
            //atomic specie
            if(Specie[I].type == 1){
                qr=1.0;
                ear=0.0;
                str=0.0;
                cpr=0.0;
            }
            else{
                //linear specie
                if(Specie[I].type == 2){
                    qr=(1.66*pow(10,-47))*((8.0*pi*pi*pia*kB*Reaction[0].temp[k])/(Specie[I].symmetry*h*h));
                    ear=r*Reaction[0].temp[k];
                    str=6.9+r*log(pia/Specie[I].symmetry)+r*log(Reaction[0].temp[k]/298.0);
                    cpr=r;
                }
                //non linear specie
                else{
                    qr=(sqrt(pi)/Specie[I].symmetry)*pow(((8.0*pi*pi*1.66*pow(10,-47)*pia*kB*Reaction[0].temp[k])/(h*h)),1.5);
                    str=11.5+0.5*r*log(pow(pia,3.0)/pow(Specie[I].symmetry,2.0))+1.5*r*log(Reaction[0].temp[k]/298.0);
                    ear=1.5*r*Reaction[0].temp[k];
                    cpr=2.5*r;
                }
            }
            //calculation of the partition function of vibration
            qv=1.0;
            stv=0.0;
            cpv=0.0;
            eav=0.0;
            fac=h*c/(kB*Reaction[0].temp[k]);
            for(j=0; j<Specie[I].nFreq; j++){
                x=fac*Specie[I].freq[j];
                qv=qv*pow((1.0-exp(-x)),(-1.0));
                ss=r*(x*exp(-x)/(1.0-exp(-x))-log(1.0-exp(-x)));
                stv=stv+ss;
                eav=eav+r*Reaction[0].temp[k]*x/(exp(x)-1.0);
                cpv=cpv+r*x*x*exp(x)/((exp(x)-1.0)*(exp(x)-1.0));
            }

            // final calculation of the partition function
            qtt[I][k]=qt[k]*qr*qv;
            //thermodynamic properties
            //entropy
            Specie[I].sto[k]=stt[k]+str+stv-2*r;
            //heat capacity
            Specie[I].cpo[k]=cpt[k]+cpr+cpv;
            //enthalpy
            Specie[I].eao[k]=(eat[k]+ear+eav-ener298)*pow(10,-3.0);
            Reaction[0].Temp[k] = Reaction[0].temp[k];
            //fprintf(f,"%.4f %.4e %.4e %.4e\n",Reaction[0].temp[k],Specie[I].sto[k],Specie[I].cpo[k],Specie[I].eao[k]);
            //fprintf(f1,"%.4f %e %e %e %e\n",Reaction[0].temp[k],qt[k],qr,qv,qtt[I][k]/(pow(10,27)));
        }

        //NASA coefficient for polynomial fitting of the thermodynamic properties
        for(i=0; i<12; i++)
            for(j=0; j<12; j++)
                tT[i][j][0]=tT[i][j][1]=0.0;

        for(i=0; i<2; i++){
            if(i==0){
                j=0;
                K=14;
            }
            else{
                j=14;
                K=60;
            }
            //tT[0][0]
            for(k=j; k<K; k++)
                tT[0][0][i]=pow(Reaction[0].temp[k],-4.0);
            tT[0][0][i]=tT[0][0][i]*9.0/4;
            //tT[0][1]
            for(k=j; k<K; k++)
                tT[0][1][i]+=(3.0/2-log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],-3.0);
            //tT[0][2]
            for(k=j; k<K; k++)
                tT[0][2][i]+=log(Reaction[0].temp[k])*pow(Reaction[0].temp[k],-2.0);
            tT[0][2][i]=tT[0][2][i]*-0.5;
            //tT[0][3]
            tT[0][3][i]=0.0;
            //tT[0][4]
            tT[0][4][i]=5.0/12*(K-j);
            //tT[0][5]
            for(k=j; k<K; k++)
                tT[0][5][i]+=Reaction[0].temp[k];
            tT[0][5][i]=tT[0][5][i]*7.0/12;
            //tT[0][6]
            for(k=j; k<K; k++)
                tT[0][6][i]+=pow(Reaction[0].temp[k],2.0);
            tT[0][6][i]=tT[0][6][i]*27.0/40;
            //tT[0][7]
            for(k=j; k<K; k++)
                tT[0][7][i]+=pow(Reaction[0].temp[k],-3.0);
            tT[0][7][i]=tT[0][7][i]*-1.0;
            //tT[0][8]
            for(k=j; k<K; k++)
                tT[0][8][i]+=pow(Reaction[0].temp[k],-2.0);
            tT[0][8][i]=tT[0][8][i]*-0.5;
            //tT[0][9]
            tT[0][9][i]=pow(Reaction[0].temp[14],-2.0);
            //tT[0][10]
            tT[0][10][i]=-1.0*tT[0][9][i];
            //tT[0][11]
            tT[0][11][i]=-0.5*tT[0][9][i];

            //tT[1][1]
            for(k=j; k<K; k++)
                tT[1][1][i]+=(2+pow(log(Reaction[0].temp[k]),2.0))*pow(Reaction[0].temp[k],-2.0);
            //tT[1][2]
            for(k=j; k<K; k++)
                tT[1][2][i]+=pow(Reaction[0].temp[k],-1.0);
            //tT[1][3]
            for(k=j; k<K; k++)
                tT[1][3][i]+=log(Reaction[0].temp[k]);
            tT[1][3][i]=tT[1][3][i]*0.5;
            //tT[1][4]
            for(k=j; k<K; k++)
                tT[1][4][i]+=(1.0/2+1.0/3*log(Reaction[0].temp[k]))*Reaction[0].temp[k];
            //tT[1][5]
            for(k=j; k<K; k++)
                tT[1][5][i]+=(2.0/3+1.0/4*log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],2.0);
            //tT[1][6]
            for(k=j; k<K; k++)
                tT[1][6][i]+=(3.0/4+1.0/5*log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],3.0);
            //tT[1][7]
            for(k=j; k<K; k++)
                tT[1][7][i]+=log(Reaction[0].temp[k])*pow(Reaction[0].temp[k],-2.0);
            //tT[1][8]
            tT[1][8][i]=-1.0*tT[1][2][i];
            //tT[1][9]
            tT[1][9][i]=pow(Reaction[0].temp[14],-1.0);
            //tT[1][10]
            tT[1][10][i]=tT[1][9][i]*log(Reaction[0].temp[14]);
            //tT[1][11]
            tT[1][11][i]=-1.0*tT[1][9][i];

            //tT[2][2]
            for(k=j; k<K; k++)
                tT[2][2][i]+=2+pow(log(Reaction[0].temp[k]),2.0);
            //tT[2][3]
            for(k=j; k<K; k++)
                tT[2][3][i]+=(3.0/2+log(Reaction[0].temp[k]))*Reaction[0].temp[k];
            //tT[2][4]
            for(k=j; k<K; k++)
                tT[2][4][i]+=(4.0/3+1.0/2*log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],2.0);
            //tT[2][5]
            for(k=j; k<K; k++)
                tT[2][5][i]+=(5.0/4+1.0/3*log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],3.0);
            //tT[2][6]
            for(k=j; k<K; k++)
                tT[2][6][i]+=(6.0/5+1.0/4*log(Reaction[0].temp[k]))*pow(Reaction[0].temp[k],4.0);
            //tT[2][7]
            for(k=j; k<K; k++)
                tT[2][7][i]+=1.0/Reaction[0].temp[k];
            //tT[2][8]
            tT[2][8][i]=2.0*tT[1][3][i];
            //tT[2][9]
            tT[2][9][i]=1.0;
            //tT[2][10]
            tT[2][10][i]=1.0;
            //tT[2][11]
            tT[2][11][i]=log(Reaction[0].temp[14]);

            //tT[3][3]
            tT[3][3][i]=40.0/27*9.0/4*tT[0][6][i];
            //tT[3][4]
            for(k=j; k<K; k++)
                tT[3][4][i]+=pow(Reaction[0].temp[k],3.0);
            tT[3][4][i]=tT[3][4][i]*5.0/3;
            //tT[3][5]
            for(k=j; k<K; k++)
                tT[3][5][i]+=pow(Reaction[0].temp[k],4.0);
            tT[3][5][i]=tT[3][5][i]*35.0/24;
            //tT[3][6]
            for(k=j; k<K; k++)
                tT[3][6][i]+=pow(Reaction[0].temp[k],5.0);
            tT[3][6][i]=tT[3][6][i]*27.0/20;
            //tT[3][7]
            tT[3][7][i]=(K-j)/2.0;
            //tT[3][8]
            tT[3][8][i]=tT[0][5][i]*12/7.0;
            //tT[3][9]
            tT[3][9][i]=Reaction[0].temp[14];
            //tT[3][10]
            tT[3][10][i]=Reaction[0].temp[14]/2.0;
            //tT[3][11]
            tT[3][11][i]=Reaction[0].temp[14];

            //tT[4][4]
            tT[4][4][i]=tT[3][5][i]*24/35.0*49/36.0;
            //tT[4][5]
            tT[4][5][i]=tT[3][6][i]*20/27.0*5/4;
            //tT[4][6]
            for(k=j; k<K; k++)
                tT[4][6][i]+=pow(Reaction[0].temp[k],6.0);
            tT[4][6][i]=tT[4][6][i]*143/120.0;
            //tT[4][7]
            tT[4][7][i]=tT[3][8][i]/3.0;
            //tT[4][8]
            tT[4][8][i]=tT[3][3][i]*2/9.0;
            //tT[4][9]
            tT[4][9][i]=pow(Reaction[0].temp[14],2.0);
            //tT[4][10]
            tT[4][10][i]=tT[4][9][i]/3.0;
            //tT[4][11]
            tT[4][11][i]=tT[4][9][i]/2.0;

            //tT[5][5]
            tT[5][5][i]=tT[4][6][i]*120/143.0*169/144.0;
            //tT[5][6]
            for(k=j; k<K; k++)
                tT[5][6][i]+=pow(Reaction[0].temp[k],7.0);
            tT[5][6][i]=tT[5][6][i]*17/14.0;
            //tT[5][7]
            tT[5][7][i]=tT[4][8][i]/2.0;
            //tT[5][8]
            tT[5][8][i]=tT[4][3][i]/5.0;
            //tT[5][9]
            tT[5][9][i]=pow(Reaction[0].temp[14],3.0);
            //tT[5][10]
            tT[5][10][i]=tT[5][9][i]/4.0;
            //tT[5][11]
            tT[5][11][i]=tT[5][9][i]/3.0;

            //tT[6][6]
            for(k=j; k<K; k++)
                tT[6][6][i]+=pow(Reaction[0].temp[k],8.0);
            tT[6][6][i]=tT[6][6][i]*441/400.0;
            //tT[6][7]
            tT[6][7][i]=tT[5][8][i]*3/5.0;
            //tT[5][8]
            tT[6][8][i]=tT[3][5][i]*24/140.0;
            //tT[5][9]
            tT[6][9][i]=pow(Reaction[0].temp[14],4.0);
            //tT[5][10]
            tT[6][10][i]=tT[6][9][i]/5.0;
            //tT[5][11]
            tT[6][11][i]=tT[6][9][i]/4.0;

            //tT[7][7]
            for(k=j; k<K; k++)
                tT[7][7][i]+=1.0/pow(Reaction[0].temp[k],2.0);
            //tT[7][10]
            tT[7][10][i]=1.0/Reaction[0].temp[14];

            //tT[8][8]
            tT[8][8][i]=K-j;
            tT[8][11][i]=1.0;
        }

        for(i=0; i<12; i++)
            for(j=0; j<12; j++)
                if(i>j){
                    tT[i][j][0]=tT[j][i][0];
                    tT[i][j][1]=tT[j][i][1];
                }
                //table OK

                for(i=0; i<12; i++)
                    vT[i][0]=vT[i][1]=0.0;
                for(i=0; i<2; i++){
                    if(i==0){
                        j=0;
                        K=14;
                    }
                    else{
                        j=14;
                        K=60;
                    }
                    for(k=j; k<K; k++)
                        vT[0][i]+=pow(Reaction[0].temp[k],-2.0)*(Specie[I].cpo[k]/r-(Specie[I].eao[k]-Specie[I].eao[4])/(r*Reaction[0].temp[k])-Specie[I].sto[k]/(2*r));
                    for(k=j; k<K; k++)
                        vT[1][i]+=pow(Reaction[0].temp[k],-1.0)*(Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(r*Reaction[0].temp[k])*log(Reaction[0].temp[k])-Specie[I].sto[k]/r);
                    for(k=j; k<K; k++)
                        vT[2][i]+=Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(r*Reaction[0].temp[k])+log(Reaction[0].temp[k])*Specie[I].sto[k]/r;
                    for(k=j; k<K; k++)
                        vT[3][i]+=Reaction[0].temp[k]*(Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(2*r*Reaction[0].temp[k])+Specie[I].sto[k]/r);
                    for(k=j; k<K; k++)
                        vT[4][i]+=pow(Reaction[0].temp[k],2.0)*(Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(3*r*Reaction[0].temp[k])+Specie[I].sto[k]/(2*r));
                    for(k=j; k<K; k++)
                        vT[5][i]+=pow(Reaction[0].temp[k],3.0)*(Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(4*r*Reaction[0].temp[k])+Specie[I].sto[k]/(3*r));
                    for(k=j; k<K; k++)
                        vT[6][i]+=pow(Reaction[0].temp[k],4.0)*(Specie[I].cpo[k]/r+(Specie[I].eao[k]-Specie[I].eao[4])/(5*r*Reaction[0].temp[k])+Specie[I].sto[k]/(4*r));
                    for(k=j; k<K; k++)
                        vT[7][i]+=((Specie[I].eao[k]-Specie[I].eao[4])/(r*Reaction[0].temp[k]))/Reaction[0].temp[k];
                    for(k=j; k<K; k++)
                        vT[8][i]+=Specie[I].sto[k]/r;
                    vT[9][i]=Specie[I].cpo[14]/r;
                    vT[10][i]=(Specie[I].eao[k]-Specie[I].eao[4])/(r*Reaction[0].temp[14]);
                    vT[11][i]=Specie[I].sto[14]/r;
                }
                for(i=0; i<10; i++){
                    vT10_1[i]=vT[i+2][0];
                    vT10_2[i]=vT[i+2][1];
                    for(j=0; j<10; j++){
                        tT10_1[i][j]=tT[i+2][j+2][0];
                        tT10_2[i][j]=tT[i+2][j+2][1];
                    }
                }
                sys_linear10(tT10_1,vT10_1,Specie[I].xT7_1);
                sys_linear10(tT10_2,vT10_2,Specie[I].xT7_2);
                for(i=0; i<12; i++){
                    vT_1[i]=vT[i][0];
                    vT_2[i]=vT[i][1];
                    for(j=0; j<12; j++){
                        tT_1[i][j]=tT[i][j][0];
                        tT_2[i][j]=tT[i][j][1];
                    }
                }
                sys_linear12(tT_1,vT_1,Specie[I].xT9_1);
                sys_linear12(tT_2,vT_2,Specie[I].xT9_2);
    }

    int r1,r2,ts;

    for(I=0; I<Nreactions ;I++){
        for(i=0; i < Nspecies; i++){
            if(strcmp(Reaction[I].react1,Specie[i].name) == 0)
                r1=i;
            if(strcmp(Reaction[I].react2,Specie[i].name) == 0)
                r2=i;
            if(strcmp(Reaction[I].ts,Specie[i].name) == 0)
                ts=i;
        }
        k=0;
        for(j=2; j<6; j++){
            tk[k]=Reaction[I].temp[j];
            Reaction[I].qt[r1][k]=qtt[r1][j];
            if(strcmp(Reaction[I].react2,"ZERO") != 0 && strcmp(Reaction[I].react2,"zero") != 0)
                Reaction[I].qt[r2][k]=qtt[r2][j];
            Reaction[I].qt[ts][k]=qtt[ts][j];
            k++;
        }
        for(j=7; j<10; j=j+2){
            tk[k]=Reaction[I].temp[j];
            Reaction[I].qt[r1][k]=qtt[r1][j];
            if(strcmp(Reaction[I].react2,"ZERO") != 0 && strcmp(Reaction[I].react2,"zero") != 0)
                Reaction[I].qt[r2][k]=qtt[r2][j];
            Reaction[I].qt[ts][k]=qtt[ts][j];
            k++;
        }
        for(j=10; j<45; j=j+2){
            tk[k]=Reaction[I].temp[j];
            Reaction[I].qt[r1][k]=qtt[r1][j];
            if(strcmp(Reaction[I].react2,"ZERO") != 0 && strcmp(Reaction[I].react2,"zero") != 0)
                Reaction[I].qt[r2][k]=qtt[r2][j];
            Reaction[I].qt[ts][k]=qtt[ts][j];
            k++;
        }
        for(j=0; j<24; j++)
            Reaction[I].temp[j]=tk[j];
    }
    //fclose(f);
    //fclose(f1);
}

//triangle bottom
void triang_btt10(double M[10][10], double b[], double Vec[]){
    int i,j;
    double sum;
    Vec[0]=b[0]/M[0][0];
    for(i=1; i<10; i++){
        sum=0.0;
        for(j=0; j<=(i-1); j++)
            sum+=M[i][j]*Vec[j];
        Vec[i]=(1.0/M[i][i])*(b[i]-sum);
    }
}

//triangle higher
void triang_high10(double M[10][10], double b[], double Vec[]){
    int i,j;
    double sum;
    Vec[9]=b[9]/M[9][9];
    for(i=8; i>=0; i=i-1){
        sum=0.0;
        for(j=i+1; j<10; j++)
            sum+=M[i][j]*Vec[j];
        Vec[i]=(1.0/M[i][i])*(b[i]-sum);
    }
}

//LU for 10
void LU10(double A[10][10], double L[10][10], double U[10][10]){
    int i,j,k;
    double sumU,sumL;
    for(i=0; i<10; i++)
        for(j=0; j<10; j++){
            L[i][j]=U[i][j]=0.0;
            if(i==j)
                L[i][i]=1.0;
        }
        for(i=0; i<10; i++)
            U[0][i]=A[0][i];
    for(i=1; i<10; i++)
        L[i][0]=A[i][0]/U[0][0];
    for(j=1; j<10; j++){
        for(i=1; i<10; i++){
            sumU=sumL=0;
            for(k=0; k<=i-1; k++)
                sumU+=L[i][k]*U[k][j];
            for(k=0; k<=j-1; k++)
                sumL+=L[i][k]*U[k][j];
            if(i<=j)
                U[i][j]=A[i][j]-sumU;
            else
                L[i][j]=(1/U[j][j])*(A[i][j]-sumL);
        }
    }
}

//linear system for 10
void sys_linear10(double M[10][10], double b[], double x[]){
    double L[10][10],U[10][10],v[10];
    LU10(M,L,U);
    triang_btt10(L,b,v);
    triang_high10(U,v,x);
}

//triangle bottom
void triang_btt12(double M[12][12], double b[], double Vec[]){
    int i,j;
    double sum;
    Vec[0]=b[0]/M[0][0];
    for(i=1; i<12; i++){
        sum=0.0;
        for(j=0; j<=(i-1); j++)
            sum+=M[i][j]*Vec[j];
        Vec[i]=(1.0/M[i][i])*(b[i]-sum);
    }
}

//triangle higher
void triang_high12(double M[12][12], double b[], double Vec[]){
    int i,j;
    double sum;
    Vec[11]=b[11]/M[11][11];
    for(i=10; i>=0; i=i-1){
        sum=0.0;
        for(j=i+1; j<12; j++)
            sum+=M[i][j]*Vec[j];
        Vec[i]=(1.0/M[i][i])*(b[i]-sum);
    }
}

//LU for 12
void LU12(double A[12][12], double L[12][12], double U[12][12]){
    int i,j,k;
    double sumU,sumL;
    for(i=0; i<12; i++)
        for(j=0; j<12; j++){
            L[i][j]=U[i][j]=0.0;
            if(i==j)
                L[i][i]=1.0;
        }
        for(i=0; i<12; i++)
            U[0][i]=A[0][i];
    for(i=1; i<12; i++)
        L[i][0]=A[i][0]/U[0][0];
    for(j=1; j<12; j++){
        for(i=1; i<12; i++){
            sumU=sumL=0;
            for(k=0; k<=i-1; k++)
                sumU+=L[i][k]*U[k][j];
            for(k=0; k<=j-1; k++)
                sumL+=L[i][k]*U[k][j];
            if(i<=j)
                U[i][j]=A[i][j]-sumU;
            else
                L[i][j]=(1/U[j][j])*(A[i][j]-sumL);
        }
    }
}

//linear system for 12
void sys_linear12(double M[12][12], double b[], double x[]){
    double L[12][12],U[12][12],v[12];
    LU12(M,L,U);
    triang_btt12(L,b,v);
    triang_high12(U,v,x);
}

//calculate reaction rate and Wigner tunneling correction
void reaction_rate_Wigner(){
    double ea,af,qttg,dzpe,ener,xb,KW,MassR2,energyR2,zepR2,MassP2,energyP2,zepP2;
    int i,I,k,r1,r2,ts,p1,p2;

    for(I=0; I < Nreactions; I++){
        for(i=0; i < Nspecies; i++){
            if(strcmp(Reaction[I].react1,Specie[i].name) == 0)
                r1=i;
            if(strcmp(Reaction[I].react2,Specie[i].name) == 0)
                r2=i;
            if(strcmp(Reaction[I].ts,Specie[i].name) == 0)
                ts=i;
            if(strcmp(Reaction[I].prod1,Specie[i].name) == 0)
                p1=i;
            if(strcmp(Reaction[I].prod2,Specie[i].name) == 0)
                p2=i;
        }
        if(strcmp(Reaction[I].react2,"ZERO") == 0 || strcmp(Reaction[I].react2,"zero") == 0){
            MassR2 = 1;
            energyR2 = 0;
            zepR2 = 0;
        }
        else{
            MassR2 = Specie[r2].Mass;
            energyR2 = Specie[r2].energy;
            zepR2 = Specie[r2].zep;
        }
        if(strcmp(Reaction[I].prod2,"ZERO") == 0 || strcmp(Reaction[I].prod2,"zero") == 0){
            MassP2 = 1;
            energyP2 = 0;
            zepP2 = 0;
        }
        else{
            MassP2 = Specie[p2].Mass;
            energyP2 = Specie[p2].energy;
            zepP2 = Specie[p2].zep;
        }
        for(k=0; k<24; k++){
            if(strcmp(Reaction[I].react2,"ZERO") == 0 || strcmp(Reaction[I].react2,"zero") == 0){
                qttg=Reaction[I].qt[ts][k]/Reaction[I].qt[r1][k];
                dzpe=Reaction[I].Psi*(Specie[ts].zep-Specie[r1].zep);
                ener=(Specie[ts].energy-Specie[r1].energy)*627.5095*pow(10,3);
                ea=ener+dzpe;
                af=(kB*Reaction[I].temp[k]/h)*qttg;
                Reaction[I].beta1=0.0;
                Reaction[I].beta2=0.0;
                //printf("%e %e %e\n",ea,af,qttg);
            }
            else{
                qttg=Reaction[I].qt[ts][k]/(Reaction[I].qt[r1][k]*Reaction[I].qt[r2][k]);
                dzpe=Reaction[I].Psi*(Specie[ts].zep-(Specie[r1].zep+zepR2));
                ener=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095*pow(10,3);
                ea=ener+dzpe;
                af=(pow(10,6)*A*kB*Reaction[I].temp[k]/h)*qttg;
                xb=(Specie[p1].Mass*MassR2)/(Specie[r1].Mass*MassP2);
                Reaction[I].beta1=(180.0/pi)*acos(sqrt(xb));
                xb=(Specie[r1].Mass-Specie[p1].Mass)*Specie[ts].Mass/(MassR2*Specie[p1].Mass);
                Reaction[I].beta2=(180.0/pi)*atan(sqrt(xb));
            }
            Reaction[I].K[k]=af*exp(-ea/(r*Reaction[I].temp[k]));
            //Wigner's correction
            KW=1.0+(1.0/30.0)*pow((h*fabs(Specie[ts].freqN)*c/(kB*Reaction[I].temp[k])),2.0);
            Reaction[I].K_Wigner[k]=Reaction[I].K[k]*KW;
        }
        Reaction[I].enthalpy=1.0*pow(10,-3)*(((Specie[p1].energy+energyP2)-(Specie[r1].energy+energyR2))*627.5095*pow(10,3)+Reaction[I].Psi*(Specie[p1].zep+zepP2)-(Specie[r1].zep+zepR2));
        Reaction[I].vg=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095;
        Reaction[I].vag=Reaction[I].vg+Reaction[I].Psi*(Specie[ts].zep-(Specie[r1].zep+zepR2))*1.0*pow(10,-3);
    }
}

//functions for Eckart tunneling correction
double xkappa(double tr1,double tr2,double delt){
    return(1.0-(cosh(2.0*pi*tr1)+cosh(2.0*pi*delt))/(cosh(2.0*pi*tr2)+cosh(2.0*pi*delt)));
}

//functions for Eckart tunneling correction
double alf(double x, double a1t){
    return(a1t*sqrt(x));
}

//functions for Eckart tunneling correction
double bt(double x,double aa1,double aa2,double b1t){
    return(b1t*sqrt((1.0+x)*aa1-aa2));
}

//functions for Eckart tunneling correction
double dtt(double aa1,double aa2){
    return((1.0/pi)*sqrt(fabs(aa1*aa2-pi*pi/4.0)));
}

//calculation of MEP and Eckart tunneling correction
void MEP_Eckart(){
    int i,I,k,r1,r2,ts,p1,p2;
    double xmp2r,xmp2p,aau,aac,vdu,vdc,bbu,bbc,fca2,alfa,xmax;
    double sou,soc,cc,a1,a2,vmax,att1,att2,vint,a1t,xpp;
    double b1t,ter1,ter2,ter3,ttest1,ttest2,f0,f2n,sm,sp;
    double xint,xhp,xmp,yu,MassR2,MassP2,energyR2,energyP2,zepR2,zepP2;

    for(I=0; I < Nreactions; I++){
        for(i=0; i < Nspecies; i++){
            if(strcmp(Reaction[I].react1,Specie[i].name) == 0)
                r1=i;
            if(strcmp(Reaction[I].react2,Specie[i].name) == 0)
                r2=i;
            if(strcmp(Reaction[I].ts,Specie[i].name) == 0)
                ts=i;
            if(strcmp(Reaction[I].prod1,Specie[i].name) == 0)
                p1=i;
            if(strcmp(Reaction[I].prod2,Specie[i].name) == 0)
                p2=i;
        }
        if(strcmp(Reaction[I].react2,"ZERO") == 0 || strcmp(Reaction[I].react2,"zero") == 0){
            MassR2 = 1;
            energyR2 = 0;
            zepR2 = 0;
        }
        else{
            MassR2 = Specie[r2].Mass;
            energyR2 = Specie[r2].energy;
            zepR2 = Specie[r2].zep;
        }
        if(strcmp(Reaction[I].prod2,"ZERO") == 0 || strcmp(Reaction[I].prod2,"zero") == 0){
            MassP2 = 1;
            energyP2 = 0;
            zepP2 = 0;
        }
        else{
            MassP2 = Specie[p2].Mass;
            energyP2 = Specie[p2].energy;
            zepP2 = Specie[p2].zep;
        }
        xmp2r=Specie[r1].energy+energyR2;
        xmp2p=Specie[p1].energy+energyP2;

        if(Reaction[I].vg > 0){
            aau=((Specie[p1].energy+energyP2)-(Specie[r1].energy+energyR2))*627.5095;
            aac=((Specie[p1].energy+energyP2)-(Specie[r1].energy+energyR2))*627.5095+Reaction[I].Psi*((Specie[p1].zep+zepP2)-(Specie[r1].zep+zepR2))*1.0*pow(10,-3);
            vdu=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095;
            vdc=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095+Reaction[I].Psi*((Specie[ts].zep-(Specie[r1].zep+zepR2)))*1.0*pow(10,-3);
            Reaction[I].veff=vdc;
        }
        else{
            if(xmp2r<=Reaction[I].energyRC){
                aau=((Specie[p1].energy+energyP2)-(Specie[r1].energy+energyR2))*627.5095;
                aac=((Specie[p1].energy+energyP2)-(Specie[r1].energy+energyR2))*627.5095+Reaction[I].Psi*((Specie[p1].zep+zepP2)-(Specie[r1].zep+zepR2))*1.0*pow(10,-3);
                vdu=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095;
                vdc=(Specie[ts].energy-(Specie[r1].energy+energyR2))*627.5095+Reaction[I].Psi*((Specie[ts].zep-(Specie[r1].zep+zepR2)))*1.0*pow(10,-3);
                Reaction[I].veff=vdc;
            }
            else{
                if(xmp2p<=Specie[ts].energy){
                    aau=((Specie[p1].energy+energyP2)-Reaction[I].energyRC)*627.5095;
                    aac=((Specie[p1].energy+energyP2)-Reaction[I].energyRC)*627.5095+Reaction[I].Psi*((Specie[p1].zep+zepP2)-(Specie[r1].zep+zepR2))*1.0*pow(10,-3);
                    vdu=(Specie[ts].energy-Reaction[I].energyRC)*627.5095;
                    vdc=(Specie[ts].energy-Reaction[I].energyRC)*627.5095+Reaction[I].Psi*((Specie[ts].zep-(Specie[r1].zep+zepR2)))*1.0*pow(10,-3);
                    Reaction[I].veff=vdc;
                }
                else{
                    aau=(Reaction[I].energyPC-Reaction[I].energyRC)*627.5095;
                    aac=(Reaction[I].energyPC-Reaction[I].energyRC)*627.5095+Reaction[I].Psi*((Specie[p1].zep+zepP2)-(Specie[r1].zep+zepR2))*1.0*pow(10,-3);
                    vdu=(Specie[ts].energy-Reaction[I].energyRC)*627.5095;
                    vdc=(Specie[ts].energy-Reaction[I].energyRC)*627.5095+Reaction[I].Psi*((Specie[ts].zep-(Specie[r1].zep+zepR2)))*1.0*pow(10,-3);
                    Reaction[I].veff=vdc;
                }
            }
        }
        bbu=(2*vdu-aau)+2*sqrt(vdu*(vdu-aau));
        bbc=(2*vdc-aac)+2*sqrt(vdc*(vdc-aac));
        Reaction[I].Mass=Specie[r1].Mass*MassR2/Specie[ts].Mass;
        Reaction[I].vref1=vdu-aau;
        Reaction[I].vref2=vdc-aac;
        fca2=627.5095*mc*pow((2.0*pi*c*0.529177249*pow(10,-10)),2.0)/(4.36*pow(10,-18));
        alfa=sqrt(fabs(fca2*Reaction[I].Mass*bbu*pow(fabs(Specie[ts].freqN),2.0)/(2.0*vdu*(vdu-aau))));
        sou=-(1.0/alfa)*(log((aau+bbu)/(fabs(bbu-aau))));
        soc=-(1.0/alfa)*(log((aac+bbc)/(fabs(bbc-aac))));
        cc=Reaction[I].Psi*(Specie[r1].zep+zepR2)*1.0*pow(10,-3);

        for(i=0; i<101; i++){
            Reaction[I].s[i]=-xmep*(1-2.0*((double)i)/100.0);
            if(Reaction[I].s[i]==-0.0)
                Reaction[I].s[i]=0.0;
            yu=exp(alfa*(Reaction[I].s[i]-sou));
            Reaction[I].yc=exp(alfa*(Reaction[I].s[i]-soc));
            Reaction[I].vmep[i]=aau*yu/(1+yu)+(bbu*yu/(pow((1+yu),2)));
            Reaction[I].vagg[i]=aac*Reaction[I].yc/(1+Reaction[I].yc)+(bbc*Reaction[I].yc/(pow((1+Reaction[I].yc),2)))+cc;
        }

        //calculate tunneling of Eckart
        a1=vdc*fca1a2/(fabs(Specie[ts].freqN));
        a2=(vdc-aac)*fca1a2/(fabs(Specie[ts].freqN));
        if(a1>a2){
            vmax=vdc;
            att1=a1;
            att2=a2;
        }
        else{
            vmax=vdc-aac;
            att1=a2;
            att2=a1;
        }
        a1=fabs(att1);
        a2=fabs(att2);
        vint=vmax*4.36*pow(10,-18)/(627.5095*kB);
        a1t=(1.0/pi)*sqrt(a1)*pow((1.0/sqrt(a1)+1.0/sqrt(a2)),-1);
        b1t=(1.0/pi)*pow((1.0/sqrt(a1)+1.0/sqrt(a2)),-1);
        xmax=(pow(a1,3)*pow(a2,2)-2*pow(a1,2.5)*pow(a2,2.5)+pow(a1,2)*pow(a2,3)-25538*pow(a1,2)*a2*pow(pi,2)+25538*a1*pow(a2,2)*pow(pi,2)+163047361*a1*pow(pi,4)+326094722*sqrt(a1)*sqrt(a2)*pow(pi,4)+163047361*a2*pow(pi,4))/(51076*pow(a1,2)*a2*pow(pi,2));
        ter1=alf(0.0,a1t);
        ter2=bt(0.0,a1,a2,b1t);
        ter3=dtt(a1,a2);
        ttest1=ter1+ter2;
        ttest2=ter1-ter2;
        f0=xkappa(ttest2,ttest1,ter3);

        for(k=0; k<24; k++){
            if(ter3 <= 113){
                ter1=alf(xmax,a1t);
                ter2=bt(xmax,a1,a2,b1t);
                ttest1=(ter1+ter2);
                ttest2=fabs(ter1-ter2);
                f2n=xkappa(ttest2,ttest1,ter3)*exp(-xmax*vint/Reaction[I].temp[k]);
                sm=0.0;
                sp=0.0;
                xhp=xmax/(2.0*5000.0);
                for(i=1; i<=2*5000; i=i+2){
                    xmp=((double)i*xhp);
                    ter1=alf(xmp,a1t);
                    ter2=bt(xmp,a1,a2,b1t);
                    ttest1=(ter1+ter2);
                    ttest2=(ter1-ter2);
                    sm=sm+xkappa(ttest2,ttest1,ter3)*exp(-xmp*vint/Reaction[I].temp[k]);
                }
                for(i=2; i<=2*5000-1; i=i+2){
                    xpp=((double)i*xhp);
                    ter1=alf(xpp,a1t);
                    ter2=bt(xpp,a1,a2,b1t);
                    ttest1=(ter1+ter2);
                    ttest2=(ter1-ter2);
                    sp=sp+xkappa(ttest2,ttest1,ter3)*exp(-xpp*vint/Reaction[I].temp[k]);
                }
                xint=(xhp/3.0)*(f0+4.0*sm+2.0*sp+f2n);
                Reaction[I].K_Eckart[k]=exp(vint/Reaction[I].temp[k])*xint*vint/Reaction[I].temp[k]*Reaction[I].K[k];
            }
            else
                Reaction[I].K_Eckart[k]=Reaction[I].K[k];
        }
    }
}

//functions for small curvature tunneling correction
double integral(double F[],int TAM, double h){
    int i;
    double S=0,aux=0;
    for(i=0; i <TAM; i++){
        if(i == 1 || i == TAM)
            aux=1;
        else{
            if(i%2 != 0)
                aux=2;
            else
                aux=4;
        }
        S=S+aux*F[i];
    }
    return (S*h/3.0);
}

//functions for small curvature tunneling correction
double max(double a, double b){
    if(a > b) return a;
    else return b;
}

//functions for small curvature tunneling correction
double Prob(double E0, double E, double vagmax,double S){
    double aux,prob=0;
    if(E0 <= E && E <= vagmax){
        prob=S;
    }
    else{
        if(vagmax <= E && E <= (2.0*vagmax-E0)){
            aux=Prob(E0,2*vagmax-E,vagmax,S);
            prob=1-aux;
        }
        else{
            if((2*vagmax-E0) <= E){
                prob=1;
            }
        }
    }
    return(prob);
}

//calculate small curvature tunneling correction
void tunneling_small_curvature(){
    int i,I,j,k;
    double E1[1000],E2[1000],S[101];
    double hs,hE1=250.0/1000,hE2,E,E0,vagmax,prob;

    hs=xmep*2.0/101;

    for(I=0; I < Nreactions; I++){
        E0=max(Reaction[I].vagg[0],Reaction[I].vagg[100]);
        vagmax=Reaction[I].vagg[0];

        //printf("\n%s %s %s %s %s\n",Reaction[I].react1,Reaction[I].react2,Reaction[I].ts,Reaction[I].prod1,Reaction[I].prod2);
        for(k=0; k<24; k++){
            E=0;
            for(i=0; i<1000; i++){
                for(j=0; j<101; j++){
                    S[j]=sqrt(2.0*Reaction[I].Mass*fabs(E-Reaction[I].vagg[j]));
                    if(vagmax < Reaction[I].vagg[j])
                        vagmax=Reaction[I].vagg[j];
                }
                prob=1.0/(1.0+exp(4*pi*integral(S,101,hs)/fca1a2));
                prob=Prob(E0,E,vagmax,prob);
                E1[i]=prob*exp(-E/(r*Reaction[I].temp[k]));
                E=E+hE1;
            }

            E=Reaction[I].vag;
            hE2=(250.0-E)/1000;
            for(i=0; i<1000; i++){
                E2[i]=exp(-E/(A*kB*Reaction[I].temp[k]));
                E=E+hE2;
            }
            Reaction[I].K_TSC[k]=integral(E1,1000,hE1)/integral(E2,1000,hE2)*Reaction[I].K[k];
            //printf("%e %e %e %e\n",log10(Reaction[I].K[k]),log10(Reaction[I].K_Wigner[k]),log10(Reaction[I].K_Eckart[k]),log10(Reaction[I].K_TSC[k]));
        }
    }
}

//calculate rate in Arrhenius form
void arrhenius(){
    int nt=24,k,i,I,inc,inc1,inc2,inc3;
    double slnt,s1t,st,stlnt,slnt2,slnt1t;
    double slnk,stlnk,slntk,serr,aux3,aux32,aux33,xc1,coefa[4][4],coefb[4],caa[4][4],cbb[4],xarh[30];

    for(I=0; I < Nreactions; I++){
        slnt=s1t=st=stlnt=slnt2=slnt1t=0.0;
        for(k=0; k<nt; k++){
            slnt=slnt+log(Reaction[I].temp[k]);
            s1t=s1t+1/Reaction[I].temp[k];
            st=st+Reaction[I].temp[k];
            stlnt=stlnt+Reaction[I].temp[k]*log(Reaction[I].temp[k]);
            slnt2=slnt2+pow((log(Reaction[I].temp[k])),2);
            slnt1t=slnt1t+(log(Reaction[I].temp[k]))/Reaction[I].temp[k];
        }
        coefa[1][1]=st;
        coefa[1][2]=stlnt;
        coefa[1][3]=-nt;
        coefa[2][1]=nt;
        coefa[2][2]=slnt;
        coefa[2][3]=-s1t;
        coefa[3][1]=slnt;
        coefa[3][2]=slnt2;
        coefa[3][3]=-slnt1t;

        //printf("\n%s %s %s %s %s\n",Reaction[I].react1,Reaction[I].react2,Reaction[I].ts,Reaction[I].prod1,Reaction[I].prod2);
        for(i=0; i<4; i++){
            slnk=stlnk=slntk=serr=0;
            for(k=0; k<nt; k++){
                if(i == 0)
                    xarh[k]=Reaction[I].K[k];
                if(i == 1)
                    xarh[k]=Reaction[I].K_Wigner[k];
                if(i == 2)
                    xarh[k]=Reaction[I].K_Eckart[k];
                if(i == 3)
                    xarh[k]=Reaction[I].K_TSC[k];

                slnk=slnk+log(xarh[k]);
                stlnk=stlnk+Reaction[I].temp[k]*log(xarh[k]);
                slntk=slntk+log(xarh[k])*log(Reaction[I].temp[k]);
            }
            coefb[1]=stlnk;
            coefb[2]=slnk;
            coefb[3]=slntk;
            caa[1][1]=coefa[1][1];
            caa[1][2]=coefa[1][2];
            caa[1][3]=coefa[1][3];
            caa[2][1]=coefa[2][1]*coefa[1][1]-coefa[2][1]*coefa[1][1];
            caa[2][2]=coefa[2][1]*coefa[1][2]-coefa[2][2]*coefa[1][1];
            caa[2][3]=coefa[2][1]*coefa[1][3]-coefa[2][3]*coefa[1][1];
            caa[3][1]=coefa[3][1]*coefa[1][1]-coefa[3][1]*coefa[1][1];
            aux32=coefa[3][1]*coefa[1][2]-coefa[3][2]*coefa[1][1];
            aux33=coefa[3][1]*coefa[1][3]-coefa[3][3]*coefa[1][1];
            caa[3][2]=0;
            caa[3][3]=caa[2][2]*aux33-aux32*caa[2][3];
            cbb[1]=coefb[1];
            cbb[2]=coefa[2][1]*coefb[1]-coefa[1][1]*coefb[2];
            aux3=coefa[3][1]*coefb[1]-coefa[1][1]*coefb[3];
            cbb[3]=caa[2][2]*aux3-aux32*cbb[2];
            inc=i*4;
            inc1=inc+1;
            inc2=inc+2;
            inc3=inc+3;
            Reaction[I].arh[inc2]=r*(cbb[3]/caa[3][3]);
            Reaction[I].arh[inc1]=(cbb[2]-caa[2][3]*Reaction[I].arh[inc2]/r)/caa[2][2];
            Reaction[I].arh[inc]=exp((cbb[1]-caa[1][2]*Reaction[I].arh[inc1]-caa[1][3]*Reaction[I].arh[inc2]/r)/caa[1][1]);
            for(k=0; k<nt; k++){
                xc1=Reaction[I].arh[inc]*(pow(Reaction[I].temp[k],Reaction[I].arh[inc1]))*exp(-Reaction[I].arh[inc2]/(r*Reaction[I].temp[k]));
                serr=serr+pow(((log(Reaction[I].K[k])-log(xc1))),2);
            }
            Reaction[I].arh[inc3]=sqrt(serr/nt);
            //printf("%e %e %e %e\n",Reaction[I].arh[inc],Reaction[I].arh[inc1],Reaction[I].arh[inc2],Reaction[I].arh[inc3]);
        }
    }
}

//calculate rovibrational levels
void rovibrational_levels(int I){
    int i,j,l,k;
    double Req,Mi,We,Be,di[5],ai[10],ca;

    //FILE *EN=fopen("En.dat","w");
    //fprintf(EN,"%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",Specie[I].di[0],Specie[I].di[1],Specie[I].di[2],Specie[I].di[3],Specie[I].di[4],Specie[I].Eref,Specie[I].De,Specie[I].Req);

    for(i=0; i<11; i++)
        for(j=0;j<11;j++)
            Specie[I].Y[i][j]=0.0;

    for(i=0; i<5; i++)
        di[i]=Specie[I].di[i];

    ca=2.99793;
    Mi=Specie[I].mass[0]*Specie[I].mass[1]*mc/(Specie[I].mass[0]+Specie[I].mass[1]);
    We=sqrt(hJ*(-Specie[I].De/627.5095*(-pow(di[0],2.0)+2.0*di[1]))/(4.0*ca*ca*pi*pi*Mi));
    Be=h/(8.0*pi*pi*Mi*Specie[I].Req*Specie[I].Req*c*pow(10,-20.0));
    Req=Specie[I].Req;

    ai[0]=We*We/(4.0*Be);
    ai[1]=Req*(2.0*pow(di[0],3.0)-6.0*di[0]*di[1]+6.0*di[2])/(3.0*(-pow(di[0],2.0)+2*di[1]));
    ai[2]=pow(Req,2.0)*(-3*pow(di[0],4.0)+12*pow(di[0],2.0)*di[1]-24*di[0]*di[2]+24*di[3])/(12*(-pow(di[0],2.0)+2*di[1]));
    ai[3]=pow(Req,3.0)*(4*pow(di[0],5.0)-20*pow(di[0],3.0)*di[1]+60*pow(di[0],2.0)*di[2]-120*di[0]*di[3]+120*di[4])/(60*(-pow(di[0],2.0)+2*di[1]));
    ai[4]=pow(Req,4.0)*(-5*pow(di[0],6.0)+30*pow(di[0],4.0)*di[1]-120*pow(di[0],3.0)*di[2]+360*pow(di[0],2.0)*di[3]-720*di[0]*di[4])/(360*(-pow(di[0],2.0)+2*di[1]));
    ai[5]=pow(Req,5.0)*(6*pow(di[0],7.0)-42*pow(di[0],5.0)*di[1]+210*pow(di[0],4.0)*di[2]-840*pow(di[0],3.0)*di[3]+2520*pow(di[0],2.0)*di[4])/(2520*(-pow(di[0],2.0)+2*di[1]));
    ai[6]=pow(Req,6.0)*(-7*pow(di[0],8.0)+56*pow(di[0],6.0)*di[1]-336*pow(di[0],5.0)*di[2]+1680*pow(di[0],4.0)*di[3]-6720*pow(di[0],3.0)*di[4])/(20160*(-pow(di[0],2.0)+2*di[1]));
    ai[7]=pow(Req,7.0)*(8*pow(di[0],9.0)-72*pow(di[0],7.0)*di[1]+504*pow(di[0],6.0)*di[2]-3024*pow(di[0],5.0)*di[3]+15120*pow(di[0],4.0)*di[4])/(181440*(-pow(di[0],2.0)+2*di[1]));
    ai[8]=pow(Req,8.0)*(-9*pow(di[0],10.0)+90*pow(di[0],8.0)*di[1]-720*pow(di[0],7.0)*di[2]+5040*pow(di[0],6.0)*di[3]-30240*pow(di[0],5.0)*di[4])/(1814400*(-pow(di[0],2.0)+2*di[1]));

    Specie[I].Y[0][0]=Be/8.0*(3.0*ai[2]-7.0*pow(ai[1],2.0)/4.0);

    Specie[I].Y[1][0]=We*(1.0+pow(Be,2.0)/(4.0*pow(We,2.0))*(25.0*ai[4]-95.0*ai[1]*ai[3]/2.0-67.0*pow(ai[2],2.0)/4.0+459.0*pow(ai[1],2.0)*ai[2]/8.0-1155.0*pow(ai[1],4.0)/64.0));

    Specie[I].Y[2][0]=Be/2.0*(3.0*(ai[2]-5.0*pow(ai[1],2.0)/4.0)+pow(Be,2.0)/(2.0*pow(We,2.0))*(245.0*ai[6]-1365.0*ai[1]*ai[5]/2.0-885.0*ai[2]*ai[4]/2.0-1085.0*pow(ai[3],2.0)/4.0+8535.0*pow(ai[1],2.0)*ai[4]/8.0+1707.0*pow(ai[2],3.0)/8.0+7335.0*ai[1]*ai[2]*ai[3]/4.0-23865.0*pow(ai[1],3.0)*ai[3]/16.0-62013.0*pow(ai[1]*ai[2],2.0)/32.0+239985.0*pow(ai[1],4.0)*ai[2]/128.0-209055.0*pow(ai[1],6.0)/512.0));

    Specie[I].Y[3][0]=pow(Be,2.0)/(2.0*We)*(10.0*ai[4]-35.0*ai[1]*ai[3]-17.0*pow(ai[2],2.0)/2.0+225.0*pow(ai[1],2.0)*ai[2]/4.0-705.0*pow(ai[1],4.0)/32.0);

    Specie[I].Y[4][0]=5*pow(Be,3.0)/pow(We,2.0)*(7*ai[6]/2.0-63/4.0*ai[1]*ai[5]-33/4.0*ai[2]*ai[4]-63/8.0*pow(ai[3],2.0)+543/16.0*pow(ai[1],2.0)*ai[4]+75/16.0*pow(ai[2],3.0)+483/8.0*ai[1]*ai[2]*ai[3]-1953/32.0*pow(ai[1],3.0)*ai[3]-4989/64.0*pow(ai[1]*ai[2],2.0)+23265/256.0*pow(ai[1],4.0)*ai[2]-23151/1024.0*pow(ai[1],6.0));

    Specie[I].Y[0][1]=Be*(1+pow(Be,2.0)/(2.0*pow(We,2.0))*(15.0+14.0*ai[1]-9.0*ai[2]+15.0*ai[3]-23.0*ai[1]*ai[2]+21/2.0*(pow(ai[1],2.0)+pow(ai[1],3.0))));

    Specie[I].Y[1][1]=pow(Be,2.0)/We*(6.0*(1.0+ai[1])+pow(Be,2.0)/pow(We,2.0)*(175.0+285.0*ai[1]-335/2.0*ai[2]+190.0*ai[3]-225/2.0*ai[4]+175.0*ai[5]+2295/8.0*pow(ai[1],2.0)-459.0*ai[1]*ai[2]+1425/4.0*ai[1]*ai[3]-795/2.0*ai[1]*ai[4]+1005/8.0*pow(ai[2],2.0)-715/2.0*ai[2]*ai[3]+1155/4.0*pow(ai[1],3.0)-9639/16.0*pow(ai[1],2.0)*ai[2]+5145/8.0*pow(ai[1],2.0)*ai[3]+4677/8.0*pow(ai[2],2.0)*ai[1]-14259/16.0*pow(ai[1],3.0)*ai[2]+31185/128.0*(pow(ai[1],4.0)+pow(ai[1],5.0))));

    Specie[I].Y[2][1]=6.0*pow(Be,3.0)/pow(We,2.0)*(5.0+10.0*ai[1]-3.0*ai[2]+5.0*ai[3]-13.0*ai[1]*ai[2]+15/2.0*(pow(ai[1],2.0)+pow(ai[1],3.0)));

    Specie[I].Y[3][1]=20.0*pow(Be,4.0)/pow(We,3.0)*(7.0+21.0*ai[1]-17/2.0*ai[2]+14.0*ai[3]-9/2.0*ai[4]+7.0*ai[5]+225/8.0*pow(ai[1],2.0)-45.0*ai[1]*ai[2]+105/4.0*ai[1]*ai[3]-51/2.0*ai[1]*ai[4]+51/8.0*pow(ai[2],2.0)-45/2.0*ai[2]*ai[3]+141/4.0*pow(ai[1],3.0)-945/16.0*pow(ai[1],2.0)*ai[2]+435/8.0*pow(ai[1],2.0)*ai[3]+411/8.0*ai[1]*pow(ai[2],2.0)-1509/16.0*pow(ai[1],3.0)*ai[2]+3807/128.0*(pow(ai[1],4.0)+pow(ai[1],5.0)));

    Specie[I].Y[0][2]=-4*pow(Be,3.0)/pow(We,2.0)*(1.0+pow(Be,2.0)/(2.0*pow(We,2.0))*(163.0+199.0*ai[1]-119.0*ai[2]+90.0*ai[3]-45.0*ai[4]-207*ai[1]*ai[2]+205/2.0*ai[1]*ai[3]-333/2.0*pow(ai[1],2.0)*ai[2]+693/4.0*pow(ai[1],2.0)+46.0*pow(ai[2],2.0)+126.0*(pow(ai[1],3.0)+pow(ai[1],4.0)/2.0)));

    Specie[I].Y[1][2]=-12.0*pow(Be,4.0)/pow(We,3.0)*(19/2.0+9.0*ai[1]+9/2.0*pow(ai[1],2.0)-4.0*ai[2]);

    Specie[I].Y[2][2]=-24.0*pow(Be,5.0)/pow(We,4.0)*(65.0+125.0*ai[1]-61.0*ai[2]+30.0*ai[3]-15.0*ai[4]+495/4.0*pow(ai[1],2.0)-117.0*ai[1]*ai[2]+26.0*pow(ai[2],2.0)+95/2.0*ai[1]*ai[3]-207/2.0*pow(ai[1],2.0)*ai[2]+90.0*(pow(ai[1],3.0)+pow(ai[1],4.0)/2.0));

    Specie[I].Y[0][3]=16.0*pow(Be,5.0)*(3.0+ai[1])/pow(We,4.0);

    Specie[I].Y[1][3]=12.0*pow(Be,6.0)/pow(We,5.0)*(233.0+279.0*ai[1]+189.0*pow(ai[1],2.0)+63.0*pow(ai[1],3.0)-88.0*ai[1]*ai[2]-120.0*ai[2]+80/3.0*ai[3]);

    Specie[I].Y[0][4]=-64.0*pow(Be,7.0)/pow(We,6.0)*(13.0+9.0*ai[1]-ai[2]+9/4.0*pow(ai[1],2.0));

    Specie[I].Y[0][10]=2048*(pow(Be,19.0)/pow(We,18.0))*(-938223*pow(ai[1],8.0)-9382230*pow(ai[1],7.0)+2918916*pow(ai[1],6.0)*ai[2]-45243198*pow(ai[1],6.0)+23351328*pow(ai[1],5.0)*ai[2]
-972972*pow(ai[1],5.0)*ai[3]-137884032*pow(ai[1],5.0)-2594592*pow(ai[1],4.0)*pow(ai[2],2.0)+86177520*pow(ai[1],4.0)*ai[2]-6949800*pow(ai[1],4.0)*ai[3]+277992*pow(ai[1],4.0)*ai[4]
-291677760*pow(ai[1],4.0)-14826240*pow(ai[1],3.0)*pow(ai[2],2.0)+1235520*pow(ai[1],3.0)*ai[2]*ai[3]+188559360*pow(ai[1],3.0)*ai[2]-22096800*pow(ai[1],3.0)*ai[3]+1710720
*pow(ai[1],3.0)*ai[4]-66528*pow(ai[1],3.0)*ai[5]-440757504*pow(ai[1],3.0)+658944*pow(ai[1],2.0)*pow(ai[2],3.0)-35354880*pow(ai[1]*ai[2],2.0)+5702400*pow(ai[1],2.0)*ai[2]*ai[3]
-228096*pow(ai[1],2.0)*ai[2]*ai[4]+259269120*pow(ai[1],2.0)*ai[2]-118800*pow(ai[1]*ai[3],2.0)-39283200*pow(ai[1],2.0)*ai[3]+4419360*pow(ai[1],2.0)*ai[4]-332640*pow(ai[1],2.0)*ai[5]
+12672*pow(ai[1],2.0)*ai[6]-467470080*pow(ai[1],2.0)+2027520*pow(ai[2],3.0)*ai[1]-253440*ai[1]*pow(ai[2],2.0)*ai[3]-41902080*pow(ai[2],2.0)*ai[1]+9820800*ai[1]*ai[2]*ai[3]
-760320*ai[1]*ai[2]*ai[4]+29568*ai[1]*ai[2]*ai[5]+213700608*ai[1]*ai[2]-39600*ai[1]*pow(ai[3],2.0)+31680*ai[1]*ai[3]*ai[4]-39283200*ai[1]*ai[3]+5713920*ai[1]*ai[4]
-624960*ai[1]*ai[5]+46080*ai[1]*ai[6]-1728*ai[1]*ai[7]-320550912*ai[1]-22528*pow(ai[2],4.0)+1745920*pow(ai[2],3.0)-422400*pow(ai[2],2.0)*ai[3]+16896*pow(ai[2],2.0)*ai[4]
-20951040*pow(ai[2],2.0)+17600*ai[2]*pow(ai[3],2.0)+6348800*ai[2]*ai[3]-714240*ai[2]*ai[4]+53760*ai[2]*ai[5]-2048*ai[2]*ai[6]+83105792*ai[2]-372000*pow(ai[3],2.0)
+57600*ai[3]*ai[4]-2240*ai[3]*ai[5]-17808384*ai[3]-1152*pow(ai[4],2.0)+3142656*ai[4]-444416*ai[5]+47616*ai[6]-2345*ai[7]+128*ai[8]-109818368);

    Specie[I].Y[1][8]=18*(pow(Be,16.0)/pow(We,15.0))*(-94616181*pow(ai[1],8.0)-756929448*pow(ai[1],7.0)+349383456*pow(ai[1],6.0)*ai[2]-2976321348*pow(ai[1],6.0)+2271758400*pow(ai[1],5.0)*ai[2]
-146214720*pow(ai[1],5.0)*ai[3]-7563816936*pow(ai[1],5.0)-369048960*pow(ai[1],4.0)*pow(ai[2],2.0)+6956738784*pow(ai[1],4.0)*ai[2]-866816640*pow(ai[1],4.0)*ai[3]+54297216*pow(ai[1],4.0)
*ai[4]-13700672910*pow(ai[1],4.0)-1740787200*pow(ai[1],3.0)*pow(ai[2],2.0)+220492800*pow(ai[1],3.0)*ai[2]*ai[3]+12936336768*pow(ai[1],3.0)*pow(ai[2],2.0)-2339418240*pow(ai[1],3.0)
*ai[3]+283627008*pow(ai[1],3.0)*ai[4]-17224704*pow(ai[1],3.0)*ai[5]-18291196440*pow(ai[1],3.0)+111716352*pow(ai[1],2.0)*pow(ai[2],3.0)-3504031488*pow(ai[1],2.0)*pow(ai[2],2.0)
+857364480*pow(ai[1],2.0)*ai[2]*ai[3]-52918272*pow(ai[1],2.0)*ai[2]*ai[4]+15544252512*pow(ai[1],2.0)*ai[2]-26496000*pow(ai[1],2.0)*pow(ai[3],2.0)-3621534720*pow(ai[1],2.0)*ai[3]
+637168896*pow(ai[1],2.0)*ai[4]-74769408*pow(ai[1],2.0)*ai[5]+4399104*pow(ai[1],2.0)*ai[6]-17776808100*pow(ai[1],2.0)+288239616*ai[1]*pow(ai[2],3.0)-54005760*ai[1]*pow(ai[2],2.0)*ai[3]
-3596524032*ai[1]*pow(ai[2],2.0)+1274496000*ai[1]*ai[2]*ai[3]-151990272*ai[1]*ai[2]*ai[4]+9117696*ai[1]*ai[2]*ai[5]+11565972288*ai[1]*ai[2]-75878400*ai[1]*pow(ai[3],2.0)+91544560*ai[1]*ai[3]*ai[4]
-3247145280*ai[1]*ai[3]+736154112*ai[1]*ai[4]-125174784*ai[1]*ai[5]+14204928*ai[1]*ai[6]-811008*ai[1]*ai[7]-11646924120*ai[1]-4575232*pow(ai[2],4.0)+213223424*pow(ai[2],3.0)
-77045760*pow(ai[2],2.0)*ai[3]+4694016*pow(ai[2],2.0)*ai[4]-1604168064*pow(ai[2],2.0)+4710400*ai[2]*pow(ai[3],2.0)+731013120*ai[2]*ai[3]-126302208*ai[2]*ai[4]+14622720*ai[2]*ai[5]
-851968*ai[2]*ai[6]+4215991968*ai[2]-62873600*pow(ai[3],2.0)+14622720*ai[3]*ai[4]-860160*ai[3]*ai[5]-1365091200*ai[3]-430080*pow(ai[4],2.0)+373189248*ai[4]-81679360*ai[5]
+13410304*ai[6]-1474560*ai[7]+81920*ai[8]-3995947365);

    Specie[I].Y[2][6]=384*(pow(Be,13.0)/pow(We,12.0))*(-933120*pow(ai[1],8.0)-5598720*pow(ai[1],7.0)+3941379*pow(ai[1],6.0)*ai[2]-17184960*pow(ai[1],6.0)+19293714*pow(ai[1],5.0)*ai[2]-1903365
*pow(ai[1],5.0)*ai[3]-35582490*pow(ai[1],5.0)-4841820*pow(ai[1],4.0)*pow(ai[2],2.0)+46211310*pow(ai[1],4.0)*ai[2]-8496900*pow(ai[1],4.0)*ai[3]+3415200*pow(ai[1],3.0)*ai[2]*ai[3]
+70093440*pow(ai[1],3.0)*ai[2]-17896680*pow(ai[1],3.0)*ai[3]+3275640*pow(ai[1],3.0)*ai[4]-317940*pow(ai[1],3.0)*ai[5]-64518840*pow(ai[1],3.0)+1730880*pow(ai[1],2.0)*pow(ai[2],3.0)
-27588384*pow(ai[1],2.0)*pow(ai[2],2.0)+10173600*pow(ai[1],2.0)*ai[2]*ai[3]-984816*pow(ai[1],2.0)*ai[2]*ai[4]+71776080*pow(ai[1],2.0)*ai[2]-500400*pow(ai[1],2.0)*pow(ai[3],2.0)
-22533600*pow(ai[1],2.0)*ai[3]+5765040*pow(ai[1],2.0)*ai[4]-1058400*pow(ai[1],2.0)*ai[5]+101920*pow(ai[1],2.0)*ai[6]-57115440*pow(ai[1],2.0)+3442176*ai[1]*pow(ai[2],3.0)-1008640
*ai[1]*pow(ai[2],2.0)*ai[3]-23339520*ai[1]*pow(ai[2],2.0)+12036000*ai[1]*ai[2]*ai[3]-2192256*ai[1]*ai[2]*ai[4]+209664*ai[1]*ai[2]*ai[5]+47536752*ai[1]*ai[2]-1118400*ai[1]*pow(ai[3],2.0)
+214960*ai[1]*ai[3]*ai[4]-17222400*ai[1]*ai[3]+5446320*ai[1]*ai[4]-1411200*ai[1]*ai[5]+258720*ai[1]*ai[6]-24528*ai[1]*ai[7]-3051520*ai[1]-84992*pow(ai[2],4.0)+2042112*pow(ai[2],3.0)
-1125120*pow(ai[2],2.0)*ai[3]+107776*pow(ai[2],2.0)*ai[4]-8991840*pow(ai[2],2.0)+109200*ai[2]*pow(ai[3],2.0)+5744000*ai[2]*ai[3]-1469136*ai[2]*ai[4]+265776*ai[2]*ai[5]-25088*ai[2]*ai[6]
+16095472*ai[2]-753400*pow(ai[3],2.0)+273840*ai[3]*ai[4]-25760*ai[3]*ai[5]-6503920*ai[3]-13040*pow(ai[4],2.0)+2378640*ai[4]-768880*ai[5]+200480*ai[6]-36288*ai[7]+3360*ai[8]-11612320);

    Specie[I].Y[3][4]=(5/32.0*pow(Be,10.0)/pow(We,9.0))*(-159014583*pow(ai[1],8.0)-636058332*pow(ai[1],7.0)+747564552*pow(ai[1],6.0)*ai[2]-1428607206*pow(ai[1],6.0)+2424873024*pow(ai[1],5.0)*ai[2]
-404593920*pow(ai[1],5.0)*ai[3]-2333591316*pow(ai[1],5.0)-1034779824*pow(ai[1],4.0)*pow(ai[2],2.0)+4229176320*pow(ai[1],4.0)*ai[2]-1192959360*pow(ai[1],4.0)*ai[3]+201000960*pow(ai[1],4.0)*ai[4]
-2984513175*pow(ai[1],4.0)-2437454016*pow(ai[1],3.0)*pow(ai[2],2.0)+827554560*pow(ai[1],3.0)*ai[2]*ai[3]+5059842720*pow(ai[1],3.0)*ai[2]-1831849920*pow(ai[1],3.0)*ai[3]+523307520*pow(ai[1],3.0)*ai[4]
-89134080*pow(ai[1],3.0)*ai[5]-3036118848*pow(ai[1],3.0)+424947456*pow(ai[1],2.0)*pow(ai[2],3.0)-2787189984*pow(ai[1]*ai[2],2.0)+1605830400*pow(ai[1],2.0)*ai[2]*ai[3]-2770790040
*pow(ai[1],2.0)*ai[2]*ai[4]+4311029880*pow(ai[1],2.0)*ai[2]-138758400*pow(ai[1],2.0)*pow(ai[3],2.0)-1834870080*pow(ai[1],2.0)*ai[3]+670510080*pow(ai[1],2.0)*ai[4]-193729536*pow(ai[1],2.0)*ai[5]
+33589248*pow(ai[1],2.0)*ai[6]-2404222272*pow(ai[1],2.0)+552752640*pow(ai[2],3.0)*ai[1]-286932480*ai[1]*pow(ai[2],2.0)*ai[3]-1839511872*ai[1]*pow(ai[2],2.0)+1351810560*ai[1]*ai[2]*ai[3]
-400997376*ai[1]*ai[2]*ai[4]+70404096*ai[1]*ai[2]*ai[5]+2446931712*ai[1]*ai[2]-198182400*ai[1]*pow(ai[3],2.0)+70471680*ai[1]*ai[3]*ai[4]-1177272480*ai[1]*ai[3]+505409280*ai[1]*ai[4]
-184859136*ai[1]*ai[5]+54878208*ai[1]*ai[6]-9805824*ai[1]*ai[7]-1357226496*ai[1]-24566784*pow(ai[2],4.0)+233680128*pow(ai[2],3.0)-208266240*pow(ai[2],2.0)*ai[3]+36584448*pow(ai[2],2.0)*ai[4]
-584617584*pow(ai[2],2.0)+36864000*ai[2]*pow(ai[3],2.0)+499001600*ai[2]*ai[3]-190373376*ai[2]*ai[4]+58189824*ai[2]*ai[5]-10379264*ai[2]*ai[6]+719257344*ai[2]-91974400*pow(ai[3],2.0)
+57108480*ai[3]*ai[4]-10465280*ai[3]*ai[5]-377807360*ai[3]-5191680*pow(ai[4],2.0)+186075840*ai[4]-78887424*ai[5]+29560832*ai[6]-9289728*ai[7]+1720320*ai[8]-412452096);

    for(i=0;i<50;i++){
        for(j=0;j<50;j++){
            Specie[I].En[i][j]=0.0;
            for(l=0;l<=5;l++){
                for(k=0;k<=10;k++)
                    Specie[I].En[i][j]+=Specie[I].Y[l][k]*pow((i+1/2.0),l)*pow(j*(j+1),k);
            }
            //fprintf(EN,"%e\t",Specie[I].En[i][j]);
        }
        //fprintf(EN,"\n");
    }

    for(i=0;i<50;i++){
        for(j=0;j<50;j++){
          Specie[I].DEn[i][j]=(Specie[I].En[i][j]-Specie[I].En[0][0])/219474.63067;
          //fprintf(EN,"%e\t",Specie[I].DEn[i][j]);
        }
      //(EN,"\n");
    }

    //fclose(EN);
}
