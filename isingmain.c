#include <math.h> 
#include <stdlib.h> 
#include <stdio.h> 
#include "cpgplot.h" 
#include "ising.h" 
#include <stdbool.h>
int i,j,k; 
int size;                       //lattice size 
double T;                       //temperature
const double minT=0.7;          //minimum temperature
double change=0.05;             //size of steps for temperature loop
long unsigned int mcs=10000;    //number of Monte Carlo steps
int transient=10000;            //number of transient steps
const int nintervals=90;
float dim[]={5, 10, 16, 25, 50, 100};
float critHCT[6];

//main program

void main()
{
    for(k=0;k<6;k++) {
        size=dim[k]; 
        printf("size=%i",size);
    int n=size*size;
    int lat[size+1][size+1];            //2d lattice for spins 
    double norm=(1.0/(double)(mcs*n));  //normalizaton for averaging

    //declaring variables to be used in calculating observables
    double E,Esq=0,etot=0,etotsq=0;
    double M=0,M_avg=0,mtot=0;
    double Mabs=0,mabstot=0;
    int de;
    lat_type pos;
    double Esq_avg[nintervals+1];
    float t[nintervals+1],C[nintervals+1],E_avg[nintervals+1],Mabs_avg[nintervals+1];

    //intitialize lattice to random configuration
    initialize(size, lat);

    //Temperature loop
    for(i=0;i<=nintervals; i++) {

    t[i]=0.7+i*change;
    T=t[i];
    printf("Temperature=%f\n",T);
    
    //Transient function
    transient_results(size, lat);

    //Observables get equilibrated lattice config values
    M=total_magnetization(size,lat);
    Mabs=abs(total_magnetization(size,lat));
    E=total_energy(size,lat);
    
    //intitialize summation variables at each temperature step
    etot=0;
    etotsq=0;
    mtot=0;
    mabstot=0;

    //Monte Carlo Loop
    for(int a=1;a<=mcs;a++)
        {
            //Metropolis Loop
            for(int b=1;b<=n;b++)
                {
                    choose_random_pos_lat(&pos); 
                    if(test_flip(pos,&de,size,lat))
                        {
                            flip(pos,size,lat);
                            
                            //adjust observables 
                            E=E+2*de;
                            M=M+2*lat[pos.x][pos.y];
                            Mabs=Mabs+abs(lat[pos.x][pos.y]);
                        }
                }
        //keep summation of observables 
        etot=etot+E/((double)2.0);
        etotsq=etotsq+(E/((double)2.0))*(E/((double)2.0));
        mtot=mtot+M;
        mabstot=mabstot+(sqrt(M*M));
        }
    //average observables
    E_avg[i]=etot*norm;            //<E> spin
    Esq_avg[i]=etotsq*norm;         //<E^2> spin
    M_avg=mtot*norm;               //<M> spin
    Mabs_avg[i]=mabstot*norm;      //<|M|> spin

    C[i]=(Esq_avg[i]-(E_avg[i]*E_avg[i]*n))/(T*T); //heat capacity/spin
    
    output(size,lat);
    
    }

    //--------------------------------------------------------------------
    // This uses pgplot to plot T vs Ms, Es, Xs, *Heat Cap, cumulant
    //--------------------------------------------------------------------

    critHCT[k]=t[find_max(nintervals,C)]; 
    printf("Curie Temp=%f\n",critHCT[k]);

    // Open a plot window
    if (!cpgopen("/XWINDOW")) return 1;

    // Heat Cap and Energy vs T curve
    cpgenv(0.8,5.,-2.,2.4,0,1);
    cpglab("T","C","Heat Cap (Green) Energy (Orange) Absolute Magnetization (Blue) vs. Temp"); 
    cpgsci(3);
    cpgpt(nintervals,t,C,7);
    cpgsci(8);
    cpgpt(nintervals,t,E_avg,6); 
    cpgsci(5);
    cpgline(nintervals,t,Mabs_avg);
    
    }

    //Size vs Curie Temp Curve
    cpgenv(5.,100.,1.9,2.3,0,1);
    cpglab("Size","Tc","Size vs. Curie Temp");
    cpgsci(1);
    cpgline(6,dim,critHCT);

    // Pause and then close plot window 
    cpgclos();
}

