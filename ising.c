#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include "ising.h"
#include <stdbool.h>
int size;               //lattice size
int n;
double T;               //starting point for temperature
const double minT;      //minimum temperature
double change;          //size of steps for temperature loop
long unsigned int mcs;  //number of Monte Carlo steps
int transient;          //number of transient steps
double norm;            //normalizaton for averaging

//function for random initialization of lattice 
void initialize(int size, int lat[size+1][size+1])
{
    int x,y;
    for(y=size;y>0;y--)
        {
            for(x=1;x<size+1;x++) 
                {
                    if(((1.0*(double)rand())/RAND_MAX) >=0.5)
                        lat[x][y]=1;
                    else
                        lat[x][y]=-1;
                }
        }
}

//output of lattice configuration to the screen 
void output (int size, int lat[size+1][size+1])
{
    for(int y=size;y>=1;y--)
        {
           for(int x=1;x<=size;x++) 
            {
                if(lat[x][y]<0)
                    printf("-");
                else
                    printf("+");
            }
          printf("\n");      //newline after each row
        }
}

//function for choosing random position on lattice 
void choose_random_pos_lat (lat_type *pos)
{
    pos->x=(int)(ceil((size-1)*(double)rand()/RAND_MAX+1)); 
    pos->y=(int)(ceil((size-1)*(double)rand()/RAND_MAX+1)); 
    if(pos->x>size||pos->y>size)
    {
        printf("error in array size!\n");
    exit(1);
    }
}

//function for calculating energy at a particular position on lattice
int energy_pos(lat_type pos, int size, int lat[size+1][size+1])
{
    //point feels effects of thing directly above, below, right & left of it 
    int up, down, left, right, e;
    if(pos.y==size)
        up=1; 
    else
        up=pos.y+1;
    
    if(pos.y==1)
        down=size;
    else
        down=pos.y-1;

    if(pos.x==1)
        left=size;
    else
        left=pos.x-1;
    
    if(pos.x==size)
        right=1;
    else
        right-pos.x+1;

    //energy for specific position.J = +1
    e=lat[pos.x][pos.y]*(lat[left][pos.y]+lat[right][pos.y]+lat[pos.x][up]+lat[pos.x][down]); 
    e*=-1;
    return e;
}

//function for testing the validity of flipping a spin at a selected position 
bool test_flip(lat_type pos, int *de, int size, int lat[size+1][size+1])
{
    *de=2*energy_pos(pos, size, lat);
    *de*=-1;
    if(*de<0)
        return true;        //flip due to lower energy
    else if(((1.0*(double)rand())/RAND_MAX)<(exp((float)(-1.*(*de))/T)))   //Actually is T*kb
        return true;        // flip due to heat bath
    else
        return false;       //no flip
}

//flip spin at given position
void flip(lat_type pos, int size, int lat [size+1][size+1])
{
    lat[pos.x][pos.y]=-lat[pos.x][pos.y];
}

//function for disreagrding transient results
void transient_results(int size, int lat[size+1][size+1])
{
    lat_type pos;
    int de=0;
    for(int a=1;a<=transient;a++)
        {
            for(int b=1;b<=n;b++)
                {
                    choose_random_pos_lat(&pos); 
                    if(test_flip(pos,&de,size,lat))
                        {
                            flip(pos, size, lat);
                        }
                }
        }
}

//function for calculating total magnetization of lattice 
int total_magnetization(int size, int lat[size+1][size+1])
{
    int m=0;
    for(int y=size;y>=1;y--)
        {
            for(int x=1; x<=size;x++)
                {
                    m=m+lat[x][y];      //b/c total m is #up-#down spins
                }
        }
    return m;
}

//function for calculating total energy of lattice 
int total energy (int size, int lat[size+1][size+1])
{
    lat_type pos;
    int ee=0;
    for(int y=size; y>0;y--)
        {
            pos.y=y;
            for(int x=1;x<size+1;x++)
                {
                    pos.x=x;
                    ee+=energy_pos(pos,size,lat);
                }
        }
    return ee;
}

int find_max(int nintervals, float C[nintervals+1])
{
    int c,index;
    float max;
    max=C[2];
    index=0;

  for(c=3;c<nintervals+1;c++)
    {
        if(C[c]>max)
            {
                index=c;
                max=C[c];
            }
    }
    return index;
}
