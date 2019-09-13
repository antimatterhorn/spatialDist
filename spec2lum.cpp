#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <vector>
#include <stddef.h>
#include <math.h>

using namespace std;

char specfile[80];
char rcfile[80];
double rc[1900][2];

double response(double wl)
{
    int i,j;
    for(i=0;i<1900;i++)
        if(rc[i][0] == wl) return rc[i][1];

    return 0.0;

}

int main()
{
	int fi,i,j,k;
    double age,ageold,lambda,flux,dl=0,lum;
	double d1,d2,d3,dummy;

    strcpy(specfile,"spectra/spec.csv");
    strcpy(rcfile,"spectra/rc.csv");
	
	ifstream readspec;
    ifstream readrc;
    
    readrc.open(rcfile);
    while(!readrc.eof()){
        readrc >> d1 >> d2 >> d3;
        rc[fi][0] = d1;
        rc[fi][1] = d3;
        fi++;
    }
    readrc.close();
    
	readspec.open(specfile);
	while(!readspec.eof()){
        //Age(yr),lambda(ang),log Flux,Flux
        readspec >> d1 >> d2 >> d3;
        //printf("%3.2e %3.2e %3.2e\n",d1,d2,d3);
		dl = d2-lambda; 
		age		= d1;
        if(age!=ageold) //new age
        {
            printf("%3.2e %3.2e\n",ageold,(19.5-2.5*log10(lum/(3.85e33))));
            lum = 0.0;
            dl = 1.0;
        }
           
		lambda	= d2;
		flux	= pow(10.0,d3);
        lum += flux*dl*response(lambda);
		
		ageold = age;
	}
	readspec.close();
	
	
	return 0;
}