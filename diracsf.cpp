#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include "spatial.h"

using namespace std;

//physical variables
double v0,r0,mass1,mass2,radius,rscale,tau,off;

//bookkeeping variables
int pix;
double ipi = 1.0/3.14159;

char csvfile[80];
char chfile[80];
char mstofile[80];
char coeffile[80];

double omegastar(double z)
{
	if (z<1.0)
	{
		return v0/1.0;
	}else {
		return v0/z;
	}
}

double omega(double z)
{
	return v0/8.0 * (1.-1./sqrt(2.)) - omegastar(z);
	//epicyclic frequency is not related to scale length
}

double phi(double z)
{
	return exp(-z/3.0);
	//hard-coding scale radius to 3.0kpc
}

void creategalaxy(double **gal)
{
	int fi,m;
	ifstream file(chfile) ;
	string line ;
	vector<string> lines ;
	while( getline( file, line ) ) lines.push_back( line ) ;
	m = lines.size();
	file.close();
	
	double age[m],band[m];
	
	ifstream readfile;
	readfile.open(chfile);
	fi=0;
	double d1,d2,dummy;
	while(!readfile.eof()){
		//Log Age,Age(yr),Age(Myr),Mbol,Mu,Mb,Mv,Mk,NUV,Mg
		readfile >> dummy >> dummy >> d1 >> dummy >> dummy >> d2 >> dummy >> dummy >> dummy >> dummy;
		age[fi] = d1;
		band[fi] = d2;
		fi++;
	}
	readfile.close();
	
	double *y2;
	y2 = (double *) malloc(m*sizeof(double));
	
	spline(age,band,m-1,0,0,y2);
	
	int imax = pix;
	int jmax = imax;
	int ctr	= pix >> 1;		//bitwise shift right by 1 equiv to div by 2
	double r,theta,L,arg;
	int i,j,k,n;
	
	for (i=0;i<imax;i++)
	{
		for (j=0;j<jmax;j++)
		{
			
			r		= i*radius/imax;
			theta	= j*2*3.14159/jmax;
			L		= 0;
			k		= fabs(ipi * tau * omega(1.0))+1.0;	//need to round off digits here
			

			for (n=-k;n<k+1;n++)
			{
				arg = (off - theta/omega(r)-3.14159*n/omega(r));
				if (arg < 0 || arg > tau) L+=0;
				else {
					L+= pow(10.0,(4.8 - splint(age,band,m,arg,y2))/2.5) * phi(r)/fabs(omega(r));
				}
			}
			
			gal[i][j] = L;
			//printf("%f\n",(double)(i*pix + j)/(double)(pix*pix));
		}
	}
	free(y2);
}

void createsnmap(double **sn)
{
	int fi,m;
	ifstream file(mstofile) ;
	string line ;
	vector<string> lines ;
	while( getline( file, line ) ) lines.push_back( line ) ;
	m = lines.size();
	file.close();
	
	double logage[m],mass[m],coef[5];
	
	ifstream readfile;
	readfile.open(mstofile);
	fi=0;
	double d1,d2,dummy;
	while(!readfile.eof()){
		//Log Age	Age(yr)	Age(Myr)	Mass	LogMass
		readfile >> d1 >> dummy >> dummy >> d2 >> dummy;
		logage[fi] = d1;
		mass[fi] = d2;
		fi++;
	}
	readfile.close();

	readfile.open(coeffile);
	fi=0;
	while(!readfile.eof()){
		readfile >> coef[fi];
		fi++;
	}	
	readfile.close();

		
	
	int imax = pix;
	int jmax = imax;
	int ctr	= pix >> 1;
	double r,theta,dth,L,arg,arg1,arg2,tauprime,m1,m2;
	int i,j,k,n;
	
	for(i=0;i<5;i++)
	{
		printf("%f\n",coef[i]);
	}
	
	tauprime = pow(10.0,linterp(mass,logage,m,mass1))*1e-6;
	
	for (i=0;i<imax;i++)
	{
		for (j=0;j<jmax;j++)
		{	

			
			r		= i*radius/imax;
			theta	= j*2*3.14159/jmax;
			dth		= 2*3.14159/jmax;
			L		= 0;
			k		= fabs(ipi * tauprime * omega(1.0))+1.0;	
			
			
			

			for (n=-k;n<k+1;n++)
			{
				arg = (off - theta/omega(r)-3.14159*n/omega(r));
				arg1 = (off - (theta-dth/2.0)/omega(r)-3.14159*n/omega(r));
				arg2 = (off - (theta+dth/2.0)/omega(r)-3.14159*n/omega(r));
				
				arg1 = log10((arg1/1.1)*1e6);
				arg2 = log10((arg2/1.1)*1e6);
				if (arg < 3 || arg > 15000) L+=0;
				else {
					m1 = pow(10.0,coef[0]+coef[1]*arg1+coef[2]*arg1*arg1+coef[3]*arg1*arg1*arg1+coef[4]*arg1*arg1*arg1*arg1);
					m2 = pow(10.0,coef[0]+coef[1]*arg2+coef[2]*arg2*arg2+coef[3]*arg2*arg2*arg2+coef[4]*arg2*arg2*arg2*arg2);
					
					if (m1>mass2) m1 = mass2;
					if (m2>mass2) m2 = mass2;
					if (m1<mass1) m1 = mass1;
					if (m2<mass1) m2 = mass1;
					
					L+= 1.1/1.3 * (pow(m2,-1.3)-pow(m1,-1.3))*phi(r)/dth;
				}				
			}
	
			

			if (L>0) sn[i][j] = L;
			else sn[i][j] = 0;
			//printf("%f\n",(double)(i*pix + j)/(double)(pix*pix));
		}
	}
}

void fruchter(double **gal, double **sn, double **fr)
{
	int i,j,n,flag;
	double totL,totS,fracL,fracS,tempr,tempt,tempv,temps;
	totL=totS=fracL=fracS=0;
	
	for (i=0;i<pix;i++)
	{
		for (j=0;j<pix;j++)
		{
			n = i*pix + j;
			fr[n][0] = i;
			fr[n][1] = j;
			fr[n][2] = gal[i][j];
			fr[n][3] = sn[i][j];
			
			//totL += gal[i][j]*(double)i*(double)radius/(double)pix;
			//totS += sn[i][j]*(double)i*(double)radius/(double)pix;
		}
	}
		
	cout << "Now Sorting...\n";
	
	n = pix*pix;
	flag = 1;
	while( flag || (n > 1))       //boolean flag (true when not equal to 0)
	{
		flag = 0;            //reset flag to 0 to check for future swaps
		n = (n+1) / 2;
		for (i = 0; i < (pix*pix - n); i++)
        {
			if (fr[i + n][2] < fr[i][2])	//sort ascending with < and descending with >
			{
				tempr = fr[i + n][0];       //swap positions i+d and i
				tempt = fr[i + n][1];
				tempv = fr[i + n][2];
				temps = fr[i + n][3];
				fr[i + n][0] = fr[i][0];
				fr[i + n][1] = fr[i][1];
				fr[i + n][2] = fr[i][2];
				fr[i + n][3] = fr[i][3];
				fr[i][0] = tempr;
				fr[i][1] = tempt;
				fr[i][2] = tempv;
				fr[i][3] = temps;
				flag = 1;                   //tells swap has occurred
			}
		}
	}
	
	for (i=0;i<(pix*pix-1);i++) {
		totL += fr[i][2]*fr[i][0]*(double)radius/(double)pix;
		totS += fr[i][3]*fr[i][0]*(double)radius/(double)pix;
		fr[i][4] = totL;
		fr[i][5] = totS;
	}
	
	
	strcpy(csvfile,"fruchter.csv");
	ofstream stream2(csvfile);
	
	stream2 << "i,j,v,s,fracL,fracS" << endl;
	
	for (i=0;i<(pix*pix-1);i++) {

		fr[i][4] = fr[i][4]/totL;
		fr[i][5] = fr[i][5]/totS;
		stream2 << fr[i][0] << "," << fr[i][1] << "," << fr[i][2] << "," << fr[i][3] << "," << fr[i][4] << "," << fr[i][5] << "," << endl;
	}
	stream2.close();
}

//int polar2xy(double **mat)
//{
//	int i,j,x,y;
//	double **mat2;
//	mat2 = (double **) malloc((2*pix+1)*sizeof(double *));
//	for (j=0;j<(2*pix+1);j++) mat2[j] = (double *) malloc((2*pix+1)*sizeof(double));
//	
//	
//	for(y=-pix+1;y<pix;y++)
//	{
//		for(x=-pix+1;x<pix;x++)
//		{
//			if(x!=0 && sqrt(x*x+y*y)<pix)
//			{
//				i = sqrt(x*x+y*y);
//				j = atan2(y,x);
//				if(j<0) j+=2*3.14159;
//				
//				mat2[y+pix][x+pix] += mat[i][j];
//			}
//			else if(x==0)
//			{
//				i = fabs(y);
//				j = pix/2;
//				mat2[y+pix][x+pix] += mat[i][j];
//			}
//		}
//	}
//	
//	for(i=0;i<pix;i++)
//	{
//		for(j=0;j<pix;j++)
//		{
//			mat[i][j] = mat2[i+pix/2][j+pix/2];
//		}
//	}
//	free(mat2);
//}

int polar2xy(double **mat)
{
	int i,j,x,y;
	double **mat2;
	mat2 = (double **) malloc((2*pix+1)*sizeof(double *));
	for (j=0;j<(2*pix+1);j++) mat2[j] = (double *) malloc((2*pix+1)*sizeof(double));
	printf("created x-y table\n");
	
	for(y=-pix+1;y<pix;y++)
	{
		for(x=-pix+1;x<pix;x++)
		{
			if(x!=0 && sqrt(x*x+y*y)<pix)
			{
				i = sqrt(x*x+y*y);
				j = (atan(y/x)+3.14159/2.0)*pix/(2.*3.14159);
				if(j<0) j=(atan(y/x)+3*3.14159/2.0)*pix/(2.*3.14159);
				
				mat2[y+pix][x+pix] += mat[i][j];
			}
			else if(x==0)
			{
				i = fabs(y);
				j = pix/2;
				mat2[y+pix][x+pix] += mat[i][j];
			}
		}
	}
	printf("table is filled\n");
	
	for(i=0;i<pix;i++)
	{
		for(j=0;j<pix;j++)
		{
			mat[i][j] = mat2[i+pix/2][j+pix/2];
		}
	}
	free(mat2);
}

int main()
{
	int j,row,col;
	strcpy(chfile,"include/ch02.dat");
	strcpy(mstofile,"include/msto02.dat");
	strcpy(coeffile,"include/coefs02.dat");
	
//	cout << "tau(Myr): ";
//	cin >> tau;
//	cout << "pixels: ";
//	cin >> pix;
//	cout << "radius(kpc): [should be about 10] ";
//	cin >> radius;
//	cout << "r0(kpc): [this is hard-coded right now] ";
//	cin >> r0;
//	cout << "v0(kpc/Myr): ";
//	cin >> v0;
//	cout << "offset(Myr): ";
//	cin >> off;
	
	tau = 20000.0;
	pix = 500.0;
	radius = 20.0;
	r0 = 8.0;
	v0 = 0.205;
	
	cout << "mass1(msun): ";
	cin >> mass1;
	cout << "mass2(msun): ";
	cin >> mass2;

	rscale = 2.0 * (radius/pix);
	
	double **galaxy;
	galaxy = (double **) malloc(pix*sizeof(double *));
	for (j=0;j<pix;j++) galaxy[j] = (double *) malloc(pix*sizeof(double));
	cout << "Creating galaxy... (" << pix << "x" << pix << ")" << endl;
	creategalaxy(galaxy);
	
		
	double **sne;
	sne = (double **) malloc(pix*sizeof(double *));
	for (j=0;j<pix;j++) sne[j] = (double *) malloc(pix*sizeof(double));
	cout << "Creating SN map... (" << pix << "x" << pix << "," << mass1 << "-" << mass2 << ")" << endl;
	createsnmap(sne);
	

	
	double **fruch;
	fruch = (double **) malloc(pix*pix*sizeof(double *));
	for (j=0;j<pix*pix;j++) fruch[j] = (double *) malloc(6*sizeof(double));
	
	fruchter(galaxy,sne,fruch);
	
	strcpy(csvfile,"gal.csv");
	ofstream stream(csvfile);
	//polar2xy(galaxy);
	
	for(row = 0; row < pix; row++)
	{
		for(col = 0; col < pix; col++)
		{
			stream << galaxy[row][col] << " ";	
		}
		stream << endl;
	}
	stream.close();
	
	strcpy(csvfile,"sne.csv");
	ofstream stream2(csvfile);
	//polar2xy(sne);
	
	for(row = 0; row < pix; row++)
	{
		for(col = 0; col < pix; col++)
		{
			stream2 << sne[row][col] << " ";	
		}
		stream2 << endl;
	}
	stream2.close();
	
	
	return 0;
}


