#include <stddef.h>
#include <math.h>

const double pi = 3.1415926;
const double onethird = 1.0 / 3.0;
const double twothirds = 2.0 / 3.0;
const double fourthirds = 4.0 / 3.0;
const double onesixth = 1.0 / 6.0;

void spline(double *xv, double *yv, int nr, double yp1, double ypn, double *ypr)
{
	double pp,qn,sig,un;
	int ii,kk;
	double uu[nr];
	
	
	if (yp1 > 0.99e99) ypr[0] = uu[0] = 0.0;
	else {
		ypr[0]	= -0.5;
		uu[0]	= (3.0/(xv[1]-xv[0])) * ((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
	}
	
	for (ii=1; ii<nr-1; ii++){
		sig		= (xv[ii]-xv[ii-1])/(xv[ii+1]-xv[ii-1]);
		pp		= sig*ypr[ii-1]+2.0;
		ypr[ii]	= (sig-1.0)/pp;
		uu[ii]	= (yv[ii+1]-yv[ii])/(xv[ii+1]-xv[ii]) - (yv[ii]-yv[ii-1])/(xv[ii]-xv[ii-1]);
		uu[ii]	= (6.0*uu[ii]/(xv[ii+1]-xv[ii-1])-sig*uu[ii-1])/pp;
	}
	
	if (ypn > 0.99e99) qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0/(xv[nr-1]-xv[nr-2])) * (ypn-(yv[nr-1]-yv[nr-2])/(xv[nr-1]-xv[nr-2]));
	}
	
	ypr[nr-1] = (un-qn*uu[nr-2])/(qn*ypr[nr-2]+1.0);
	
	for (kk=nr-1; kk>=0; kk--)
	{
		ypr[kk] = ypr[kk]*ypr[kk+1]+uu[kk];
	}
}

double splint(double *xx, double *yy, int nr, double z, double *y2a)
{
	int klo,khi,kkk;
	double yo,hh,b,a;
	
	klo = 0;
	khi = nr-1;
	
	do {
		kkk = (khi+klo) >> 1;
		if (xx[kkk] > z) khi = kkk;
		else klo = kkk;
	} while ((khi-klo)>1);
	
	hh = xx[khi]-xx[klo];
	if (hh == 0.0) throw("Bad input to routine splint\n");
	a = (xx[khi]-z)/hh;
	b = (z-xx[klo])/hh;
	
	yo = a*yy[klo] + b*yy[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi])*(hh*hh)*onesixth;
	return yo;		
}

double linterp(double *xx, double *yy, int nr, double z)
{
	int jj;
	double slope,yo;
	if ((xx[0]-xx[1])<0) {
		//normal mode
		for (jj=0;jj<nr;jj++){
			if (xx[jj]>z) {
				slope	= (yy[jj]-yy[jj-1])/(xx[jj]-xx[jj-1]);
				yo		= slope*z + yy[jj] - xx[jj]*slope;
				break;
			}
		}
	}else {
		//reverse mode
		for (jj=0;jj<nr;jj++){
			if (xx[jj]<z) {
				slope	= (yy[jj]-yy[jj-1])/(xx[jj]-xx[jj-1]);
				yo		= slope*z + yy[jj] - xx[jj]*slope;
				break;
			}
		}
	}
	return yo;			
}