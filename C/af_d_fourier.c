#include <stdio.h>
#include <math.h>
#include "af_fourier.c"
#define MaxProduct 50

/*	fourier-calculated atomic functions (derivatives)		*/
/*	version 0.2							*/
/*				*/
/*	structure: 					*/
/*	f_XXXX_prime(double x, ..., int k) - prime of atomic function "XXXX"	*/
/*	F_XXXX_prime(double t, ..., int k) - image of prime of "XXXX"*/
/*	o_XXXX(...) - any other function, used further	*/
/*							*/

double f_up_prime(double x, int k)
{
	if(fabs(x)>1.)
		return 0.;
	if(k==1)
	{
		if(x>0.)
			return -2.*f_up(2.*x-1.);
		else
			return 2.*f_up(2.*x+1.);
	}
	else
	{
		if(x>0.)
			return -2.*f_up_prime(2.*x-1.,k-1);
		else
			return 2.*f_up_prime(2.*x+1.,k-1);
	}
}

double f_ha_prime(double x, double a, int k)
{
	if(fabs(x)>1./(a-1.))
		return 0.;
	if(k==1)
	{
		if(x>0.)
			return -a*a*0.5*f_ha(a*x-1.,a);
		else
			return a*a*0.5*f_ha(a*x+1.,a);
	}
	else
	{
		if(x>0.)
			return -a*a*0.5*f_ha_prime(a*x-1.,a,k-1);
		else
			return a*a*0.5*f_ha_prime(a*x+1.,a,k-1);
	}
}

int o_binCoeff(int k, int n)
{
    int num, den, i, j ;
    if ( n < k ) 
    {
       return(0) ; 
    }
    else 
    {
	den = 1;
	num = 1 ; 
	for (i =  1  ; i <= k   ; i = i+1)
	    den =    den * i;
	for (j = n-k+1; j<=n; j=j+1)	
	    num = num * j;
	return(num/den);
    } 
}

double f_fup_prime(double x, int n, int k)
{
	int i;
	double y,z=0.;
	if(k==1)
	{
		for(i=0;i<n+2;i++)
		{
			z+=(double)(o_binCoeff(i,n+1)-o_binCoeff(i-1,n+1))*
				f_fup(2.*x+(double)n/2.+1.-(double)i,n);
		}
		return pow(2.,-n+1)*z;
	}
	if(k>1)
	{
		for(i=0;i<n+2;i++)
		{
			z+=(double)(o_binCoeff(i,n+1)-o_binCoeff(i-1,n+1))*2.*
				f_fup_prime(2.*x+(double)n/2.+1.-(double)i,n,k-1);
		}
		return pow(2.,-n+1)*z;
	}
}

double f_cup_prime(double x, int k)
{
	if(fabs(x)>2.)
		return 0.;
	double 	pik,
		y = 0.0;
		
	for(pik = M_PI; pik<=M_PI*MaxProduct; pik+=M_PI)
	{
		y+=pow(F_up(pik*0.5),2)*pow(pik*0.5,k)*cos(pik*x*0.5+M_PI*0.5*k);
	}
	return y*0.5;
}

double f_chan_prime(double x, double a, double n, int k)
{
	if(fabs(x)>(double)n/(double)(a-1))
		return 0.;
	double 	pik,
		y = 0.0,
		h = M_PI*(double)(a-1)/(double)n;
		
	for(pik = h; pik<=h*MaxProduct; pik+=h)
	{
		y+=pow(F_ha(pik,a),n)*pow(pik,k)*cos(pik*x+0.5*M_PI*k);
		//y+=pow(F_ha(pik,a),n)*pow(pik,k)*cos(pik*x*0.5+M_PI*0.5*k);
	}

	return y*(a-1.)/(double)n;
}

/*
int main()
{
	FILE *graph;
	double i;
	graph = fopen("fup.dat", "w");
	for(i=-2.; i<=2.; i+=0.01)
		fprintf(graph,"%f %f\n", i, f_chan_prime(i,2.,2.,1));
	fclose(graph);
	return 0;
}
*/
