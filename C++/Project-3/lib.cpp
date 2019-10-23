#include "lib.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <time.h>
using namespace std;



/*
** The function
**              gaussLegendre()
** takes the lower and upper limits of integration x1, x2, calculates
** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
** of length n of the Gauss--Legendre n--point quadrature formulae.
*/
void gaussLegendre(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/(j);
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > pow(10.,-6.));

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


//  Function to set up Gauss-Laguerre integration points and weights given
//  as examples https://github.com/CompPhysics/ComputationalPhysics/blob/master/doc/Projects/2019/Project3/CodeExamples/gauss-laguerre.cpp

#define EPS 3.0e-14
#define MAXIT 10

double gammln(double);

//  Note that you need to call it with a given value of alpha,
// called alf here. This comes from x^{alpha} exp(-x)

void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaussian-Laguerre" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
#undef EPS
#undef MAXIT


/* This function is from lib.h belonging to prof. Morten Hjorth-Jensen and found through
 * available github repository.
double ran0(long *idum)
    ** is an "Minimal" random number generator of Park and Miller
    ** (see Numerical recipe page 279). Set or reset the input value
    ** idum to any integer value (except the unlikely value MASK)
    ** to initialize the sequence; idum must not be altered between
    ** calls for sucessive deviates in a sequence.
    ** The function returns a uniform deviate between 0.0 and 1.0.
    */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
{
   long     k;
   double   ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
double int_MC_brute(double *x){
    double alpha = 2.0;
    double cons = -2*alpha;
    double exp1 = (cons*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
    double exp2 = (cons*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]));
    double denom = sqrt(pow((x[0]-x[3]),2) +pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2));

    if(denom < pow(10.,-6.)){ // cheking if denominator is zero
        return 0;
    }else{
        return exp(exp1+exp2)/denom;
    }
}
void  MCBrute(int dim, double a, double b, int samples, ofstream &file){
    clock_t start, end;

    start = clock();
    double x[dim];
    long idum = -1;
    double int_mc = 0.;
    double variance = 0.;
    double sum_sigma = 0.;

    double fx;
    // evaluating integral
    for(int i = 0;i<=samples;i++){
        //creates points to integrate over
        for(int j = 0;j<dim;j++){
            x[j] = a + (b-a)*ran0(&idum);  // z = a +(b-a)x -> x [0,1] p. 273
        }
        fx = int_MC_brute(x);
        int_mc +=fx;
        sum_sigma += fx*fx;
    }
    double jacobidet = pow((b-a),dim);
    int_mc = jacobidet*int_mc/((double) samples);
    sum_sigma = jacobidet*sum_sigma/((double) samples);
    variance = fabs(sum_sigma-int_mc*int_mc);
    end = clock();
    double time = ((double)(end-start)/CLOCKS_PER_SEC);
    cout<<"Brute forcce Monte Carlo = "<<int_mc<<", N = "<<samples<<endl;
    file.open("/Users/andreas/Computational Physics/Fys4150/Project 3.2/C++/Project-3/MCBruteData.txt", ios::app);//std::ios_base::app
    file<<setw(5)<<setprecision(6)<<int_mc<<setw(20)<<setprecision(4)<<variance<<setw(15)<<samples<<setprecision(3)<<setw(10)<<time<<endl;
    file.close();

}

double int_MC_spherical(double *x){
    double alpha = 2;
    double cosbeta = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
    double r12 = sqrt(x[0]*x[0] +x[1]*x[1] - 2*x[0]*x[1]*cosbeta);
    if(r12 < pow(10.0, -6)||isnan(r12)){ //checking if demoniator is zero
        return 0;
    }else{
        return sin(x[2])*sin(x[3])*x[0]*x[0]*x[1]*x[1]/r12;  //1.0/(32*pow(alpha,5))*
    }


}
void MonteCarlo(int n, ofstream &file){
    clock_t start, end;
    start = clock();
    double int_mc = 0.;
    double variance = 0.;
    double sum_sigma = 0.;
    double fx;
    long idum = -1;
    double x[6];
    double pi = 4*atan(1);
    double jacobidet = 4*pow(pi,4.)*1.0/16.0;// Volume contains 4 jacobideterminants(pi,pi,2pi,2pi) and a scaling factor (1/16)
    // evaluating integral
    for(int i = 1;i<n;i++){

        for(int j = 0;j<2;j++){         //radial points
            double y = ran0(&idum);
            x[j] = -0.25*log(1.-y);     //p.367
        }
        for(int j = 2;j<4;j++){         //angular (theta)
            x[j] = pi*ran0(&idum);
        }
        for(int j = 4;j<6;j++){        //angular (phi)
            x[j] = 2*pi*ran0(&idum);
        }

        fx = int_MC_spherical(x);
        int_mc +=fx;
        sum_sigma +=fx*fx;
    }

    int_mc = jacobidet*int_mc/n;
    sum_sigma = jacobidet*sum_sigma/n;
    variance = fabs(sum_sigma-int_mc*int_mc);
    end = clock();
    double time = ((double)(end-start)/CLOCKS_PER_SEC);
    cout<<"Monte Carlo i.s = "<<int_mc<<", N = "<<n<<endl;
    file.open("/Users/andreas/Computational Physics/Fys4150/Project 3.2/C++/Project-3/MCData.txt", ios::app);//std::ios_base::app
    file<<setw(5)<<setprecision(6)<<int_mc<<setw(20)<<setprecision(4)<<variance<<setw(15)<<setprecision(3)<<setw(10)<<time<<endl;
    file.close();
}
// first try at setting up integral function to evaluate but gives data that are slightly off
// dont know whats causing this but have rewritten the function ti untegrate in integrandSPherical2
double integrandSPherical(double r1, double theta1, double phi1,double r2, double theta2, double phi2){
    double alpha = 2;
    double cosbeta = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double r12 = sqrt(r1*r1 +r2*r2 - 2*r1*r2*cosbeta);
    double exponential = -(2*alpha-1)*(r1+r2);
    if(r12 < pow(10.0, -6)||isnan(r12)){ //checking if demoniator is zero
        return 0;
    }else{
        //alternative return
        //return sin(theta1)*sin(theta2)/r12;
        return exp(exponential)*sin(theta1)*sin(theta2)/r12;
    }
}
//integrand for laguerre to evaluate. This gives right values.
double integrandSPherical2(double r1, double theta1, double phi1,double r2, double theta2, double phi2){
    double alpha = 2;
    double cosbeta = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double r12 = sqrt(r1*r1 +r2*r2 - 2*r1*r2*cosbeta);
    double exponential = -(2*alpha-1)*(r1+r2);
    if(r12 < pow(10.0, -6)||isnan(r12)){ //checking if demoniator is zero
        return 0;
    }else{
        return 1.0/(32*pow(alpha,5))*sin(theta1)*sin(theta2)/r12;
    }
}
