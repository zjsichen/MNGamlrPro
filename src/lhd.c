// likelihood functionals

#include <stdlib.h>
#include <stdio.h>
#include <Rmath.h>
#include "lhd.h"
#include "vec.h"


double intercept(int n, double *e,
                 double *v, double *z, double vsum,double *tocounts){
  int i;

  double rsum = 0.0;
  for(i=0; i<n; i++)
    rsum += v[i]*(z[i]-e[i]-log(tocounts[i]));
  return rsum/vsum;
}

// Negative Log LHD
double sse(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts){
  double l = 0.0;
  double r;
  double vi = 1.0;
  int i;
  for(i=0; i<n; i++){
    if(v[0]!=0) vi = v[i];
    r = (y[i] - a - e[i]);
    l += 0.5*r*r*vi;
  }
  return l;
}

double bin_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts){
  double l = 0.0;
  double f;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    l += -y[i]*f + log(1 + exp(f));
  }
  return l;
}

double po_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts){
  double l = 0.0;
  double f;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    l += exp(f) - y[i]*f;
  }
  return l;
}

double nb_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts){
  double l = 0.0;
  double f;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    l += (1/alphanb)*log(1+alphanb*exp(f))-y[i]*log(alphanb*exp(f)/(1+alphanb*exp(f)));
  }
  return l;
}
// Re-weightings
double bin_reweight(int n, double a, double alphanb, double *e,
                    double *y, double *v, double *z, int *vzf, double *tocounts){
  double q, ee, vs, f;
  vs = 0.0;
  if(!isfinite(a)) return 0.0;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    ee = exp(f);
    q = ee/(1.0+ee);
    v[i] = 1.0/(2.0 + 1.0/ee + ee);
    if(v[i]<1e-12){ // perfect separation
      v[i] = 0.0;
      z[i] = y[i];
      *vzf = 1;
    }else{
      z[i] = f + (y[i]-q)/v[i];
      vs += v[i];
    }
  }
  return vs;
}

double po_reweight(int n, double a, double alphanb, double *e,
                   double *y, double *v, double *z, int *vzf, double *tocounts){
  double vs, f;
  vs = 0.0;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    v[i] = exp(f);
    if(v[i]<1e-12){ // perfect fit
      v[i] = 0.0;
      z[i] = y[i];
      *vzf = 1;
    }else{
      vs += v[i];
      z[i] = f + y[i]/v[i] - 1.0;
    }
  }
  return vs;
}


double nb_reweight(int n, double a, double alphanb, double *e,
                   double *y, double *v, double *z, int *vzf, double *tocounts){
  double vs, f;
  vs = 0.0;
  int i;
  for(i=0; i<n; i++){
    f = a + e[i]+log(tocounts[i]);
    //v[i] = -exp(f)*(1+alphanb*y[i])/pow((1+alphanb*exp(f)),2);
    v[i] = pow(exp(f),2)/(exp(f)+alphanb*pow(exp(f),2));
    if(v[i]<1e-12){ // perfect fit
      v[i] = 0.0;
      z[i] = y[i];
      *vzf = 1;
    }else{
      vs += v[i];
      z[i] = f + y[i]/exp(f) - 1.0;
    }
  }
  return vs;
}
