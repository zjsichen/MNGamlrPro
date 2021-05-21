// likelihood functionals

#ifndef __LHDMOD_H__
#define __LHDMOD_H__

double intercept(int n, double *e, double *v, double *z, double vsum, double *tocounts);
double sse(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts);
double bin_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts);
double po_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts);
double nb_nllhd(int n, double a, double alphanb, double *e, double *y, double *v, double *tocounts);

double bin_reweight(int n, double a, double alphanb, double *e,
                    double *y, double *v, double *z, int *vzf, double *tocounts);
double po_reweight(int n, double a, double alphanb, double *e,
                   double *y, double *v, double *z, int *vzf, double *tocounts);
double nb_reweight(int n, double a, double alphanb, double *e,
                   double *y, double *v, double *z, int *vzf, double *tocounts);

#endif
