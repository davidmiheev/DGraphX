#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "DGraphX.h"

void interpH(double* x, int n, double *m);
void interpS(double* x, int n, double *m);
double f(double x);
double df(double x);
double r(double y, double x);
void dtoa(double x, char s[], int ndeg);
void itoa(int n, char s[]);
void scat(const char* s, char* str, unsigned long len,char* msg);
void solveEq(int mode, int n, Point *pt);
double accuracy(double (*sol)(int, int, int), double* a, int n);
void RichardsonIter(double *a, int n, double min, double max, int *niter);
void SeidelGauss(double *a, int n, int *niter);
void SetParameter(double par);
