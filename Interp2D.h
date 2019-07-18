#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double* interpH(double* x, int n);
double* interpS(double* x, int n);
double f(double x);
double df(double x);
double r(double y, double x);
void dtoa(double x, char s[]);
void itoa(int n, char s[]);
void scat(char* s, char* str, unsigned long len,char* msg);
