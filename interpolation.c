#include "DGraphX.h"
#include "Interp2D.h"
#define m(i,j) m[4*i+j]

/*double f(double x) {
    return sin(x);
}*/

double r(double y, double x) {
    return (f(y)-f(x))/(y-x);
}

double* interpH(double* x, int n) {
    double *m = (double*)malloc(4*(n-1)*sizeof(double));
    for (int i = 0; i < n-1; i++) {
        m(i,0)=f(x[i]);
        m(i,1)=df(x[i]);
        m(i,2)=(3*r(x[i+1],x[i])-2*df(x[i])-df(x[i+1]))/(x[i+1]-x[i]);
        m(i,3)=(df(x[i+1])+df(x[i])-2*r(x[i+1],x[i]))/(x[i+1]-x[i])/(x[i+1]-x[i]);
    }
    return m;
}

double* interpS(double* x, int n) {
    double *m,*d,*a,*b,*c;
    m = (double*)malloc(4*(n-1)*sizeof(double));
    d = (double*)malloc(n*sizeof(double));
    a = (double*)malloc((n-1)*sizeof(double));
    b = (double*)malloc(n*sizeof(double));
    c = (double*)malloc((n-1)*sizeof(double));
    b[0] = 1;
    for (int i=0; i < n-1; i++) {
        if(i < n-2) a[i] = x[i+2]-x[i+1]; else a[i] = 0;
        if(i > 0) c[i] = (x[i+1]-x[i])/b[i]; else c[i] = 0;
        if(i < n-2) b[i+1] = 2*(x[i+2]-x[i])-a[i]*c[i]; else b[i+1] = 1 - a[i]*c[i];
    } d[0] = df(x[0]);
    for (int i=1; i < n; i++) {
        if(i < n-1) d[i] = (3*(r(x[i],x[i-1])*(x[i+1]-x[i])+r(x[i+1],x[i])*(x[i]-x[i-1]))-a[i-1]*d[i-1])/b[i];
        else d[i] = (df(x[i]) - a[i-1]*d[i-1])/b[i];
    }
    for (int i = n - 2; i >= 0; i--) d[i] = d[i] - c[i]*d[i+1];
    for (int i = 0; i < n-1; i++) {
        m(i,0)=f(x[i]);
        m(i,1)=d[i];
        m(i,2)=(3*r(x[i+1],x[i])-2*d[i]-d[i+1])/(x[i+1]-x[i]);
        m(i,3)=(d[i+1]+d[i]-2*r(x[i+1],x[i]))/(x[i+1]-x[i])/(x[i+1]-x[i]);
    }
    free(a); free(b); free(c); free(d);
    return m;
}
