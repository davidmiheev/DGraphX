#include "DGraphX.h"
#include "Interp2D.h"
#define m(i,j) m[4*i+j]
#define N_ITER 100
#define min(a,b) (a) < (b) ? (a) : (b)

double shoot(int mode, double par, int n, Point *pt);

/*double f(double x) {
    return sin(x);
}*/

double r(double y, double x) {
    return (f(y)-f(x))/(y-x);
}

void interpH(double* x, int n, double *m) {
    for (int i = 0; i < n-1; i++) {
        m(i,0)=f(x[i]);
        m(i,1)=df(x[i]);
        m(i,2)=(3*r(x[i+1],x[i])-2*df(x[i])-df(x[i+1]))/(x[i+1]-x[i]);
        m(i,3)=(df(x[i+1])+df(x[i])-2*r(x[i+1],x[i]))/(x[i+1]-x[i])/(x[i+1]-x[i]);
    }
}

void interpS(double* x, int n, double *m) {
    double *d,*a,*b,*c;
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
}




void solveEq(int mode, int n, Point *pt) {
    FILE *f;
    double h = 1./(n-1), a, b, c = 5,//exp(h)/5. - exp(sqrt(6)*h)/20. - 3*exp(-sqrt(6)*h)/20.,-(exp(h)/5. - exp(sqrt(6)*h)/20. - 3*exp(-sqrt(6)*h)/20.)
    d = -1, x;
    f = fopen("res.txt","w");
    x = shoot(mode, c, n, pt);
    for(int j = 0; fabs(x) > 1e-5 && j < N_ITER ; j++) {
        a = shoot(mode, c, n, pt);
        b = shoot(mode, d, n, pt);
        if(a*b < 0) {
            x = shoot(mode, c/2.+d/2., n, pt);
            if(x != 0) {
            switch (x < 0) {
                case 1:
                    if(min(a,b) == a) c = c/2.+d/2.;
                    else d = c/2. + d/2.;
                    break;
                case 0:
                    if(min(a,b) == a) d = c/2. + d/2.;
                    else c = c/2. + d/2.;
                    break;
            }
          }
        } else { if(j == 0) printf("%f %f\n", a, b); }
    } printf("%f\n", x);
    //for(int i = 0; i < n - 1; i++) fprintf(f,"(%f,%f)--", pt[i].x, pt[i].y);
    //fprintf(f,"(%f,%f)", pt[n-1].x, pt[n-1].y);
    fclose(f);
}

double shoot(int mode, double par, int n, Point *pt) {
    double h = 1./(n-1), tmp[3], res = 0, a;
    tmp[0] = 0; tmp[1] = par;
    for (int i = 0; i < n - 2; i++) {
        res += h*(tmp[0]*(sin(i*h) + tmp[0]*tmp[0])/2. + tmp[1]*(sin((i+1)*h) + tmp[1]*tmp[1])/2.);
        pt[i].x = ((double) i)/(n-1); pt[i].y = tmp[0];
        a = tmp[1];
        switch (mode) {
            case 1:
                tmp[2] = (6*h*h + 2.)*tmp[1] - 100*sin(20*(i+1)*h)*h*h - tmp[0];
                //printf("%f\n", tmp[2]);
                break;
            case 2:
                tmp[2] = ((6*h*h + 2.)*tmp[1] - exp((i+1)*h/2.)*h*h - tmp[0]);
                //printf("%f\n", tmp[2]);
                break;
            case 3:
                tmp[2] = ((6*h*h + 2.)*tmp[1] - exp((i+1)*h*(i+1)*h)*h*h - tmp[0]);
                //printf("%f\n", tmp[2]);
                break;
        } tmp[0] = a; tmp[1] = tmp[2];
    } res += h*((tmp[0]*(sin((n-2)*h) + tmp[0]*tmp[0]))/2. + (tmp[1]*(sin(1) + tmp[1]*tmp[1]))/2.);
     pt[n-2].x = ((double) (n-2))/(n-1); pt[n-2].y = tmp[0];
    pt[n-1].x = 1; pt[n-1].y = tmp[1];
    return res;
}

/*void solveEq(int mode, int n, Point *pt) {
    double *f, *a, *b, *c, *d, *x, h = 1./(n-1), tmp[3], tmp1;
    a = (double*)malloc((n-1)*sizeof(double));
    b = (double*)malloc(n*sizeof(double));
    c = (double*)malloc((n-1)*sizeof(double));
    d = (double*)malloc(n*sizeof(double));
    f = (double*)malloc((n-1)*sizeof(double));
    x = (double*)malloc(n*sizeof(double));
    switch (mode) {
        case 1:
            for (int i = 1; i < n; i++) x[i] = exp(i*h)/5. - exp(sqrt(6)*i*h)/20. - 3*exp(-sqrt(6)*i*h)/20.;
            break;
        case 2:
            for (int i = 1; i < n; i++) x[i] = 4*exp(i*h/2.)/23. - exp(sqrt(6)*i*h)/23. - 3*exp(-sqrt(6)*i*h)/23.;
            break;
        case 3:
            for (int i = 0; i < n; i++) x[i] = 0;
            break;
    } x[0] = 0; f[n - 2] = 0;
    for(int j = 0; j < N_ITER; j++) {
        for(int i = 1; i < n; i++) f[n - 2] += x[i]*(sin(i*h) + x[i]*x[i]);
        f[n - 2] -= x[n - 1]*(sin(1) + x[n - 1]*x[n - 1])/2.;
        tmp1 = f[n - 2];
        switch (mode) {
            case 1:
                for (int i = 1; i < n - 1; i++) f[i - 1] = (6 + 2./h/h)*x[i] - x[i + 1]/h/h - x[i - 1]/h/h - exp(i*h);
                break;
            case 2:
                for (int i = 1; i < n - 1; i++) f[i - 1] = (6 + 2./h/h)*x[i] - x[i + 1]/h/h - x[i - 1]/h/h - exp(i*h/2.);
                break;
            case 3:
                for (int i = 1; i < n - 1; i++) f[i - 1] = (6 + 2./h/h)*x[i] - x[i + 1]/h/h - x[i - 1]/h/h - exp(i*h*i*h);
                break;
        }
        tmp[0] = (sin(h) + 3*x[1]*x[1]); tmp[1] = (sin(2*h) + 3*x[2]*x[2]); tmp[2] = (sin(3*h) + 3*x[3]*x[3]);
        for (int i = 1; i < n - 2; i++) {
            a[0] = tmp[1] + 6*tmp[0]*h*h + 2*tmp[0];
            a[1] = tmp[2] - tmp[0];
            if(i < n - 3) tmp[2] = (sin((i + 3)*h) + 3*x[i + 3]*x[i + 3]);
            if(i == n - 4) tmp[2] /= 2.;
            f[n - 2] += f[i]*tmp[0]*h*h;
            tmp[0] = a[0]; tmp[1] = a[1];
        }
        b[0] = 1;
        for (int i = 0; i < n-1; i++) {
            if(i < n-2) a[i] = -1./h/h; else a[i] = tmp[0];
            if(i > 0) c[i] = (-1./h/h)/b[i]; else c[i] = 0;
            if(i < n-2) b[i+1] = 6 + 2./h/h - a[i]*c[i]; else b[i+1] = tmp[1] - a[i]*c[i];
        } d[0] = 0;
        for (int i = 1; i < n; i++) d[i] = (f[i - 1] - a[i-1]*d[i-1])/b[i];
        for (int i = n - 2; i >= 0; i--) d[i] = d[i] - c[i]*d[i+1];
        for (int i = 1; i < n; i++) x[i] -= d[i];
    }
    for(int i = 0; i < n; i++) { pt[i].x = ((double) i)/(n-1); pt[i].y = x[i]; }
    //for(int i = 0; i < n - 2; i++) printf("%f ", f[i]);
    printf("%f\n", tmp1*h);
    free(a); free(b); free(c); free(d); free(f); free(x);
}*/
