#include "DGraphX.h"
#include "Interp2D.h"
#define m(i,j) m[4*i+j]
#define N_ITER 100
#define PI 3.14159265358
#define min(a,b) (a) < (b) ? (a) : (b)

static int permut[5] = {5, 2, 4, 1, 3};
static int tmp[5];
static double coeff;

double shoot(int mode, double par, int n, Point *pt);


double func(double x, double y) {
    return 20*sin(20*x*y);//10;//(coeff*8*PI*PI + 1 + sin(x*x))*sin(2*PI*x)*sin(2*PI*y)/2;//20*sin(20*x*y);
}

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

void SetParameter(double eps) {
    coeff = eps;
}

double diff(double * a, double * b, int i, int j, int n) {
    double res, h = 1./(n-1);
    res = b[(i-1)*(n-2)+j-1] - (4 + h*h*(1+sin((i%(n-1))*(i%(n-1))*h*h))/coeff)*a[i*n+j] + a[i*n+j-1] + a[i*n+j+1] + a[(i-1)*n + j] + a[(i+1)*n + j];
    return  res;
}

int err(double * a, double * b, double eps, int n) {
    double res = 0.;
    for(int i = 1; i < n - 1; i++) {
        for(int j = 1; j < n - 1; j++) {
            res += diff(a,b,i,j,n)*diff(a,b,i,j,n)/(n-1);
        }
    } //printf("%f\n", res);
    return (res > eps);
}

double accuracy(double (*sol)(int, int, int), double* a, int n) {
    double res = fabs(sol(0, 0, n) - a[0]);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(res < fabs(sol(i, j, n) - a[i*n +j])) res = fabs(sol(i, j, n) - a[i*n + j]);
        }
    }
    return res;
}

void permutation(int *a) {
    for (int i = 0; i < 5; i++) tmp[i] = a[i];
    for (int i = 0; i < 5; i++) a[i] = tmp[permut[i]-1];
}

void RichardsonIter(double *a, int n, double min, double max, int *niter ) {
    double h = 1./(n-1); int k = 0, d[5] = {1, 9, 3, 7, 5};//{3, 4, 5, 1, 2};;
    double* b = (double*) malloc((n-2)*(n-2)*sizeof(double));
    double c[5];
    //permutation(d);
    for (int i = 0; i < 5; i++) {
        c[i] = 2./(max + min + (max-min)*cos(PI*(d[i])/10.));
    }
    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < n - 2; j++) {
            b[k] = h*h*func((j+1)*h,(i+1)*h)/coeff;
            k++;
        }
    } k = 0;
    while(err(a,b, 1e-10, n) && k < 1000) {
        for(int l = 0; l < 5; l++) {
            for(int i = 1; i < n - 1; i++) {
                for(int j = 1; j < n - 1; j++) {
                    a[i*n + j] += c[l%5]*diff(a,b,i,j,n);
                }
            }
        } k++;
    } *niter = k*5;//printf("%d\n", k);
    /*for (int i = 0; i < n*n; i++) {
        printf("%f\n", a[i]);
    }*/
    free(b);
}

void SeidelGauss(double *a, int n, int *niter) {
    double h = 1./(n-1); int k = 0;
    double* b = (double*) malloc((n-2)*(n-2)*sizeof(double));
    double* d = (double*) malloc((n-2)*(n-2)*sizeof(double));
    for (int i = 0; i < n - 2; i++) {
        for (int j = 0; j < n - 2; j++) {
            b[k] = h*h*func((j+1)*h,(i+1)*h)/coeff;
            k++;
        }
    } k = 0;
    while(err(a,b,1e-10,n) && k < 10000) { k++;
        for(int i = 1; i < n - 1; i++) {
            for(int j = 1; j < n - 1; j++) {
                d[(i-1)*(n-2)+j-1] = diff(a,b,i,j,n);
            }
        }
        for(int i = 0; i < (n - 2)*(n - 2); i++) {
            if(i < n - 2) {
                if(i == 0) d[i] /= (4 + (1+sin(h*h))*h*h/coeff);
                else d[i] = (d[i]+d[i-1])/(4 + (1 + sin((i%(n-2) + 1)*(i%(n-2) + 1)*h*h))*h*h/coeff);
            }
            else { d[i] = (d[i] + d[i-1] + d[i - (n-2)])/(4 + (1 + sin((i%(n-2) + 1)*(i%(n-2) + 1)*h*h))*h*h/coeff); }
        }
        for (int i = 1; i < (n-1); i++) {
            for(int j = 1; j < n-1; j++) {
                a[i*n + j] += d[(i-1)*(n-2)+j-1];
            }
        }
    } *niter = k;//printf("%d\n", k);
    /*for (int i = 0; i < n*n; i++) {
        printf("%f\n", a[i]);
    }*/
    free(b); free(d);
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
