#include "DGraphX.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

#define a(ix,jx,nx)   A[(ix)*nx+jx]
#define b(ix,jx,nx)   B[(ix)*nx+jx]

#define m(i,j) V[3*i+j]
#define n(i,j) N[3*i+j]
#define t(i,j) T[3*i+j]

#define MAX_OUTPUT_SIZE 5

static double tmp[3];
static double t;

void PrintMatrix(int n, double* A) {
    int i, j, nPrint;
    
    nPrint = (n > MAX_OUTPUT_SIZE) ? MAX_OUTPUT_SIZE : n;
    
    for (i = 0; i < nPrint; ++i)
    {
        printf("\t[ ");
        for (j = 0; j < nPrint; ++j) {
            if(fabs(a(i,j,n)) > 1e-10)
                printf("%10.3g ", a(i,j,n));
            else printf("%10.3g ", 0.);
        }
        printf("\t]\n");
    }
}

int InvMatrix(int n, double* A, double* B) {
    int i;
    int j;
    int k;
    double tmp1,tmp2;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            b(i, j, n) = (double)(i == j);
    
    for (i = 0; i < n; ++i) {
        tmp1 = 0.0;
        for (j = i + 1; j < n; j++)
            tmp1 += a(j, i,n) * a(j, i,n);
        
        tmp2 = sqrt(tmp1 + a(i ,i,n) * a(i ,i,n));
        
        if (tmp2 < 1e-100)
            return -1;
        
        a(i ,i,n) -= tmp2;
        
        tmp1 = sqrt(tmp1 + a(i ,i,n) * a(i ,i,n));
        
        if (tmp1 < 1e-100) {
            a(i ,i,n) += tmp2;
            continue;
        }
        
        tmp1 = 1.0 / tmp1;
        for (j = i; j < n; ++j)
            a(j ,i,n) *= tmp1;
        
        for (k = i + 1; k < n; ++k) {
            tmp1 = 0.0;
            for (j = i; j < n; ++j)
                tmp1 += a(j ,i,n) * a(j ,k,n);
            
            tmp1 *= 2.0;
            for (j = i; j < n; ++j)
                a(j ,k,n) -= tmp1 * a(j ,i,n);
        }
        
        for (k = 0; k < n; ++k) {
            tmp1 = 0.0;
            for (j = i; j < n; ++j)
                tmp1 += a(j ,i,n) * b(k ,j,n);
            
            tmp1 *= 2.0;
            for (j = i; j < n; ++j)
                b(k ,j,n) -= tmp1 * a(j ,i,n);
        }
        a(i,i,n) = tmp2;
    }
    
    for (i = 0; i < n; ++i)
        for (j = i + 1; j < n; ++j) {
            tmp1 = b(i,j,n);
            b(i,j,n) = b(j,i,n);
            b(j,i,n) = tmp1;
        }
    
    for (k = 0; k < n; ++k)
        for (i = n - 1; i >= 0; --i) {
            tmp1 = b(i,k,n);
            for (j = i + 1; j < n; ++j)
                tmp1 -= a(i ,j,n) * b(j,k,n);
            b(i,k,n) = tmp1 / a(i,i,n);
        }
    return 0;
}

double inner(double *v, double* u, int n) {
    //double res = 0.;
    t = 0.;
    for (int i = 0; i < n; i++) {
        t += v[i]*u[i];
    }
    return t;
}

void InitV(double *v, double *u) {
    for (int i = 0; i < 3; i++) v[i] = u[i];
}

void InitM(double *V, double *N) {
    for (int i = 0; i < 9; i++) V[i] = N[i];
}

void MxV ( double* V, double* v) {
    InitV(tmp,v);
    for (int i = 0; i < 3 ; i++) {
        v[i] = m(i,0)*tmp[0] + m(i,1)*tmp[1] + m(i,2)*tmp[2];
    }
}

void MxM ( double* V, double* N ) {
    //double t[3];
    for (int i = 0; i < 3 ; i++) {
        tmp[0] = m(i,0); tmp[1] = m(i,1); tmp[2] = m(i,2);
        for(int j = 0; j < 3; j++)
            m(i,j) = tmp[0]*n(0,j) + tmp[1]*n(1,j) + tmp[2]*n(2,j);//m(i,j) = t(i,0)*n(0,j) + t(i,1)*n(1,j) + t(i,2)*n(2,j);
    }
}

void intxvec(double c, double* v, double* res) {
    for (int i = 0; i < 3; i++)
    res[i] = c*v[i];
    //return res;
}

void plus(double *v, double *u, double* res) {
    
    for (int i = 0; i < 3; i++)
        res[i] = v[i] + u[i];
}

void minus(double *v, double *u, double* res) {
    
    for (int i = 0; i < 3; i++)
        res[i] = v[i] - u[i];
    
}

void Init(double *v, double x, double y, double z) {
    v[0] = x; v[1] = y; v[2] = z;
}

void MatrixId ( double *N ) {
    for (int i = 0; i < 9; i++) N[i] = 0.;
    n(0,0) = 1; n(1,1) = 1; n(2,2) = 1;
}

void cross(double* a, double* b, double *c) {
     Init(c, (-1)*(a[1]*b[2] - a[2]*b[1]), (-1)*(a[2]*b[0] - a[0]*b[2]), (-1)*(a[0]*b[1] - a[1]*b[0]));
}

void normalize(double* a, int n) {
    t = sqrt(inner(a,a,n));
    for (int i = 0; i < n; i++) {
        a[i] /= t;
    }
}

void reflect(double *a, double *b, double* res) {
    normalize(b,3); //printf("%f\n",inner(a, b, 3));
    intxvec(inner(a, b, 3), b, res);
    minus(res, a, res);
    intxvec(2, res, res);
    plus(a, res, res);
}

double signed_area(Point a, Point b, Point c)  {
    return ((b.x-a.x) * (c.y-a.y) - (b.y-a.y) * (c.x-a.x));
}
