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
    double res = 0.;
    for (int i = 0; i < n; i++) {
        res += v[i]*u[i];
    }
    return res;
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
    double t[3];
    for (int i = 0; i < 3 ; i++) {
        t[0] = m(i,0); t[1] = m(i,1); t[2] = m(i,2);
        for(int j = 0; j < 3; j++)
            m(i,j) = t[0]*n(0,j) + t[1]*n(1,j) + t[2]*n(2,j);//m(i,j) = t(i,0)*n(0,j) + t(i,1)*n(1,j) + t(i,2)*n(2,j);
    }
}

void Init(double *v, double x, double y, double z) {
    v[0] = x; v[1] = y; v[2] = z;
}

void MatrixId ( double *N ) {
    for (int i =0; i<9; i++) N[i] = 0.;
    n(0,0) = 1; n(1,1) = 1; n(2,2) = 1;
}
