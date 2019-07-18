// File for creating 2D and 3D graphics
//
// David Mikheev 7/18/19

#include "DGraphX.h"
//#include "Interp3D.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define m(i,j) V[3*i+j]
#define n(i,j) N[3*i+j]
#define t(i,j) T[3*i+j]

#define ALPHA 64.
#define THETA 45.
#define PHI 35.
#define DPHI 1.7
#define X 1
#define Y 2
#define Z 3
#define CW ( -1 )
#define CCW 1
#define K  0.01745329251
#define PI 3.14159265358

#define min(a,b) fabs(a) < fabs(b) ? fabs(a) : fabs(b)

// If there can be an X-axis and a Y-axis, why not a Z-axis?

double xc, yc;

static double M[9];

Point map(int x, int y) {
    Point P; double kx = (-xmin)/(xmax - xmin), ky = (ymax)/(ymax - ymin);
    P.x = (x - width * kx) * xc;
    P.y = (height * ky - y) * yc;
    return P;
}

XPoint invmap(double x, double y) {
    XPoint pt; double kx = (-xmin)/(xmax - xmin), ky = (ymax)/(ymax - ymin);
    pt.x = x/xc + width * kx;
    pt.y = -y/yc + height * ky;
    return pt;
}

void scale(int mode, double a, double b, double lambda) {
    switch (mode) {
        case 1:
            xmin = (b + a)/2 - lambda * ((b + a)/2 - xmin);
            xmax = (b + a)/2 + lambda * (-(b + a)/2 + xmax);
            ymin *= lambda; ymax *= lambda;
            xc = (xmax - xmin)/width; yc = (ymax - ymin)/height;
            break;
        case -1:
            xmin = (b + a)/2 - (1./lambda)*((b + a)/2 - xmin);
            xmax = (b + a)/2 + (1./lambda)*(-(b + a)/2 + xmax);
            ymin /= lambda; ymax /= lambda;
            xc = (xmax - xmin)/width; yc = (ymax - ymin)/height;
            break;
    }
}

void xshift(int mode, double valueOfShift) {
    switch (mode) {
        case 1:
            xmin += valueOfShift * (xmax - xmin); xmax += valueOfShift * (xmax - xmin);
            break;
        case -1:
            xmin -= valueOfShift * (xmax - xmin); xmax -= valueOfShift * (xmax - xmin);
            break;
    }
}

void yshift(int mode, double valueOfShift) {
    switch (mode) {
        case 1:
            ymin += valueOfShift*(ymax-ymin); ymax += valueOfShift*(ymax-ymin);
            break;
        case -1:
            ymin -= valueOfShift*(ymax-ymin); ymax -= valueOfShift*(ymax-ymin);
            break;
    }
}

void DrawGraph2D(double f(double), double a, double b) {
    double dx = 1e-3;
    WSetColor (RED);
    SetLineWidth(4);
    for (double t = a; t < b; t += dx) {
        if(t > xmin && t < xmax && f(t) > ymin && f(t) < ymax &&
           t + dx > xmin && t + dx < xmax && f(t + dx) > ymin && f(t + dx) < ymax)
            DrawLine(t, f(t), t + dx, f(t + dx));
    }
}

void DrawParametric2D(double (*parf[]) (double), double it, double ft) {
    double dt = 1e-3;
    WSetColor (MAGENTA);
    SetLineWidth(4);
    for (double t = it; t < ft + dt; t += dt) {
        if(parf[0](t) > xmin && parf[0](t) < xmax && parf[1](t) > ymin && parf[1](t) < ymax &&
           parf[0](t + dt) > xmin && parf[0](t + dt) < xmax &&
           parf[1](t + dt) > ymin && parf[1](t + dt) < ymax)
            DrawLine(parf[0](t), parf[1](t), parf[0](t + dt), parf[1](t + dt));
    }
}


void DrawLine(double x1, double y1, double x2, double y2) {
    WDrawLine(invmap(x1,y1).x, invmap(x1,y1).y,
              invmap(x2,y2).x, invmap(x2,y2).y);
}

void DrawGrid() {
    SetLineWidth(0);
    for (int i = (int) xmin; i <= (int) xmax ; i++)
        if(i != 0) DrawLine((double) i, ymin, (double) i, ymax);
    for (int i = (int) ymin; i <= (int) ymax ; i++)
        if(i != 0) DrawLine(xmin, (double) i, xmax, (double) i);
}

void DrawAxes2D() {
    WSetColor (LIGHTGREEN);
    SetLineWidth(3);
    WDrawLine (0, invmap(0,0).y, width-29, invmap(0,0).y);
    WFillTriangle (invmap(0,0).x, 5, invmap(0,0).x-5, 30, invmap(0,0).x+5, 30);
    WDrawLine (invmap(0,0).x, 29, invmap(0,0).x, height);
    WFillTriangle (width-5, invmap(0,0).y, width-30, invmap(0,0).y - 5, width-30, invmap(0,0).y + 5);
    WDrawString ("X", width - 20, invmap(0,0).y + 20);
    WDrawString ("Y", invmap(0,0).x - 20, 20);
    DrawPoint(invmap(0,0).x, invmap(0,0).y,4);
}

void invmapX(Point * pts, XPoint *pt, int n) {
    for (int i = 0; i < n; i++)
        pt[i] = invmap(pts[i].x,pts[i].y);
}


void DrawLineX(Point * pts, XPoint *pt) {
    invmapX(pts, pt, 2);
    WDrawLine(pt[0].x, pt[0].y, pt[1].x, pt[1].y);
}

void pallette(unsigned long* cx) {
    unsigned long tmpx;
    for (int i = 0; i < 100; i++) {
       WGetColor(0,(i/100.)*65535, (1-i/100.)*65535, &tmpx);
       cx[i] = tmpx;
    }
}

void gradient( double x, double y, unsigned long color, double (*z[]) (double,double), int i) {
    double t = 0.5 + atan(0.9 * z[i](x,y))/PI;
     WSetColor(t*color/* + ((1-t)*color2)*/);//cx[(int)(t*100)]
}

void Rotation(int Axis, double angle, int mode) {
    double N[9];
    MatrixId(N);
    switch (Axis) {
        case 1:
            n(1,1) = n(2,2) = cos(angle);
            n(2,1) = n(1,2) = mode*sin(angle);
            n(2,1) *= (-1.);
            break;
        case 2:
            n(0,0) = n(2,2) = cos(angle);
            n(2,0) = n(0,2) = mode*sin(angle);
            n(2,0) *= (-1.);
            break;
        case 3:
            n(0,0) = n(1,1) = cos(angle);
            n(1,0) = n(0,1) = mode*sin(angle);
            n(1,0) *= (-1.);
            break;
        default:
            break;
    } //printf("%f\n", n(1,0));
    MxM(M,N);
}

void RotateX(int mode) {
    Rotation(X,K*DPHI,mode);
}

void RotateY(int mode) {
    Rotation(Y,K*DPHI,mode);
}

void RotateZ(int mode) {
    Rotation(Z,K*DPHI,mode);
}

void IdMatrix() {
    MatrixId(M);
}

void Projection(Point *p, double* v) {
    MxV(M,v);
    p->x = v[0]*20./(20-v[1]);
    p->y = v[2]*20./(20-v[1]);
}

Point ObliqueProjection ( /*PointR3 p*/ double x, double y, double z ) {
    Point P;
    P.x = x - y*cos(K*ALPHA)/2;
    P.y = z - y*sin(K*ALPHA)/2;
    return P;
}

void OrthoProjection () {
    Rotation(X,K*PHI,CCW);
    Rotation(Z,K*THETA,CW);
}

void DrawAxes() {
    Point p1,p2; double u[3];
    SetLineWidth(0);
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Init(u,5-t,t*cos(K*phi)/2,t*sin(K*phi)/2);
            Projection(&p1,u);
            Init(u,5-t,t*cos(K*(phi+0.5))/2,t*sin(K*(phi+0.5))/2);
            Projection(&p2,u);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Init(u,t*cos(K*phi)/2,5-t,t*sin(K*phi)/2);
            Projection(&p1,u);
            Init(u,t*cos(K*(phi+0.5))/2,5-t,t*sin(K*(phi+0.5))/2);
            Projection(&p2,u);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Init(u,t*cos(K*phi)/2,t*sin(K*phi)/2,5-t);
            Projection(&p1,u);
            Init(u,t*cos(K*(phi+0.5))/2,t*sin(K*(phi+0.5))/2,5-t);
            Projection(&p2,u);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    } SetLineWidth(3);
    Init(u,-5,0,0); Projection(&p1,u);
    Init(u,5,0,0); Projection(&p2,u);
    DrawLine(p1.x,p1.y,p2.x,p2.y); Init(u,5.3,0,0); Projection(&p1,u);
    WDrawString ("X", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
    Init(u,0,-5,0); Projection(&p1,u);
    Init(u,0,5,0); Projection(&p2,u);
    DrawLine(p1.x,p1.y,p2.x,p2.y); Init(u,0,5.3,0); Projection(&p1,u);
    WDrawString ("Y", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
    Init(u,0,0,-5); Projection(&p1,u);
    Init(u,0,0,5); Projection(&p2,u);
    DrawLine(p1.x,p1.y,p2.x,p2.y); Init(u,0,0,5.3); Projection(&p1,u);
    WDrawString ("Z", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
}

void VectorSight(double* v) {
    double N[9],T[9];
    Init(v,0,1,0);
    InitM(T,M);
    InvMatrix(3,T,N);
    MxV(N,v); //printf("(%f,%f,%f)\n", v[0],v[1],v[2]);
}

int LineSightControlAxes(double* v,double x, double y, double z) { // to improve
    double tmp[3]; InitV(tmp,v); v[0] = (20*v[0] - x)/20.; v[1] = (20*v[1] - y)/20.; v[2] = (20*v[2] - z)/20.;
    if((fabs(y*v[2]-z*v[1]) < 0.03 && y/v[1] < 0 && fabs(x-(y/v[1])*v[0]) < 5) ||
       (fabs(x*v[2]-z*v[0]) < 0.03 && z/v[2] < 0 && fabs(y-(z/v[2])*v[1]) < 5) ||
       (fabs(x*v[1]-y*v[0]) < 0.03 && x/v[0] < 0 && fabs(z-(x/v[0])*v[2]) < 5)) {
       InitV(v,tmp);
        return 1;
    } InitV(v,tmp);
    return 0;
}

void elemDraw(double x, double y, double (*z[]) (double, double), double dx, double dy, int i,
              double *u, double *w, Point *pts, XPoint *pt) {
    Init(u,x,y,z[i](x,y)); Init(w,x+dx,y,z[i](x+dx,y));
    Projection(&pts[0],u); Projection(&pts[1],w);
    Init(u,x,y+dy,z[i](x,y+dy)); Init(w,x+dx,y+dy,z[i](x+dx,y+dy));
    Projection(&pts[3],u); Projection(&pts[2],w);
    gradient(x+dx/2.,y+dy/2., BLUE , z, i);
    invmapX(pts, pt, 4); //DrawLineX(p1.x,p1.y,p2.x,p2.y,p,q);
    WFillPolygon (pt, 4);
}

int comp (const void *i, const void *j)
{
    if(((const Pair*) i)->z > ((const Pair*) j)->z) return 1;
    else return -1;
}
    
void DrawGraph3DX(double a, double b, double (*z[]) (double, double), int mode, int n) {
    double dx = 5*(1e-2), dy = 5*(1e-2), v[3]; Pair m[32]; 
    double u[3],w[3]; Point pts[4]; XPoint pt[4];
    VectorSight(v); WSetColor(BLUE);
        if( 20.*v[1] - (a + b)/2. > 1e-20 && fabs(20.*v[1] - (a + b)/2.) > fabs(20.*v[0] - (a + b)/2.)) {
            for (double y = a; y < b + 1e-2; y+=dy) {
                for (double x = a; x < b + 1e-2; x+=dx) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(v[2] > 1e-20) {
                            if((fabs(z[m[i].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[i].num](x, y)))) {
                                elemDraw(x, y, z, dx, dy, m[i].num, u, w , pts, pt);
                            }
                        } else {
                            if((fabs(z[m[n-i-1].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[n-i-1].num](x, y)))) {
                                elemDraw(x, y, z, dx, dy, m[n-i-1].num, u, w , pts, pt);
                            }
                        }
                    }
                }
            }
        }
        if( 20.*v[1] - (a + b)/2.< 1e-20 && fabs(20.*v[1] - (a + b)/2.) > fabs(20.*v[0] - (a + b)/2.)) {
            for (double y = b; y > a - 1e-2; y -= dy) {
                for (double x = a; x < b + 1e-2; x += dx) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(v[2] > 1e-20) {
                            if((fabs(z[m[i].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[i].num](x, y)))) {
                                elemDraw(x, y, z, dx, -dy, m[i].num, u, w , pts, pt);
                            }
                        } else {
                            if((fabs(z[m[n-i-1].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[n-i-1].num](x, y)))) {
                                elemDraw(x, y, z, dx, -dy, m[n-i-1].num, u, w , pts, pt);
                            }
                        }
                    }
                }
            }
        }
        if( 20.*v[0] - (a + b)/2. > 1e-20 && fabs(20.*v[1] - (a + b)/2.) < fabs(20.*v[0] - (a + b)/2.)) {
            for (double x = a; x < b + 1e-2; x+=dx) {
                for (double y = a; y < b + 1e-2; y+=dy) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(v[2] > 1e-20) {
                            if((fabs(z[m[i].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[i].num](x, y)))) {
                                elemDraw(x, y, z, dx, dy, m[i].num, u, w , pts, pt);
                            }
                        } else {
                            if((fabs(z[m[n-i-1].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[n-i-1].num](x, y)))) {
                                elemDraw(x, y, z, dx, dy, m[n-i-1].num, u, w , pts, pt);
                            }
                        }
                    }
                }
            }
        }
        if( 20.*v[0] - (a + b)/2. < 1e-20 && fabs(20.*v[1] - (a + b)/2.) < fabs(20.*v[0] - (a + b)/2.)) {
            for (double x = b; x > a - 1e-2; x-=dx) {
                for (double y = a; y < b + 1e-2; y+=dy) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(v[2] > 1e-20) {
                            if((fabs(z[m[i].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[i].num](x, y)))) {
                                elemDraw(x, y, z, -dx, dy, m[i].num, u, w , pts, pt);
                            }
                        } else {
                            if((fabs(z[m[n-i-1].num](x, y)) < 5) &&
                               !(mode && LineSightControlAxes(v, x, y, z[m[n-i-1].num](x, y)))) {
                                elemDraw(x, y, z, -dx, dy, m[n-i-1].num, u, w , pts, pt);
                            }
                        }
                    }
                }
            }
        }
}

void ParametricCurve3D(double (*curvefun[]) (double), double it, double s, int mode) {
    double dt = 1e-4; double u[3],w[3],v[3]; Point pts[2]; XPoint pt[2];
    VectorSight(v);
    WSetColor(RED); SetLineWidth(3);
    for(double t = it; t < s; t += dt) {
        if(!(mode && LineSightControlAxes(v, curvefun[0](t),curvefun[1](t),curvefun[2](t)))) {
            Init(u,curvefun[0](t),curvefun[1](t),curvefun[2](t)); Init(w,curvefun[0](t+dt),curvefun[1](t+dt),curvefun[2](t+dt));
            Projection(&pts[0],u); Projection(&pts[1],w);
            DrawLineX(pts,pt);
        }
    }
}

void ParametricGraph3D(double (*parfun[]) (double, double), double it, double ft, double is, double fs, int mode) {
    double dt = 5e-2, ds = 5e-2; double u[3],w[3],v[3],tmpx[3]; Point pts[4], m[0x8000]; XPoint pt[4]; Pair z[0x8000]; int j = 0;
    VectorSight(v);
    WSetColor(BLUE);
    for(double t = it; t < ft; t += dt) {
        for (double s = is; s < fs; s += ds) {
            m[j].x = t; m[j].y = s;
            Init(tmpx, 20*v[0] - parfun[0](t,s), 20*v[1] - parfun[1](t,s), 20*v[2] - parfun[2](t,s));
            z[j].z = inner(tmpx, tmpx, 3); z[j].num = j; j++;
        }
    } qsort(z, j, sizeof(*z), &comp);
    for (int i = j - 1; i >= 0; i--) {
        if((fabs(parfun[2](m[z[i].num].x, m[z[i].num].y)) < 5) &&
           !(mode && LineSightControlAxes(v, parfun[0](m[z[i].num].x, m[z[i].num].y),
                    parfun[1](m[z[i].num].x, m[z[i].num].y), parfun[2](m[z[i].num].x, m[z[i].num].y)))) {
            Init(u, parfun[0](m[z[i].num].x, m[z[i].num].y),  parfun[1](m[z[i].num].x, m[z[i].num].y),  parfun[2](m[z[i].num].x, m[z[i].num].y));
            Projection(&pts[0],u);
            Init(w, parfun[0](m[z[i].num].x + dt, m[z[i].num].y),
                  parfun[1](m[z[i].num].x + dt, m[z[i].num].y),  parfun[2](m[z[i].num].x + dt, m[z[i].num].y));
            Projection(&pts[1],w);
            Init(u, parfun[0](m[z[i].num].x, m[z[i].num].y + ds),
                  parfun[1](m[z[i].num].x, m[z[i].num].y + ds),  parfun[2](m[z[i].num].x, m[z[i].num].y + ds));
            Projection(&pts[3],u);
            Init(w, parfun[0](m[z[i].num].x + dt, m[z[i].num].y + ds),
                  parfun[1](m[z[i].num].x + dt, m[z[i].num].y + ds),  parfun[2](m[z[i].num].x + dt, m[z[i].num].y + ds));
            Projection(&pts[2],w);
            WSetColor((0.5 + atan(0.9 * parfun[2](m[z[i].num].x, m[z[i].num].y))/PI)*BLUE);
            invmapX(pts, pt, 4);
            WFillPolygon (pt, 4);
        }
    }
    
}
