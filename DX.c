// Creating 2D and 3D graphics
//
// David Mikheev 7/18/19

#include "DGraphX.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define n(i,j) N[3*i+j]

#define ALPHA 64.
#define THETA 45.
#define PHI 35.
#define DPHI 1.7
#define COEF 1.
#define X 1
#define Y 2
#define Z 3
//#define EYE 40.
#define CW ( -1 )
#define CCW 1
#define K  0.01745329251
#define PI 3.14159265358
#define EPS_X 5e-2
#define EPS_Y 5e-2

//#define min(a,b) fabs(a) < fabs(b) ? fabs(a) : fabs(b)

// If there can be an X-axis and a Y-axis, why not a Z-axis?

double xc, yc;

static double M[9];

static double view[3];

static double light[3];

static double normal[3];

static unsigned long cx[128];

static double cameraPos = 20.;

static double nearY, farY;

static double tmp[3];

static unsigned long SHADINGCOLOR;

static Point pt[4];

static XPoint xpt[4];

static int key = 0;

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

void InitCameraPosition(double pos) {
    cameraPos = pos;
}

void ChangeCameraPosition(double valueOfShift) {
    cameraPos += valueOfShift;
    key = 1;
}

void SetShadingColor(unsigned long color) {
    SHADINGCOLOR = color;
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
//============ colors (to improve) ==========
void pallette(unsigned long firstColor, unsigned long secondColor, int modeColor) {
    unsigned long tmpx;
    XColor fc, sc;
    if(modeColor == 0) fc.pixel = firstColor; else fc.pixel = BLACK;
    QueryColor(&fc);
    if(modeColor == 0) sc.pixel = secondColor; else sc.pixel = SHADINGCOLOR;
    QueryColor(&sc);
    for (int i = 0; i < 100; i++) {
       WGetColor(fc.red + (sc.red - fc.red)*i/100,
                 fc.green + (sc.green - fc.green)*i/100,
                 fc.blue + (sc.blue - fc.blue)*i/100, &tmpx);
       cx[i] = tmpx;
    }
}
//=============================
void elemRotation(int Axis, double angle, int mode) {
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
    elemRotation(X,K*DPHI,mode);
}

void RotateY(int mode) {
    elemRotation(Y,K*DPHI,mode);
}

void RotateZ(int mode) {
    elemRotation(Z,K*DPHI,mode);
}

void IdMatrix() {
    MatrixId(M);
}

void Projection(Point *p, double x, double y, double z) {
    Init(tmp, x, y, z);
    MxV(M, tmp);
    p->x = 2*tmp[0]*(nearY - farY)/(cameraPos - tmp[1]);
    p->y = 2*tmp[2]*(nearY - farY)/(cameraPos - tmp[1]);
}

Point ObliqueProjection ( /*PointR3 p*/ double x, double y, double z ) {
    Point P;
    P.x = x - y*cos(K*ALPHA)/2;
    P.y = z - y*sin(K*ALPHA)/2;
    return P;
}

void InitialPosition () {
    MatrixId(M);
    elemRotation(X,K*PHI,CCW);
    elemRotation(Z,K*THETA,CW);
}

void DrawAxes() {
    Point p1,p2;
    SetLineWidth(0);
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Projection(&p1,5-t,t*cos(K*phi)/2,t*sin(K*phi)/2);
            Projection(&p2,5-t,t*cos(K*(phi+0.5))/2,t*sin(K*(phi+0.5))/2);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Projection(&p1,t*cos(K*phi)/2,5-t,t*sin(K*phi)/2);
            Projection(&p2,t*cos(K*(phi+0.5))/2,5-t,t*sin(K*(phi+0.5))/2);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360.5; phi+=0.5) {
            Projection(&p1,t*cos(K*phi)/2,t*sin(K*phi)/2,5-t);
            Projection(&p2,t*cos(K*(phi+0.5))/2,t*sin(K*(phi+0.5))/2,5-t);
            DrawLine(p1.x,p1.y,p2.x,p2.y);
        }
    } SetLineWidth(3);
    Projection(&p1,-5,0,0);
    Projection(&p2,5,0,0);
    DrawLine(p1.x,p1.y,p2.x,p2.y);
    Projection(&p1,5.3,0,0);
    WDrawString ("X", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
    Projection(&p1,0,-5,0);
    Projection(&p2,0,5,0);
    DrawLine(p1.x,p1.y,p2.x,p2.y);
    Projection(&p1,0,5.3,0);
    WDrawString ("Y", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
    Projection(&p1,0,0,-5);
    Projection(&p2,0,0,5);
    DrawLine(p1.x,p1.y,p2.x,p2.y);
    Projection(&p1,0,0,5.3);
    WDrawString ("Z", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
}

void VectorSight() {
    double N[9], T[9];
    Init(view, 0, 1, 0);
    InitM(T, M);
    InvMatrix(3, T, N);
    MxV(N, view); //printf("(%f,%f,%f)\n", view[0],view[1],view[2]);
}

void VectorLight(double x, double y, double z) {
    double N[9], T[9];
    Init(light, x, y, z);
    normalize(light, 3);
    InitM(T, M);
    InvMatrix(3, T, N);
    MxV(N, light);
}

int LineSightControlAxes(double x, double y, double z) {
    double tmpx[3]; InitV(tmpx, view);
    view[0] = (cameraPos*view[0] - x)/cameraPos; view[1] = (cameraPos*view[1] - y)/cameraPos;
    view[2] = (cameraPos*view[2] - z)/cameraPos;
    if((fabs(y*view[2]-z*view[1]) < 0.03 && y/view[1] < 0 && fabs(x-(y/view[1])*view[0]) < 5) ||
       (fabs(x*view[2]-z*view[0]) < 0.03 && z/view[2] < 0 && fabs(y-(z/view[2])*view[1]) < 5) ||
       (fabs(x*view[1]-y*view[0]) < 0.03 && x/view[0] < 0 && fabs(z-(x/view[0])*view[2]) < 5)) {
       InitV(view, tmpx);
        return 1;
    } InitV(view, tmpx);
    return 0;
}

void elemDraw(double x, double y, double (*z[]) (double, double), int i, double *u, double *w, int modeColor) {
    Init(tmp, x, y, z[i](x,y)); MxV(M, tmp);
    if(cameraPos - tmp[1] > 5e-1) {
    Projection(&pt[0], x, y, z[i](x, y));
    Projection(&pt[1], x + EPS_X, y, z[i](x + EPS_X, y));
    if(EPS_X*EPS_Y > 1e-2) Projection(&pt[2], x, y + EPS_Y, z[i](x, y + EPS_Y));
    else Projection(&pt[2], x + EPS_X, y + EPS_Y, z[i](x + EPS_X, y + EPS_Y));
    if(modeColor == 0)  WSetColor(cx[(int)((0.5 + atan(COEF * z[i](x+EPS_X/2., y + EPS_Y/2.))/PI)*100)]);
    else {
        Init(u, 0, EPS_Y, (z[i](x, y + EPS_Y) - z[i](x, y)));
        Init(w, EPS_X, 0, (z[i](x + EPS_X, y) -  z[i](x, y)));
        cross(u, w, normal);
        normalize(normal, 3);
        if(normal[2] < 1e-20) Init(normal, (-1)*normal[0], (-1)*normal[1], (-1)*normal[2]);
        WSetColor(cx[(int)((inner(normal, light, 3)/2. + 1./2.)*99)]);
    }
    if(EPS_X*EPS_Y < 1e-2)  {
        Projection(&pt[3], x, y + EPS_Y, z[i](x, y + EPS_Y));
        invmapX(pt, xpt, 4);
        WFillPolygon (xpt, 4);
    } else {
    invmapX(pt, xpt, 3); //DrawLineX(p1.x,p1.y,p2.x,p2.y,p,q);
    WFillPolygon (xpt, 3);
    
    Projection(&pt[0], x + EPS_X, y + EPS_Y, z[i](x + EPS_X, y + EPS_Y));
    
        invmapX(pt, xpt, 3);
        WFillPolygon (xpt, 3);
    }
    }
}


int comp (const void *i, const void *j)
{
    if(((const Pair*) i) -> z > ((const Pair*) j) -> z) return 1;
    else return -1;
}

int ufun(double a, double b, double t) {
    if (a > b) return (t < a);
    else return (t > a);
}
    
void DrawGraph3DX(double a, double b, double (*z[]) (double, double), unsigned long firstColor,
                  unsigned long secondColor, int mode, int modeColor, int n) {
    double dx = EPS_X, dy = EPS_Y, miny = a, maxy = b, minx = a, maxx = b; Pair m[8]; 
    double u[3], w[3];
    pallette(firstColor, secondColor, modeColor);
    VectorSight();
    if(key == 0) InitCameraPosition(b + 10);
    farY = a; nearY = b;
        if(fabs(cameraPos*view[1] - (a + b)/2.) > fabs(cameraPos*view[0] - (a + b)/2.)) {
            if( cameraPos*view[1] - (a + b)/2.< 1e-20 ) { dy *= -1; miny = b; maxy = a; }
            for (double y = miny; ufun(maxy, miny, y); y+=dy) {
                for (double x = a; x < b; x+=dx) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(view[2] > 1e-20) {
                            if(!(mode && LineSightControlAxes(x, y, z[m[i].num](x, y))))
                                elemDraw(x, y, z, m[i].num, u, w, modeColor);
                        } else {
                            if(!(mode && LineSightControlAxes(x, y, z[m[n-i-1].num](x, y))))
                                elemDraw(x, y, z, m[n-i-1].num, u, w, modeColor);
                        }
                    }
                }
            }
        }
        if(fabs(cameraPos*view[1] - (a + b)/2.) < fabs(cameraPos*view[0] - (a + b)/2.)) {
            if( cameraPos*view[0] - (a + b)/2. < 1e-20 ) { dx *= -1; maxx = a; minx = b;}
            for (double x = minx; ufun(maxx, minx, x); x+=dx) {
                for (double y = a; y < b; y+=dy) {
                    for (int i = 0; i < n; i++)  { m[i].z = z[i](x,y); m[i].num = i; }
                    qsort(m, n, sizeof(*m), &comp);
                    for (int i = 0; i < n; i++)  {
                        if(view[2] > 1e-20) {
                            if(!(mode && LineSightControlAxes(x, y, z[m[i].num](x, y))))
                                elemDraw(x, y, z, m[i].num, u, w, modeColor);
                        } else {
                            if(!(mode && LineSightControlAxes(x, y, z[m[n-i-1].num](x, y))))
                                elemDraw(x, y, z, m[n-i-1].num, u, w, modeColor);
                        }
                    }
                }
            }
        }
}

void ParametricCurve3D(double (*curvefun[]) (double), double it, double s, int mode) {
    double dt = 1e-4; Point pts[2]; XPoint ptx[2];
    VectorSight();
    nearY = farY = curvefun[1](it);
    SetLineWidth(3);
    pallette(RED, BLUE, 0);
    for(double t = it; t < s; t += dt) {
        if(curvefun[1](t) > nearY) nearY = curvefun[1](t);
        if(curvefun[1](t) < farY) farY = curvefun[1](t);
    }
    for(double t = it; t < s; t += dt) {
        if(!(mode && LineSightControlAxes(curvefun[0](t), curvefun[1](t), curvefun[2](t)))) {
            Projection(&pts[0], curvefun[0](t), curvefun[1](t), curvefun[2](t));
            Projection(&pts[1], curvefun[0](t+dt), curvefun[1](t+dt), curvefun[2](t+dt));
            WSetColor(cx[(int)((0.5 + atan(COEF * curvefun[2](t + dt/2.))/PI)*100)]);
            DrawLineX(pts, ptx);
        }
    }
}

void ParametricGraph3D(double (*parfun[]) (double, double), unsigned long firstColor, unsigned long secondColor,
                       double it, double ft, double is, double fs, int mode, int modeColor) {
    double dt = (ft-it)/100., ds = (fs-is)/100.; double u[3],w[3],tmpx[3];
    Pair z[0x8000]; Point m[0x8000];
    int j = 0;
    VectorSight();
    pallette(firstColor, secondColor, modeColor);
    nearY = farY = parfun[1](it, is);
    for(double t = it; t < ft; t += dt) {
        for (double s = is; s < fs; s += ds) {
            if(parfun[1](t,s) > nearY) nearY = parfun[1](t,s);
            if(parfun[1](t,s) < farY) farY = parfun[1](t,s);
        }
    } if(key == 0) InitCameraPosition(nearY + 10);
    for(double t = it; t < ft; t += dt) {
        for (double s = is; s < fs; s += ds) {
            m[j].x = t; m[j].y = s;
            Init(tmpx, cameraPos*view[0] - parfun[0](t,s),
                 cameraPos*view[1] - parfun[1](t,s),
                 cameraPos*view[2] - parfun[2](t,s));
            z[j].z = inner(tmpx, tmpx, 3); z[j].num = j; j++;
        }
    } qsort(z, j, sizeof(*z), &comp);
    for (int i = j - 1; i >= 0; i--) {
        Init(tmp, parfun[0](m[z[i].num].x, m[z[i].num].y), parfun[1](m[z[i].num].x, m[z[i].num].y), parfun[2](m[z[i].num].x, m[z[i].num].y)); MxV(M, tmp);
        if(cameraPos -  tmp[1] > 5e-1) {
        if(!(mode && LineSightControlAxes(parfun[0](m[z[i].num].x, m[z[i].num].y),
                    parfun[1](m[z[i].num].x, m[z[i].num].y), parfun[2](m[z[i].num].x, m[z[i].num].y)))) {
            Projection(&pt[0], parfun[0](m[z[i].num].x, m[z[i].num].y),
                       parfun[1](m[z[i].num].x, m[z[i].num].y),
                       parfun[2](m[z[i].num].x, m[z[i].num].y));
            Projection(&pt[1], parfun[0](m[z[i].num].x + dt, m[z[i].num].y),
                       parfun[1](m[z[i].num].x + dt, m[z[i].num].y),
                       parfun[2](m[z[i].num].x + dt, m[z[i].num].y));
            Projection(&pt[2], parfun[0](m[z[i].num].x, m[z[i].num].y + ds),
                       parfun[1](m[z[i].num].x, m[z[i].num].y + ds),
                       parfun[2](m[z[i].num].x, m[z[i].num].y + ds));
            if(modeColor == 0) WSetColor(cx[(int)((0.5 + atan(COEF * parfun[2](m[z[i].num].x, m[z[i].num].y))/PI)*100)]);
            else {
                Init(u, parfun[0](m[z[i].num].x + dt, m[z[i].num].y) - parfun[0](m[z[i].num].x, m[z[i].num].y),
                     parfun[1](m[z[i].num].x + dt, m[z[i].num].y) - parfun[1](m[z[i].num].x, m[z[i].num].y),
                     parfun[2](m[z[i].num].x + dt, m[z[i].num].y) - parfun[2](m[z[i].num].x, m[z[i].num].y));
                Init(w, parfun[0](m[z[i].num].x, m[z[i].num].y + ds) - parfun[0](m[z[i].num].x, m[z[i].num].y),
                      parfun[1](m[z[i].num].x, m[z[i].num].y + ds) - parfun[1](m[z[i].num].x, m[z[i].num].y),
                     parfun[2](m[z[i].num].x, m[z[i].num].y + ds) - parfun[2](m[z[i].num].x, m[z[i].num].y));
                cross(u, w, normal);
                normalize(normal, 3); //printf("%f\n", inner(normal, normal, 3));
                WSetColor(cx[(int)((inner(normal, light, 3)/2. + 1./2.)*99)]);
            }
                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
            
            Projection(&pt[0], parfun[0](m[z[i].num].x + dt, m[z[i].num].y + ds),
                       parfun[1](m[z[i].num].x + dt, m[z[i].num].y + ds),
                       parfun[2](m[z[i].num].x + dt, m[z[i].num].y + ds));

                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
        }
        }
    }
}


/*if(modeColor != 0) {
 Init(u, 0, (-1)*dy, (z[i](x + dx, y) - z[i](x + dx, y + dy)));
 Init(w, (-1)*dx, 0, (z[i](x, y + dy) -  z[i](x + dx, y + dy)));
 cross(u, w, normal);
 normalize(normal, 3);
 if(normal[2] < 1e-20) Init(normal, (-1)*normal[0], (-1)*normal[1], (-1)*normal[2]);
 WSetColor(cx[(int)((inner(normal, light, 3)/2. + 1./2.)*99)]);
 }*/
