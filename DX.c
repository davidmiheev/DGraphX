// Creating 2D and 3D graphics
//
// David Miheev 7/18/19
// If there can be an X-axis and a Y-axis, why not a Z-axis?

#include "DGraphX.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define n(i,j) N[3*i+j]
#define cx(i,j) CX[BW_COLOR*i+j]

#define ALPHA 64.
#define THETA 45.
#define PHI 35.
#define DPHI 1.7
#define COEF 1.
#define X 1
#define Y 2
#define Z 3
#define N_COLOR 64
#define BW_COLOR 32


#define CW ( -1 )
#define CCW 1

#define K  0.01745329251
#define PI 3.14159265358
#define EPS_X 5e-2
#define EPS_Y 5e-2

#define BUF_SIZE 1000000

double xc = INIT_SCALE, yc =  INIT_SCALE;

double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

static double M[9];// matrix: local object coord -> global static coord

static double TMP_1[9];

static double TMP_2[9];

static double view[3];

static double light[3];

static double normal[3];

//static double m_normal[1000][3];

static unsigned long CX[N_COLOR*BW_COLOR];//[128];

static double centre[3];

static double cameraPos = 20.;

static double screenPos = 10.;

static double nearY, farY;

static double tmp[3];

static unsigned long SHADINGCOLOR;

static Point pt[4];

static XPoint xpt[4];

static double gridleft, gridright, griddown, gridup;

static int key = 0;

static int key_1 = 0;

static Pair z[BUF_SIZE];

static vertex zm[BUF_SIZE];
static polygon zb[BUF_SIZE];

void SetSize() {
    double xcent = (xmax+xmin)/2., ycent = (ymax+ymin)/2.;
    xmin = xcent-width*xc/2.; 
    xmax = xcent + width*xc/2.; 
    ymin = ycent - height*yc/2.; 
    ymax = ycent + height*yc/2.;
    gridleft = xmin; gridright = xmax; griddown = ymin; gridup = ymax;
}

Point map(int x, int y) {
    Point P; double kx = (-xmin)/(xmax - xmin), ky = (ymax)/(ymax - ymin);
    P.x = (x - width * kx) * xc;
    P.y = (height * ky - y) * yc;
    return P;
}

XPoint invmap(double x, double y) {
    XPoint pt;
    double kx = (-xmin)/(xmax - xmin), ky = (ymax)/(ymax - ymin);
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

void InitCameraPosition(double pos, double focus) {
    cameraPos = pos + focus;
    screenPos = pos;
}

void InitKey() {
    key = 0;
}

/*void InitScreenPosition(double pos) {
    screenPos = pos;
}*/

void ChangeScreenPosition(double valueOfShift) {
    screenPos += valueOfShift;
    key = 1;
}

void ChangeCameraPosition(double valueOfShift) {
    cameraPos += valueOfShift;
    screenPos += valueOfShift;
    key = 1;
}

void SetShadingColor(unsigned long color) {
    SHADINGCOLOR = color;
}

void DrawLinear(Point *pt, int n) {
    pallette(RED, BLUE, 0, n - 1);
    SetLineWidth(4);
    for (int i = 0; i < n - 1; i++) {
        WSetColor(cx(i, 0));
        DrawLine(pt[i].x, pt[i].y, pt[i + 1].x, pt[i + 1].y);
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

void DrawGrid(int mode) {
    double a = gridright, b = gridleft, c = gridup, d = griddown;
    SetLineWidth(0);
    if(gridright - gridleft > 2*(xmax - xmin) ||
       gridup - griddown > 2*(ymax - ymin))  {
        key_1++; a = xmax; b = xmin; c = ymax; d = ymin;
    }
    if(gridright - gridleft < (xmax - xmin)/2. ||
       gridup - griddown < (ymax - ymin)/2.) {
        key_1--; a = xmax; b = xmin; c = ymax; d = ymin;
    }
    gridright = a; gridleft = b; gridup = c; griddown = d;
    if(key_1 != 0) {
    switch (key_1 > 0) {
        case 1:
            for (int i = 0; i < key_1; i++) {
            xmin *= 2; xmax *= 2; ymin *= 2; ymax *= 2; xc *= 2; yc *= 2;
            }
            break;
        case 0:
            for (int i = 0; i > key_1; i--) {
            xmin /= 2; xmax /= 2; ymin /= 2; ymax /= 2; xc /= 2; yc /= 2;
            }
            break;
    }
    }
    if(mode == 1) {
    for (int i = (int) xmin; i <= (int) xmax ; i++)
    { if(i != 0) DrawLine((double) i, ymin, (double) i, ymax); }
    for (int i = (int) ymin; i <= (int) ymax; i++)
    {   if(i != 0) DrawLine(xmin, (double) i, xmax, (double) i); }
    }
     if(key_1 != 0) {
         switch (key_1 > 0) {
        case 1:
            for (int i = 0; i < key_1; i++) {
            xmin /= 2; xmax /= 2; ymin /= 2; ymax /= 2; xc /= 2; yc /= 2;
            }
            break;
        case 0:
            for (int i = 0; i > key_1; i--) {
            xmin *= 2; xmax *= 2; ymin *= 2; ymax *= 2; xc *= 2; yc *= 2;
            }
            break;
         }
     }
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
    DrawPoint(invmap(0,0).x, invmap(0,0).y, 4);
}

void invmapX(Point * pts, XPoint *pt, int n) {
    for (int i = 0; i < n; i++)
        pt[i] = invmap(pts[i].x,pts[i].y);
}


void DrawLineX(Point * pts, XPoint *pt) {
    invmapX(pts, pt, 2);
    WDrawLine(pt[0].x, pt[0].y, pt[1].x, pt[1].y);
}
//============ colors ==========
void pallette(unsigned long firstColor, unsigned long secondColor, int modeColor, int size) {
    unsigned long tmpx;
    XColor fc, sc, w;
    w.pixel = WHITE; QueryColor(&w);
    if(size > N_COLOR) size = N_COLOR;
    if(modeColor == 0) fc.pixel = firstColor; else fc.pixel = BLACK;
    QueryColor(&fc);
    if(modeColor == 0) sc.pixel = secondColor; else sc.pixel = SHADINGCOLOR;
    QueryColor(&sc);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < BW_COLOR; j++) {
            WGetColor(0.65*(fc.red + (sc.red - fc.red)*i/size + modeColor*(w.red - (fc.red + (sc.red - fc.red)*i/size))*j/BW_COLOR),
                      0.65*(fc.green + (sc.green - fc.green)*i/size + modeColor*(w.green - (fc.green + (sc.green - fc.green)*i/size))*j/BW_COLOR),
                      0.65*(fc.blue + (sc.blue - fc.blue)*i/size + modeColor*(w.blue - (fc.blue + (sc.blue - fc.blue)*i/size))*j/BW_COLOR), &tmpx);
            cx(i,j) = tmpx;
        }
    }
}

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

void SetCentre(double x, double y, double z) {
    centre[0] = x;
    centre[1] = y;
    centre[2] = z;
}

void translation(double *p, double x, double y, double z) {
    p[0] += x;
    p[1] += y;
    p[2] += z;
}

void Projection(Point *p, double x, double y, double z) {
    Init(tmp, x - centre[0], y - centre[1], z - centre[2]);
    MxV(M, tmp);
    p->x = tmp[0]*(cameraPos-screenPos)/(cameraPos - tmp[1]);
    p->y = tmp[2]*(cameraPos-screenPos)/(cameraPos - tmp[1]);
}

Point ObliqueProjection ( /*PointR3 p*/ double x, double y, double z ) {
    Point P;
    P.x = x - y*cos(K*ALPHA)/2;
    P.y = z - y*sin(K*ALPHA)/2;
    return P;
}

void InitialPosition () {
    for(int i = 0; i < BUF_SIZE; i++) zm[i].data = 0;
    MatrixId(M);
    elemRotation(X,K*PHI,CCW);
    elemRotation(Z,K*THETA,CW);
}

void DrawAxes(int *n) {
    Point p1,p2;
    double tmpx[3];
    SetLineWidth(0);
    *n = 0;
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            
            Init(tmpx, (5-t), t*cos(K*phi)/2., t*sin(K*phi)/2.);
             
            z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
            
            Projection(&p1,5-t,t*cos(K*phi)/2.,t*sin(K*phi)/2.);
             
            Projection(&p2,5-t,t*cos(K*(phi+2))/2.,t*sin(K*(phi+2))/2.);
            zm[*n].x = p1.x; zm[*n].y = p1.y;//DrawLine(p1.x,p1.y,p2.x,p2.y);
            zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1;
            (*n)++;
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            Init(tmpx, t*cos(K*phi)/2., (5-t), t*sin(K*phi)/2.);
            z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
            Projection(&p1,t*cos(K*phi)/2., 5-t, t*sin(K*phi)/2.);
            Projection(&p2,t*cos(K*(phi+2))/2.,5-t,t*sin(K*(phi+2))/2.);
            zm[*n].x = p1.x; zm[*n].y = p1.y;//DrawLine(p1.x,p1.y,p2.x,p2.y);
            zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1;
            (*n)++;
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            Init(tmpx, t*cos(K*phi)/2., t*sin(K*phi)/2.,(5-t));
            z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
            Projection(&p1,t*cos(K*phi)/2, t*sin(K*phi)/2, 5-t);
            Projection(&p2,t*cos(K*(phi+2))/2., t*sin(K*(phi+2))/2., 5-t);
            zm[*n].x = p1.x; zm[*n].y = p1.y;
            zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1;
            (*n)++;
        }
    } //SetLineWidth(3);
    for(double t = 0; t < 1-0.03; t += 0.03) {
        Init(tmpx, (-5 + 10*t),0,0);
        z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
        Projection(&p1,-5 + 10*t,0,0);
        Projection(&p2,-5 + 10*(t+0.03),0,0);
        zm[*n].x = p1.x; zm[*n].y = p1.y;
        zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1; (*n)++;
    }
    for(double t = 0; t < 1-0.03; t += 0.03) {
        Init(tmpx, 0, (-5 + 10*t),0);
        z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
        Projection(&p1,0,-5 + 10*t,0);
        Projection(&p2,0,-5 + 10*(t+0.03),0);
        zm[*n].x = p1.x; zm[*n].y = p1.y;//DrawLine(p1.x,p1.y,p2.x,p2.y);
        zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1; (*n)++;
    }
     for(double t = 0; t < 1-0.03; t += 0.03) {
        Init(tmpx, 0, 0,(-5 + 10*t));
        z[*n].z = inner(view, tmpx, 3); z[*n].num = *n;
        Projection(&p1,0,0,-5 + 10*t);
        Projection(&p2,0,0,-5 + 10*(t+0.03));
        zm[*n].x = p1.x; zm[*n].y = p1.y;//DrawLine(p1.x,p1.y,p2.x,p2.y);
        zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1; (*n)++;
    }
    //WSetColor(LIGHTGREEN);
        Projection(&p1,5.3,0,0);
        WDrawString ("X", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
        Projection(&p1,0,5.3,0);
        WDrawString ("Y", invmap(p1.x,p1.y).x, invmap(p1.x,p1.y).y);
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

int SightControl(double x, double y, double z, double eps) {
    Init(tmp, x, y, z);
    MxV(M, tmp);
    return (screenPos - tmp[1] > eps);
}

void DrawVector(double x1, double y1, double z1, double x2, double y2, double z2) {
    Point p1,p2;
    Projection(&p1, x1, y1, z1);
    Projection(&p2, x2, y2, z2);
    DrawLine(p1.x, p1.y, p2.x, p2.y);
}

int comp (const void *i, const void *j)
{
    if(((const Pair*) i) -> z > ((const Pair*) j) -> z) return 1;
    else return -1;
}

int comp_1 (const void *i, const void *j)
{
    if(((const polygon*) i) -> z > ((const polygon*) j) -> z) return 1;
    else return -1;
}

int ufun(double a, double b, double t) {
    if (a > b) return (t < a);
    else return (t > a);
}

void normalMap(polygon *m, int n) {
    for(int i = 0; i < n; i++) {
        normalize(m[i].normal, 3);
    }
}

/*int Inside(polygon p, Point a) {
    int j = 0;
    for(int i = 0; i < p.n; i++) {
        if((signed_area(p.vert[i%(p.n)],p.vert[(i+1)%(p.n)],p.vert[(i+2)%(p.n)])*signed_area(p.vert[i%(p.n)], p.vert[(i+1)%(p.n)], a)) > 0.) ++j;
    }
    if(j == p.n) return 1;
    else return 0;
}*/

void InvProj(double x, double y, double *res) {
    res[0] = x*(nearY - farY)*cameraPos/500.;
    res[1] = 0;
    res[2] = y*(nearY - farY)*cameraPos/500.;
    InitM(TMP_1, M);
    InvMatrix(3, TMP_1, TMP_2);
    MxV(TMP_2, res);
    res[0] = cameraPos*view[0] - res[0];
    res[1] = cameraPos*view[1] - res[1];
    res[2] = cameraPos*view[2] - res[2];
    normalize(res, 3);
}


void DrawLinear3D(double* a, int n, unsigned long Color) {
    double vec[3], res[3], h = 1./(n-1); int k = 0;
    pallette(Color, BLACK, 1, N_COLOR);
    VectorSight();
    if(key == 0)  InitCameraPosition(2, 1);
    printf("[DGraphX]:\tcameraPos = %f\n", cameraPos);
    for(int i = 0; i < n - 1; i++) {
        for(int j = 0; j < n - 1; j++) { //zb[k].n = 3;
            Init(zb[k].vert[0], j*h, i*h, a[i*n + j]);
            Init(zb[k].vert[1], (j+1)*h, i*h, a[i*n + j+1]);
            Init(zb[k].vert[2], j*h, (i+1)*h, a[(i+1)*n + j]);
            plus(zb[k].vert[0], zb[k].vert[1], vec);
            plus(vec, zb[k].vert[2], vec);
            intxvec(1./3., vec, vec);
            minus(vec, centre, vec);
            zb[k].z = inner(view, vec, 3);
            Init(vec, h, 0, a[i*n + j + 1] - a[i*n + j]);
            Init(res, 0, h, a[(i+1)*n + j] - a[i*n + j]);
            cross(vec, res, zb[k].normal); //zb[k + 1].n = 3;
            Init(zb[k+1].vert[0], j*h, (i+1)*h, a[(i+1)*n + j]);
            Init(zb[k+1].vert[1], (j+1)*h, i*h, a[i*n + j+1]);
            Init(zb[k+1].vert[2], (j+1)*h, (i+1)*h, a[(i+1)*n + j+1]);
            plus(zb[k+1].vert[0], zb[k+1].vert[1], vec);
            plus(vec, zb[k+1].vert[2], vec);
            intxvec(1./3., vec, vec);
            minus(vec, centre, vec);
            zb[k+1].z = inner(view, vec,3);
            Init(vec, h, 0, a[(i+1)*n + j + 1] - a[(i+1)*n + j]);
            Init(res, 0, h, a[(i+1)*n + j + 1] - a[i*n + j + 1]);
            cross(vec, res, zb[k+1].normal);
            k+=2;
        }
    }
    qsort(zb, k, sizeof(*zb), &comp_1);
    normalMap(zb, k); //WSetColor(RED);
    for (int i = 0; i < k; i++) {
        if(SightControl(zb[i].vert[0][0], zb[i].vert[0][1], zb[i].vert[0][2], 0.1)) {
        Projection(&pt[0], zb[i].vert[0][0], zb[i].vert[0][1], zb[i].vert[0][2]);
        Projection(&pt[1], zb[i].vert[1][0], zb[i].vert[1][1], zb[i].vert[1][2]);
        Projection(&pt[2], zb[i].vert[2][0], zb[i].vert[2][1], zb[i].vert[2][2]);
        //InvProj(zb[i].vert[0].x, zb[i].vert[0].y, res);
        Init(res, (cameraPos*view[0] -  zb[i].vert[0][0] +centre[0]), (cameraPos*view[1] - zb[i].vert[0][1] + centre[1]), (cameraPos*view[2] - zb[i].vert[0][2])+centre[2]);
        normalize(res,3);
        //WSetColor(RED); DrawVector(zb[i].vert[0][0], zb[i].vert[0][1], zb[i].vert[0][2], zb[i].vert[0][0] + res[0], zb[i].vert[0][1] + res[1], zb[i].vert[0][2] + res[2]);
        if(inner(res, zb[i].normal, 3) < -1e-10) Init(zb[i].normal, (-1)*zb[i].normal[0], (-1)*zb[i].normal[1], (-1)*zb[i].normal[2]);
        reflect(light, zb[i].normal, vec);
        WSetColor(cx((int)((inner(zb[i].normal, light, 3)/2. + 1./2.)*N_COLOR),(int)((inner(vec, res, 3)/2. + 1./2.)*BW_COLOR)));
        
        invmapX(pt, xpt, 3);
        WFillPolygon (xpt, 3);
        }
    }
}
    
void DrawGraph3DX(double a, double b, double (*f[]) (double, double), unsigned long firstColor,
                  unsigned long secondColor, int mode, int modeColor, int n) {
    double dx = EPS_X, dy = EPS_Y;
    double u[3], w[3], tmpx[3], res[3];
    int j = 0;
    pallette(firstColor, secondColor, modeColor, N_COLOR);
    //DrawLinear3D();
    VectorSight();
    farY = a; nearY = b;
    if(key == 0)  InitCameraPosition(nearY + 5, 2*(nearY-farY));
    printf("[DGraphX]:\tcameraPos = %f\n", cameraPos);
    if(mode) DrawAxes(&j); //printf("%d:(\n", j);
    for(int i = 0; i < n; i++) {
        for (double y = a; y < b; y += dy) {
            for (double x = a; x < b; x += dx) {
                zm[j].x = x; zm[j].y = y; zm[j].z = f[i](x, y); zm[j].data = i;
                Init(tmpx, /*cameraPos*view[0] - */x,
                 /*cameraPos*view[1] - */y,
                 /*cameraPos*view[2] - */f[i](x, y));
                 z[j].z = inner(tmpx, view, 3); z[j].num = j; j++;
            }
        }
    } qsort(z, j, sizeof(*z), &comp);
    for (int i = 0; i < j; i++) {
        if(SightControl(zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z, 0.1)) {
            if(zm[z[i].num].data != -1) {
            Projection(&pt[0], zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z);
            Projection(&pt[1], zm[z[i].num].x + EPS_X, zm[z[i].num].y, f[zm[z[i].num].data](zm[z[i].num].x + EPS_X, zm[z[i].num].y));
            
            if(EPS_X*EPS_Y > 1e-2) Projection(&pt[2], zm[z[i].num].x, zm[z[i].num].y + EPS_Y, f[zm[z[i].num].data](zm[z[i].num].x, zm[z[i].num].y + EPS_Y));
            
            else Projection(&pt[2], zm[z[i].num].x + EPS_X, zm[z[i].num].y + EPS_Y, f[zm[z[i].num].data](zm[z[i].num].x + EPS_X, zm[z[i].num].y + EPS_Y));
            
            if(modeColor == 0)  WSetColor(cx((int)((0.5 + atan(COEF * f[zm[z[i].num].data](zm[z[i].num].x + EPS_X/2., zm[z[i].num].y + EPS_Y/2.))/PI)*N_COLOR),0));//[0]
            else {
                Init(u, 0, EPS_Y, (f[zm[z[i].num].data](zm[z[i].num].x, zm[z[i].num].y + EPS_Y) - f[zm[z[i].num].data](zm[z[i].num].x, zm[z[i].num].y)));
                Init(w, EPS_X, 0, (f[zm[z[i].num].data](zm[z[i].num].x + EPS_X, zm[z[i].num].y) - f[zm[z[i].num].data](zm[z[i].num].x, zm[z[i].num].y)));
                cross(u, w, normal);
                normalize(normal, 3);
                InitV(tmp, view);
                Init(view, (cameraPos*view[0] - zm[z[i].num].x), (cameraPos*view[1] - zm[z[i].num].y), (cameraPos*view[2] - zm[z[i].num].z));
                normalize(view, 3);
                if(inner(view, normal, 3) < 0.) Init(normal, (-1)*normal[0], (-1)*normal[1], (-1)*normal[2]);
                reflect(light, normal, res);
                WSetColor(cx((int)((inner(normal, light, 3)/2. + 1./2.)*N_COLOR), (int)((inner(view, res, 3)/2. + 1./2.)*BW_COLOR)));
                InitV(view, tmp);
            }
            if(EPS_X*EPS_Y < 1e-2)  {
                Projection(&pt[3], zm[z[i].num].x, zm[z[i].num].y + EPS_Y, f[zm[z[i].num].data](zm[z[i].num].x, zm[z[i].num].y + EPS_Y));
                invmapX(pt, xpt, 4);
                WFillPolygon (xpt, 4);
            } else {
            invmapX(pt, xpt, 3); //DrawLineX(p1.x,p1.y,p2.x,p2.y,p,q);
            WFillPolygon (xpt, 3);
            
            Projection(&pt[0], zm[z[i].num].x + EPS_X, zm[z[i].num].y + EPS_Y, f[zm[z[i].num].data](zm[z[i].num].x + EPS_X, zm[z[i].num].y + EPS_Y));
            
                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
            }
            }
            else { WSetColor(LIGHTGREEN);  DrawLine(zm[z[i].num].x,zm[z[i].num].y,zm[z[i].num].z,zm[z[i].num].t); }
        }
    }
}

void ParametricCurve3D(double (*curvefun[]) (double), double it, double s, int mode) {
    double dt = 1e-4; 
    Point pts[2]; 
    XPoint ptx[2];
    VectorSight();
    nearY = farY = curvefun[1](it);
    SetLineWidth(3);
    pallette(RED, BLUE, 0, N_COLOR);
    if(key == 0) InitCameraPosition(nearY + 5, 2*(nearY-farY));
    for(double t = it; t < s; t += dt) {
        if(curvefun[1](t) > nearY) nearY = curvefun[1](t);
        if(curvefun[1](t) < farY) farY = curvefun[1](t);
    }
    for(double t = it; t < s; t += dt) {
        if(SightControl(curvefun[0](t), curvefun[1](t), curvefun[2](t), 5)) {
        if(!(mode && LineSightControlAxes(curvefun[0](t), curvefun[1](t), curvefun[2](t)))) {
            Projection(&pts[0], curvefun[0](t), curvefun[1](t), curvefun[2](t));
            Projection(&pts[1], curvefun[0](t+dt), curvefun[1](t+dt), curvefun[2](t+dt));
            WSetColor(cx((int)((0.5 + atan(COEF * curvefun[2](t + dt/2.))/PI)*N_COLOR), 0));
            DrawLineX(pts, ptx);
        }
        }
    }
}

void ParametricGraph3D(double (*parfun[]) (double, double), unsigned long firstColor, unsigned long secondColor,
                       double it, double ft, double is, double fs, int mode, int modeColor) {
    double dt = (ft-it)/150., ds = (fs-is)/100.; 
    double u[3],w[3],tmpx[3], res[3]; 
    int j = 0;
    SetShadingColor(MAGENTA);
    pallette(firstColor, secondColor, modeColor, N_COLOR);
    VectorSight();
    nearY = farY = parfun[1](it, is);
    for(double t = it; t < ft; t += dt) {
        for (double s = is; s < fs; s += ds) {
            if(parfun[1](t,s) > nearY) nearY = parfun[1](t,s);
            if(parfun[1](t,s) < farY) farY = parfun[1](t,s);
        }
    } if(key == 0) InitCameraPosition(nearY + 5, 2*(nearY-farY));
    printf("[DGraphX]:\tcameraPos = %f\n", cameraPos);
    if(mode) DrawAxes(&j);
    for(double t = it; t < ft; t += dt) {
        for (double s = is; s < fs; s += ds) {
            zm[j].x = t; zm[j].y = s;
            Init(tmpx,parfun[0](t,s), parfun[1](t,s), parfun[2](t,s));
            z[j].z = inner(view, tmpx, 3); z[j].num = j; j++;
        }
    } qsort(z, j, sizeof(*z), &comp);
    for (int i = 0; i < j; i++) {
        if(SightControl(parfun[0](zm[z[i].num].x, zm[z[i].num].y), parfun[1](zm[z[i].num].x, zm[z[i].num].y), parfun[2](zm[z[i].num].x, zm[z[i].num].y), 0.1)) {
            if(zm[z[i].num].data != -1) {
            Projection(&pt[0], parfun[0](zm[z[i].num].x, zm[z[i].num].y),
                       parfun[1](zm[z[i].num].x, zm[z[i].num].y),
                       parfun[2](zm[z[i].num].x, zm[z[i].num].y));
            Projection(&pt[1], parfun[0](zm[z[i].num].x + dt, zm[z[i].num].y),
                       parfun[1](zm[z[i].num].x + dt, zm[z[i].num].y),
                       parfun[2](zm[z[i].num].x + dt, zm[z[i].num].y));
            Projection(&pt[2], parfun[0](zm[z[i].num].x, zm[z[i].num].y + ds),
                       parfun[1](zm[z[i].num].x, zm[z[i].num].y + ds),
                       parfun[2](zm[z[i].num].x, zm[z[i].num].y + ds));
            
            if(modeColor == 0) WSetColor(cx((int)((0.5 + atan(COEF * parfun[2](zm[z[i].num].x, zm[z[i].num].y))/PI)*N_COLOR), 0));
            else {
                Init(u, parfun[0](zm[z[i].num].x + dt, zm[z[i].num].y) - parfun[0](zm[z[i].num].x, zm[z[i].num].y),
                     parfun[1](zm[z[i].num].x + dt, zm[z[i].num].y) - parfun[1](zm[z[i].num].x, zm[z[i].num].y),
                     parfun[2](zm[z[i].num].x + dt, zm[z[i].num].y) - parfun[2](zm[z[i].num].x, zm[z[i].num].y));
                Init(w, parfun[0](zm[z[i].num].x, zm[z[i].num].y + ds) - parfun[0](zm[z[i].num].x, zm[z[i].num].y),
                      parfun[1](zm[z[i].num].x, zm[z[i].num].y + ds) - parfun[1](zm[z[i].num].x, zm[z[i].num].y),
                     parfun[2](zm[z[i].num].x, zm[z[i].num].y + ds) - parfun[2](zm[z[i].num].x, zm[z[i].num].y));
                cross(u, w, normal);
                normalize(normal, 3);
                InitV(tmpx, view);
                Init(view, (cameraPos*view[0] - parfun[0](zm[z[i].num].x, zm[z[i].num].y)), (cameraPos*view[1] - parfun[1](zm[z[i].num].x, zm[z[i].num].y)), (cameraPos*view[2] - parfun[2](zm[z[i].num].x, zm[z[i].num].y)));
                normalize(view, 3);
                if(inner(view, normal, 3) < 0.) { intxvec((-1), normal, res); InitV(normal, res); }
                reflect(light, normal, res);
                
                  WSetColor(cx((int)((inner(normal, light, 3)/2. + 1./2.)*N_COLOR), (int)((inner(view, res, 3)/2. + 1./2.)*BW_COLOR)));
                
                InitV(view, tmpx);
            }
                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
            
            Projection(&pt[0], parfun[0](zm[z[i].num].x + dt, zm[z[i].num].y + ds),
                       parfun[1](zm[z[i].num].x + dt, zm[z[i].num].y + ds),
                       parfun[2](zm[z[i].num].x + dt, zm[z[i].num].y + ds));

                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
           }
            else { WSetColor(LIGHTGREEN);  DrawLine(zm[z[i].num].x,zm[z[i].num].y,zm[z[i].num].z,zm[z[i].num].t); }
        }
    }
}


void DrawPolytope(Face *f, int n) {
    double tmpx[3], a = 0, b = 0, c = 0, u[3], w[3]; Point vert[10]; XPoint vertx[10];
    VectorSight();
    SetShadingColor(BLUE);
    pallette(RED, BLUE, 1, N_COLOR);
    if(key == 0) InitCameraPosition(5, 2);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < f[i].n; j++) {
            a += f[i].vertex[0][j]/(f[i].n);
            b += f[i].vertex[1][j]/(f[i].n);
            c += f[i].vertex[2][j]/(f[i].n);
        }
        Init(tmpx, cameraPos*view[0] - a,
             cameraPos*view[1] - b,
             cameraPos*view[2] - c);
        z[i].z = inner(tmpx, tmpx, 3); z[i].num = i;
        a = 0; b = 0; c = 0;
    } a = 0; b = 0; c = 0;
    qsort(z, n, sizeof(*z), &comp);
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < f[i].n; j++) {
            a += f[z[i].num].vertex[0][j]/(f[z[i].num].n);
            b += f[z[i].num].vertex[1][j]/(f[z[i].num].n);
            c += f[z[i].num].vertex[2][j]/(f[z[i].num].n);
        } Init(tmpx, a, b, c);
        if(SightControl(a, b, c, 0.01)) {
        Init(u, f[z[i].num].vertex[0][1] - f[z[i].num].vertex[0][0],
            f[z[i].num].vertex[1][1]-f[z[i].num].vertex[1][0],
                 f[z[i].num].vertex[2][1]-f[z[i].num].vertex[2][0]);
        Init(w, f[z[i].num].vertex[0][f[z[i].num].n-1] - f[z[i].num].vertex[0][0],
            f[z[i].num].vertex[1][f[z[i].num].n-1]-f[z[i].num].vertex[1][0],
                 f[z[i].num].vertex[2][f[z[i].num].n-1]-f[z[i].num].vertex[2][0]);
        cross(u,w,normal);
        if(inner(normal, tmpx, 3) < 0) Init(normal, (-1)*normal[0], (-1)*normal[1], (-1)*normal[2]);
        normalize(normal, 3);
        Init(u, (cameraPos*view[0] -  f[z[i].num].vertex[0][0] +centre[0]), (cameraPos*view[1] - f[z[i].num].vertex[0][1] + centre[1]), (cameraPos*view[2] - f[z[i].num].vertex[0][2])+centre[2]);
        normalize(u,3);
        reflect(light, normal, w);
        WSetColor(cx((int)((inner(normal, light, 3)/2. + 1./2.)*N_COLOR), (int)((inner(w, u, 3)/2. + 1./2.)*BW_COLOR)));
        for (int j = 0; j < f[z[i].num].n; j++)
                Projection(&vert[j], f[z[i].num].vertex[0][j],
                       f[z[i].num].vertex[1][j],
                       f[z[i].num].vertex[2][j]);
        
        invmapX(vert, vertx, f[z[i].num].n);
        WFillPolygon (vertx, f[z[i].num].n);
        } a = 0; b = 0; c = 0;
  }
}
    

//////////////////////////////////////////////////////

/*
   for(int i = 0; i < width; i++) {
       for(int j = 0; j < height; j++) {
           p1 = map(i,j); back = 1;
           for(l = k - 1; l >= 0; l--) {
               if(Inside(zb[l], p1)) {
                   back = 0; break;
               }
           }
           if(back == 0) {
               InvProj(p1.x, p1.y, res);
               if(inner(res, zb[l].normal, 3) < 0.) Init(zb[l].normal, (-1)*zb[l].normal[0], (-1)*zb[l].normal[1], (-1)*zb[l].normal[2]);
               n1 = (int)((inner(zb[l].normal, light, 3)/2. + 1./2.)*N_COLOR);
               reflect(light, zb[l].normal, vec);
               n2 = (int)((inner(vec, res, 3)/2. + 1./2.)*BW_COLOR);
               if(curcolor[0] != n1 || curcolor[1] != n2) {
                   WSetColor(cx(n1, n2));
                   curcolor[0] = n1; curcolor[1] = n2;
               }
           }
           else {
               if(curcolor[0] != 0 || curcolor[1] != 0) {
                    WSetColor(BLACK);
                   curcolor[0] = 0; curcolor[1] = 0;
               }
           }
           WDrawPoint(i,j);
     }
  }*/
