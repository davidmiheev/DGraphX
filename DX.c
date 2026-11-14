// Creating 2D and 3D graphics
//
// David Miheev 7/18/19
// If there can be an X-axis and a Y-axis, why not a Z-axis?

#include "DGraphX.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
#define TUBE_NR 10


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
static Pair zt[BUF_SIZE];

// 4-pass LSD radix sort on .z; depth keys are quantized to uint32 in-place.
static void radixSortPairs(Pair *arr, int n) {
    double mn = arr[0].z, mx = arr[0].z;
    for (int i = 1; i < n; i++) {
        if (arr[i].z < mn) mn = arr[i].z;
        if (arr[i].z > mx) mx = arr[i].z;
    }
    double sc = (mx - mn) > 1e-12 ? 4294967295./(mx - mn) : 0.;
    for (int i = 0; i < n; i++) arr[i].z = (double)(unsigned)((arr[i].z - mn)*sc);
    Pair *src = arr, *dst = zt, *sw;
    for (int p = 0; p < 4; p++) {
        unsigned cnt[256] = {0}, pos[256], acc = 0;
        int sh = 8*p;
        for (int i = 0; i < n; i++) cnt[((unsigned)src[i].z >> sh) & 255]++;
        for (int c = 0; c < 256; c++) { pos[c] = acc; acc += cnt[c]; }
        for (int i = 0; i < n; i++) dst[pos[((unsigned)src[i].z >> sh) & 255]++] = src[i];
        sw = src; src = dst; dst = sw;
    }
    if (src != arr) memcpy(arr, src, n*sizeof(*arr));
}

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
    double cy = (ymin + ymax)/2;
    switch (mode) {
        case 1:
            xmin = (b + a)/2 - lambda * ((b + a)/2 - xmin);
            xmax = (b + a)/2 + lambda * (-(b + a)/2 + xmax);
            ymin = cy - lambda * (cy - ymin); ymax = cy + lambda * (ymax - cy);
            xc = (xmax - xmin)/width; yc = (ymax - ymin)/height;
            break;
        case -1:
            xmin = (b + a)/2 - (1./lambda)*((b + a)/2 - xmin);
            xmax = (b + a)/2 + (1./lambda)*(-(b + a)/2 + xmax);
            ymin = cy - (1./lambda) * (cy - ymin); ymax = cy + (1./lambda) * (ymax - cy);
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
    if(modeColor == 1) fc.pixel = BLACK; else fc.pixel = firstColor;
    QueryColor(&fc);
    if(modeColor == 1) sc.pixel = SHADINGCOLOR; else sc.pixel = secondColor;
    QueryColor(&sc);
    for (int i = 0; i < size; i++) {
        double br = fc.red + (sc.red - fc.red)*i/size;
        double bg = fc.green + (sc.green - fc.green)*i/size;
        double bb = fc.blue + (sc.blue - fc.blue)*i/size;
        for (int j = 0; j < BW_COLOR; j++) {
            double t = j/(double)BW_COLOR;
            double kr = modeColor == 1 ? (w.red - br)*t : (modeColor == 2 ? -0.85*br*t : 0.);
            double kg = modeColor == 1 ? (w.green - bg)*t : (modeColor == 2 ? -0.85*bg*t : 0.);
            double kb = modeColor == 1 ? (w.blue - bb)*t : (modeColor == 2 ? -0.85*bb*t : 0.);
            WGetColor(0.65*(br + kr), 0.65*(bg + kg), 0.65*(bb + kb), &tmpx);
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

void RotateX_Angle(double angle) {
    elemRotation(X,angle,1);
}

void RotateY_Angle(double angle) {
    elemRotation(Y,angle,1);
}

void RotateZ_Angle(double angle) {
    elemRotation(Z,angle,1);
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


void DrawLinear3D(double* a, int n, unsigned long Color, int mode) {
    double vec[3], res[3], h = 10./(n-1); int k = 0;
    double zmin = a[0], zmax = a[0];
    for (int i = 0; i < n*n; i++) {
        if (a[i] < zmin) zmin = a[i];
        if (a[i] > zmax) zmax = a[i];
    }
    double zmid = (zmin + zmax)/2, zscale = (zmax - zmin) > 1e-12 ? 5./(zmax - zmin) : 1;
#define ZS(v) (((v) - zmid)*zscale)
    pallette(Color, BLACK, 1, N_COLOR);
    VectorSight();
    if(key == 0)  InitCameraPosition(10, 20);
    int ax = 0;
    if(mode) DrawAxes(&ax);
    int st = interactive ? 2 : 1;
    for(int i = 0; i < n - 1; i += st) {
        int sy = i + st <= n - 1 ? st : n - 1 - i;
        for(int j = 0; j < n - 1; j += st) { //zb[k].n = 3;
            int sx = j + st <= n - 1 ? st : n - 1 - j;
            Init(zb[k].vert[0], j*h - 5, i*h - 5, ZS(a[i*n + j]));
            Init(zb[k].vert[1], (j+sx)*h - 5, i*h - 5, ZS(a[i*n + j+sx]));
            Init(zb[k].vert[2], j*h - 5, (i+sy)*h - 5, ZS(a[(i+sy)*n + j]));
            plus(zb[k].vert[0], zb[k].vert[1], vec);
            plus(vec, zb[k].vert[2], vec);
            intxvec(1./3., vec, vec);
            minus(vec, centre, vec);
            z[ax + k].z = inner(view, vec, 3); z[ax + k].num = ~k;
            Init(vec, sx*h, 0, (a[i*n + j + sx] - a[i*n + j])*zscale);
            Init(res, 0, sy*h, (a[(i+sy)*n + j] - a[i*n + j])*zscale);
            cross(vec, res, zb[k].normal); //zb[k + 1].n = 3;
            Init(zb[k+1].vert[0], j*h - 5, (i+sy)*h - 5, ZS(a[(i+sy)*n + j]));
            Init(zb[k+1].vert[1], (j+sx)*h - 5, i*h - 5, ZS(a[i*n + j+sx]));
            Init(zb[k+1].vert[2], (j+sx)*h - 5, (i+sy)*h - 5, ZS(a[(i+sy)*n + j+sx]));
            plus(zb[k+1].vert[0], zb[k+1].vert[1], vec);
            plus(vec, zb[k+1].vert[2], vec);
            intxvec(1./3., vec, vec);
            minus(vec, centre, vec);
            z[ax + k + 1].z = inner(view, vec,3); z[ax + k + 1].num = ~(k + 1);
            Init(vec, sx*h, 0, (a[(i+sy)*n + j + sx] - a[(i+sy)*n + j])*zscale);
            Init(res, 0, sy*h, (a[(i+sy)*n + j + sx] - a[i*n + j + sx])*zscale);
            cross(vec, res, zb[k+1].normal);
            k+=2;
        }
    }
    normalMap(zb, k); //WSetColor(RED);
    radixSortPairs(z, ax + k);
    for (int i = 0; i < ax + k; i++) {
        if (z[i].num >= 0) {
            WSetColor(LIGHTGREEN);
            DrawLine(zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z, zm[z[i].num].t);
            continue;
        }
        polygon *pl = &zb[~z[i].num];
        if(SightControl(pl->vert[0][0], pl->vert[0][1], pl->vert[0][2], 0.1)) {
        Projection(&pt[0], pl->vert[0][0], pl->vert[0][1], pl->vert[0][2]);
        Projection(&pt[1], pl->vert[1][0], pl->vert[1][1], pl->vert[1][2]);
        Projection(&pt[2], pl->vert[2][0], pl->vert[2][1], pl->vert[2][2]);
        Init(res, (cameraPos*view[0] -  pl->vert[0][0] +centre[0]), (cameraPos*view[1] - pl->vert[0][1] + centre[1]), (cameraPos*view[2] - pl->vert[0][2])+centre[2]);
        normalize(res,3);
        if(inner(res, pl->normal, 3) < -1e-10) Init(pl->normal, (-1)*pl->normal[0], (-1)*pl->normal[1], (-1)*pl->normal[2]);
        reflect(light, pl->normal, vec);
        WSetColor(cx((int)((inner(pl->normal, light, 3)/2. + 1./2.)*N_COLOR),(int)((inner(vec, res, 3)/2. + 1./2.)*BW_COLOR)));

        invmapX(pt, xpt, 3);
        WFillPolygon (xpt, 3);
        }
    }
#undef ZS
}

void DrawGraph3DX(double a, double b, double (*f[]) (double, double), unsigned long firstColor,
                  unsigned long secondColor, int mode, int modeColor, int n) {
    double dx = EPS_X, dy = EPS_Y;
    double u[3], w[3], tmpx[3], res[3];
    int j = 0;
    static double *gv = NULL, *gm = NULL;
    static size_t gcap = 0;
    static int ge = -1, gn = -1, gnxc = -1;
    static double ga = 1e30, gb = 1e30;
    int nxc = (int)((b - a)/dx + 0.5), nyc = (int)((b - a)/dy + 0.5), nxv = nxc + 1, nf = n*nxv*(nyc + 1);
    int st = interactive ? 2 : 1;
    if(ge != renderEpoch || gn != n || gnxc != nxc || ga != a || gb != b) {
        if((size_t)nf > gcap) {
            free(gv); free(gm); gcap = 0;
            gv = malloc(nf*sizeof(*gv));
            gm = malloc(nf*sizeof(*gm));
            if(gv && gm) gcap = nf;
        }
        if(gv && gm) {
            for (int i = 0; i < n; i++)
                for (int iy = 0; iy <= nyc; iy++)
                    for (int ix = 0; ix <= nxc; ix++)
                        gv[i*nxv*(nyc + 1) + iy*nxv + ix] = f[i](a + ix*dx, a + iy*dy);
            for (int i = 0; i < n; i++)
                for (int iy = 0; iy < nyc; iy++)
                    for (int ix = 0; ix < nxc; ix++)
                        gm[i*nxv*(nyc + 1) + iy*nxv + ix] = f[i](a + (ix + .5)*dx, a + (iy + .5)*dy);
        }
        ge = renderEpoch; gn = n; gnxc = nxc; ga = a; gb = b;
    }
    pallette(firstColor, secondColor, modeColor, N_COLOR);
    VectorSight();
    farY = a; nearY = b;
    if(key == 0)  InitCameraPosition(nearY + 5, 2*(nearY-farY));
    if(mode) DrawAxes(&j); //printf("%d:(\n", j);
    if(gv && gm)
    for(int i = 0; i < n; i++) {
        for (int iy = 0; iy < nyc; iy += st) {
            for (int ix = 0; ix < nxc; ix += st) {
                int vi = iy*nxv + ix;
                double x = a + ix*dx, y = a + iy*dy;
                zm[j].x = x; zm[j].y = y; zm[j].z = gv[i*nf/n + vi];
                zm[j].data = i; zm[j].t = vi;
                Init(tmpx, x, y, zm[j].z);
                 z[j].z = inner(tmpx, view, 3); z[j].num = j; j++;
            }
        }
    } radixSortPairs(z, j);
    for (int i = 0; i < j; i++) {
        if(zm[z[i].num].data == -1) {
            WSetColor(LIGHTGREEN);  DrawLine(zm[z[i].num].x,zm[z[i].num].y,zm[z[i].num].z,zm[z[i].num].t);
            continue;
        }
        int v = (int)zm[z[i].num].t, d = zm[z[i].num].data;
        int cix = v % nxv, ciy = v / nxv;
        int sx = cix + st <= nxc ? st : nxc - cix;
        int sy = ciy + st <= nyc ? st : nyc - ciy;
        double *gv_ = &gv[d*nf/n], *gm_ = &gm[d*nf/n];
        if(SightControl(zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z, 0.1)) {
            Projection(&pt[0], zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z);
            Projection(&pt[1], zm[z[i].num].x + sx*dx, zm[z[i].num].y, gv_[v + sx]);

            if(EPS_X*EPS_Y > 1e-2) Projection(&pt[2], zm[z[i].num].x, zm[z[i].num].y + sy*dy, gv_[v + sy*nxv]);

            else Projection(&pt[2], zm[z[i].num].x + sx*dx, zm[z[i].num].y + sy*dy, gv_[v + sx + sy*nxv]);

            if(modeColor == 0)  WSetColor(cx((int)((0.5 + atan(COEF * gm_[v])/PI)*N_COLOR),0));//[0]
            else {
                Init(u, 0, sy*EPS_Y, (gv_[v + sy*nxv] - zm[z[i].num].z));
                Init(w, sx*EPS_X, 0, (gv_[v + sx] - zm[z[i].num].z));
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
                Projection(&pt[3], zm[z[i].num].x, zm[z[i].num].y + sy*dy, gv_[v + sy*nxv]);
                invmapX(pt, xpt, 4);
                WFillPolygon (xpt, 4);
            } else {
            invmapX(pt, xpt, 3); //DrawLineX(p1.x,p1.y,p2.x,p2.y,p,q);
            WFillPolygon (xpt, 3);

            Projection(&pt[0], zm[z[i].num].x + sx*dx, zm[z[i].num].y + sy*dy, gv_[v + sx + sy*nxv]);

                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
            }
        }
    }
}

static void tubeTangent(double (*curvefun[]) (double), double t, double ftol, double P[3], double T[3]) {
    Init(P, curvefun[0](t), curvefun[1](t), curvefun[2](t));
    Init(T, curvefun[0](t+ftol) - curvefun[0](t-ftol), curvefun[1](t+ftol) - curvefun[1](t-ftol), curvefun[2](t+ftol) - curvefun[2](t-ftol));
    if(inner(T, T, 3) < 1e-18) Init(T, 1., 0., 0.);
    normalize(T, 3);
}

static void tubeRingAt(double P[3], double n1[3], double n2[3], double r, double ring[][3]) {
    for (int j = 0; j < TUBE_NR; j++)
        Init(ring[j], P[0] + r*(cos(2.*PI*j/TUBE_NR)*n1[0] + sin(2.*PI*j/TUBE_NR)*n2[0]),
             P[1] + r*(cos(2.*PI*j/TUBE_NR)*n1[1] + sin(2.*PI*j/TUBE_NR)*n2[1]),
             P[2] + r*(cos(2.*PI*j/TUBE_NR)*n1[2] + sin(2.*PI*j/TUBE_NR)*n2[2]));
}

// Rotate v by the rotation taking unit tangent a to unit tangent b (Rodrigues).
static void transportFrame(double a[3], double b[3], double v[3]) {
    double ax[3], cr[3], res[3];
    cross(a, b, ax);
    double s2 = inner(ax, ax, 3);
    if(s2 < 1e-18) return;
    double sn = sqrt(s2), c = inner(a, b, 3), d;
    intxvec(1./sn, ax, ax);
    cross(ax, v, cr);
    d = inner(ax, v, 3);
    Init(res, v[0]*c + cr[0]*sn + ax[0]*d*(1. - c),
         v[1]*c + cr[1]*sn + ax[1]*d*(1. - c),
         v[2]*c + cr[2]*sn + ax[2]*d*(1. - c));
    InitV(v, res);
}

static void tubeTri(int k, double v0[3], double v1[3], double v2[3]) {
    double vec[3], e1[3], e2[3];
    Init(zb[k].vert[0], v0[0], v0[1], v0[2]);
    Init(zb[k].vert[1], v1[0], v1[1], v1[2]);
    Init(zb[k].vert[2], v2[0], v2[1], v2[2]);
    plus(v0, v1, vec); plus(vec, v2, vec); intxvec(1./3., vec, vec);
    minus(vec, centre, vec);
    zb[k].z = inner(view, vec, 3);
    minus(v1, v0, e1); minus(v2, v0, e2);
    cross(e1, e2, zb[k].normal);
}

void ParametricCurve3D(double (*curvefun[]) (double), double it, double s, int mode) {
    const int NT = 400;
    static double rings[401][TUBE_NR][3];
    static int tre = -1;
    static double r = 0.05;
    double dt = (s - it)/NT, ftol = dt*0.5;
    double vec[3], res[3];
    int k = 0, ax = 0;
    SetShadingColor(BLUE);
    SetCentre(0, 0, 0);
    pallette(BLUE, BLACK, 1, N_COLOR);
    VectorSight();
    if(tre != renderEpoch) {
        nearY = farY = curvefun[1](it);
        for(double t = it; t < s; t += (s - it)/2000.) {
            if(curvefun[1](t) > nearY) nearY = curvefun[1](t);
            if(curvefun[1](t) < farY) farY = curvefun[1](t);
        }
        r = (nearY - farY)/40.;
        if(r < 0.05) r = 0.05;
        double P[3], T[3], Tp[3], ref[3], n1[3], n2[3];
        tubeTangent(curvefun, it, ftol, P, Tp);
        Init(ref, 0., 1., 0.);
        cross(Tp, ref, n1);
        if(inner(n1, n1, 3) < 1e-12) { Init(ref, 1., 0., 0.); cross(Tp, ref, n1); }
        normalize(n1, 3);
        cross(Tp, n1, n2);
        normalize(n2, 3);
        tubeRingAt(P, n1, n2, r, rings[0]);
        for (int i = 1; i <= NT; i++) {
            tubeTangent(curvefun, it + i*dt, ftol, P, T);
            transportFrame(Tp, T, n1);
            normalize(n1, 3);
            cross(T, n1, n2);
            normalize(n2, 3);
            tubeRingAt(P, n1, n2, r, rings[i]);
            InitV(Tp, T);
        }
        tre = renderEpoch;
    }
    if(key == 0) InitCameraPosition(nearY + 5, 2*(nearY-farY));
    if(mode) DrawAxes(&ax);
    int st = interactive ? 2 : 1;
    for (int i = st; i <= NT; i += st) {
        for (int j = 0; j < TUBE_NR; j++) {
            int jn = (j + 1) % TUBE_NR;
            tubeTri(k, rings[i-st][j], rings[i-st][jn], rings[i][j]);
            tubeTri(k + 1, rings[i][j], rings[i-st][jn], rings[i][jn]);
            k += 2;
        }
    }
    for (int i = 0; i < k; i++) { z[ax + i].z = zb[i].z; z[ax + i].num = ~i; }
    normalMap(zb, k);
    radixSortPairs(z, ax + k);
    for (int i = 0; i < ax + k; i++) {
        if (z[i].num >= 0) {
            WSetColor(LIGHTGREEN);
            DrawLine(zm[z[i].num].x, zm[z[i].num].y, zm[z[i].num].z, zm[z[i].num].t);
            continue;
        }
        polygon *pl = &zb[~z[i].num];
        if(SightControl(pl->vert[0][0], pl->vert[0][1], pl->vert[0][2], 0.1)) {
        Projection(&pt[0], pl->vert[0][0], pl->vert[0][1], pl->vert[0][2]);
        Projection(&pt[1], pl->vert[1][0], pl->vert[1][1], pl->vert[1][2]);
        Projection(&pt[2], pl->vert[2][0], pl->vert[2][1], pl->vert[2][2]);
        Init(res, (cameraPos*view[0] -  pl->vert[0][0] +centre[0]), (cameraPos*view[1] - pl->vert[0][1] + centre[1]), (cameraPos*view[2] - pl->vert[0][2])+centre[2]);
        normalize(res,3);
        if(inner(res, pl->normal, 3) < -1e-10) Init(pl->normal, (-1)*pl->normal[0], (-1)*pl->normal[1], (-1)*pl->normal[2]);
        reflect(light, pl->normal, vec);
        WSetColor(cx((int)((inner(pl->normal, light, 3)/2. + 1./2.)*N_COLOR),(int)((inner(vec, res, 3)/2. + 1./2.)*BW_COLOR)));
        invmapX(pt, xpt, 3);
        WFillPolygon (xpt, 3);
        }
    }
}

#define PGT 150
#define PGS 100
static double pgv[(PGT + 1)*(PGS + 1)*3];
static int pge = -1;

void ParametricGraph3D(double (*parfun[]) (double, double), unsigned long firstColor, unsigned long secondColor,
                       double it, double ft, double is, double fs, int mode, int modeColor) {
    double dt = (ft-it)/PGT, ds = (fs-is)/PGS;
    double u[3],w[3],tmpx[3], res[3];
    int j = 0;
    SetShadingColor(MAGENTA);
    pallette(firstColor, secondColor, modeColor, N_COLOR);
    VectorSight();
    if(pge != renderEpoch) {
        int idx = 0;
        nearY = farY = parfun[1](it, is);
        for (int i = 0; i <= PGT; i++) {
            double t = it + i*dt;
            for (int jj = 0; jj <= PGS; jj++) {
                double s = is + jj*ds;
                pgv[idx] = parfun[0](t, s);
                pgv[idx+1] = parfun[1](t, s);
                pgv[idx+2] = parfun[2](t, s);
                if(pgv[idx+1] > nearY) nearY = pgv[idx+1];
                if(pgv[idx+1] < farY) farY = pgv[idx+1];
                idx += 3;
            }
        }
        pge = renderEpoch;
    }
    if(key == 0) InitCameraPosition(nearY + 5, 2*(nearY-farY));
    if(mode) DrawAxes(&j);
    int st = interactive ? 2 : 1;
    for (int i = 0; i < PGT; i += st) {
        for (int jj = 0; jj < PGS; jj += st) {
            int v = i*(PGS + 1) + jj;
            zm[j].x = it + i*dt; zm[j].y = is + jj*ds; zm[j].data = 0; zm[j].t = v;
            Init(tmpx, pgv[3*v], pgv[3*v+1], pgv[3*v+2]);
            z[j].z = inner(view, tmpx, 3); z[j].num = j; j++;
        }
    } radixSortPairs(z, j);
    for (int i = 0; i < j; i++) {
        if(zm[z[i].num].data == -1) {
            WSetColor(LIGHTGREEN);  DrawLine(zm[z[i].num].x,zm[z[i].num].y,zm[z[i].num].z,zm[z[i].num].t);
            continue;
        }
        int v = (int)zm[z[i].num].t;
        int ci = v/(PGS + 1), cj = v%(PGS + 1);
        int sx = ci + st <= PGT ? st : PGT - ci;
        int sy = cj + st <= PGS ? st : PGS - cj;
        double *p00 = &pgv[3*v], *p10 = &pgv[3*(v + sx*(PGS + 1))], *p01 = &pgv[3*(v + sy)], *p11 = &pgv[3*(v + sx*(PGS + 1) + sy)];
        if(SightControl(p00[0], p00[1], p00[2], 0.1)) {
            Projection(&pt[0], p00[0], p00[1], p00[2]);
            Projection(&pt[1], p10[0], p10[1], p10[2]);
            Projection(&pt[2], p01[0], p01[1], p01[2]);

            if(modeColor == 0) WSetColor(cx((int)((0.5 + atan(COEF * p00[2])/PI)*N_COLOR), 0));
            else {
                Init(u, p10[0] - p00[0], p10[1] - p00[1], p10[2] - p00[2]);
                Init(w, p01[0] - p00[0], p01[1] - p00[1], p01[2] - p00[2]);
                cross(u, w, normal);
                normalize(normal, 3);
                InitV(tmpx, view);
                Init(view, (cameraPos*view[0] - p00[0]), (cameraPos*view[1] - p00[1]), (cameraPos*view[2] - p00[2]));
                normalize(view, 3);
                if(inner(view, normal, 3) < 0.) { intxvec((-1), normal, res); InitV(normal, res); }
                reflect(light, normal, res);

                  WSetColor(cx((int)((inner(normal, light, 3)/2. + 1./2.)*N_COLOR), (int)((inner(view, res, 3)/2. + 1./2.)*BW_COLOR)));

                InitV(view, tmpx);
            }
                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);

            Projection(&pt[0], p11[0], p11[1], p11[2]);

                invmapX(pt, xpt, 3);
                WFillPolygon (xpt, 3);
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
    radixSortPairs(z, n);
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
