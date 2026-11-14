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

// Left-multiply M by a camera-space axis rotation: M = N*M.
// Camera space: x = screen right, y = view depth, z = screen down.
static void camRotation(int Axis, double angle) {
    double N[9];
    MatrixId(N);
    switch (Axis) {
        case 1:
            n(1,1) = n(2,2) = cos(angle);
            n(1,2) = sin(angle);
            n(2,1) = -sin(angle);
            break;
        case 3:
            n(0,0) = n(1,1) = cos(angle);
            n(0,1) = sin(angle);
            n(1,0) = -sin(angle);
            break;
        default:
            return;
    }
    double c0, c1, c2;
    for (int j = 0; j < 3; j++) {
        c0 = M[j]; c1 = M[3 + j]; c2 = M[6 + j];
        for (int i = 0; i < 3; i++)
            M[3*i + j] = N[3*i]*c0 + N[3*i + 1]*c1 + N[3*i + 2]*c2;
    }
}

// Trackball drag: horizontal -> rotate about screen-up axis,
// vertical -> rotate about screen-right axis.
void RotateDrag(double dx, double dy) {
    camRotation(X, dy*0.005);
    camRotation(Z, -dx*0.005);
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

#define NEAR_EPS 0.1
static double nearW(void) { return cameraPos - screenPos + NEAR_EPS; }
static double mapX(double x) { return x/xc + width*((-xmin)/(xmax - xmin)); }
static double mapY(double y) { return -y/yc + height*(ymax/(ymax - ymin)); }

// World point -> projected scene coords; returns camera-space w (< NEAR_EPS if behind).
static double ProjectW(double x, double y, double z, Point *p) {
    double t[3];
    Init(t, x - centre[0], y - centre[1], z - centre[2]);
    MxV(M, t);
    double w = cameraPos - t[1];
    p->x = t[0]*(cameraPos - screenPos)/w;
    p->y = t[2]*(cameraPos - screenPos)/w;
    return w;
}

// Per-pixel shading context: to-camera direction of the triangle.
typedef struct { double res[3]; } ShadeCtx;

static int shadeIdx(double t, int n) {
    int i = (int)((t/2. + .5)*n);
    if (i < 0) i = 0;
    if (i >= n) i = n - 1;
    return i;
}

// Phong-ish: interpolated normal -> normalize -> diffuse+specular palette lookup.
static unsigned long shadePhong(const float *a, void *ctx_) {
    ShadeCtx *c = ctx_;
    double n[3] = {a[0], a[1], a[2]}, r[3];
    normalize(n, 3);
    if (inner(n, c->res, 3) < 0) intxvec(-1., n, n);
    reflect(light, n, r);
    return cx(shadeIdx(inner(n, light, 3), N_COLOR),
              shadeIdx(inner(r, c->res, 3), BW_COLOR));
}

// Interpolated height -> color-ramp lookup.
static unsigned long shadeHeight(const float *a, void *ctx_) {
    (void)ctx_;
    return cx((int)((0.5 + atan(COEF*a[0])/PI)*N_COLOR), 0);
}

// Flat per-face shading (fallback when phongEnabled is off).
static unsigned long flatShade(double tri[][3], ShadeCtx *sc) {
    double e1[3], e2[3], n[3], r[3];
    minus(tri[1], tri[0], e1); minus(tri[2], tri[0], e2);
    cross(e1, e2, n);
    normalize(n, 3);
    if (inner(n, sc->res, 3) < 0) intxvec(-1., n, n);
    reflect(light, n, r);
    return cx(shadeIdx(inner(n, light, 3), N_COLOR),
              shadeIdx(inner(r, sc->res, 3), BW_COLOR));
}

// Near-plane clip (at the screen plane) + project + rasterize.
// attr: per-vertex float attributes lerped through the clip; NULL = flat color.
#define FB_NA 4
static void rasterTriA(double v[][3], float attr[][FB_NA], int na,
                       unsigned long (*shade)(const float *, void *), void *ctx,
                       unsigned long flatcol) {
    double cv[3][3], w[3], nw = nearW();
    for (int i = 0; i < 3; i++) {
        Init(cv[i], v[i][0] - centre[0], v[i][1] - centre[1], v[i][2] - centre[2]);
        MxV(M, cv[i]);
        w[i] = cameraPos - cv[i][1];
    }
    double px[4], py[4], pw[4];
    float pa[4][FB_NA];
    int np = 0;
    for (int i = 0; i < 3; i++) {
        int j = (i + 1)%3;
        double fi = w[i] - nw, fj = w[j] - nw;
        if (fi > 0) {
            px[np] = cv[i][0]; py[np] = cv[i][2]; pw[np] = w[i];
            for (int k = 0; k < na; k++) pa[np][k] = attr[i][k];
            np++;
        }
        if ((fi > 0) != (fj > 0)) {
            double t = fi/(fi - fj);
            px[np] = cv[i][0] + (cv[j][0] - cv[i][0])*t;
            py[np] = cv[i][2] + (cv[j][2] - cv[i][2])*t;
            pw[np] = w[i] + (w[j] - w[i])*t;
            for (int k = 0; k < na; k++) pa[np][k] = attr[i][k] + (attr[j][k] - attr[i][k])*t;
            np++;
        }
    }
    if (np < 3) return;
    double sx[4], sy[4], sq[4];
    for (int i = 0; i < np; i++) {
        sx[i] = mapX(px[i]*(cameraPos - screenPos)/pw[i]);
        sy[i] = mapY(py[i]*(cameraPos - screenPos)/pw[i]);
        sq[i] = 1./pw[i];
    }
    for (int i = 1; i + 1 < np; i++) {
        if (na)
            FB_TriA(sx[0], sy[0], sq[0], pa[0], sx[i], sy[i], sq[i], pa[i],
                    sx[i + 1], sy[i + 1], sq[i + 1], pa[i + 1], na, shade, ctx);
        else
            FB_Tri(sx[0], sy[0], sq[0], sx[i], sy[i], sq[i], sx[i + 1], sy[i + 1], sq[i + 1], flatcol);
    }
}

static void rasterTri(double v[][3], unsigned long col) {
    rasterTriA(v, NULL, 0, NULL, NULL, col);
}

static float axq[BUF_SIZE][2];
static Point axlbl[3];
static int axlbl_n = 0;

static void drawAxisSegs(int ax) {
    for (int i = 0; i < ax; i++) {
        if (axq[i][0] < 0 || axq[i][1] < 0) continue;
        FB_Seg(mapX(zm[i].x), mapY(zm[i].y), axq[i][0],
               mapX(zm[i].z), mapY(zm[i].t), axq[i][1], LIGHTGREEN);
    }
}

void DrawAxisLabels(void) {
    static const char *names[3] = {"X", "Y", "Z"};
    WSetColor(LIGHTGREEN);
    for (int i = 0; i < axlbl_n; i++)
        WDrawString(names[i], invmap(axlbl[i].x, axlbl[i].y).x, invmap(axlbl[i].x, axlbl[i].y).y);
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
    Point p1, p2;
    double tmpx[3];
    SetLineWidth(0);
    *n = 0;
    axlbl_n = 0;
#define AXSEG(X1, Y1, Z1, X2, Y2, Z2) do { \
    double wq0_ = ProjectW(X1, Y1, Z1, &p1), wq1_ = ProjectW(X2, Y2, Z2, &p2); \
    Init(tmpx, X1, Y1, Z1); \
    z[*n].z = inner(view, tmpx, 3); z[*n].num = *n; \
    zm[*n].x = p1.x; zm[*n].y = p1.y; \
    zm[*n].z = p2.x; zm[*n].t = p2.y; zm[*n].data = -1; \
    axq[*n][0] = wq0_ > nearW() ? 1./wq0_ : -1.f; \
    axq[*n][1] = wq1_ > nearW() ? 1./wq1_ : -1.f; \
    (*n)++; } while(0)
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            AXSEG(5-t, t*cos(K*phi)/2., t*sin(K*phi)/2.,
                  5-t, t*cos(K*(phi+2))/2., t*sin(K*(phi+2))/2.);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            AXSEG(t*cos(K*phi)/2., 5-t, t*sin(K*phi)/2.,
                  t*cos(K*(phi+2))/2., 5-t, t*sin(K*(phi+2))/2.);
        }
    }
    for(double t = 0.3; t > 0; t-=0.05) {
        for (double phi = 0; phi < 360; phi+=2) {
            AXSEG(t*cos(K*phi)/2., t*sin(K*phi)/2., 5-t,
                  t*cos(K*(phi+2))/2., t*sin(K*(phi+2))/2., 5-t);
        }
    }
    for(double t = 0; t < 1-0.03; t += 0.03)
        AXSEG(-5 + 10*t, 0, 0, -5 + 10*(t+0.03), 0, 0);
    for(double t = 0; t < 1-0.03; t += 0.03)
        AXSEG(0, -5 + 10*t, 0, 0, -5 + 10*(t+0.03), 0);
    for(double t = 0; t < 1-0.03; t += 0.03)
        AXSEG(0, 0, -5 + 10*t, 0, 0, -5 + 10*(t+0.03));
#undef AXSEG
    Projection(&axlbl[0], 5.3, 0, 0);
    Projection(&axlbl[1], 0, 5.3, 0);
    Projection(&axlbl[2], 0, 0, 5.3);
    axlbl_n = 3;
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
    double h = 10./(n-1);
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
    if (FB_Begin(BLACK)) {
        ShadeCtx sc;
        double tri[3][3];
        float attr[4][FB_NA];
        for(int i = 0; i < n - 1; i += st) {
            int sy = i + st <= n - 1 ? st : n - 1 - i;
            for(int j = 0; j < n - 1; j += st) {
                int sx = j + st <= n - 1 ? st : n - 1 - j;
                for (int vi = 0; vi < 4; vi++) {
                    int ii = vi >= 2 ? i + sy : i, jj = (vi == 1 || vi == 3) ? j + sx : j;
                    double dfx, dfy;
                    if (jj > 0 && jj < n - 1) dfx = (a[ii*n + jj + 1] - a[ii*n + jj - 1])/(2.*h);
                    else if (jj < n - 1) dfx = (a[ii*n + jj + 1] - a[ii*n + jj])/h;
                    else dfx = (a[ii*n + jj] - a[ii*n + jj - 1])/h;
                    if (ii > 0 && ii < n - 1) dfy = (a[(ii + 1)*n + jj] - a[(ii - 1)*n + jj])/(2.*h);
                    else if (ii < n - 1) dfy = (a[(ii + 1)*n + jj] - a[ii*n + jj])/h;
                    else dfy = (a[ii*n + jj] - a[(ii - 1)*n + jj])/h;
                    attr[vi][0] = -zscale*dfx; attr[vi][1] = -zscale*dfy; attr[vi][2] = 1.;
                }
                Init(tri[0], j*h - 5, i*h - 5, ZS(a[i*n + j]));
                Init(tri[1], (j+sx)*h - 5, i*h - 5, ZS(a[i*n + j + sx]));
                Init(tri[2], (j+sx)*h - 5, (i+sy)*h - 5, ZS(a[(i+sy)*n + j + sx]));
                Init(sc.res, cameraPos*view[0] - tri[0][0] + centre[0], cameraPos*view[1] - tri[0][1] + centre[1], cameraPos*view[2] - tri[0][2] + centre[2]);
                normalize(sc.res, 3);
                float at0[3][FB_NA] = {{attr[0][0], attr[0][1], attr[0][2], 0}, {attr[1][0], attr[1][1], attr[1][2], 0}, {attr[3][0], attr[3][1], attr[3][2], 0}};
                rasterTriA(tri, at0, phongEnabled ? 3 : 0, shadePhong, &sc, flatShade(tri, &sc));
                Init(tri[1], j*h - 5, (i+sy)*h - 5, ZS(a[(i+sy)*n + j]));
                float at1[3][FB_NA] = {{attr[0][0], attr[0][1], attr[0][2], 0}, {attr[2][0], attr[2][1], attr[2][2], 0}, {attr[3][0], attr[3][1], attr[3][2], 0}};
                rasterTriA(tri, at1, phongEnabled ? 3 : 0, shadePhong, &sc, flatShade(tri, &sc));
            }
        }
        drawAxisSegs(ax);
        FB_Flush();
    }
    DrawAxisLabels();
#undef ZS
}

void DrawGraph3DX(double a, double b, double (*f[]) (double, double), unsigned long firstColor,
                  unsigned long secondColor, int mode, int modeColor, int n) {
    double dx = EPS_X, dy = EPS_Y;
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
    if(gv && gm && FB_Begin(BLACK)) {
        ShadeCtx sc;
        double tri[3][3];
        float attr[4][FB_NA], at[3][FB_NA];
        for(int i = 0; i < n; i++) {
            double *gv_ = &gv[i*nf/n], *gm_ = &gm[i*nf/n];
            for (int iy = 0; iy < nyc; iy += st) {
                for (int ix = 0; ix < nxc; ix += st) {
                    int v = iy*nxv + ix;
                    int sx = ix + st <= nxc ? st : nxc - ix;
                    int sy = iy + st <= nyc ? st : nyc - iy;
                    double x = a + ix*dx, y = a + iy*dy;
                    int na;
                    if(modeColor == 0) {
                        na = 1;
                        attr[0][0] = gv_[v]; attr[1][0] = gv_[v + sx];
                        attr[2][0] = gv_[v + sy*nxv]; attr[3][0] = gv_[v + sx + sy*nxv];
                    } else {
                        na = 3;
                        for (int vi = 0; vi < 4; vi++) {
                            int iiy = vi >= 2 ? iy + sy : iy, iix = (vi == 1 || vi == 3) ? ix + sx : ix;
                            int vv = iiy*nxv + iix;
                            double dfx, dfy;
                            if (iix > 0 && iix < nxc) dfx = (gv_[vv + 1] - gv_[vv - 1])/(2.*dx);
                            else if (iix < nxc) dfx = (gv_[vv + 1] - gv_[vv])/dx;
                            else dfx = (gv_[vv] - gv_[vv - 1])/dx;
                            if (iiy > 0 && iiy < nyc) dfy = (gv_[vv + nxv] - gv_[vv - nxv])/(2.*dy);
                            else if (iiy < nyc) dfy = (gv_[vv + nxv] - gv_[vv])/dy;
                            else dfy = (gv_[vv] - gv_[vv - nxv])/dy;
                            attr[vi][0] = -dfx; attr[vi][1] = -dfy; attr[vi][2] = 1.;
                        }
                    }
                    Init(tri[0], x, y, gv_[v]);
                    Init(tri[1], x + sx*dx, y, gv_[v + sx]);
                    Init(tri[2], x + sx*dx, y + sy*dy, gv_[v + sx + sy*nxv]);
                    Init(sc.res, cameraPos*view[0] - x + centre[0], cameraPos*view[1] - y + centre[1], cameraPos*view[2] - gv_[v] + centre[2]);
                    normalize(sc.res, 3);
                    memcpy(at[0], attr[0], sizeof(at[0])); memcpy(at[1], attr[1], sizeof(at[1])); memcpy(at[2], attr[3], sizeof(at[2]));
                    rasterTriA(tri, at, phongEnabled ? na : 0, modeColor == 0 ? shadeHeight : shadePhong, &sc,
                               modeColor == 0 ? cx((int)((0.5 + atan(COEF * gm_[v])/PI)*N_COLOR),0) : flatShade(tri, &sc));
                    Init(tri[1], x, y + sy*dy, gv_[v + sy*nxv]);
                    memcpy(at[1], attr[2], sizeof(at[1]));
                    rasterTriA(tri, at, phongEnabled ? na : 0, modeColor == 0 ? shadeHeight : shadePhong, &sc,
                               modeColor == 0 ? cx((int)((0.5 + atan(COEF * gm_[v])/PI)*N_COLOR),0) : flatShade(tri, &sc));
                }
            }
        }
        drawAxisSegs(j);
        FB_Flush();
    }
    DrawAxisLabels();
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
    static double gc[3] = {0., 0., 0.};
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
        Init(gc, 0., 0., 0.);
        for (int i = 0; i <= NT; i++)
            for (int j = 0; j < TUBE_NR; j++)
                plus(gc, rings[i][j], gc);
        intxvec(1./((NT + 1)*TUBE_NR), gc, gc);
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
    normalMap(zb, k);
    // Signed volume detects winding: >0 = outward normals, <0 = inward.
    double sv = 0, cr[3], a0[3], a1[3], a2[3];
    for (int i = 0; i < k; i++) {
        minus(zb[i].vert[0], gc, a0);
        minus(zb[i].vert[1], gc, a1);
        minus(zb[i].vert[2], gc, a2);
        cross(a1, a2, cr);
        sv += inner(a0, cr, 3);
    }
    int k2 = 0;
    for (int i = 0; i < k; i++) {
        Init(res, cameraPos*view[0] - zb[i].vert[0][0] + centre[0],
                  cameraPos*view[1] - zb[i].vert[0][1] + centre[1],
                  cameraPos*view[2] - zb[i].vert[0][2] + centre[2]);
        normalize(res, 3);
        double d = inner(res, zb[i].normal, 3);
        if((sv > 0) == (d < 0)) continue; // backface
        zb[k2] = zb[i]; k2++;
    }
    if (FB_Begin(BLACK)) {
        for (int i = 0; i < k2; i++) {
            polygon *pl = &zb[i];
            Init(res, (cameraPos*view[0] -  pl->vert[0][0] +centre[0]), (cameraPos*view[1] - pl->vert[0][1] + centre[1]), (cameraPos*view[2] - pl->vert[0][2])+centre[2]);
            normalize(res,3);
            if(inner(res, pl->normal, 3) < -1e-10) Init(pl->normal, (-1)*pl->normal[0], (-1)*pl->normal[1], (-1)*pl->normal[2]);
            reflect(light, pl->normal, vec);
            rasterTri(pl->vert, cx((int)((inner(pl->normal, light, 3)/2. + 1./2.)*N_COLOR),(int)((inner(vec, res, 3)/2. + 1./2.)*BW_COLOR)));
        }
        drawAxisSegs(ax);
        FB_Flush();
    }
    DrawAxisLabels();
}

#define PGT 150
#define PGS 100
static double pgv[(PGT + 1)*(PGS + 1)*3];
static int pge = -1;

void ParametricGraph3D(double (*parfun[]) (double, double), unsigned long firstColor, unsigned long secondColor,
                       double it, double ft, double is, double fs, int mode, int modeColor) {
    double dt = (ft-it)/PGT, ds = (fs-is)/PGS;
    double u[3],w[3],tmpx[3];
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
    if (FB_Begin(BLACK)) {
        ShadeCtx sc;
        double tri[3][3];
        float attr[4][FB_NA], at[3][FB_NA];
        for (int i = 0; i < PGT; i += st) {
            for (int jj = 0; jj < PGS; jj += st) {
                int v = i*(PGS + 1) + jj;
                int sx = i + st <= PGT ? st : PGT - i;
                int sy = jj + st <= PGS ? st : PGS - jj;
                double *p00 = &pgv[3*v];
                int na;
                if(modeColor == 0) {
                    na = 1;
                    attr[0][0] = pgv[3*v + 2]; attr[1][0] = pgv[3*(v + sx*(PGS + 1)) + 2];
                    attr[2][0] = pgv[3*(v + sy) + 2]; attr[3][0] = pgv[3*(v + sx*(PGS + 1) + sy) + 2];
                } else {
                    na = 3;
                    for (int vi = 0; vi < 4; vi++) {
                        int ti = vi == 1 || vi == 3 ? i + sx : i, sj = vi >= 2 ? jj + sy : jj;
                        int tp = ti > 0 ? ti - 1 : ti, tn = ti < PGT ? ti + 1 : ti;
                        int sp = sj > 0 ? sj - 1 : sj, sn = sj < PGS ? sj + 1 : sj;
                        double *pa = &pgv[3*(tn*(PGS + 1) + sj)], *pb = &pgv[3*(tp*(PGS + 1) + sj)];
                        double *pc = &pgv[3*(ti*(PGS + 1) + sn)], *pd = &pgv[3*(ti*(PGS + 1) + sp)];
                        Init(u, pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]);
                        Init(w, pc[0] - pd[0], pc[1] - pd[1], pc[2] - pd[2]);
                        cross(u, w, tmpx);
                        attr[vi][0] = tmpx[0]; attr[vi][1] = tmpx[1]; attr[vi][2] = tmpx[2];
                    }
                }
                double *p10 = &pgv[3*(v + sx*(PGS + 1))], *p01 = &pgv[3*(v + sy)], *p11 = &pgv[3*(v + sx*(PGS + 1) + sy)];
                Init(tri[0], p00[0], p00[1], p00[2]);
                Init(tri[1], p10[0], p10[1], p10[2]);
                Init(tri[2], p11[0], p11[1], p11[2]);
                Init(sc.res, cameraPos*view[0] - p00[0] + centre[0], cameraPos*view[1] - p00[1] + centre[1], cameraPos*view[2] - p00[2] + centre[2]);
                normalize(sc.res, 3);
                memcpy(at[0], attr[0], sizeof(at[0])); memcpy(at[1], attr[1], sizeof(at[1])); memcpy(at[2], attr[3], sizeof(at[2]));
                rasterTriA(tri, at, phongEnabled ? na : 0, modeColor == 0 ? shadeHeight : shadePhong, &sc,
                           modeColor == 0 ? cx((int)((0.5 + atan(COEF * p00[2])/PI)*N_COLOR), 0) : flatShade(tri, &sc));
                Init(tri[1], p01[0], p01[1], p01[2]);
                memcpy(at[1], attr[2], sizeof(at[1]));
                rasterTriA(tri, at, phongEnabled ? na : 0, modeColor == 0 ? shadeHeight : shadePhong, &sc,
                           modeColor == 0 ? cx((int)((0.5 + atan(COEF * p00[2])/PI)*N_COLOR), 0) : flatShade(tri, &sc));
            }
        }
        drawAxisSegs(j);
        FB_Flush();
    }
    DrawAxisLabels();
}

void DrawPolytope(Face *f, int n) {
    double a = 0, b = 0, c = 0, u[3], w[3];
    VectorSight();
    SetShadingColor(BLUE);
    pallette(RED, BLUE, 1, N_COLOR);
    if(key == 0) InitCameraPosition(5, 2);
    if (!FB_Begin(BLACK)) return;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < f[i].n; j++) {
            a += f[i].vertex[0][j]/(f[i].n);
            b += f[i].vertex[1][j]/(f[i].n);
            c += f[i].vertex[2][j]/(f[i].n);
        }
        Init(u, f[i].vertex[0][1] - f[i].vertex[0][0],
            f[i].vertex[1][1]-f[i].vertex[1][0],
                 f[i].vertex[2][1]-f[i].vertex[2][0]);
        Init(w, f[i].vertex[0][f[i].n-1] - f[i].vertex[0][0],
            f[i].vertex[1][f[i].n-1]-f[i].vertex[1][0],
                 f[i].vertex[2][f[i].n-1]-f[i].vertex[2][0]);
        cross(u, w, normal);
        Init(u, cameraPos*view[0] - a, cameraPos*view[1] - b, cameraPos*view[2] - c);
        if(inner(normal, u, 3) < 0)
            Init(normal, (-1)*normal[0], (-1)*normal[1], (-1)*normal[2]);
        normalize(normal, 3);
        Init(u, (cameraPos*view[0] -  f[i].vertex[0][0] +centre[0]), (cameraPos*view[1] - f[i].vertex[0][1] + centre[1]), (cameraPos*view[2] - f[i].vertex[0][2])+centre[2]);
        normalize(u,3);
        reflect(light, normal, w);
        unsigned long col = cx((int)((inner(normal, light, 3)/2. + 1./2.)*N_COLOR), (int)((inner(w, u, 3)/2. + 1./2.)*BW_COLOR));
        for (int j = 1; j + 1 < f[i].n; j++) {
            double tri[3][3];
            Init(tri[0], f[i].vertex[0][0], f[i].vertex[1][0], f[i].vertex[2][0]);
            Init(tri[1], f[i].vertex[0][j], f[i].vertex[1][j], f[i].vertex[2][j]);
            Init(tri[2], f[i].vertex[0][j+1], f[i].vertex[1][j+1], f[i].vertex[2][j+1]);
            rasterTri(tri, col);
        }
        a = 0; b = 0; c = 0;
    }
    FB_Flush();
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
