// DGraphX -- free grapher, version 0.1
//
// David Miheev 7/18/19

#include "DGraphX.h"
#include "Interp2D.h"
#include "matrix.h"
#include "tinyexpr.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define m(ix,jx) m[4*(ix)+jx]

#define NMAX 100 
#define STEP 0.05
#define LAMBDA 1.1
#define SHIFT 0.01
#define EPS 1e-3
#define N 6
#define K  0.01745329251
#define PI 3.14159265358

static int variant = 0;
static int key_1 = 0;
static int key_2 = 0;
static int key_3 = 0;
static int key_4 = 0;
static int key_5 = 0;
static int key_6 = 0;
static int key_7 = 3;
static int key_8 = 1;
static int key_9 = 0;
static int key_10 = 1;
static int key_11 = 1;
static int key_12 = 1;
static int ngrid = 10;
static int iter = 0;
static int niter;
static int max_iter = 1000;
static int mode = 1;
static int dim = 3;
static int count = 0;
static int count1 = 0;
static double m[4*NMAX];
//static double d[NMAX];
static Point pt[NMAX];
static double fval[100000];
static clock_t t;

// TinyExpr integration
static char expr_f[256] = "10./x/x-9./x";
static char expr_u[256] = "3*cos(3*sqrt(x*x+y*y))*exp(-sqrt(x*x+y*y)/3.)";
static char expr_v[256] = "3*exp(-x*x/2.-2*y*y)";
static char expr_w[256] = "sin(x*y)";
static char expr_cx[256] = "20*cos(t)^2";
static char expr_cy[256] = "(10./2.)*sin(t)";
static char expr_cz[256] = "10*sin(2*t)";
static char expr_px[256] = "3*(3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*cos(t)";
static char expr_py[256] = "3*(3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*sin(t)";
static char expr_pz[256] = "3*sin(t/2.)*sin(s)+3*cos(t/2.)*sin(2*s)";

static double cur_x, cur_y, cur_t, cur_s;
static te_expr *te_f = NULL, *te_u = NULL, *te_v = NULL, *te_w = NULL;
static te_expr *te_cx = NULL, *te_cy = NULL, *te_cz = NULL;
static te_expr *te_px = NULL, *te_py = NULL, *te_pz = NULL;

static void trim_newline(char *s) {
    int l = strlen(s);
    if (l > 0 && s[l-1] == '\n') s[l-1] = '\0';
}

int compile_functions() {
    int err;
    te_variable v1[] = {{"x", &cur_x, 0, 0}};
    te_variable v2[] = {{"x", &cur_x, 0, 0}, {"y", &cur_y, 0, 0}};
    te_variable v_t[] = {{"t", &cur_t, 0, 0}};
    te_variable v_ts[] = {{"t", &cur_t, 0, 0}, {"s", &cur_s, 0, 0}};

    int total_err = 0;
#define COMPILE(expr, vars, count, target) \
    { te_expr *tmp = te_compile(expr, vars, count, &err); \
      if (tmp) { if (target) te_free(target); target = tmp; } \
      else { printf("[DGraphX] TinyExpr error in '%s' at position %d\n", expr, err); total_err++; } }

    COMPILE(expr_f, v1, 1, te_f);
    COMPILE(expr_u, v2, 2, te_u);
    COMPILE(expr_v, v2, 2, te_v);
    COMPILE(expr_w, v2, 2, te_w);
    COMPILE(expr_cx, v_t, 1, te_cx);
    COMPILE(expr_cy, v_t, 1, te_cy);
    COMPILE(expr_cz, v_t, 1, te_cz);
    COMPILE(expr_px, v_ts, 2, te_px);
    COMPILE(expr_py, v_ts, 2, te_py);
    COMPILE(expr_pz, v_ts, 2, te_pz);
#undef COMPILE
    return total_err == 0;
}

// ================ functions for drawing (input data) ===================

double f(double x) {
    if (!te_f) compile_functions();
    cur_x = x;
    double res = te_f ? te_eval(te_f) : 0;
    return !isfinite(res) ? 0 : res;
}

double df(double x) {
    const double h = 1e-6;
    return (f(x + h) - f(x - h)) / (2*h);
}

double parx2D(double t) {
    if (!te_cx) compile_functions();
    cur_t = t;
    double res = te_cx ? te_eval(te_cx) : 0;
    return !isfinite(res) ? 0 : res;
}

double pary2D(double t) {
    if (!te_cy) compile_functions();
    cur_t = t;
    double res = te_cy ? te_eval(te_cy) : 0;
    return !isfinite(res) ? 0 : res;
}

double u(double x, double y) {
    if (!te_u) compile_functions();
    cur_x = x; cur_y = y;
    double res = te_u ? te_eval(te_u) : 0;
    return !isfinite(res) ? 0 : res;
}

double w(double x, double y) {
    if (!te_w) compile_functions();
    cur_x = x; cur_y = y;
    double res = te_w ? te_eval(te_w) : 0;
    return !isfinite(res) ? 0 : res;
}

double v(double x, double y) {
    if (!te_v) compile_functions();
    cur_x = x; cur_y = y;
    double res = te_v ? te_eval(te_v) : 0;
    return !isfinite(res) ? 0 : res;
}

double curvex(double t) {
    if (!te_cx) compile_functions();
    cur_t = t;
    double res = te_cx ? te_eval(te_cx) : 0;
    return !isfinite(res) ? 0 : res;
}

double curvey(double t) {
    if (!te_cy) compile_functions();
    cur_t = t;
    double res = te_cy ? te_eval(te_cy) : 0;
    return !isfinite(res) ? 0 : res;
}

double curvez(double t) {
    if (!te_cz) compile_functions();
    cur_t = t;
    double res = te_cz ? te_eval(te_cz) : 0;
    return !isfinite(res) ? 0 : res;
}

double parx(double t, double s) {
    if (!te_px) compile_functions();
    cur_t = t; cur_s = s;
    double res = te_px ? te_eval(te_px) : 0;
    return !isfinite(res) ? 0 : res;
}

double pary(double t, double s) {
    if (!te_py) compile_functions();
    cur_t = t; cur_s = s;
    double res = te_py ? te_eval(te_py) : 0;
    return !isfinite(res) ? 0 : res;
}

double parz(double t, double s) {
    if (!te_pz) compile_functions();
    cur_t = t; cur_s = s;
    double res = te_pz ? te_eval(te_pz) : 0;
    return !isfinite(res) ? 0 : res;
}

// Mouse interaction state
static int is_dragging = 0;
static int last_mouse_x, last_mouse_y;
static int active_button;

int MouseEventHandler(int type, int x, int y, int button) {
    if (type == ButtonPress) {
        if (button == Button4 || button == Button5) {
            if (dim == 2) {
                double cx = (xmin + xmax)/2;
                scale(button == Button4 ? -1 : 1, cx, cx, 1.1);
            } else {
                if (count == 0) { InitialPosition(); count = 1; }
                ChangeCameraPosition(button == Button4 ? -0.5 : 0.5);
            }
            return KEY_PRESS_EXPOSE;
        }
        is_dragging = 1;
        last_mouse_x = x;
        last_mouse_y = y;
        active_button = button;
        return KEY_PRESS_NOTHING;
    } else if (type == ButtonRelease) {
        is_dragging = 0;
        return KEY_PRESS_NOTHING;
    } else if (type == MotionNotify) {
        if (!is_dragging) return KEY_PRESS_NOTHING;
        int dx = x - last_mouse_x;
        int dy = y - last_mouse_y;
        last_mouse_x = x;
        last_mouse_y = y;

        if (dim == 3) {
            if (count == 0) {
                InitialPosition();
                count = 1;
            }
            if (active_button == Button1) {
                RotateY_Angle(dx * 0.005);
                RotateX_Angle(dy * 0.005);
                return KEY_PRESS_EXPOSE;
            } else if (active_button == Button3) {
                ChangeCameraPosition(dy * 0.1);
                return KEY_PRESS_EXPOSE;
            }
        } else { // dim == 2
            if (active_button == Button1) {
                if (dx != 0) xshift(dx > 0 ? -1 : 1, 0.002 * abs(dx));
                if (dy != 0) yshift(dy > 0 ? 1 : -1, 0.002 * abs(dy));
                return KEY_PRESS_EXPOSE;
            } else if (active_button == Button3) {
                scale(dy > 0 ? -1 : 1, (xmin+xmax)/2, (xmin+xmax)/2, 1.0 + abs(dy)*0.01);
                return KEY_PRESS_EXPOSE;
            }
        }
    }
    return KEY_PRESS_NOTHING;
}

double polynom(double y,double *x,double *m,int i) { // for interpolation
    return m(i,0) + m(i,1)*(y-x[i]) + m(i,2)*(y-x[i])*(y-x[i]) + m(i,3)*(y-x[i])*(y-x[i])*(y-x[i]);
}

double solution(int i, int j, int n) {
    return sin(2*PI*j/(n-1))*sin(2*PI*i/(n-1))/2;
}
void init_val() {
    for(int i = 0; i < ngrid; i++) {
        for(int j = 0; j < ngrid; j++) {
            if(i == 0 || i == ngrid-1 || j == 0 || j == ngrid-1)
                fval[i*ngrid + j] = sin(2*PI*j/(ngrid-1))*sin(2*PI*i/(ngrid-1))/2;
            else
                fval[i*ngrid + j] = 0;
        }
    }
}

//==========================================================================

static void DrawWindowContent (double* x, int* n, double a, double b) {
    
    char s[128], str[128];
    int j;
    double tmp;
    static Face faces[100];
    static int nFaces;
    FILE *fi;
    double stp, (*fun[NMAX])(double, double), (*curvefun[3])(double), (*parfun[3])(double, double), (*parf[2])(double);
    
    double dx = EPS;
    
    SetFont(HELVETICA16);
    WSetColor (BLACK);
    WFillRectangle (0, 0, width, height);
    if (dim == 2) {
    if(variant == 0) {
        SetLineWidth(1);
        WSetColor(MAGENTA); 
        WDrawString ("DGraphX", 300, 169);
        WSetColor(LIGHTGREEN);
        for(int i = 0; i < 9; i++) DrawArc(193 + 3*i, 170, 50 + i*(20 + 2*i), 50 + i*(20 + 2*i), 0, 360*64 );
        WDrawString("This is 2D mode", 100, 500);
        SetFont(HELVETICA12);
        WDrawString("' Ctrl ' -- change dimension; ' 1 ' -- draw interpolation with Hermite polynomial; ' 2 ' -- draw interpolation with spline",100, 550);
        WDrawString("' 3 ' -- draw graph; ' 4 ' -- draw parametric curve (x(t), y(t))", 100, 575);
        WDrawString("' i ' -- enter new function expression in terminal", 100, 600);
    } else {
        switch(mode) {
            case 1:
                for (int i = 0; i < *n; i++) x[i] = a + i*(b-a)/(*n-1);
                break;
            case 3:
                for (int i = *n - 1; i >= 0; i--)
                    x[*n - i - 1] = a/2 + b/2 + (b-a)*cos((2*i + 1)*acos(0)/(*n))/2;
                break;
            default:
                break;
        }
    if(key_1 == 1) { if(*n<NMAX) {
        switch(mode) {
            case 1:
                (*n)+=1; for (int i = 0; i < *n; i++) x[i] = a + i*(b-a)/(*n-1);
                break;
            case 2:
                (*n)+=1; tmp = EPS + (double)(rand()%((int)(100*(b-a))))/100 + a;
                for(j = 0; j < *n - 1 && fabs(x[j] - tmp) > EPS; j++);
                if(j<*n-2 || fabs(x[*n-2]-tmp) < EPS) {
                    for (int i = *n - 2; i > j; i--) x[i + 1] = x[i];
                    if(fabs(x[*n-2]-tmp) > EPS) x[j + 1] = x[j]/2 + x[j + 1]/2;
                    else { x[j] = x[j-1]/2 + x[j]/2; x[j + 1] = b; }
                }
                else {
                    for(j = 0; j < *n - 1 && x[j+1] < tmp; j++);
                    for (int i = *n - 2; i > j; i--) x[i + 1] = x[i];
                    x[j + 1] = tmp;
                }
                break;
            case 3:
                (*n)++;
                for (int i = *n - 1; i >= 0; i--)
                    x[*n - i - 1] = a/2 + b/2 + (b-a)*cos((2*i + 1)*acos(0)/(*n))/2;
                break;
            }
        } key_1 = 0;
    }
    if(key_1 == -1) { if(*n>2) {
        switch(mode) {
            case 1:
                (*n)-=1; for (int i = 0; i < *n; i++) x[i] = a + i*(b-a)/(*n-1);
                break;
            case 2:
                (*n)-=1; j = 1 + rand()%(*n - 1);
                for (int i = j; i < *n; i++) x[i] = x[i + 1];
                break;
            case 3:
               (*n)-=1;
                for (int i = *n - 1; i >= 0; i--)
                    x[*n - i - 1] = a/2 + b/2 + (b-a)*cos((2*i + 1)*acos(0)/(*n))/2;
                break;
        }
      } key_1 = 0;
    }
    if(key_2 != 0) { scale(key_2, a, b, LAMBDA); key_2 = 0; }
    if(key_3 != 0) { xshift(key_3, SHIFT);  key_3 = 0; }
    if(key_4 != 0) { yshift(key_4, SHIFT); key_4 = 0; }
    DrawAxes2D();
        
    stp = STEP * (ymax - ymin)/2;
    DrawGrid(key_7);
    if(variant != 4 && variant != 5) {
    SetLineWidth(0);
        if(a > xmin && a < xmax) {
            DrawLine (a, stp, a, -stp);
            DrawPoint(invmap(a,0).x, invmap(a,0).y, 3);
            if(key_5 == 1&&x[0] < xmax) { SetLineStyle(DASH);
                DrawPoint(invmap(x[0],0).x, invmap(x[0],0).y,3);
                DrawPoint(invmap(x[0],f(x[0])).x, invmap(x[0],f(x[0])).y,3);
                DrawLine (x[0],0,x[0],f(x[0]));
                SetLineStyle(SOLID);
            }
        }
        if(b > xmin && b < xmax) {
            DrawLine (b, stp, b, -stp);
            DrawPoint(invmap(b, 0).x, invmap(b, 0).y, 3);
            if(key_5 == 1 && x[*n - 1] > xmin) { SetLineStyle(DASH);
                DrawPoint(invmap(x[*n - 1], 0).x, invmap(x[*n - 1], 0).y, 3);
                DrawPoint(invmap(x[*n - 1], f(x[*n - 1])).x, invmap(x[*n - 1], f(x[*n - 1])).y, 3);
                DrawLine (x[*n - 1], 0, x[*n - 1], f(x[*n - 1]));
                SetLineStyle(SOLID);
            }
        } dtoa(a,s,2); dtoa(b,str,2);
        if(a > xmin && a < xmax) WDrawString (s, invmap(a - 1.1*stp, -1.9*stp).x , invmap(a - 1.1*stp, -1.9*stp).y);
        if(b > xmin && b < xmax) WDrawString (str, invmap(b - 1.1*stp, -1.9*stp).x , invmap(b - 1.1*stp, -1.9*stp).y);
        WSetColor (LIGHTGREEN); SetLineWidth(0);
        for (int i = 0; i < *n - 1; i++) {
                if(i > 0 && x[i] > xmin && x[i] < xmax && key_5 == 1) {
                    DrawLine (x[i], stp/2, x[i], -stp/2);
                    SetLineStyle(DASH);
                    DrawPoint(invmap(x[i], 0).x, invmap(x[i], 0).y, 3);
                    DrawPoint(invmap(x[i], f(x[i])).x, invmap(x[i], f(x[i])).y, 3);
                    DrawLine(x[i], 0, x[i], f(x[i]));
                    SetLineStyle(SOLID);
                }
            }
        }
        SetFont(HELVETICA16);
        itoa(*n,s); scat("Number of interpolation nodes = ", s, 32, str); WDrawString(str, 20, 20);
        SetLineWidth(4);
        if(key_6 == 0 && variant != 4 && variant != 5) DrawGraph2D(f, a, b);
        switch (variant) {
            case 1:
                interpH(x, *n, m);
                WSetColor (BLUE);
                break;
            case 2:
                interpS(x, *n, m);
                WSetColor (MAGENTA);
                break;
            case 4:
                parf[0] = parx2D; parf[1] = pary2D;
                if(key_6 == 0) DrawParametric2D(parf, 0, 2*PI);
                break;
            case 5:
                DrawPoint(invmap(1,0).x, invmap(1,0).y, 4);
                switch (key_10) {
                    case 1:
                        WDrawString("Equation: -u'' + 6u = exp(x)", 20, 40);
                        break;
                    case 2:
                        WDrawString("Equation: -u'' + 6u = exp(x/2)", 20, 40);
                        break;
                    case 3:
                        WDrawString("Equation: -u'' + 6u = exp(x^2)", 20, 40);
                        break;
                }
                solveEq(key_10, *n, pt);
                DrawLinear(pt, *n);
                break;
        }  if(variant != 3 && variant != 4 && variant != 5) {
            for (int i = 0; i < *n - 1; i++) {
                for(double t = x[i]; t < x[i + 1]; t += dx) {
                    if(t > xmin && t < xmax && polynom(t, x, m, i) > ymin && polynom(t, x, m, i) < ymax &&
                       t + dx > xmin && t + dx < xmax && polynom(t + dx, x, m, i) > ymin && polynom(t + dx, x, m, i) < ymax)
                        DrawLine(t, polynom(t, x, m, i), t + dx, polynom(t + dx, x, m, i));
                }
            }
        } 
  }
    } else {
        if (count == 0) {
            InitialPosition();
            count++;
        }
        if(variant == 0) {
            SetLineWidth(1);
            WSetColor(MAGENTA); 
            WDrawString ("DGraphX", 300, 169);
            WSetColor(LIGHTGREEN);
            for(int i = 0; i < 9; i++) DrawArc(193+3*i, 170, 50+i*(20+2*i), 50+i*(20+2*i), 0,360*64);
            WDrawString("This is 3D mode", 100, 500);  SetFont(HELVETICA12);
            WDrawString("' Ctrl ' -- change dimension; ' 1 ' -- draw 3D graph z(x,y); ' 2 ' -- draw 3D parametric curve (x(t), y(t), z(t))", 100, 550);
            WDrawString("' 3 ' -- draw 3D parametric surface (x(t,s), y(t,s), z(t,s)); ' 4 ' -- draw polytope (from data.txt)", 100, 575);
            WDrawString("' i ' -- enter new function expression in terminal", 100, 600);
        } else {
            if(key_2 != 0) { ChangeCameraPosition(key_2 * -0.5); key_2 = 0; }
            if(key_3 != 0) {
                switch (key_7) {
                    case 1: RotateX(key_3); break;
                    case 2: RotateY(key_3); break;
                    case 3: RotateZ(key_3); break;
                } key_3 = 0;
            }
            if(key_4 != 0) {
                if(key_12 == 1) ChangeCameraPosition(key_4 * -0.5);
                if(key_12 == 2) ChangeScreenPosition(key_4 * -0.5);
                key_4 = 0;
            }
            VectorLight(5.*cos(K*iter), 0, 5.*sin(K*iter));
            switch (key_10) {
                case 1:
                    SetParameter(1);
                    break;
                case 2:
                    SetParameter(0.1);
                    break;
                case 3:
                    SetParameter(0.01);
                    break;
            }
            switch (variant) {
                case 1:
                    SetShadingColor(GREEN);
                    fun[0] = u; fun[1] = v; fun[2] = w;
                    DrawGraph3DX(-5, 5, fun, BLUE, GREEN, key_5, key_8, 1);
                    break;
                case 2:
                    curvefun[0] = curvex; curvefun[1] = curvey; curvefun[2] = curvez;
                    ParametricCurve3D(curvefun, 0, 10*PI, key_5);
                    break;
                case 3:
                    parfun[0] = parx; parfun[1] = pary; parfun[2] = parz;
                    ParametricGraph3D(parfun, RED, BLUE, -PI, PI, -PI, PI, key_5, key_8);
                    break;
                case 4:
                    if(key_11 == 1) {
                    fi = fopen("data.txt", "r");
                    if(fi == NULL) { printf("[DGraphX]: data.txt not found :(\n"); variant = 0; break; }
                    int bad = 0;
                    if (fscanf(fi, "%d", &nFaces) != 1 || nFaces < 0 || nFaces > 100) bad = 1;
                    for (int i = 0; i < nFaces && !bad; i++) {
                        if (fscanf(fi, "%d", &faces[i].n) != 1 || faces[i].n < 0 || faces[i].n > 10) { bad = 1; break; }
                        for (int j = 0; j < faces[i].n; j++) {
                            if (fscanf(fi, "%lf %lf %lf", &faces[i].vertex[0][j], &faces[i].vertex[1][j], &faces[i].vertex[2][j]) != 3) { bad = 1; break; }
                        }
                    }
                    fclose(fi); key_11 = 0;
                    if (bad) { printf("[DGraphX]: bad data.txt format\n"); nFaces = 0; variant = 0; break; }
                    } DrawPolytope(faces, nFaces);
                    break;
                case 5:
                    SetShadingColor(GREEN);
                    SetCentre(0, 0, 0);
                    for(int i = 0; i < ngrid; i++) {
                        for(int jj = 0; jj < ngrid; jj++) {
                            fval[i*ngrid + jj] = solution(i, jj, ngrid);
                        }
                    }
                    DrawLinear3D(fval, ngrid, BLUE, key_5);
                    WSetColor(LIGHTGREEN);
                    itoa(ngrid,s); scat("[Solution] Size of grid: ", s, 25, str); scat(str,"x", strlen(str), s); itoa(ngrid, str); scat(s, str, strlen(s), s);
                    WDrawString(s, 20, 20);
                    break;
                case 6:
                    SetShadingColor(BLUE);
                    SetCentre(0, 0, 0);
                    if(key_11 == 1) { init_val(); t = clock(); SeidelGauss(fval, ngrid, &niter, max_iter); t = clock() - t; key_11 = 0; }
                    DrawLinear3D(fval, ngrid, BLUE, key_5);  WSetColor(LIGHTGREEN);
                    itoa(ngrid,s); scat("[Seidel Gauss] Size of grid: ", s, 29, str); scat(str,"x", strlen(str), s); itoa(ngrid, str); scat(s, str, strlen(s), s);
                    WDrawString(s, 20, 20);
                    dtoa(accuracy(solution, fval, ngrid), s, 5); scat("error: ", s, 7, str); WDrawString(str, 20, 40);
                    itoa(niter, s); scat("number of iterations: ", s, 22, str); WDrawString(str, 20, 60);
                    dtoa((double) t/CLOCKS_PER_SEC, s, 5); scat("time: ", s, 6, str); scat(str, " sec", strlen(str), s); WDrawString(s, 20, 80);
                    break;
                case 7:
                    SetShadingColor(RED);
                    SetCentre(0, 0, 0);
                    if(key_11 == 1) { init_val(); t = clock(); RichardsonIter(fval, ngrid, 2, 8, &niter, max_iter); t = clock() - t; key_11 = 0; }
                    DrawLinear3D(fval, ngrid, BLUE, key_5); WSetColor(LIGHTGREEN);
                    itoa(ngrid,s); scat("[Richardson iter] Size of grid: ", s, 32, str); scat(str,"x", strlen(str), s); itoa(ngrid, str); scat(s, str, strlen(s), s);
                    WDrawString(s, 20, 20);
                    dtoa(accuracy(solution, fval, ngrid), s, 5); scat("error: ", s, 7, str); WDrawString(str, 20, 40);
                    itoa(niter, s); scat("number of iterations: ", s, 22, str); WDrawString(str, 20, 60);
                    dtoa((double) t/CLOCKS_PER_SEC, s, 5); scat("time: ", s, 6, str); scat(str, " sec", strlen(str), s); WDrawString(s, 20, 80);
                    break;
            }
        }
    }
}

static int KeyPressFunction (int nKeySym) {
    char input[256];
    switch (nKeySym) {
        case XK_Q:
        case XK_q:
            return KEY_PRESS_QUIT;
            
        case XK_F1:
        case XK_1:
            variant = 1;
            InitKey();
            break;
        case XK_F2:
        case XK_2:
            variant = 2;
            InitKey();
            break;
        case XK_F3:
        case XK_3:
            variant = 3;
            InitKey();
            break;
        case XK_F4:
        case XK_4:
            variant = 4;
            key_11 = 1;
            InitKey();
            break;
        case XK_F5:
        case XK_5:
            variant = 5;
            key_11 = 1;
            InitKey();
            break;
        case XK_F6:
        case XK_6:
            variant = 6;
            key_11 = 1;
            InitKey();
            break;
        case XK_F7:
        case XK_7:
            variant = 7;
            key_11 = 1;
            InitKey();
            break;
        case XK_plus :
        case XK_equal :
            if(ngrid < 100) ++ngrid;
            key_11 = 1;
            renderEpoch++;
            break;
        case XK_minus :
            if(ngrid > 10) --ngrid;
            key_11 = 1;
            renderEpoch++;
            break;
        case XK_B:
        case XK_b:
            key_1 = -1;
            break;
        case XK_N:
        case XK_n:
            key_1 = 1;
            break;
        case XK_Right:
            if(key_9 == 0) key_3 = 1;
            else iter++;
            break;
        case XK_Left:
            if(key_9 == 0) key_3 = -1;
            else iter--;
            break;
        case XK_Up:
            key_4 = 1;
            break;
        case XK_Down:
            key_4 = -1;
            break;
        case XK_L:
            if(key_8 == 0) key_8 = 1;
            else key_8 = 0;
            break;
        case XK_l:
            if (key_8 == 1) {
                if(key_9 == 0) key_9 = 1;
                else key_9 = 0;
            }
            break;
        case XK_C:
        case XK_c:
            key_2 = -1;
            break;
        case XK_V:
        case XK_v:
            key_2 = 1;
            break;
        case XK_S:
        case XK_s:
            key_12 = 1;
            break;
        case XK_F:
        case XK_f:
            key_12 = 2;
            break;
        case XK_X:
        case XK_x:
            key_7 = 1;
            break;
        case XK_Y:
        case XK_y:
            key_7 = 2;
            break;
        case XK_Z:
        case XK_z:
            key_7 = 3;
            break;
        case XK_Return:
            if(key_5 == 0) key_5 = 1;
            else key_5 = 0;
            break;
        case XK_BackSpace:
            if(key_6 == 0) key_6 = 1;
            else key_6 = 0;
            break;
        case XK_i:
        case XK_I:
            printf("\n[DGraphX] Runtime Function Input:\n");
            printf("Available functions to change:\n");
            printf("  f: 2D function f(x)\n");
            printf("  u, v, w: 3D surfaces u(x,y), v(x,y), w(x,y)\n");
            printf("  cx, cy, cz: 3D parametric curve x(t), y(t), z(t)\n");
            printf("  px, py, pz: 3D parametric surface x(t,s), y(t,s), z(t,s)\n");
            printf("  iters: Max iterations for solvers\n");
            
            printf("Enter function name (e.g. 'f'): ");
            fflush(stdout);
            
            if (fgets(input, sizeof(input), stdin)) {
                trim_newline(input);
                if (strlen(input) == 0) {
                    if (!fgets(input, sizeof(input), stdin)) break;
                    trim_newline(input);
                }
                
                char func_name[16];
                strncpy(func_name, input, 15);
                func_name[15] = '\0';
                
                if (strcmp(func_name, "iters") == 0) {
                    printf("Enter new value for max iterations: ");
                } else {
                    printf("Enter new expression for %s: ", func_name);
                }
                fflush(stdout);
                
                if (fgets(input, sizeof(input), stdin)) {
                    trim_newline(input);
                    if (strlen(input) > 0) {
                        int found = 1;
                        char old_expr[256];
                        
                        if (strcmp(func_name, "f") == 0) { strcpy(old_expr, expr_f); strcpy(expr_f, input); }
                        else if (strcmp(func_name, "u") == 0) { strcpy(old_expr, expr_u); strcpy(expr_u, input); }
                        else if (strcmp(func_name, "v") == 0) { strcpy(old_expr, expr_v); strcpy(expr_v, input); }
                        else if (strcmp(func_name, "w") == 0) { strcpy(old_expr, expr_w); strcpy(expr_w, input); }
                        else if (strcmp(func_name, "cx") == 0) { strcpy(old_expr, expr_cx); strcpy(expr_cx, input); }
                        else if (strcmp(func_name, "cy") == 0) { strcpy(old_expr, expr_cy); strcpy(expr_cy, input); }
                        else if (strcmp(func_name, "cz") == 0) { strcpy(old_expr, expr_cz); strcpy(expr_cz, input); }
                        else if (strcmp(func_name, "px") == 0) { strcpy(old_expr, expr_px); strcpy(expr_px, input); }
                        else if (strcmp(func_name, "py") == 0) { strcpy(old_expr, expr_py); strcpy(expr_py, input); }
                        else if (strcmp(func_name, "pz") == 0) { strcpy(old_expr, expr_pz); strcpy(expr_pz, input); }
                        else if (strcmp(func_name, "iters") == 0) {
                            int new_val = atoi(input);
                            if (new_val > 0) {
                                max_iter = new_val;
                                printf("[DGraphX] max_iter updated to %d\n", max_iter);
                                key_11 = 1;
                            } else {
                                printf("[DGraphX] Invalid value for iters: %s\n", input);
                            }
                            found = 0; // Skip compilation
                        }
                        else {
                            printf("[DGraphX] Unknown function name: %s\n", func_name);
                            found = 0;
                        }
                        
                        if (found) {
                            if (compile_functions()) {
                                printf("[DGraphX] Function %s updated successfully to: %s\n", func_name, input);
                                renderEpoch++;
                            } else {
                                printf("[DGraphX] Failed to compile expression. Reverting to previous state.\n");
                                // Simple revert
                                if (strcmp(func_name, "f") == 0) strcpy(expr_f, old_expr);
                                else if (strcmp(func_name, "u") == 0) strcpy(expr_u, old_expr);
                                else if (strcmp(func_name, "v") == 0) strcpy(expr_v, old_expr);
                                else if (strcmp(func_name, "w") == 0) strcpy(expr_w, old_expr);
                                else if (strcmp(func_name, "cx") == 0) strcpy(expr_cx, old_expr);
                                else if (strcmp(func_name, "cy") == 0) strcpy(expr_cy, old_expr);
                                else if (strcmp(func_name, "cz") == 0) strcpy(expr_cz, old_expr);
                                else if (strcmp(func_name, "px") == 0) strcpy(expr_px, old_expr);
                                else if (strcmp(func_name, "py") == 0) strcpy(expr_py, old_expr);
                                else if (strcmp(func_name, "pz") == 0) strcpy(expr_pz, old_expr);
                                compile_functions();
                            }
                        }
                    } else {
                        printf("[DGraphX] Empty expression, update cancelled.\n");
                    }
                }
            }
            return KEY_PRESS_EXPOSE;
        case XK_Control_L:
        case XK_Control_R:
            if(dim == 2) { dim = 3;
                printf("[DGraphX] Initializing 3D mode...\n");
                variant = 0; key_1 = 0; key_2 = 0; key_3 = 0; key_4 = 0;
                key_5 = 0; key_6 = 0; key_7 = 3; mode = 1; count = 0; count1 = 0;
                key_8 = 0; key_9 = 0;
                xc = yc = INIT_SCALE;
                InitialPosition();
                count = 1;
                SetSize();
            }
            else { dim = 2;
                printf("[DGraphX] Initializing 2D mode...\n");
                variant = 0; key_1 = 0; key_2 = 0; key_3 = 0; key_4 = 0;
                key_5 = 0; key_6 = 0; key_7 = 3; mode = 1; count = 0; count1 = 0;
                key_8 = 0; key_9 = 0;
                xc = yc = INIT_SCALE;
                SetSize();
            }
            break;
        case XK_8:
            mode = 1;
            break;
        case XK_9:
            mode = 3;
            break;
        case XK_Tab:
            if(key_10 < 3) key_10++;
            else key_10 = 1;
            key_11 = 1;
            renderEpoch++;
            break;
        case XK_0:
            variant = 0;
            break;
        default:
            return KEY_PRESS_NOTHING;
    }

    //
    return KEY_PRESS_EXPOSE;
}

int main (void) {
    int ret_code; int n = N;
    double x[NMAX], a = -5, b = 5;
    
    IdMatrix();
    compile_functions();
    SetParameter(1.0);

    if(dim == 2) {
    printf("\n[DGraphX]:\t+----------------------- DGraphX --------------------------+\n");
    printf("\t\t| Еnter the endpoints of the segment of interpolation:     |\n");
    printf("\t\t| Choose segment split mode (uniform -- 1, manually -- 2): |\n");
    printf("\t\t  mode: ");
    if (scanf("%d",&mode) != 1) { printf("[DGraphX]: invalid input\n"); return 1; }
    switch (mode) {
        case 1:
            printf("\t\t| === === === === === === === === === === === === === ===  |\n");
            printf("\t\t  Segment of interpolation: [%.2f,%.2f]; n = %d\n",a,b,n);
            printf("\t\t| Segment split mode -- uniform                            |\n");
            printf("\t\t+----------------------------------------------------------+\n\n");
            break;
        case 2:
            x[0] = a;
            printf("\t\t| Enter number of interpolation nodes:                     |\n");
            printf("\t\t   ");
            if (scanf("%d",&n) != 1 || n < 2 || n > NMAX) { printf("[DGraphX]: invalid node count\n"); return 1; }
            x[n-1] = b;
            printf("\t\t| Enter interpolation nodes:                               |\n");
            for (int i = 1; i < n - 1; i++) {
                printf("\t\t   x[%d] = ",i);
                if (scanf("%lf",&x[i]) != 1) { printf("[DGraphX]: invalid node\n"); return 1; }
            }
            printf("\t\t+----------------------------------------------------------+\n\n");
            break;
        case 3:
            printf("\t\t| === === === === === === === === === === === === === ===  |\n");
            printf("\t\t  Segment of interpolation: [%.2f,%.2f]; n = %d\n",a,b,n);
            printf("\t\t| Segment split mode -- optimal (Chebyshev nodes)          |\n");
            printf("\t\t+----------------------------------------------------------+\n\n");
            break;
        default:
            printf("[DGraphX]: incorrect mode\n");
            return 1;
    }
    }
    
    ret_code = DrawWindow (DrawWindowContent, KeyPressFunction, MouseEventHandler, x, &n,a,b);
    
    if (te_f) te_free(te_f);
    if (te_u) te_free(te_u);
    if (te_v) te_free(te_v);
    if (te_w) te_free(te_w);
    if (te_cx) te_free(te_cx);
    if (te_cy) te_free(te_cy);
    if (te_cz) te_free(te_cz);
    if (te_px) te_free(te_px);
    if (te_py) te_free(te_py);
    if (te_pz) te_free(te_pz);

    switch (ret_code) {
        case X11_ERR_1:
            printf ("[DGraphX]:\t%s\n", X11_ERR_MSG_1);
            break;
        case X11_ERR_2:
            printf ("[DGraphX]:\t%s\n", X11_ERR_MSG_2);
            break;
        case X11_ERR_3:
            printf ("[DGraphX]:\t%s\n", X11_ERR_MSG_3);
            break;
        case X11_ERR_4:
            printf ("[DGraphX]:\t%s\n", X11_ERR_MSG_4);
            break;
        case X11_ERR_5:
            printf ("[DGraphX]:\t%s\n", X11_ERR_MSG_5);
            break;
        default:
            break;
    }
    
    return 0;
}
