// DGraphX -- free grapher, version 0.1
//
// David Miheev 7/18/19

#include "DGraphX.h"
#include "Interp2D.h"
#include "matrix.h"
#include <time.h>

#define m(ix,jx) m[4*(ix)+jx]

#define NMAX 100 
#define STEP 0.05
#define LAMBDA 1.1
#define SHIFT 0.01
#define EPS 1e-3
#define N 6
#define ALPHA 64
#define DALPHA 1.7
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
static int mode = 1;
static int dim = 3;
static int count = 0;
static int count1 = 0;
static double m[4*NMAX];
//static double d[NMAX];
static Point pt[NMAX];
static double fval[100000];
static clock_t t;

double mypow(double x, int n) {
    double res = 1.;
    for (int i = 0; i < n; i++) res *= x;
    return res;
}

// ================ functions for drawing (input data) ===================

double f(double x) { // function for draw by DrawGraph2D
    return 10./x/x-9./x;//sin(x);
}

double df(double x) { // derivative of f(x)
    return cos(x);
}

double parx2D(double t) { // x-component for DrawParametric2D
    return cos(t);
}

double pary2D(double t) { // y-component for DrawParametric2D
    return sin(t);
}

double u(double x, double y) { // function for DrawGraph3DX
    return 3*cos(3*sqrt(x*x+y*y))*exp(-sqrt(x*x+y*y)/3.);
}

double w(double x, double y) { // function for DrawGraph3DX
    return sin(x*y);
}

double v(double x, double y) { // function for DrawGraph3DX
    return 3*exp(-x*x/2.-2*y*y);
}

double curvex(double t) { // x-component for ParametricCurve3D
    return 20*mypow(cos(t),2);
}

double curvey(double t) { // y-component for ParametricCurve3D
    return (10./2.)*sin(t);
}

double curvez(double t) { // z-component for ParametricCurve3D
    return 10*sin(2*t);
}
double parx(double t, double s) { // x-component for ParametricGraph3D
    return  3*(3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*cos(t);
}

double pary(double t, double s) { // y-component for ParametricGraph3D
    return 3*(3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*sin(t);
}

double parz(double t, double s) { // z-component for ParametricGraph3D
    return 3*sin(t/2.)*sin(s)+3*cos(t/2.)*sin(2*s);
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
        if(i != 10 && j!= 10) fval[i*ngrid + j] = 0;
        else fval[i*ngrid + j] = sin(2*PI*j/(ngrid-1))*sin(2*PI*i/(ngrid-1))/2;
    }
 }
}

//==========================================================================

static void DrawWindowContent (double* x, int* n, double a, double b) {
    
    char s[128], str[128];
    int j, nFaces; 
    double tmp; Face faces[100]; 
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
        if(variant == 0) {
            SetLineWidth(1);
            WSetColor(MAGENTA); 
            WDrawString ("DGraphX", 300, 169);
            WSetColor(LIGHTGREEN);
            for(int i = 0; i < 9; i++) DrawArc(193+3*i, 170, 50+i*(20+2*i), 50+i*(20+2*i), 0,360*64);
            WDrawString("This is 3D mode", 100, 500);  SetFont(HELVETICA12);
            WDrawString("' Ctrl ' -- change dimension; ' 1 ' -- draw 3D graph z(x,y); ' 2 ' -- draw 3D parametric curve (x(t), y(t), z(t))", 100, 550);
            WDrawString("' 3 ' -- draw 3D parametric graph (x(t, s), y(t, s), z(t, s))", 100, 575);
            WDrawString("' 4 ' -- draw polyhedra", 100, 600);
            WDrawString("' 5,6,7 ' -- draw interpolation plots", 230, 600);
        } else {
            if(key_2 != 0) {
                scale(key_2, -1, 1, LAMBDA);  key_2 = 0;;
            }
            if(key_4 != 0) {

                if(key_12 == 1) ChangeCameraPosition(key_4 * -0.5);
                if(key_12 == 2) ChangeScreenPosition(key_4 * -0.5);
                key_4 = 0;
            } WSetColor (LIGHTGREEN);
            if(count == 0) {
                InitialPosition(); //IdMatrix();
            }
            if(key_3 == 1) {
                count++;
                switch (key_7) {
                    case 1:
                        RotateX(1);
                        break;
                    case 2:
                        RotateY(1);
                        break;
                    case 3:
                        RotateZ(1);
                        break;
                    default:
                        break;
                } key_3 = 0;
            }
            if(key_3 == -1) {
                count++;
                switch (key_7) {
                    case 1:
                        RotateX(-1);
                        break;
                    case 2:
                        RotateY(-1);
                        break;
                    case 3:
                        RotateZ(-1);
                        break;
                    default:
                        break;
                } key_3 = 0;
            } VectorLight(cos(K*DALPHA*iter), 0, sin(K*DALPHA*iter));
            //if(key_5 == 1) DrawAxes();
            SetShadingColor(GREEN);
            SetCentre(0,0,0);
            fun[0] = u; fun[1] = w; fun[2] = v;
            curvefun[0] = curvex; curvefun[1] = curvey; curvefun[2] = curvez;
            parfun[0] = parx; parfun[1] = pary; parfun[2] = parz;
            fi = fopen("data.txt", "r");
            fscanf(fi,"%d",&nFaces);
            for (int i = 0; i < nFaces; i++) {
                fscanf(fi,"%d",&faces[i].n);
                for (int k = 0; k < faces[i].n; k++) {
                    fscanf(fi,"%lf",&faces[i].vertex[0][k]);
                    fscanf(fi,"%lf",&faces[i].vertex[1][k]);
                    fscanf(fi,"%lf",&faces[i].vertex[2][k]);
                    //faces[i].vertex[0][k] *= 3;
                    //faces[i].vertex[1][k] *= 1./3;
                    //faces[i].vertex[2][k] *= 1./3;
                }
            } /*printf("%f, %f, %f\n", faces[0].vertex[0][2], faces[0].vertex[1][2], faces[0].vertex[2][2]);*/ 
            fclose(fi);
            if(key_6 == 0 && variant == 1) DrawGraph3DX(a, b, fun, BLUE, GREEN, key_5, key_8, 1);
            if(key_6 == 0 && variant == 2) ParametricCurve3D(curvefun, -PI, PI, key_5);
            if(key_6 == 0 && variant == 3) ParametricGraph3D(parfun, RED, BLUE, /*-2.5*/ 0, 2*PI, 0, 2*PI, key_5, key_8);
            if(key_6 == 0 && variant == 4) DrawPolytope(faces, nFaces);
            SetCentre(0.5,0.5,0);
            
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
            
            if(variant == 5) { 
                SetShadingColor(GREEN);
                for(int i = 0; i < ngrid; i++) {
                    for(int j = 0; j < ngrid; j++) {
                        fval[i*ngrid + j] = solution(i, j, ngrid);
                    }
                } 
                
                DrawLinear3D(fval, ngrid, BLUE);  
                WSetColor(LIGHTGREEN);
                itoa(ngrid,s); 
                scat("[Solution] Size of grid: ", s, 25, str); 
                scat(str,"x", strlen(str), s); 
                itoa(ngrid, str); 
                scat(s, str, strlen(s), s);
                WDrawString(s, 20, 20);
            }
            
            if(variant == 6) { 
                SetShadingColor(BLUE);
                if(key_11 == 1) { init_val(); t = clock(); SeidelGauss(fval, ngrid, &niter); t = clock() - t; key_11 = 0; }
                DrawLinear3D(fval, ngrid, BLUE);  WSetColor(LIGHTGREEN);
                //printf("%f\n", accuracy(solution, fval, ngrid));
                itoa(ngrid,s); scat("[Seidel Gauss] Size of grid: ", s, 29, str); scat(str,"x", strlen(str), s); itoa(ngrid, str); scat(s, str, strlen(s), s);
                WDrawString(s, 20, 20);
                dtoa(accuracy(solution, fval, ngrid), s, 5); scat("error: ", s, 7, str); WDrawString(str, 20, 40);
                itoa(niter, s); scat("number of iterations: ", s, 22, str); WDrawString(str, 20, 60);
                dtoa((double) t/CLOCKS_PER_SEC, s, 5); scat("time: ", s, 6, str); scat(str, " sec", strlen(str), s); WDrawString(s, 20, 80);
            }
            
            if(variant == 7) { 
                SetShadingColor(RED);
                if(key_11 == 1) { init_val(); t = clock(); RichardsonIter(fval, ngrid, 2, 8, &niter); t = clock() - t; key_11 = 0; }
                DrawLinear3D(fval, ngrid, BLUE); WSetColor(LIGHTGREEN);
                itoa(ngrid,s); scat("[Richardson iter] Size of grid: ", s, 32, str); scat(str,"x", strlen(str), s); itoa(ngrid, str); scat(s, str, strlen(s), s);
                WDrawString(s, 20, 20);
                dtoa(accuracy(solution, fval, ngrid), s, 5); scat("error: ", s, 7, str); WDrawString(str, 20, 40);
                itoa(niter, s); scat("number of iterations: ", s, 22, str); WDrawString(str, 20, 60);
                dtoa((double) t/CLOCKS_PER_SEC, s, 5); scat("time: ", s, 6, str); scat(str, " sec", strlen(str), s); WDrawString(s, 20, 80);
            }
            
        }
    }
}

static int KeyPressFunction (int nKeySym) {
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
            InitKey();
            break;
        case XK_F5:
        case XK_5:
            variant = 5;
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
            break;
        case XK_minus :
            if(ngrid > 10) --ngrid;
            key_11 = 1;
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
        case XK_Control_L:
        case XK_Control_R:
            if(dim == 2) { dim = 3;
                variant = 0; key_1 = 0; key_2 = 0; key_3 = 0; key_4 = 0;
                key_5 = 0; key_6 = 0; key_7 = 3; mode = 1; count = 0; count1 = 0;
                key_8 = 0; key_9 = 0;
                xc = yc = INIT_SCALE;
                SetSize();//xmin = -8., xmax = 8., ymin = -6., ymax = 6.;
            }
            else { dim = 2;
                variant = 0; key_1 = 0; key_2 = 0; key_3 = 0; key_4 = 0;
                key_5 = 0; key_6 = 0; key_7 = 3; mode = 1; count = 0; count1 = 0;
                key_8 = 0; key_9 = 0;
                xc = yc = INIT_SCALE;
                SetSize();//xmin = -8., xmax = 8., ymin = -6., ymax = 6.;
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
    // printf("[DGraphX]: Choose dimension:\n");
    // printf("[DGraphX]: "); scanf("%d", &dim);
    if(dim == 2) {
    printf("\n[DGraphX]:\t+----------------------- DGraphX --------------------------+\n");
    printf("\t\t| Еnter the endpoints of the segment of interpolation:     |\n");
    printf("\t\t| Choose segment split mode (uniform -- 1, manually -- 2): |\n");
    printf("\t\t  mode: "); scanf("%d",&mode);
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
            printf("\t\t   "); scanf("%d",&n); x[n-1] = b;
            printf("\t\t| Enter interpolation nodes:                               |\n");
            printf("\t\t    "); for(int i = 1; i < n - 1; i++) scanf("%lf",&x[i]);
            printf("\t\t| === === === === === === === === === === === === === ===  |\n");
            printf("\t\t  Segment of interpolation: [%.2f,%.2f]; n = %d\n",a,b,n);
            printf("\t\t| Segment split mode -- manually                           |\n");
            printf("\t\t+----------------------------------------------------------+\n\n");
            break;
        case 3:
            printf("\t\t| === === === === === === === === === === === === === ===  |\n");
            printf("\t\t  Segment of interpolation: [%.2f,%.2f]; n = %d\n",a,b,n);
            printf("\t\t| Segment split mode -- Chebyshev nodes                    |\n");
            printf("\t\t+----------------------------------------------------------+\n\n");
            break;
        default:
            printf("[DGraphX]: uncorrect mode\n");
            return 0;
            break;
    }
    } else {
        printf("\n[DGraphX]:\t+----------------------- DGraphX --------------------------+\n");
        printf("\t\t| Еnter the endpoints of the square of interpolation:      |\n");
        printf("\t\t| === === === === === === === === === === === === === ===  |\n");
        printf("\t\t  Square of interpolation: [%.1f,%.1f]x[%.1f,%.1f]; n = %d\n",a, b, a, b, n*n);
        printf("\t\t| This is 3D mode!                                         |\n");
        printf("\t\t+----------------------------------------------------------+\n\n");
    }
    
    ret_code = DrawWindow (DrawWindowContent, KeyPressFunction, x, &n,a,b);
    
    if (ret_code)
    {
        switch (ret_code) {
            case X11_ERR_1:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_1);
                break;
            case X11_ERR_2:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_2);
                break;
            case X11_ERR_3:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_3);
                break;
            case X11_ERR_4:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_4);
                break;
            case X11_ERR_5:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_5);
                break;
            default:
                printf ("[DGraphX]: %s\n", X11_ERR_MSG_DEF);
                break;
        }
        return ret_code;
    }
    return 0;
}
