// DGraphX -- free grapher, version 1.0
//
// David Mikheev 7/18/19

#include "DGraphX.h"
#include "Interp2D.h"
#include "matrix.h"

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
static int key_8 = 0;
static int key_9 = 0;
static int iter = 0;
static int mode = 1;
static int dim = 3;
static int count = 0;
static int count1 = 0;

double xmin = -8., xmax = 8., ymin = -6., ymax = 6.;

double mypow(double x, int n) {
    double res = 1.;
    for (int i = 0; i < n; i++) res *= x;
    return res;
}

// ================ functions for drawing (input data) ===================

double f(double x) { // function for draw by DrawGraph2D
    return sin(x);
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
    return 3*cos(3*sqrt(x*x+y*y))*exp(-sqrt(x*x+y*y)/3.);//sin(3*sqrt(x*x+y*y))/sqrt(x*x+y*y);
}

double w(double x, double y) { // function for DrawGraph3DX
    return sin(x*y);
}

double z(double x, double y) { // function for DrawGraph3DX
    return 1-(x*x+y*y)/6.;
}

double curvex(double t) { // x-component for ParametricCurve3D
    return /*t*sin(5*t)*/20*mypow(cos(t),2);
}

double curvey(double t) { // y-component for ParametricCurve3D
    return /*t*cos(5*t)*/(10./2.)*sin(t);
}

double curvez(double t) { // z-component for ParametricCurve3D
    return /*t*/10*sin(2*t);
}
double parx(double t, double s) { // x-component for ParametricGraph3D
    return /*3*cos(t)*cos(s);//t*(1-t*t/3.+s*s)/3.;//2*cosh(t)*cos(s);//s*cos(t);//2*(4+cos(t))*cos(s)// (3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*cos(t)*/-(2./15.)*cos(t)*(3*cos(s) - 30*sin(t) + 90*mypow(cos(t), 4)*sin(t)-60*mypow(cos(t), 6)*sin(t)+5*cos(t)*sin(t)*cos(s));
}

double pary(double t, double s) { // y-component for ParametricGraph3D
    return /*3*cos(t)*sin(s);//-s*(1-s*s/3.+t*t)/3.;//2*cosh(t)*sin(s);//s*sin(t);//2*(4+cos(t))*sin(s)//(3+cos(t/2.)*sin(s) - sin(t/2.)*sin(2*s))*sin(t)*/-(1./15.)*sin(t)*(3*cos(s) - 3*mypow(cos(t),2)*cos(s)-48*mypow(cos(t),4)*cos(s)+48*mypow(cos(t),6)*cos(s)-60*sin(t)+5*cos(t)*cos(s)*sin(t)-5*mypow(cos(t),3)*cos(s)*sin(t)-80*mypow(cos(t),5)*cos(s)*sin(t)+80*mypow(cos(t),7)*cos(s)*sin(t));
}

double parz(double t, double s) { // z-component for ParametricGraph3D
    return /*3*sin(t);//(t*t - s*s)/3.;//2*t//sin(t/2.)*sin(s)+cos(t/2.)*sin(2*s)*/(2./15.)*(3+5*cos(t)*sin(t))*sin(s);
}

double polynom(double y,double *x,double *m,int i) { // for interpolation
    return m(i,0) + m(i,1)*(y-x[i]) + m(i,2)*(y-x[i])*(y-x[i]) + m(i,3)*(y-x[i])*(y-x[i])*(y-x[i]);
}

//==========================================================================

static void DrawWindowContent (double* x, int* n, double a, double b) {
    
    char s[128],str[128];
    int j; double tmp;
    double *m, stp, (*fun[NMAX])(double, double), (*curvefun[3])(double), (*parfun[3])(double, double), (*parf[2])(double);
    //unsigned long color[NMAX];//double X[NMAX],Y[NMAX];
    double dx = EPS;
    
    SetFont(HELVETICA16);
    WSetColor (BLACK);
    WFillRectangle (0, 0, width, height);
    if (dim == 2) {
    if(variant == 0) {
        SetLineWidth(1);
        WSetColor(MAGENTA); WDrawString ("DGraphX", 300, 169);
        WSetColor(LIGHTGREEN);
        for(int i = 0; i < 9; i++) DrawArc( 193 + 3*i, 170, 50 + i*(20 + 2*i), 50 + i*(20 + 2*i), 0, 360*64 );
        WDrawString("This is 2D mode", 100, 500);
        SetFont(HELVETICA12);
        WDrawString("' Ctrl ' -- change dimension; ' 1 ' -- draw interpolation with Hermite polynomial; ' 2 ' -- draw interpolation with spline", 100, 525);
        WDrawString("' 3 ' -- draw graph; ' 4 ' -- draw parametric curve (x(t), y(t))", 100, 550);
        xc = (xmax - xmin)/width;
        yc = (ymax - ymin)/height;
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
        } dtoa(a,s); dtoa(b,str);
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
    if(key_2 != 0)  { scale(key_2, a, b, LAMBDA); key_2 = 0; }
    if(key_3 != 0) { xshift(key_3, SHIFT);  key_3 = 0; }
    if(key_4 != 0) { yshift(key_4, SHIFT); key_4 = 0; }
    DrawAxes2D();
        
    stp = STEP * (ymax - ymin)/2;
    if(key_7 == 1) DrawGrid();
    if(variant != 4) {
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
        }
        if(a > xmin && a < xmax) WDrawString (s, invmap(a - 1.1*stp, -1.9*stp).x , invmap(a - 1.1*stp, -1.9*stp).y);
        if(b > xmin && b < xmax) WDrawString (str, invmap(b - 1.1*stp, -1.9*stp).x , invmap(b - 1.1*stp, -1.9*stp).y);
        SetFont(HELVETICA12);
        itoa(*n,s); scat("Number of interpolation nodes = ", s, 32, str); WDrawString(str, 20, 20);
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
        } SetLineWidth(4);
        if(key_6 == 0 && variant != 4) DrawGraph2D(f, a, b);
        switch (variant) {
            case 1:
                m = interpH(x, *n);
                WSetColor (BLUE);
                break;
            case 2:
                m = interpS(x, *n);
                WSetColor (MAGENTA);
                break;
            case 4:
                parf[0] = parx2D; parf[1] = pary2D;
                if(key_6 == 0) DrawParametric2D(parf, 0, 2*PI);
                break;
        }  if(variant != 3 && variant != 4) {
            for (int i = 0; i < *n - 1; i++) {
                for(double t = x[i]; t < x[i + 1]; t += dx) {
                    if(t > xmin && t < xmax && polynom(t, x, m, i) > ymin && polynom(t, x, m, i) < ymax &&
                       t + dx > xmin && t + dx < xmax && polynom(t + dx, x, m, i) > ymin && polynom(t + dx, x, m, i) < ymax)
                        DrawLine(t, polynom(t, x, m, i), t + dx, polynom(t + dx, x, m, i));
                }
            } free(m);
        }
  }
    } else {
        if(variant == 0 || variant == 4) {
            SetLineWidth(1);
            WSetColor(MAGENTA); WDrawString ("DGraphX", 300, 169);
            WSetColor(LIGHTGREEN);
            for(int i = 0; i < 9; i++) DrawArc( 193+3*i, 170, 50+i*(20+2*i), 50+i*(20+2*i), 0,360*64);
            WDrawString("This is 3D mode", 100, 500);  SetFont(HELVETICA12);
            WDrawString("' Ctrl ' -- change dimension; ' 1 ' -- draw 3D graph z(x,y); ' 2 ' -- draw 3D parametric curve (x(t), y(t), z(t))", 100, 525);
            WDrawString("' 3 ' -- draw 3D parametric graph (x(t, s), y(t, s), z(t, s))", 100, 550);
            xc = (xmax-xmin)/width;
            yc = (ymax-ymin)/height;
        } else {
            if(key_2 != 0) {
                xshift(key_2, SHIFT);  key_2 = 0;;
            }
            if(key_4 != 0) {
                //if(count1 + key_2 < 9 && count1 + key_2 > -9) {
                ChangeCameraPosition(key_4 * -0.5);/*scale(key_2, a, b, LAMBDA);  count1 += key_2; }*/   key_4 = 0;
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
            if(key_5 == 1) DrawAxes();
            SetShadingColor(GREEN);
            fun[0] = u; fun[1] = w; fun[2] = z;
            curvefun[0] = curvex; curvefun[1] = curvey; curvefun[2] = curvez;
            parfun[0] = parx; parfun[1] = pary; parfun[2] = parz;
            if(key_6 == 0 && variant == 1) DrawGraph3DX(a, b, fun, BLUE, GREEN, key_5, key_8, 1);
            if(key_6 == 0 && variant == 2) ParametricCurve3D(curvefun, -PI, PI, key_5);
            if(key_6 == 0 && variant == 3) ParametricGraph3D(parfun, RED, BLUE, /*-2.5*/ 0, PI, 0, 2*PI, key_5, key_8);
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
            break;
        case XK_F2:
        case XK_2:
            variant = 2;
            break;
        case XK_F3:
        case XK_3:
            variant = 3;
            break;
        case XK_F4:
        case XK_4:
            variant = 4;
            break;
        case XK_Shift_L:
            key_1 = -1;
            break;
        case XK_Shift_R:
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
                
                xmin = -8., xmax = 8., ymin = -6., ymax = 6.;
            }
            else { dim = 2;
                variant = 0; key_1 = 0; key_2 = 0; key_3 = 0; key_4 = 0;
                key_5 = 0; key_6 = 0; key_7 = 3; mode = 1; count = 0; count1 = 0;
                key_8 = 0; key_9 = 0;
                
                xmin = -8., xmax = 8., ymin = -6., ymax = 6.;
            }
            break;
        case XK_8:
            mode = 1;
            break;
        case XK_9:
            mode = 3;
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
    double x[NMAX], a, b;
    printf("[DGraphX]: Choose dimension:\n");
    printf("[DGraphX]: "); scanf("%d", &dim);
    if(dim == 2) {
    printf("\n[DGraphX]:\t+----------------------- DGraphX --------------------------+\n");
    printf("\t\t| Еnter the endpoints of the segment of interpolation:     |\n");
    printf("\t\t   "); scanf("%lf%lf",&a,&b); //print//if(n<=1) printf("n = %d < 2 Impossible!
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
        printf("\t\t   "); scanf("%lf%lf",&a,&b);
        printf("\t\t| === === === === === === === === === === === === === ===  |\n");
        printf("\t\t  Square of interpolation: [%.1f,%.1f]x[%.1f,%.1f]; n = %d\n",a,b,a,b,n*n);
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
