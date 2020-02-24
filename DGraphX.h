// DGraphX -- free grapher, version 1.0 (also can draw interpolation of function of one variable)
//
// David Mikheev 7/18/19


#ifndef DGRAPHX_H
#define DGRAPHX_H

#include <X11/keysym.h>
#include <X11/Xlib.h>
 
typedef struct {
  unsigned long black;
  unsigned long blue;
  unsigned long green; 
  unsigned long cyan;
  unsigned long red;
  unsigned long magenta;
  unsigned long brown;
  unsigned long light_gray;
  unsigned long dark_gray;
  unsigned long light_blue;
  unsigned long light_green;
  unsigned long light_cyan;
  unsigned long yellow;
  unsigned long white;
} COLORS;

//======= Colors ===================================
#define BLACK           (rColors.black)
#define BLUE            (rColors.blue)
#define GREEN           (rColors.green)
#define CYAN            (rColors.cyan)
#define RED             (rColors.red)
#define MAGENTA         (rColors.magenta)
#define BROWN           (rColors.brown)
#define DARKGRAY        (rColors.dark_gray)
#define LIGHTGRAY       (rColors.light_gray)
#define LIGHTBLUE       (rColors.light_blue)
#define LIGHTGREEN      (rColors.light_green)
#define LIGHTCYAN       (rColors.light_cyan)
#define YELLOW          (rColors.yellow)
#define WHITE           (rColors.white)
//================================================

//======= Fonts ====================================================
#define HELVETICA16 "-adobe-helvetica-medium-o-normal-*-16-*-*-*-*-*-*-*"
#define HELVETICA12 "-adobe-helvetica-medium-o-normal-*-12-*-*-*-*-*-*-*"
//========================================================================

//======= Styles ==============
#define SOLID LineSolid
#define DASH LineOnOffDash
//=============================

#define INIT_SCALE 0.02

//======= Global data =========
extern COLORS rColors;
extern int width;
extern int height;
extern double xc, yc;
extern double xmin, xmax, ymin, ymax;
extern int x_0, y_0;
//=============================

#define KEY_PRESS_NOTHING	0
#define KEY_PRESS_EXPOSE	1	
#define KEY_PRESS_QUIT		( -1 )


typedef struct {
    double x;
    double y;
} Point;

typedef  struct {
    double z;
    int num;
} Pair;

typedef struct {
    double vertex[3][10];
    int n;
} Face;

typedef struct {
    double x;
    double y;
    double z;
    double t;
    int data;
    int data_1;
} vertex;

int DrawWindow (void (*DrawWindowContent) (double*,int*,double,double),
		int (*KeyPressFunction) (int), double*, int*,double,double);


void WDrawString (const char *string, int x, int y);
void WDrawPoint (int x, int y);
void DrawPoint(int x, int y,int d);
void WDrawLine (int x_start, int y_start, int x_end, int y_end);
void WDrawRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right);
void WFillRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right);
void WFillTriangle (int x_1, int y_1, int x_2, int y_2, int x_3, int y_3);
void WFillPolygon (XPoint * points, int num);
void DrawArc(int ,int , unsigned int, unsigned int, int, int);
void DrawLine(double,double,double,double);
void DrawAxes2D();
void DrawAxes(int *n);
void DrawGrid(int mode);
void scale(int mode, double a, double b, double lambda);
void xshift(int mode, double valueOfShift);
void yshift(int mode, double valueOfShift);

void DrawLinear(Point *pt, int n);

void DrawGraph2D(double f(double), double a, double b);

void DrawParametric2D(double (*parf[]) (double), double it, double ft);

void DrawGraph3DX(double a, double b, double (*z[]) (double, double), unsigned long, unsigned long,
                  int, int, int);

void ParametricCurve3D(double (*curvefun[]) (double), double it, double s, int mode);

void ParametricGraph3D(double (*parfun[]) (double, double),
                       unsigned long, unsigned long, double, double, double, double, int, int);

void DrawPolytope(Face *f, int n);

void SetSize();
void SetLineWidth(int);

void SetLineStyle(int style);

void WSetFillStyle (int style);

void SetFont(const char* name);

void WSetColor (unsigned long color);

void WSetTitle (const char *s);

XPoint invmap(double x, double y);

Point map(int x, int y);
Point ObliqueProjection ( double x, double y, double z );
void ChangeCameraPosition(double valueOfShift);
void SetShadingColor(unsigned long color);

void VectorLight(double x, double y, double z);
void InitialPosition ();
void RotateX(int mode);
void RotateY(int mode);
void RotateZ(int mode);
void IdMatrix();

void pallette(unsigned long firstColor, unsigned long secondColor, int modeColor, int size);
int WGetColor (unsigned int, unsigned int, unsigned int, unsigned long*);
void QueryColor(XColor *color);

#define X11_ERR_1	1
#define X11_ERR_2	2
#define X11_ERR_3	3
#define X11_ERR_4	4
#define X11_ERR_5	5


#define X11_ERR_MSG_1	"Cannot connect to the X server :("
#define X11_ERR_MSG_2	"Not enough memory to create X Windows data structures :("
#define X11_ERR_MSG_3	"Cannot allocate enough colors in the X Windows pallete :("
#define X11_ERR_MSG_4	"Cannot allocate enough X bitmaps :("
#define X11_ERR_MSG_5  "Bad Window :("

#define X11_ERR_MSG_DEF	"Unknown error in an X Windows interface code :("

#endif
