#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include "DGraphX.h"
#include "Interp2D.h"

#define WND_X           0
#define WND_Y           0 
#define WND_WDT         800 
#define WND_HGT         600
#define WND_MIN_WDT     50
#define WND_MIN_HGT     50
#define WND_BORDER_WDT  5

#define FOREVER for(;;)

#define WND_TITLE       "DGraphX"
#define WND_ICON_TITLE  "DGraphX"
#define PRG_CLASS       "Grapher"

static int SetWindowManagerHints ( Display * prDisplay,
		       const char *psPrgClass,	
		       Window nWnd,	
		       int nMinWidth,	
		       int nMinHeight,	
		       const char *psTitle,	
		       const char *psIconTitle,	
		       Pixmap nIconPixmap ) {
  XSizeHints rSizeHints;	
  XWMHints rWMHints;
  XClassHint rClassHint;
  XTextProperty prWindowName, prIconName;

  if (!XStringListToTextProperty ((char**)&psTitle, 1, &prWindowName)
      || !XStringListToTextProperty ((char**)&psIconTitle, 1, &prIconName))
    return 1;			// No memory! 

  rSizeHints.flags = PPosition | PSize | PMinSize;
  rSizeHints.min_width = nMinWidth;
  rSizeHints.min_height = nMinHeight;

  rWMHints.flags = StateHint | IconPixmapHint | InputHint;
  rWMHints.initial_state = NormalState;
  rWMHints.input = True;
  rWMHints.icon_pixmap = nIconPixmap;

  rClassHint.res_name = 0;
  rClassHint.res_class = (char*) psPrgClass;

  XSetWMProperties (prDisplay, nWnd, &prWindowName, &prIconName, 0, 0, &rSizeHints, &rWMHints, &rClassHint);

  return 0;
}

static int AllocColors (Display * prDisplay, Colormap colormap, COLORS * colors) {
  XColor rColor;
  XColor rColorBase;

  if (!XAllocNamedColor (prDisplay, colormap, "black", &rColor, &rColorBase)) return 1;
  colors->black = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "blue", &rColor, &rColorBase)) return 2;
  colors->blue = rColor.pixel;
   
  if (!XAllocNamedColor (prDisplay, colormap, "green", &rColor, &rColorBase)) return 3;
    
  colors->green = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "cyan", &rColor, &rColorBase)) return 4;
  colors->cyan = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "red", &rColor, &rColorBase))  return 5;
  colors->red = rColor.pixel;
   
  if (!XAllocNamedColor(prDisplay, colormap, "magenta", &rColor, &rColorBase))  return 6;
  colors->magenta = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "brown", &rColor, &rColorBase))  return 7;
  colors->brown = rColor.pixel;

  if (!XAllocNamedColor(prDisplay, colormap, "light gray", &rColor, &rColorBase)) return 8;
  colors->light_gray = rColor.pixel;

  if (!XAllocNamedColor(prDisplay, colormap, "dark gray", &rColor, &rColorBase)) return 9;
  colors->dark_gray = rColor.pixel;

  if (!XAllocNamedColor(prDisplay, colormap, "light blue", &rColor, &rColorBase))  return 10;
  colors->light_blue = rColor.pixel;

  if (!XAllocNamedColor(prDisplay, colormap, "light green", &rColor, &rColorBase))  return 11;
  colors->light_green = rColor.pixel;

  if (!XAllocNamedColor(prDisplay, colormap, "light cyan", &rColor, &rColorBase))  return 12;
  colors->light_cyan = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "yellow", &rColor, &rColorBase))   return 13;
  colors->yellow = rColor.pixel;

  if (!XAllocNamedColor (prDisplay, colormap, "white", &rColor, &rColorBase)) return 14;
  colors->white = rColor.pixel;

  return 0;//printf("%u %u %u\n",rColor.blue, rColor.green, rColor.red);
}

//============== XDraw =====================
static Display *prDisplay;
static int nScreenNum;		
static GC prGC; // graphic context
static Window nWnd;		
static Drawable draw;
static Font font;
static XGCValues val;
static int key = 0;
unsigned long tmp;
COLORS rColors;

int width = WND_WDT;
int height = WND_HGT;

int DrawWindow (void (*DrawWindowContent) (double*, int*, double, double),
	    int (*KeyPressFunction) (int), double* x, int* n, double a, double b) {
  XEvent rEvent;
  XWindowAttributes rAttributes;
  KeySym nKeySym;
  int nDepth;
  
  if ((prDisplay = XOpenDisplay (NULL)) == NULL)
    return PLOT_X11_ERR_1;
  
  nScreenNum = DefaultScreen (prDisplay);

  nWnd = XCreateSimpleWindow (prDisplay, RootWindow (prDisplay, nScreenNum),
			      WND_X, WND_Y, WND_WDT, WND_HGT, WND_BORDER_WDT,
			      BlackPixel (prDisplay, nScreenNum),
			      WhitePixel (prDisplay, nScreenNum));

  if (SetWindowManagerHints (prDisplay, PRG_CLASS, nWnd,
			     WND_MIN_WDT, WND_MIN_HGT, WND_TITLE, WND_ICON_TITLE, 0))
    return PLOT_X11_ERR_2;
  
  XSelectInput (prDisplay, nWnd, ExposureMask | KeyPressMask);

  XMapWindow (prDisplay, nWnd);
    
  prGC = XCreateGC (prDisplay, nWnd, GCLineWidth, &val);
    
  if (AllocColors(prDisplay, DefaultColormap (prDisplay, nScreenNum), &rColors))
    return PLOT_X11_ERR_3;

  if (XGetWindowAttributes (prDisplay, nWnd, &rAttributes)) {
      width = rAttributes.width;
      height = rAttributes.height;
      nDepth = rAttributes.depth;
    } else return PLOT_X11_ERR_5;

  draw = XCreatePixmap (prDisplay, nWnd, width, height, nDepth);
    
  if (!draw)
    return PLOT_X11_ERR_4;
    
  //=========== XLoop ==============
  FOREVER {
      
      XNextEvent (prDisplay, &rEvent);

      switch (rEvent.type) {
          case Expose:
              if (rEvent.xexpose.count != 0)
                  break;

              if (XGetWindowAttributes (prDisplay, nWnd, &rAttributes)) {
                  if (width != rAttributes.width
                      || height != rAttributes.height
                      || nDepth != rAttributes.depth) {
                      XFreePixmap (prDisplay, draw);
                      width = rAttributes.width;
                      height = rAttributes.height;
                      nDepth = rAttributes.depth;
                      draw = XCreatePixmap (prDisplay, nWnd, width, height, nDepth);
                      if (!draw) return PLOT_X11_ERR_4;
                  }
              } else return PLOT_X11_ERR_5;

              DrawWindowContent (x,n,a,b);
              XCopyArea (prDisplay, draw, nWnd, prGC, 0, 0, width, height, 0, 0);
              break;
              
          case KeyPress:
              XLookupString (&rEvent.xkey, NULL, 0, &nKeySym, NULL);
              switch (KeyPressFunction (nKeySym)) {
                  case KEY_PRESS_NOTHING:
                      break;

                  case KEY_PRESS_EXPOSE:
                      rEvent.type = Expose;
                      rEvent.xexpose.count = 0;
                      XSendEvent (prDisplay, nWnd, True, 0, &rEvent);
                      break;

                  default:
                  case KEY_PRESS_QUIT:
                      XFreePixmap (prDisplay, draw);
                      if(key != 0) XUnloadFont(prDisplay,font);
                      XFreeGC (prDisplay, prGC);
                      XCloseDisplay (prDisplay);
                      return 0;
              }
              break;
      } //while (XCheckMaskEvent(prDisplay, KeyPressMask | KeyReleaseMask, &rEvent));
    }
   //=====================================
    
  return 0;
}
//============================

//========== XInstruments =================
void WSetColor (unsigned long color) {
  XSetForeground (prDisplay, prGC, color);
}

void WDrawString (const char *string, int x, int y) {
  XDrawString (prDisplay, draw, prGC, x, y, string, strlen (string));
}

void WDrawPoint (int x, int y) {
  XDrawPoint (prDisplay, draw, prGC, x, y);
}

void WDrawLine (int x_start, int y_start, int x_end, int y_end) {
  XDrawLine (prDisplay, draw, prGC, x_start, y_start, x_end, y_end);
}

void WDrawRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right) {
  XDrawRectangle (prDisplay, draw, prGC, x_top_left, y_top_left,
		  x_bottom_right - x_top_left + 1,
		  y_bottom_right - y_top_left + 1);
}

void WFillRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right) {
  XFillRectangle (prDisplay, draw, prGC, x_top_left, y_top_left,
		  x_bottom_right - x_top_left + 2,
		  y_bottom_right - y_top_left + 2);
}

void WFillTriangle (int x_1, int y_1, int x_2, int y_2, int x_3, int y_3) {
  XPoint points[3] = { {x_1, y_1}, {x_2, y_2}, {x_3, y_3} };
  XFillPolygon (prDisplay, draw, prGC, points, 3, Convex, CoordModeOrigin);
}

void WFillPolygon (XPoint * points, int num) {
  XFillPolygon (prDisplay, draw, prGC, points, num, Convex, CoordModeOrigin);
}
   /* FillSolid, FillTiled, FillStippled, FillOpaeueStippled. */
void WSetFillStyle (int Style) {
  XSetFillStyle (prDisplay, prGC, Style);
}

int WGetColor (unsigned int Red, unsigned int Green, unsigned int Blue, unsigned long *pColor) {
  XColor rColor;

  rColor.red = Red << 8;
  rColor.green = Green << 8;
  rColor.blue = Blue << 8;
  rColor.flags = DoRed | DoGreen | DoBlue;

  if (XAllocColor (prDisplay, DefaultColormap (prDisplay, nScreenNum), &rColor) == 0)	return 1;
  *pColor = rColor.pixel;
  return 0;
}

void WSetTitle (const char *s) {
  XTextProperty tProp;
  XStringListToTextProperty ((char**)&s, 1, &tProp);
  XSetWMName (prDisplay, nWnd, &tProp);
}

void DrawArc(int x,int y, unsigned int width, unsigned int height, int a1, int a2) {
    XDrawArc(prDisplay,draw,prGC,x,y,width,height,a1,a2);
}

void DrawPoint(int x, int y, int d) {
    XPoint points[4] = { {x-d, y-d}, {x-d, y+d}, {x+d, y+d}, {x+d, y-d} };
    XFillPolygon (prDisplay, draw, prGC, points, 4, Convex, CoordModeOrigin); //XFillArc(prDisplay, draw, prGC,x,y,d,d,0,360*64);
}

void SetFont(const char* name) {
    if(key != 0) XUnloadFont(prDisplay, font);
    font = XLoadFont(prDisplay,name);
    XSetFont(prDisplay,prGC,font);
    key = 1;
}

void SetLineWidth(int line_width) {
    val . line_width = line_width;
    XChangeGC(prDisplay, prGC, GCLineWidth, &val);
}

void SetLineStyle(int style) {
    val . line_style = style;
    XChangeGC(prDisplay, prGC, GCLineStyle, &val);
}
//=========



