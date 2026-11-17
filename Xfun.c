#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/select.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#include "DGraphX.h"
#include "Interp2D.h"

#define WND_X           0
#define WND_Y           0 
#define WND_WDT         800 
#define WND_HGT         700
#define WND_MIN_WDT     50
#define WND_MIN_HGT     50
#define WND_BORDER_WDT  5
#define PIXMAP_SIZE     3000


#define FOREVER for(;;)

#define WND_TITLE       "DGraphX"
#define WND_ICON_TITLE  "DGraphX"
#define PRG_CLASS       "Grapher"

static int SetWindowManagerHints ( Display * display,
		       const char *psPrgClass,	
		       Window wnd,
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

  XSetWMProperties (display, wnd, &prWindowName, &prIconName, 0, 0, &rSizeHints, &rWMHints, &rClassHint);

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
static void FreeFonts(void);
static XGCValues val;
int renderEpoch = 0;
int interactive = 0;
int refineEnabled = 1;
int phongEnabled = 1;
int gridEnabled = 1;
#define FB_MAXATTR 8
static XImage *fbImg = NULL;
static uint32_t *fbPix = NULL;
static float *fbZ = NULL;
static int fbW = 0, fbH = 0;
static Visual *fbVisual = NULL;
static int fbVisualDepth = 0;

COLORS rColors;

int width = WND_WDT;
int height = WND_HGT;
int x_0 = PIXMAP_SIZE/2 - WND_WDT/2;
int y_0 = PIXMAP_SIZE/2 - WND_HGT/2;

int DrawWindow (void (*DrawWindowContent) (double*, int*, double, double),
	    int (*KeyPressFunction) (int), MouseEventFunc MouseFunc, double* x, int* n, double a, double b) {
  XEvent rEvent;
  XWindowAttributes rAttributes;
  KeySym nKeySym;
  int nDepth;
  
  if ((prDisplay = XOpenDisplay (NULL)) == NULL)
    return X11_ERR_1;
  
  nScreenNum = DefaultScreen (prDisplay);

  nWnd = XCreateSimpleWindow (prDisplay, RootWindow (prDisplay, nScreenNum),
			      WND_X, WND_Y, WND_WDT, WND_HGT, WND_BORDER_WDT,
			      BlackPixel (prDisplay, nScreenNum),
			      BlackPixel (prDisplay, nScreenNum));

  if (SetWindowManagerHints (prDisplay, PRG_CLASS, nWnd,
			     WND_MIN_WDT, WND_MIN_HGT, WND_TITLE, WND_ICON_TITLE, 0))
    return X11_ERR_2;
  
  XSelectInput (prDisplay, nWnd, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask);

  XMapWindow (prDisplay, nWnd);
    
    prGC = XCreateGC (prDisplay, nWnd, GCLineWidth, &val);
    // XCopyArea from the off-screen pixmap never needs exposure events; without
    // this every copy queues a NoExpose that wakes the event loop for nothing.
    XSetGraphicsExposures (prDisplay, prGC, False);

  if (AllocColors(prDisplay, DefaultColormap (prDisplay, nScreenNum), &rColors))
    return X11_ERR_3;

  if (XGetWindowAttributes (prDisplay, nWnd, &rAttributes)) {
      width = rAttributes.width;
      height = rAttributes.height;
      nDepth = rAttributes.depth;
      fbVisual = rAttributes.visual;
      fbVisualDepth = rAttributes.depth;
    } else return X11_ERR_5;

  draw = XCreatePixmap (prDisplay, nWnd, PIXMAP_SIZE, PIXMAP_SIZE, nDepth);
    SetSize();
    
  if (!draw)
    return X11_ERR_4;
    
  //=========== XLoop ==============
  int coarseShown = 0, pendingRedraw = 0;
  FOREVER {
      int redraw = 0, quit = 0, ev = 0, hadInput = 0;

      XNextEvent (prDisplay, &rEvent);

      // Process the current event plus everything already queued,
      // then redraw once per batch (coalesces mouse/scroll/key floods).
      for (;;) {
          switch (rEvent.type) {
              case Expose:
                  if (rEvent.xexpose.count != 0)
                      break;

                  if (XGetWindowAttributes (prDisplay, nWnd, &rAttributes)) {
                      if (width != rAttributes.width
                          || height != rAttributes.height
                          || nDepth != rAttributes.depth) {
                          width = rAttributes.width;
                          height = rAttributes.height;
                          nDepth = rAttributes.depth;
                          x_0 = PIXMAP_SIZE/2 - width/2;
                          y_0 = PIXMAP_SIZE/2 - height/2;
                          SetSize();
                          if (!draw) return X11_ERR_4;
                      }
                  } else return X11_ERR_5;

                  redraw = 1;
                  break;

              case KeyPress:
                  hadInput = 1;
                  XLookupString (&rEvent.xkey, NULL, 0, &nKeySym, NULL);
                  switch (KeyPressFunction (nKeySym)) {
                      case KEY_PRESS_NOTHING:
                          break;
                      case KEY_PRESS_EXPOSE:
                          redraw = 1;
                          break;
                      default:
                          quit = 1;
                          break;
                  }
                  break;

              case ButtonPress:
              case ButtonRelease:
              case MotionNotify:
                  hadInput = 1;
                  if (rEvent.type == MotionNotify)
                      while (XCheckTypedWindowEvent (prDisplay, nWnd, MotionNotify, &rEvent));
                  if (MouseFunc) {
                      int res = 0;
                      if (rEvent.type == MotionNotify) {
                          res = MouseFunc(rEvent.type, rEvent.xmotion.x, rEvent.xmotion.y, 0);
                      } else {
                          res = MouseFunc(rEvent.type, rEvent.xbutton.x, rEvent.xbutton.y, rEvent.xbutton.button);
                      }
                      if (res == KEY_PRESS_EXPOSE) redraw = 1;
                      else if (res == KEY_PRESS_QUIT) quit = 1;
                  }
                  break;
          }
          if (quit || !XPending (prDisplay) || ++ev > 64) break;
          XNextEvent (prDisplay, &rEvent);
      }

      if (quit) {
          XFreePixmap (prDisplay, draw);
          FreeFonts();
          XFreeGC (prDisplay, prGC);
          XCloseDisplay (prDisplay);
          return 0;
      }

      // A redraw skipped by the frame-rate cap is carried to the next batch,
      // so it is not lost when the queued events don't request one themselves.
      redraw |= pendingRedraw;
      if (redraw) {
          // Frame-rate cap: during input floods, skip redraws more frequent
          // than ~30fps; the skipped frame is drawn after the next batch.
          struct timespec now;
          clock_gettime (CLOCK_MONOTONIC, &now);
          static struct timespec lastDraw;
          long elapsed = (now.tv_sec - lastDraw.tv_sec)*1000 +
                         (now.tv_nsec - lastDraw.tv_nsec)/1000000;
          if (elapsed < 33 && XPending (prDisplay) > 0) { redraw = 0; pendingRedraw = 1; }
          else { clock_gettime (CLOCK_MONOTONIC, &lastDraw); pendingRedraw = 0; }
      }
      if (redraw) {
          // Mid-interaction (progressive refinement): render a coarse frame.
          interactive = hadInput && refineEnabled;
          DrawWindowContent (x,n,a,b);
          XCopyArea (prDisplay, draw, nWnd, prGC, x_0, y_0, width, height, 0, 0);
          coarseShown = interactive;
          interactive = 0;
      }
      // A coarse frame is on screen: once input has been quiet for 40ms,
      // replace it with a full-resolution frame. This runs after every batch,
      // including one whose events requested no redraw (e.g. the ButtonRelease
      // that ends a drag), so a coarse frame is never left on screen.
      // XPending also flushes the coarse frame's XCopyArea to the server.
      if (coarseShown && !XPending (prDisplay)) {
          fd_set fds; struct timeval tv;
          FD_ZERO(&fds); FD_SET(ConnectionNumber (prDisplay), &fds);
          tv.tv_sec = 0; tv.tv_usec = 40000;
          if (select (ConnectionNumber (prDisplay) + 1, &fds, NULL, NULL, &tv) == 0) {
              DrawWindowContent (x,n,a,b);
              XCopyArea (prDisplay, draw, nWnd, prGC, x_0, y_0, width, height, 0, 0);
              coarseShown = 0;
          }
      }
    }
   //======================================
    
  return 0;
}
//=========================================

//========== XInstruments =================
void WSetColor (unsigned long color) {
  XSetForeground (prDisplay, prGC, color);
}

void WDrawString (const char *string, int x, int y) {
  XDrawString (prDisplay, draw, prGC, x_0 + x, y_0 + y, string, strlen (string));
}

void WDrawPoint (int x, int y) {
  XDrawPoint (prDisplay, draw, prGC, x_0 + x, y_0 + y);
}

void WDrawLine (int x_start, int y_start, int x_end, int y_end) {
  XDrawLine (prDisplay, draw, prGC, x_0 + x_start, y_0 + y_start, x_0 + x_end, y_0 + y_end);
}

void WDrawRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right) {
  XDrawRectangle (prDisplay, draw, prGC, x_0 + x_top_left, y_0 + y_top_left,
		  x_0 + x_bottom_right - x_top_left + 1,
		  y_0 + y_bottom_right - y_top_left + 1);
}

void WFillRectangle (int x_top_left, int y_top_left, int x_bottom_right, int y_bottom_right) {
  XFillRectangle (prDisplay, draw, prGC, x_0 + x_top_left, y_0 + y_top_left,
		  x_0 + x_bottom_right - x_top_left + 2,
		  y_0 + y_bottom_right - y_top_left + 2);
}

void WFillTriangle (int x_1, int y_1, int x_2, int y_2, int x_3, int y_3) {
  XPoint points[3] = { {x_0 + x_1, y_0 + y_1}, {x_0 + x_2, y_0 + y_2}, {x_0 + x_3, y_0 + y_3} };
  XFillPolygon (prDisplay, draw, prGC, points, 3, Convex, CoordModeOrigin);
}

void WFillPolygon (XPoint * points, int num) {
  if (num <= 0) return;
  
  // Check if all points are within the pixmap bounds after translation
  for (int i = 0; i < num; i++) {
    int tx = points[i].x + x_0;
    int ty = points[i].y + y_0;
    if (tx < 0 || tx >= PIXMAP_SIZE || ty < 0 || ty >= PIXMAP_SIZE) return;
  }

  // Temporarily translate points for drawing
  for (int i = 0; i < num; i++) {
    points[i].x += x_0;
    points[i].y += y_0;
  }

  XFillPolygon (prDisplay, draw, prGC, points, num, Convex, CoordModeOrigin);

  // Restore points to their original state
  for (int i = 0; i < num; i++) {
    points[i].x -= x_0;
    points[i].y -= y_0;
  }
}
   /* FillSolid, FillTiled, FillStippled, FillOpaeueStippled. */
void WSetFillStyle (int Style) {
  XSetFillStyle (prDisplay, prGC, Style);
}

//======= software z-buffered framebuffer ========
static int fbEnsure (void) {
  if (width <= 0 || height <= 0 || !fbVisual) return 0;
  if (fbW == width && fbH == height && fbImg && fbPix && fbZ) return 1;
  if (fbImg) { XDestroyImage (fbImg); fbImg = NULL; }
  free (fbPix); free (fbZ); fbPix = NULL; fbZ = NULL;
  fbW = width; fbH = height;
  fbPix = malloc ((size_t)fbW*fbH*sizeof (*fbPix));
  fbZ = malloc ((size_t)fbW*fbH*sizeof (*fbZ));
  fbImg = XCreateImage (prDisplay, fbVisual, fbVisualDepth, ZPixmap, 0,
			malloc ((size_t)fbW*fbH*4), fbW, fbH, 32, fbW*4);
  if (!fbPix || !fbZ || !fbImg || !fbImg->data) return 0;
  return 1;
}

// Draw commands are recorded during FB_Begin..FB_Flush and replayed
// across per-row bands in worker threads at flush time.
typedef struct {
  int type;           // 0 flat tri, 1 attr tri, 2 seg
  double x[3], y[3], q[3];
  float at[3][FB_MAXATTR];
  int na;
  unsigned long col;
  unsigned long (*shade) (const float *, void *);
  unsigned char ctx[24];   // ShadeCtx snapshot (callers reuse the ctx struct per-tri)
} FBCmd;
#define FB_MAXTHREADS 8
static FBCmd *fbCmds = NULL;
static size_t fbCmdN = 0, fbCmdCap = 0;
static int fbThreads = 0;

static FBCmd *fbCmdNew (void) {
  if (fbCmdN == fbCmdCap) {
    fbCmdCap = fbCmdCap ? fbCmdCap*2 : 4096;
    FBCmd *p = realloc (fbCmds, fbCmdCap*sizeof (*fbCmds));
    if (!p) { fbCmdCap = fbCmdN = 0; return NULL; }
    fbCmds = p;
  }
  return &fbCmds[fbCmdN++];
}

int FB_Begin (unsigned long bg) {
  if (!fbEnsure ()) return 0;
  fbCmdN = 0;
  for (long i = 0; i < (long)fbW*fbH; i++) {
    fbPix[i] = (uint32_t)bg; fbZ[i] = 0.f;
  }
  return 1;
}

void FB_Tri (double ax, double ay, double aq,
	     double bx, double by, double bq,
	     double cx2, double cy2, double cq,
	     unsigned long col) {
  FBCmd *c = fbCmdNew ();
  if (!c) return;
  c->type = 0;
  c->x[0] = ax; c->y[0] = ay; c->q[0] = aq;
  c->x[1] = bx; c->y[1] = by; c->q[1] = bq;
  c->x[2] = cx2; c->y[2] = cy2; c->q[2] = cq;
  c->col = col;
}

void FB_TriA (double ax, double ay, double aq, const float *a0,
	      double bx, double by, double bq, const float *a1,
	      double cx2, double cy2, double cq, const float *a2,
	      int na, unsigned long (*shade) (const float *, void *), void *ctx) {
  FBCmd *c = fbCmdNew ();
  if (!c) return;
  c->type = 1;
  c->x[0] = ax; c->y[0] = ay; c->q[0] = aq;
  c->x[1] = bx; c->y[1] = by; c->q[1] = bq;
  c->x[2] = cx2; c->y[2] = cy2; c->q[2] = cq;
  c->na = na;
  for (int k = 0; k < na; k++) {
    c->at[0][k] = a0[k]; c->at[1][k] = a1[k]; c->at[2][k] = a2[k];
  }
  c->shade = shade;
  memset (c->ctx, 0, sizeof c->ctx);
  if (ctx) memcpy (c->ctx, ctx, sizeof c->ctx);
}

void FB_Seg (double ax, double ay, double aq,
	     double bx, double by, double bq,
	     unsigned long col) {
  FBCmd *c = fbCmdNew ();
  if (!c) return;
  c->type = 2;
  c->x[0] = ax; c->y[0] = ay; c->q[0] = aq;
  c->x[1] = bx; c->y[1] = by; c->q[1] = bq;
  c->col = col;
}

// Per-band raster cores: only pixels in rows [y0b, y1b] are written.
static void fbTriRast (const FBCmd *c, int y0b, int y1b) {
  double ax = c->x[0], ay = c->y[0], aq = c->q[0];
  double bx = c->x[1], by = c->y[1], bq = c->q[1];
  double cx2 = c->x[2], cy2 = c->y[2], cq = c->q[2];
  double denom = (bx - ax)*(cy2 - ay) - (cx2 - ax)*(by - ay);
  if (denom == 0) return;
  int xlo = (int)fmin (ax, fmin (bx, cx2)), xhi = (int)ceil (fmax (ax, fmax (bx, cx2)));
  int ylo = (int)fmin (ay, fmin (by, cy2)), yhi = (int)ceil (fmax (ay, fmax (by, cy2)));
  if (xlo < 0) xlo = 0;
  if (ylo < 0) ylo = 0;
  if (xhi >= fbW) xhi = fbW - 1;
  if (yhi >= fbH) yhi = fbH - 1;
  if (ylo < y0b) ylo = y0b;
  if (yhi > y1b) yhi = y1b;
  if (xlo > xhi || ylo > yhi) return;
  double inv = 1./denom;
  for (int y = ylo; y <= yhi; y++) {
    double py = y + 0.5;
    long row = (long)y*fbW;
    for (int x = xlo; x <= xhi; x++) {
      double px = x + 0.5;
      double w0 = ((bx - px)*(cy2 - py) - (cx2 - px)*(by - py))*inv;
      double w1 = ((cx2 - px)*(ay - py) - (ax - px)*(cy2 - py))*inv;
      double w2 = 1. - w0 - w1;
      if (w0 < 0 || w1 < 0 || w2 < 0) continue;
      float q = (float)(w0*aq + w1*bq + w2*cq);
      long idx = row + x;
      if (q > fbZ[idx]) { fbZ[idx] = q; fbPix[idx] = (uint32_t)c->col; }
    }
  }
}

static void fbTriARast (const FBCmd *c, int y0b, int y1b) {
  double ax = c->x[0], ay = c->y[0], aq = c->q[0];
  double bx = c->x[1], by = c->y[1], bq = c->q[1];
  double cx2 = c->x[2], cy2 = c->y[2], cq = c->q[2];
  const float *a0 = c->at[0], *a1 = c->at[1], *a2 = c->at[2];
  int na = c->na;
  double denom = (bx - ax)*(cy2 - ay) - (cx2 - ax)*(by - ay);
  if (denom == 0) return;
  int xlo = (int)fmin (ax, fmin (bx, cx2)), xhi = (int)ceil (fmax (ax, fmax (bx, cx2)));
  int ylo = (int)fmin (ay, fmin (by, cy2)), yhi = (int)ceil (fmax (ay, fmax (by, cy2)));
  if (xlo < 0) xlo = 0;
  if (ylo < 0) ylo = 0;
  if (xhi >= fbW) xhi = fbW - 1;
  if (yhi >= fbH) yhi = fbH - 1;
  if (ylo < y0b) ylo = y0b;
  if (yhi > y1b) yhi = y1b;
  if (xlo > xhi || ylo > yhi) return;
  double inv = 1./denom;
  float at[FB_MAXATTR];
  for (int y = ylo; y <= yhi; y++) {
    double py = y + 0.5;
    long row = (long)y*fbW;
    for (int x = xlo; x <= xhi; x++) {
      double px = x + 0.5;
      double w0 = ((bx - px)*(cy2 - py) - (cx2 - px)*(by - py))*inv;
      double w1 = ((cx2 - px)*(ay - py) - (ax - px)*(cy2 - py))*inv;
      double w2 = 1. - w0 - w1;
      if (w0 < 0 || w1 < 0 || w2 < 0) continue;
      double qd = w0*aq + w1*bq + w2*cq;
      float q = (float)qd;
      long idx = row + x;
      if (q <= fbZ[idx]) continue;
      // Perspective-correct: q = 1/w is affine in screen space but the
      // attributes are not, so interpolate a/w and divide by the interpolated 1/w.
      double p0 = w0*aq/qd, p1 = w1*bq/qd, p2 = w2*cq/qd;
      for (int k = 0; k < na; k++) at[k] = p0*a0[k] + p1*a1[k] + p2*a2[k];
      fbZ[idx] = q; fbPix[idx] = (uint32_t)c->shade (at, (void *)c->ctx);
    }
  }
}

static void fbSegRast (const FBCmd *c, int y0b, int y1b) {
  double ax = c->x[0], ay = c->y[0], aq = c->q[0];
  double bx = c->x[1], by = c->y[1], bq = c->q[1];
  if (fmax (ay, by) < y0b || fmin (ay, by) > y1b) return;
  int steps = (int)fmax (fabs (bx - ax), fabs (by - ay));
  for (int i = 0; i <= steps; i++) {
    double t = steps ? (double)i/steps : 0.;
    int x = (int)(ax + (bx - ax)*t + 0.5), y = (int)(ay + (by - ay)*t + 0.5);
    if (x < 0 || x >= fbW || y < y0b || y > y1b) continue;
    float q = (float)(aq + (bq - aq)*t);
    long idx = (long)y*fbW + x;
    if (q > fbZ[idx]) { fbZ[idx] = q; fbPix[idx] = (uint32_t)c->col; }
  }
}

static void fbBandRun (int y0b, int y1b) {
  for (size_t i = 0; i < fbCmdN; i++) {
    const FBCmd *c = &fbCmds[i];
    switch (c->type) {
      case 0: fbTriRast (c, y0b, y1b); break;
      case 1: fbTriARast (c, y0b, y1b); break;
      default: fbSegRast (c, y0b, y1b); break;
    }
  }
}

static void *fbWorker (void *arg) {
  long b = (long)arg;
  int y0b = (int)((long)fbH*b/fbThreads), y1b = (int)((long)fbH*(b + 1)/fbThreads) - 1;
  fbBandRun (y0b, y1b);
  return NULL;
}

void FB_Flush (void) {
  if (!fbImg || !fbPix) return;
  if (fbThreads == 0) {
    long np = sysconf (_SC_NPROCESSORS_ONLN);
    fbThreads = np > 1 ? (np > FB_MAXTHREADS ? FB_MAXTHREADS : (int)np) : 1;
  }
  if (fbThreads > 1 && fbCmdN) {
    pthread_t th[FB_MAXTHREADS];
    int spawned[FB_MAXTHREADS] = {0};
    for (long b = 1; b < fbThreads; b++)
      spawned[b] = pthread_create (&th[b], NULL, fbWorker, (void *)b) == 0;
    fbBandRun (0, (int)((long)fbH/fbThreads) - 1);
    for (long b = 1; b < fbThreads; b++) {
      if (spawned[b]) pthread_join (th[b], NULL);
      else fbBandRun ((int)((long)fbH*b/fbThreads), (int)((long)fbH*(b + 1)/fbThreads) - 1);
    }
  } else {
    fbBandRun (0, fbH - 1);
  }
  fbCmdN = 0;
  if (fbImg->bits_per_pixel == 32) {
    for (int y = 0; y < fbH; y++) {
      uint32_t *dst = (uint32_t *)(fbImg->data + (long)y*fbImg->bytes_per_line);
      for (int x = 0; x < fbW; x++) dst[x] = fbPix[(long)y*fbW + x];
    }
  } else {
    for (int y = 0; y < fbH; y++)
      for (int x = 0; x < fbW; x++)
	XPutPixel (fbImg, x, y, fbPix[(long)y*fbW + x]);
  }
  XPutImage (prDisplay, draw, prGC, fbImg, 0, 0, x_0, y_0, fbW, fbH);
}

// Pixel value of a 16-bit colour channel in a TrueColor visual: keep the top
// bits that fit the channel mask, as the X server does for XAllocColor.
static unsigned long TrueColorChannel (unsigned short c, unsigned long mask) {
  int shift = 0, bits = 0;
  if (!mask) return 0;
  while (!((mask >> shift) & 1)) shift++;
  while ((mask >> (shift + bits)) & 1) bits++;
  unsigned long v = bits >= 16 ? c : (unsigned long)(c >> (16 - bits));
  return (v << shift) & mask;
}

int WGetColor (unsigned int red, unsigned int green, unsigned int blue, unsigned long *myColor) {
  XColor Color;

  Color.red = red;
  Color.green = green;
  Color.blue = blue;

  // TrueColor visuals (practically all modern displays) map RGB to pixels
  // arithmetically, so compute the pixel locally instead of waiting for an
  // XAllocColor reply from the server for every colour.
  if (fbVisual && fbVisual->class == TrueColor) {
    *myColor = TrueColorChannel (Color.red, fbVisual->red_mask)
             | TrueColorChannel (Color.green, fbVisual->green_mask)
             | TrueColorChannel (Color.blue, fbVisual->blue_mask);
    return 0;
  }

  if (XAllocColor(prDisplay, DefaultColormap (prDisplay, nScreenNum), &Color) == 0) return 1;
  *myColor = Color.pixel;
    
  return 0;
}

void QueryColor(XColor *color) {
    XQueryColor(prDisplay, DefaultColormap (prDisplay, nScreenNum), color);
}

void WSetTitle (const char *s) {
  XTextProperty tProp;
  XStringListToTextProperty ((char**)&s, 1, &tProp);
  XSetWMName (prDisplay, nWnd, &tProp);
}

void DrawArc(int x,int y, unsigned int width, unsigned int height, int a1, int a2) {
    XDrawArc(prDisplay,draw,prGC,x_0+x,y_0+y,width,height,a1,a2);
}

void DrawPoint(int x, int y, int d) {
    XPoint points[4] = { {x_0+x-d, y_0+y-d}, {x_0+x-d, y_0+y+d}, {x_0+x+d, y_0+y+d}, {x_0+x+d, y_0+y-d} };
    XFillPolygon (prDisplay, draw, prGC, points, 4, Convex, CoordModeOrigin); //XFillArc(prDisplay, draw, prGC,x,y,d,d,0,360*64);
}

// Fonts are loaded once per name and cached: SetFont runs several times per
// frame, and reloading leaked an XFontStruct (and a server font) every call.
#define FONT_CACHE 8
static struct { char name[128]; XFontStruct *fs; } fontCache[FONT_CACHE];
static int fontCacheN = 0;

static XFontStruct *LoadFontWithFallback(const char *name) {
    XFontStruct *fs = XLoadQueryFont(prDisplay, name);
    if (fs) return fs;
    int count = 0;
    char **names = XListFonts(prDisplay, "-adobe-*-*-*-*-*-*-*-*-*-*-*-*-*", 1, &count);
    if (names && count > 0) fs = XLoadQueryFont(prDisplay, names[0]);
    if (names) XFreeFontNames(names);
    if (!fs) fs = XLoadQueryFont(prDisplay, "fixed");   // core X font, always present
    fprintf(stderr, "[DGraphX]: font %s not found, using %s\n", name,
            fs ? "a fallback font" : "the server default font");
    return fs;
}

void SetFont(const char* name) {
    int i;
    for (i = 0; i < fontCacheN; i++)
        if (strcmp(fontCache[i].name, name) == 0) break;
    if (i == fontCacheN) {
        XFontStruct *fs = LoadFontWithFallback(name);
        if (fontCacheN == FONT_CACHE) {   // cache full: use without caching
            if (fs) { XSetFont(prDisplay, prGC, fs->fid); XFreeFont(prDisplay, fs); }
            return;
        }
        snprintf(fontCache[i].name, sizeof fontCache[i].name, "%s", name);
        fontCache[i].fs = fs;   // NULL = keep the GC's default font
        fontCacheN++;
    }
    if (fontCache[i].fs) XSetFont(prDisplay, prGC, fontCache[i].fs->fid);
}

static void FreeFonts(void) {
    for (int i = 0; i < fontCacheN; i++)
        if (fontCache[i].fs) XFreeFont(prDisplay, fontCache[i].fs);
    fontCacheN = 0;
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



