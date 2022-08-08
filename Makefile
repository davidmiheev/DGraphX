NAME            = DGraphX

CC              = clang
LINK            = clang

DEBUG           = 
OPTIMAZE_COMMON   = -O2 
OPTIMAZE_SPECIFIC =
OPTIMAZE        = $(OPTIMAZE_COMMON) $(OPTIMAZE_SPECIFIC)

WARNINGS        = -W -Wall
INCLUDE_DIR     = -I. -I/opt/X11/include 
CFLAGS          = -c $(DEBUG) $(OPTIMAZE) $(WARNINGS) $(INCLUDE_DIR) 

LIB_DIR         = -L. -L/opt/X11/lib
LIB             = -lX11 
LDFLAGS         = $(DEBUG)

OBJS            = Xfun.o DGraphX.o interpolation.o dtoa.o DX.o matrix.o

all:    $(NAME)

$(NAME)         : $(OBJS)
	$(LINK) $(LDFLAGS) $(OBJS) $(LIB_DIR) $(LIB) -lm -o $(NAME)

clean:
	rm -f $(OBJS) $(NAME)

.c.o:
	$(CC) $(CFLAGS) $<

#-------------------------------------------------------------------------------
Xfun.o      : Xfun.c DGraphX.h Interp2D.h

DGraphX.o        : DGraphX.c DGraphX.h Interp2D.h matrix.h

interpolation.o		: interpolation.c Interp2D.h DGraphX.h 

dtoa.o		: dtoa.c Interp2D.h DGraphX.h

DX.o		: DX.c DGraphX.h matrix.h 

matrix.o	: matrix.c DGraphX.h matrix.h

#-------------------------------------------------------------------------------
