PROG =  grd_interp
CSRCS = grd_interp.c bilinear.c


NETCDF = $(NETCDFHOME)/lib
NETINC = $(NETCDFHOME)/include
GMTLIB = $(GMTHOME)/lib
GMTINC = $(GMTHOME)/include
FSRCS =
FC = 
CC = gcc

OBJS =  $(CSRCS:.c=.o)	$(FSRCS:.f=.o)
INC = ../../include

CLIBS = -L../../lib -lflex -lm
FLIBS = 
CFLAGS = -O3 -Wall 

$(PROG): $(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(CLIBS) $(FLIBS) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o -I$(INC) -I$(GMTINC) -I$(NETINC)

install:
	\mv $(PROG) ../../bin
clean:
	-rm -f $(OBJS) tags core
