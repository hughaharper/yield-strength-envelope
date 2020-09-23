#
# creates fftfault library from C and fortran code
#
INC = -I/usr/local/include -I/usr/local/include/gmt -I/usr/local/include -I../include
CFLAGS = $(INC)  -O2 -m64
FFLAGS =  -O3 -m64
FC = gfortran
#
OBJECTS=	fourt.o		readgrd.o	writegrd.o	merc.o\
		moment.o	set_litho_defaults.o	temp_plt.o  temp_sleep.o\
		depth_sflr.o	pressure.o	byerlee.o	ductile.o\
		elbendstress.o	yse_moment.o	mechthk.o	yse_moment_old.o

liba:	$(OBJECTS)
	ar rv libflex.a $(OBJECTS)
	rm *.o
	chmod 644 libflex.a
	ranlib libflex.a
