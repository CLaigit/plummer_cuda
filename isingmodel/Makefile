
COMPILER1 =  gcc
COMPILER2 = /usr/local/cuda-5.5/bin/nvcc
OBJS1 = src/ising.c
OBJS2 = src/ising3d.c
OBJS3 = src/ising2d.cu
LIB = -lm
# CFLAG = -arch=sm_20
CFLAG += -Wall

ifdef gflag
 CFLAG += -g
endif

ifdef picture
 CFLAG += -DPICORE=1
endif

default: clean ising ising3d
all: clean ising ising3d cuda

ising:	${OBJS}
	${COMPILER1}	${CFLAG}	$(LIB)	${OBJS1}	-O3	-o	src/ising.o

ising3d:	${OBJS}
	${COMPILER1}	${CFLAG}	$(LIB)	${OBJS2}	-O3	-o	src/ising3d.o

cuda:	${OBJS}
	${COMPILER2}	${CFLAG}	$(LIB)	${OBJS3}	-O3	-o	src/ising2d.o
clean:
	-rm	-f	src/*.o
cleandata:	clean2ddata	clean3ddata
clean2ddata:
	-rm	-f	data/2d/*.dat
clean3ddata:
	-rm	-f	data/3d/*.dat
