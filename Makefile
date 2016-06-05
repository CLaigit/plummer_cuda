
COMPILER =  /usr/local/cuda-5.5/bin/nvcc
OBJS1 = stable.cu
OBJS2 = unstable.cu
OBJS3 = reduction.cu

CFLAG = -arch=sm_20
LIB = -lm
ifdef gflag
 CFLAG += -g
endif

default: clean stable
all: clean stable unstable
reduction: clean reduction

stable:${OBJS}
${COMPILER} ${CFLAG} $(LIB) ${OBJS1} -o stable.o
unstable:${OBJS}
${COMPILER} ${CFLAG} $(LIB) ${OBJS2} -o unstable.o
reduction:${OBJS}
${COMPILER} ${CFLAG} $(LIB) ${OBJS3} -o reduction.o

clean:
    -rm -f *.o
