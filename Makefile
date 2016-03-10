
COMPILER =  /usr/local/cuda-5.5/bin/nvcc
OBJS = stable.cu
CFLAG = -arch=sm_20
LIB = -lm
ifdef gflag
 CFLAG += -g
endif

default: clean stable
all: clean stable unstable

stable:${OBJS}
    ${COMPILER} ${CFLAG} $(LIB) ${OBJS} -o stable
unstable:${OBJS}
    ${COMPILER} ${CFLAG} $(LIB) ${OBJS} -o unstable

clean:
    -rm -f *.o
