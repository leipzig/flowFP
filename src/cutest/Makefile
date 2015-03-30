########################################################################
#                                                                      #
#   Copyright (C) 2008 by University of Pennsylvania,                  #
#   Philadelphia, PA USA. All rights reserved.                         #
#                                                                      #
########################################################################
LIBS = -L$$R_HOME/Libraries -lR

.PHONY: utest clean add

CC=gcc
IDIR = ..
DEPS = ${IDIR}/flowFP.h

ALL_TEST_FUNCS = test_example.c test_idx_lut.c test_median.c test_variance.c
SRCS = CuTest.c ${ALL_TEST_FUNCS}
OBJS = $(SRCS:%.c=%.o)

APP_SRCS = ../flowFP.c ../init.c ../split_utils.c
APP_OBJS = flowFP.o init.o split_utils.o

CFLAGS = -g3 -D_CUTEST_ -std=gnu99 -I$$R_HOME/Headers  -I/usr/local/include

Alltests utest: ${APP_OBJS} ${OBJS}
	./make-tests.sh ${ALL_TEST_FUNCS} >Alltests.c
	${CC} -c -o Alltests.o Alltests.c ${CFLAGS} -I ${IDIR}
	@echo "Compiling and linking Alltests"
	$(CC) $(CFLAGS) $(LIBS) -o Alltests *.o
	@echo "Running Alltests"
	./Alltests

%.o: ../%.c ${DEPS}
	${CC} -c -o $@ $< ${CFLAGS} -I ${IDIR} ${INC_FLAGS}


clean:
	rm -f core.* *.o Alltests.c Alltests

add:
	-@if [ -d CSV ] ; then cvs add test_*.c; fi

# END

