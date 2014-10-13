########################################################################
#                                                                      #
#   Copyright (C) 2008 by University of Pennsylvania,                  #
#   Philadelphia, PA USA. All rights reserved.                         #
#                                                                      #
########################################################################
LIBS = -L/netusr/lib64/R/lib -lR

.PHONY: utest clean add

CC=gcc
IDIR =..
DEPS = $(wildcard ${IDIR}/*.h) 

SRCS =$(wildcard *.c)
ALL_TEST_FUNCS = $(wildcard test_*.c)
OBJS = $(SRCS:%.c=%.o)

APP_SRCS = $(wildcard ../*.c)
_APP_SRCS = $(subst ../,,$(APP_SRCS))    # subsitutes ../ with an empty string, and then changes .c with .o
APP_OBJS = $(_APP_SRCS:%.c=%.o)

CFLAGS = -g3 -D_CUTEST_ -std=gnu99 -I/netusr/lib64/R/include  -I/usr/local/include

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
