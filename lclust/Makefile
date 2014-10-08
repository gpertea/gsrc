# Useful directories

THISCODEDIR := .
TGICLASSDIR := ../tgi_cl

# Directories to search for header files
SEARCHDIRS := -I- -I${THISCODEDIR} -I${TGICLASSDIR}


SYSTYPE :=     $(shell uname)

# C compiler

CC      := g++
CFLAGS  = -O2 -Wall ${SEARCHDIRS} -fno-exceptions -fno-rtti -D_REENTRANT

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++
LDFLAGS    =
LOADLIBES := 

.PHONY : all
all:    lclust

lclust:  ./lclust.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBES}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core lclust *.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core lclust *.o ${TGICLASSDIR}/gcl/GBase.o ${TGICLASSDIR}/gcl/GString.o ${TGICLASSDIR}/gcl/GArgs.o


