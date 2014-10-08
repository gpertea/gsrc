# Useful directories

THISCODEDIR := .
GCLDIR := ../gclib

# Directories to search for header files
SEARCHDIRS := -I${THISCODEDIR} -I${GCLDIR}


SYSTYPE :=     $(shell uname)

# C compiler

CC      := g++
CFLAGS  = -O2 -Wall ${SEARCHDIRS} -fno-exceptions -fno-rtti -fno-strict-aliasing -D_REENTRANT

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
all:    nrcl

nrcl:  ./nrcl.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBES}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core nrcl *.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core nrcl *.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o


