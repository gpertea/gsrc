# Useful directories

GCLDIR := ../gclib

# Directories to search for header files
SEARCHDIRS :=  -I${GCLDIR}

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

# C compiler

CC      := g++
LINKER := g++

BASEFLAGS  = -Wall ${SEARCHDIRS} -fno-exceptions -fno-rtti -fno-strict-aliasing -D_REENTRANT

ifeq ($(findstring debug,$(MAKECMDGOALS)),)
  CFLAGS = -O2 $(MARCH) $(BASEFLAGS)
  LDFLAGS =
else
  CFLAGS = -g $(BASEFLAGS)
  LDFLAGS = -g
endif

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

LOADLIBS := -lz

.PHONY : all
all:    mgmerge
debug: mgmerge
mgmerge:  ./mgmerge.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBS}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core mgmerge *.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core mgmerge *.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o


