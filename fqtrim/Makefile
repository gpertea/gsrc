GDIR := ../gclib

SEARCHDIRS := -I${GDIR}

#SYSTYPE :=     $(shell uname)

CC      := g++

NOTHREADS :=
TRIMDEBUG :=

ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
 NOTHREADS=1
endif

ifneq (,$(findstring fulldebug,$(MAKECMDGOALS)))
 TRIMDEBUG=1
endif

ifneq (,$(findstring trimdebug,$(MAKECMDGOALS)))
 TRIMDEBUG=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CC} -dumpmachine)))
 WINDOWS=1
endif

# MinGW32 GCC 4.5 link problem fix
ifdef WINDOWS
ifneq (,$(findstring 4.5.,$(shell g++ -dumpversion)))
 LFLAGS += -static-libstdc++ -static-libgcc
endif
endif

# Misc. system commands
ifdef WINDOWS
RM = del /Q
else
RM = rm -f
endif

# File endings
ifdef WINDOWS
EXE = .exe
else
EXE =
endif


BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -D_REENTRANT -fno-strict-aliasing -fno-exceptions -fno-rtti

# C/C++ linker

LINKER  := g++

LIBS :=
#LIBS := -lbam -lz

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS += -lpthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

ifneq (,$(findstring release,$(MAKECMDGOALS)))
  CFLAGS := -O2 -DNDEBUG -D_NDEBUG -DNODEBUG $(BASEFLAGS) $(CFLAGS)
  LDFLAGS := $(LDFLAGS)
  #-L${BAM} 
else
  CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(CFLAGS)
  ifdef TRIMDEBUG
     CFLAGS += -DTRIMDEBUG
  endif
  CFLAGS += $(BASEFLAGS)
  LDFLAGS := -g $(LDFLAGS)
  #-L${BAM}
endif

###----- generic build rule

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

OBJS = ${GDIR}/GBase.o ${GDIR}/gdna.o ${GDIR}/GArgs.o ${GDIR}/GStr.o \
 ${GDIR}/GAlnExtend.o

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

.PHONY : all release trimdebug fulldebug nothreads
all: fqtrim
debug:  fqtrim
nothreads: fqtrim
release: fqtrim
fulldebug:  fqtrim
trimdebug:  fqtrim

fqtrim.o ${GDIR}/gdna.o ${GDIR}/GAlnExtend.o: ${GDIR}/GAlnExtend.h ${GDIR}/gdna.h

fqtrim: ${OBJS} ./fqtrim.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
# target for removing all object files

.PHONY : clean release debug nothreads
clean:
	${RM} core core.* fqtrim.exe fqtrim ${OBJS} *.o* *.~*


