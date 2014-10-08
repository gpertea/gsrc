# Useful directories

#this is the path where the samtools source package was built
#so all the samtools headers (*.h) and (*.a) library files are in there
SAM := /mylocal/src/samtools-0.1.18

GDIR := ../gclib
# Directory where libz.a can be found
# (please build that first before making this package) 
#ZDIR := ../zlib
# Directories to search for header files
#SEARCHDIRS := -I${ZDIR} -I${SAM} -I${GDIR}

SEARCHDIRS := -I${SAM} -I${GDIR}

SYSTYPE :=     $(shell uname)

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686
else
    MARCH = 
endif    

CC      := g++
BASEFLAGS  = -Wall ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-exceptions -fno-rtti -fno-strict-aliasing \
-D_REENTRANT

ifneq (,$(findstring release,$(MAKECMDGOALS)))
  CFLAGS = -O2 -msse2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS = -L${SAM}
else
  CFLAGS = -g -DGDEBUG -DDEBUG $(BASEFLAGS)
  LDFLAGS = -g -L${SAM}
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++
#LDFLAGS    = -L${ZDIR} -lz
#LDFLAGS    = ${LDFLAGS}
LIBS := -lz

OBJS := ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o ${GDIR}/gdna.o ./GapAssem.o

#ifndef NOTHREADS
# OBJS += ${GDIR}/GThreads.o 
#endif

#ifdef GDEBUG
# OBJS += ${GDIR}/proc_mem.o
#endif



.PHONY : all debug release
#all:    mblasm mblaor
#debug:  all
debug : all
release : all
all : bamcons
bamcons :  ./bamcons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/GBam.o ${OBJS} 
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS) -lbam ${LIBS}

mblasm :  ./mblasm.o ${GDIR}/GCdbYank.o ${GDIR}/gcdb.o ${OBJS}
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS) ${LIBS}

./GapAssem.o: GapAssem.h

mblaor :  ./mblaor.o ./GapAssem.o ${GDIR}/GCdbYank.o ${GDIR}/gcdb.o ${OBJS}
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS) ${LIBS}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core mblasm mblasm.exe *.o ${GDIR}/GCdbYank.o ${GDIR}/gdna.o ${GDIR}/gcdb.o ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core mblasm mblasm.exe *.o ${GDIR}/GCdbYank.o ${GDIR}/gdna.o ${GDIR}/gcdb.o ${GDIR}/GBase.o ${GDIR}/GStr.o ${GDIR}/GArgs.o


