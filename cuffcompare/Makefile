# Useful directories

THISCODEDIR := .
GCLDIR := ../gclib
# Directory where libz.a can be found
# (please build that first before making this package) 
# ZDIR := ../zlib
# Directories to search for header files
#SEARCHDIRS := -I${ZDIR} -I${THISCODEDIR} -I${GCLDIR}
SEARCHDIRS := -I${THISCODEDIR} -I${GCLDIR}

SYSTYPE :=     $(shell uname)

# C compiler

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
else
    MARCH = 
endif    

# CVS checked in
CC      := g++
BASEFLAGS  = -Wall -Wextra ${SEARCHDIRS} $(MARCH) \
 -fno-exceptions -fno-rtti -D_REENTRANT 


ifeq ($(findstring debug,$(MAKECMDGOALS)),)
  CFLAGS = -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS = 
  LIBS = 
else
  CFLAGS = -g -DDEBUG $(BASEFLAGS)
  LDFLAGS = -g -Wl,--eh-frame-hdr -L/opt/geo/lib
  LIBS = 
  # use these instead, for HEAP Profiling with gproftools 
  #CFLAGS = -g -DDEBUG -DHEAPROFILE -I/opt/geo/include $(BASEFLAGS)
  #LDFLAGS = -g -Wl,--eh-frame-hdr -L/opt/geo/lib
  #LIBS = -ltcmalloc
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

LINKER    := g++

OBJS = ${GCLDIR}/GFastaIndex.o ${GCLDIR}/GFaSeqGet.o ${GCLDIR}/gff.o \
 ./gtf_tracking.o ${GCLDIR}/gdna.o ${GCLDIR}/codons.o ${GCLDIR}/GBase.o \
 ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o

.PHONY : all
all:    cuffcompare test_load
#all: test_load
debug:  cuffcompare test_load
${GCLDIR}/gff.o  : ${GCLDIR}/gff.h
./gtf_tracking.o : ./gtf_tracking.h
./cuffcompare.o : ./gtf_tracking.h
cuffcompare: ${OBJS} ./cuffcompare.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
test_load: ${OBJS} ./test_load.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

testfai: ${OBJS} ./testfai.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
# target for removing all object files

.PHONY : clean
clean:: 
	@${RM} core core.* testfai cuffcompare cuffcompare.exe ${OBJS} *.o* 


