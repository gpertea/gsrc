#	makefile for mgblast programs
#  
#>>>> Directory where NCBI toolkit was unpacked & built:
#>>>> (modify this as needed!)
NCBIDIR = /mylocal/src/ncbi

#----------------
#     (based on NCBI's makedemo.unx)
#	$Id: Makefile,v 1.1.1.1 2009/12/05 04:48:57 gpertea Exp $
#
#	Sun with unbundled ANSI compiler [ make CC=acc ]
#	Sun with Gnu C [ make CC=gcc ]
#   Sun with Solaris 2.x OS [ make CC="cc -Xa" OTHERLIBS="-lthread -lm" ]
#	Silicon Graphics [ make OTHERLIBS="-lm -lPW -lsun" ]
#	IBM 3090 with AIX [ make  ]
#	ULTRIX [ make ]
#	Apple AUX [ make SHELL=\usr\bin ]
#	NeXt [ make CFLAGS1="-c -ansi"  ]
#   DEC Alpha under OSF/1  [ make CC=cc ]
#   BeOS [  make OTHERLIBS="" ]
#	
#	Subsequent makes should make ... nocopy instead of make ... all
#	to avoid recopying the include files
#
#   NOTE: It is now possible to build using shared libraries, which 
#     consumes significantly less disk space.  To date, this has only
#     been tested using Solaris.  To use shared libraries, try:
#       make ... NCBI_LINKINGLIBDIR=../shlib
#            NCBI_OTHERLINKOPTIONS=-R`pwd`/../shlib
#

SUFFIXLCL = $(LCL)
CFLAGS1 = -c -O2
LDFLAGS1 = -O2
CC = gcc
RAN = ls -l
OTHERLIBS = -lm
THREAD_OTHERLIBS = -lthread
THREAD_OTHERLIBS =
THREAD_OBJ = ncbithr.o
THREAD_OBJ =

SUFFIXENV = unx
ASNDIR = asnstat
LIB1 = -lncbi
LIB2 = -lncbiobj
LIB3 = -lncbicdr
LIB4 = -lvibrant
LIB4 =
LIB5 = -lncbiacc
LIB20 = -lncbidesk
LIB20 =
LIB22 = -lncbimmdb
LIB23 = -lncbitool
LIBCOMPADJ = -lblastcompadj
LIB50 = -lregexp
LIB60 = -lblast
LIB61 = -lblastapi
NCBI_INCDIR = $(NCBIDIR)/include
NCBI_LIBDIR = $(NCBIDIR)/lib
NCBI_LINKINGLIBDIR = $(NCBI_LIBDIR)
# A kludge which tells the linker not to try to use shared libraries
# when linking programs which use main() instead of Nlm_Main()
USES_MAIN_KLUDGE = -L$(NCBI_LIBDIR)
NCBI_OTHERLINKOPTIONS =
CFLAGS = $(CFLAGS1) -Wno-deprecated-declarations -I$(NCBI_INCDIR)
LDFLAGS = $(LDFLAGS1) -Wno-deprecated-declarations -I$(NCBI_INCDIR) $(NCBI_OTHERLINKOPTIONS) -L$(NCBI_LINKINGLIBDIR)

# ENTREZLIBS = $(LIB5) $(LIB3)

## All things to make
##
all : mgblast

## To clean out the directory without removing make
##
clean :
	- rm -f *.[ao] mgblast mgblast.exe

## Implicit actions
##

.c.o :
	$(CC) $(CFLAGS) $<

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@
MG_LIBS = $(LIB61) $(LIB60) $(LIB23) $(LIBCOMPADJ) $(LIB2) $(LIB1) $(OTHERLIBS) $(THREAD_OTHERLIBS)
MG_OBJ =  ./mgblast.o ./mgblast_util.o $(THREAD_OBJ)
# mgblast
mgblast:  $(MG_OBJ)
	$(CC) $(LDFLAGS) -o mgblast $(MG_OBJ) $(MG_LIBS)
#mgblast : mgblast.c $(THREAD_OBJ)
#	$(CC) -o mgblast $(LDFLAGS) mgblast.c $(THREAD_OBJ) $(LIB61) \
#		$(LIB60) $(LIB23) $(LIBCOMPADJ) $(LIB2) $(LIB1) $(OTHERLIBS) \
#		$(THREAD_OTHERLIBS)
