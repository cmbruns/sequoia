VERSION = 0.9.9a
CFLAGS1 = -DVERSION=\"$(VERSION)\"

# Mac OS X Tiger
# compile on koobi
RM = /bin/rm -f
CC = gcc-3.3
C++ = g++-3.3
#CFLAGS = $(CFLAGS1) -b i486-linuxaout -g -Wall
CFLAGS = $(CFLAGS1) -O3 -g -Wall -DMAC_OS_X=1
LIBS = -lstaticreadline -lncurses
XLIBS = -L /usr/X11/lib -lXm -lXt -lX11 -lXext -lstdc++ -lreadline -lcurses -lm
LD = gcc-3.3

# linux
# compile on baxter
# RM = /bin/rm -f
# CC = gcc
# C++ = g++
# #CFLAGS = $(CFLAGS1) -b i486-linuxaout -g -Wall
# CFLAGS = $(CFLAGS1) -g -O3 -Wall -static
# LIBS = -lreadline -lncurses
# XLIBS = -L /usr/X11/lib -lXm -lXt -lX11 -lXext -lstdc++ -lreadline -lcurses -lm
# LD = gcc

# OSF1 V4.0
# compile on castle
# # [compile on alf1.mrc-lmb.cam.ac.uk]
# RM = /bin/rm -f
# CC = gcc
# C++ = g++
# CFLAGS = $(CFLAGS1) -g -Wall -static -I/u/binf/cbruns/devel/include -L/u/binf/cbruns/devel/lib
# LIBS = -lreadline -lcurses
# LIBS = -lXm -lXt -lX11 -lXext -lstdc++ -lreadline -lcurses -lm
# LD = gcc

## SunOS5.6
## compile on spicerack
#RM = /bin/rm -f
#CFLAGS = $(CFLAGS1) -s -O2 -fguiding-decls -I /tools/gnu/include/ -L /tools/gnu/lib/ -DNO_READLINE_HISTORY
#XLIBS = -lsocket -lnsl -lw -ldl -lstdc++ -lXm -lXt -lX11 -lXext -lm
#LIBS = -lstdc++ -lreadline -lcurses -lm
#CC = gcc
#C++ = gcc

# IRIX6.5
# compile on pilgrim
#RM = /bin/rm -f
#CC = gcc
#C++ = g++
##CFLAGS = $(CFLAGS1) -g -fguiding-decls
#CFLAGS = $(CFLAGS1) -O3 -fguiding-decls
#LIBS = -lstdc++ -lreadline -lcurses -lm
#XLIBS = -lXm -lXt -lX11 -lXext -lstdc++ -lm
#LD = gcc

#IRIX5.3
# compile on si
#RM = /bin/rm -f
##CFLAGS = $(CFLAGS1) -g -fguiding-decls -Wall
#CFLAGS = $(CFLAGS1) -O2 -fguiding-decls -I /asd/d3p81/gnu/irix5/include/ -L /asd/d3p81/gnu/irix5/lib/
#XLIBS = -lXm -lXt -lX11 -lXext -lstdc++ -lm
#LIBS = -lstdc++ -lreadline -lcurses -lm
#CC = gcc
#C++ = g++
#LD = gcc

# HPUX
#RM = /bin/rm -f
#CC = gcc
#C++ = g++
##CFLAGS = $(CFLAGS1) -g -Wall
#CFLAGS = $(CFLAGS1) -O4
#LIBS = -lstdc++ -lm
#LD = gcc

# ms-dos/djgpp
# compile on atreus
#RM = del
#CC = gcc
#C++ = gxx
##CFLAGS = $(CFLAGS1) -s -O5 -Wall -fguiding-decls
#CFLAGS = $(CFLAGS1) -s -O5 -Wall -static
##LIBS = -lgpp -lm
#LIBS = -lstdcx -lreadline -lm
#LD = gcc

# object files only for X interface
WITH_X_OBJS = xoshook.o xsequoia.o mtwin.o seqwin.o matwin.o strwin.o cvector.o

# object files only for tty interface
NO_X_OBJS = oshook.o sequoia.o

# object files for matrix data types
MATRIX_OBJS = array.o vector.o matrix.o rmat.o sequence.o mystring.o \
              brookhav.o path_mat.o parse.o overlay.o comp_mat.o ltable.o parse2.o \
	      instanti.o instant2.o

#object files in common
BOTH_OBJS =  $(MATRIX_OBJS) consens.o readline.o \
	geoscore.o timelice.o sort.o \
	variable.o

CK_OBJS = $(BOTH_OBJS) $(NO_X_OBJS)

CK_HDRS = array.hxx array2d.hxx matrix.hxx sequence.hxx \
	comp_mat.hxx path_mat.hxx parse.hxx mystring.hxx \
	vector.hxx brookhav.hxx overlay.hxx rmat.hxx \
	rvec.hxx cmbmacro.hxx consens.hxx oshook.hxx variable.hxx

XCK_OBJS = $(BOTH_OBJS) $(WITH_X_OBJS)

XCK_HDRS = array.hxx array2d.hxx matrix.hxx sequence.hxx \
	comp_mat.hxx path_mat.hxx parse.hxx mystring.hxx \
	vector.hxx brookhav.hxx overlay.hxx rmat.hxx \
	rvec.hxx oshook.hxx cmbmacro.hxx consens.hxx xsequoia.hxx \
	mtwin.h

default: sequoia

all: sequoia xsequoia

2Dgeom.o: 2Dgeom.cxx 2Dgeom.hxx vector.hxx array.hxx
	$(C++) $(CFLAGS) -c 2Dgeom.cxx

array.o: array.cxx array.hxx cmbmacro.hxx
	$(C++) $(CFLAGS) -c array.cxx

array2d.o: array.o array2d.hxx array2d.cxx
	$(C++) $(CFLAGS) -c array2d.cxx

brookhav.o: brookhav.hxx brookhav.cxx vector.cxx rvec.cxx rmat.cxx mystring.cxx sequence.hxx
	$(C++) $(CFLAGS) -c brookhav.cxx

consens.o: $(CK_HDRS) consens.cxx
	$(C++) $(CFLAGS) -c consens.cxx

clean:
	$(RM) *.o 
	$(RM) *.obj
	$(RM) *.bak
	$(RM) core
	$(RM) *~

cmbmacro.o: cmbmacro.hxx cmbmacro.cxx
	$(C++) $(CFLAGS) -c cmbmacro.cxx

comp_mat.o: comp_mat.cxx comp_mat.hxx matrix.hxx
	$(C++) $(CFLAGS) -c comp_mat.cxx

cvector.o: cvector.c cvector.h
	$(CC) $(CFLAGS) -c cvector.c

equivpdb: $(EQ_OBJS) $(EQ_HDRS) equivpdb.cxx
	$(C++) $(CFLAGS) $(EQ_OBJS) equivpdb.cxx -o equivpdb -lm

geoscore.o: geoscore.cxx geoscore.hxx ltable.hxx
	$(C++) $(CFLAGS) -c geoscore.cxx

instanti.o: array.cxx array2d.cxx vector.cxx matrix.cxx brookhav.hxx instanti.cxx rvec.hxx parse.hxx
	$(C++) $(CFLAGS) -c instanti.cxx

instant2.o: cmbmacro.hxx cmbmacro.cxx instant2.cxx
	$(C++) $(CFLAGS) -c instant2.cxx

matrix.o: vector.o matrix.hxx matrix.cxx
	$(C++) $(CFLAGS) -c matrix.cxx

ltable.o: ltable.cxx ltable.hxx
	$(CC) $(CFLAGS) -c ltable.cxx

matwin.o: matwin.c
	$(CC) $(CFLAGS) -c matwin.c

mtwin.o: mtwin.c mtwin.h
	$(CC) $(CFLAGS) -c mtwin.c

mystring.o: mystring.hxx mystring.cxx
	$(C++) $(CFLAGS) -c mystring.cxx

overlay.o: overlay.hxx overlay.cxx rvec.hxx rmat.hxx brookhav.hxx
	$(C++) $(CFLAGS) -c overlay.cxx

oshook.o: oshook.hxx oshook.cxx
	$(C++) $(CFLAGS) -c oshook.cxx

parse.o: parse.cxx parse.hxx oshook.hxx mystring.hxx matrix.hxx \
 sequence.hxx comp_mat.hxx path_mat.hxx brookhav.hxx overlay.hxx \
 consens.hxx
	$(C++) $(CFLAGS) -c parse.cxx

parse2.o: parse2.cxx mystring.hxx array.hxx parse.hxx
	$(C++) $(CFLAGS) -c parse2.cxx

path_mat.o: path_mat.cxx path_mat.hxx matrix.hxx sequence.hxx \
              brookhav.hxx geoscore.hxx
	$(C++) $(CFLAGS) -c path_mat.cxx

rmat.o: rmat.cxx rmat.hxx matrix.cxx ltable.hxx
	$(C++) $(CFLAGS) -c rmat.cxx

readline.o: readline.c readline.h
	$(CC) $(CFLAGS) -c readline.c

rvec.o: rvec.cxx rvec.hxx vector.cxx
	$(C++) $(CFLAGS) -c rvec.cxx

sequence.o: sequence.cxx sequence.hxx mystring.hxx comp_mat.hxx
	$(C++) $(CFLAGS) -c sequence.cxx

sequoia: $(CK_OBJS) $(CK_HDRS)
	$(C++) $(CFLAGS) $(CK_OBJS) -o sequoia $(LIBS)

sequoia.exe: sequoia
	coff2exe -s d:/djgpp/bin/go32.exe sequoia

sequoia.o: $(CK_HDRS) sequoia.cxx parse.hxx timestamp
	$(C++) $(CFLAGS) -c sequoia.cxx

seqwin.o: seqwin.c
	$(CC) $(CFLAGS) -c seqwin.c

strwin.o: strwin.c strctyp.h cvector.h
	$(CC) $(CFLAGS) -c strwin.c

sort.o: sort.cc sort.h
	$(C++) $(CFLAGS) -c sort.cc

timelice.o: timelice.cxx timelice.hxx
	$(C++) $(CFLAGS) -c timelice.cxx

timestamp:
	 echo 'char * timestamp = "'`date`' on host '`uname -a`'";' > timestmp.h

variable.o: variable.cxx variable.hxx
	$(C++) $(CFLAGS) -c variable.cxx

vector.o: vector.cxx vector.hxx array.cxx
	$(C++) $(CFLAGS) -c vector.cxx

xoshook.o: oshook.hxx xoshook.cxx strctyp.h cvector.h parse.hxx
	$(C++) $(CFLAGS) -c xoshook.cxx

xsequoia: $(XCK_OBJS) $(XCK_HDRS)
	$(CC) $(CFLAGS) $(XCK_OBJS) -o xsequoia $(XLIBS)

xsequoia.o: $(XCK_HDRS) xsequoia.cxx
	$(C++) $(CFLAGS) -c xsequoia.cxx

#

vartest: vartest.cxx $(BOTH_OBJS) oshook.o
	$(CC) $(CFLAGS) $(BOTH_OBJS) oshook.o vartest.cxx -o vartest $(LIBS)

test: test.cxx array.cxx mystring.cxx sequence.cxx
	$(CC) $(CFLAGS) -g test.cxx -o test $(LIBS)

testread: testread.cc readline.o
	$(CC) $(CFLAGS) testread.cc readline.o -o testread -lreadline -lstdc++ -lcurses
