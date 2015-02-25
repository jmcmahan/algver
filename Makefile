# Simple make file that compiles any .c file in the directory as a 
# self-contained executable.

#CC=g++
VPATH=/home/jerry/files/Dropbox/math/research/ncsu_postdoc/queso_algver_allcase
CC=mpic++
LIBS=-lm -lgsl -lgslcblas -L/usr/local/lib -lqueso 
CFLAGS=-I/usr/local/include -O2 
#OBJS=main.o algver.o linalgwrapper.o rngwrapper.o
OBJS=main.o algver.o linalgwrapper.o rngwrapper.o JeffreysJointPdf.o JeffreysVectorRealizer.o JeffreysVectorRV.o InformJointPdf.o InformVectorRealizer.o InformVectorRV.o




.cpp:
	$(CC) $(CFLAGS) $(LIBS) $< -o $@

all: vertest

vertest: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o vertest

#main.o: main.cpp main.hpp algver.hpp linalgwrapper.hpp rngwrapper.hpp
main.o: main.cpp main.hpp algver.hpp linalgwrapper.hpp rngwrapper.hpp JeffreysJointPdf.h JeffreysVectorRealizer.h JeffreysVectorRV.h InformJointPdf.h InformVectorRealizer.h InformVectorRV.h
	$(CC) $(CFLAGS) -c $(VPATH)/main.cpp

algver.o: algver.cpp algver.hpp linalgwrapper.hpp
	$(CC) $(CFLAGS) -c $(VPATH)/algver.cpp

linalgwrapper.o: linalgwrapper_gsl.cpp linalgwrapper.hpp
	$(CC) $(CFLAGS) -c $(VPATH)/linalgwrapper_gsl.cpp -o linalgwrapper.o

rngwrapper.o: rngwrapper.cpp rngwrapper.hpp
	$(CC) $(CFLAGS) -c $(VPATH)/rngwrapper.cpp -o rngwrapper.o

JeffreysJointPdf.o: JeffreysJointPdf.C JeffreysJointPdf.h
	$(CC) $(CFLAGS) -c $(VPATH)/JeffreysJointPdf.C -o JeffreysJointPdf.o

JeffreysVectorRealizer.o: JeffreysVectorRealizer.C JeffreysVectorRealizer.h
	$(CC) $(CFLAGS) -c $(VPATH)/JeffreysVectorRealizer.C -o JeffreysVectorRealizer.o

JeffreysVectorRV.o: JeffreysVectorRV.C JeffreysVectorRV.h
	$(CC) $(CFLAGS) -c $(VPATH)/JeffreysVectorRV.C -o JeffreysVectorRV.o

InformJointPdf.o: InformJointPdf.C InformJointPdf.h
	$(CC) $(CFLAGS) -c $(VPATH)/InformJointPdf.C -o InformJointPdf.o

InformVectorRealizer.o: InformVectorRealizer.C InformVectorRealizer.h
	$(CC) $(CFLAGS) -c $(VPATH)/InformVectorRealizer.C -o InformVectorRealizer.o

InformVectorRV.o: InformVectorRV.C InformVectorRV.h
	$(CC) $(CFLAGS) -c $(VPATH)/InformVectorRV.C -o InformVectorRV.o


.PHONY: clean

clean:
	rm $(OBJS) vertest

