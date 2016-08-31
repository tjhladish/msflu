SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)
CPP = g++
CFLAGS = -std=c++11 -O2 -Wno-deprecated-declarations
#CFLAGS = -std=c++11 -Wall -O2
ABCDIR = $(HOME)/work/AbcSmc
SQLDIR = $(ABCDIR)/sqdb
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L $(HPC_GSL_LIB) $(TACC_GSL_LIB) -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABCDIR) -I$(SQLDIR) -I$(EPIFIRE) -I $(HPC_GSL_INC) $(TACC_GSL_INC)

EPIFIRE = /home/tjhladish/work/EpiFire/src/

all: compiler libabc libepifire msmsc_abc

compiler:
ifdef G++VER
CPP = g++-4.9
endif

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

libepifire:
	$(MAKE) -C $(EPIFIRE) -f Makefile

msmsc_abc: libabc libepifire msmsc_abc.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) MSMS_Cluster_Sim.h msmsc_abc.cpp -o msmsc_abc $(INCLUDE) $(EPIFIRE)*.o $(ABC_LIB) $(GSL_LIB)

msmsc_sim: msmsc_sim.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) $(INCLUDE) $(EPIFIRE)*.o MSMS_Cluster_Sim.h msmsc_sim.cpp -o msmsc_sim
