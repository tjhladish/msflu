SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)
CPP = g++
CFLAGS = -std=c++11 -O2
#CFLAGS = -std=c++11 -Wall -O2
ABCDIR = $(HOME)/work/AbcSmc
SQLDIR = $(ABCDIR)/sqdb
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABCDIR) -I$(SQLDIR) -I$(EPIFIRE)

EPIFIRE = "/home/tjhladish/work/EpiFire/src/"

all: compiler libabc msmsc_abc

compiler:
ifdef G++VER
CPP = g++-4.9
endif

libabc:  
	$(MAKE) -C $(ABCDIR) -f Makefile

msmsc_abc: msmsc_abc.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) $(INCLUDE) $(EPIFIRE)*.o MSMS_Cluster_Sim.h msmsc_abc.cpp -o msmsc_abc $(ABC_LIB) $(GSL_LIB)

msmsc_sim: msmsc_sim.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) $(INCLUDE) $(EPIFIRE)*.o MSMS_Cluster_Sim.h msmsc_sim.cpp -o msmsc_sim
