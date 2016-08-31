CPP = g++
CFLAGS = -std=c++11 -O2 -Wno-deprecated-declarations
#CFLAGS = -std=c++11 -Wall -O2
ABCDIR = $(HOME)/work/AbcSmc
SQLDIR = $(ABCDIR)/sqdb
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABCDIR) -I$(SQLDIR) -I$(EPIFIRE)

ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

EPIFIRE = /home/tjhladish/work/EpiFire/src/

all: libabc libepifire msmsc_abc

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

libepifire:
	$(MAKE) -C $(EPIFIRE) -f Makefile

msmsc_abc: libabc libepifire msmsc_abc.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) $(INCLUDE) $(EPIFIRE)*.o $(ABCDIR)/sqlite3.o MSMS_Cluster_Sim.h msmsc_abc.cpp -o msmsc_abc $(ABC_LIB) $(GSL_LIB)

msmsc_sim: msmsc_sim.cpp MSMS_Cluster_Sim.h
	$(CPP) $(CFLAGS) $(INCLUDE) $(EPIFIRE)*.o MSMS_Cluster_Sim.h msmsc_sim.cpp -o msmsc_sim
