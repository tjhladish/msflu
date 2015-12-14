EPIFIRE = "/home/tjhladish/work/EpiFire/src/"

msmsc_sim: msmsc_sim.cpp MSMS_Cluster_Sim.h
	g++ --std=c++11 -O2 -I $(EPIFIRE) $(EPIFIRE)*.o MSMS_Cluster_Sim.h msmsc_sim.cpp -o msmsc_sim
