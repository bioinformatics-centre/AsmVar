BAMTOOLS_ROOT = $(CURDIR)/bamtools
CC = g++
CFLAGS = -Wl,-rpath,$(BAMTOOLS_ROOT)/lib  -Wall -O3  -I $(BAMTOOLS_ROOT)/include -L $(BAMTOOLS_ROOT)/lib
#No Warning meseger#CFLAGS = -Wl,-rpath,$(BAMTOOLS_ROOT)/lib  -O3  -I $(BAMTOOLS_ROOT)/include -L $(BAMTOOLS_ROOT)/lib
TABIX = tabix/tabix.o  -L./tabix -ltabix -lz

#Stats For Bases : 
#	$(CC) $(CFLAGS) StatsMain.cpp Stats.cpp Coverage.cpp BamInsert.cpp utility.cpp -lbamtools -o BaseStats $(TABIX)

Stats For Bases Samplely Version :
	$(CC) $(CFLAGS) StatsMain.cpp Stats.cpp Coverage.cpp utility.cpp -lbamtools -o ASV_BaseStatsSimply $(TABIX) 
	#$(CC) $(CFLAGS) StatsMain.cpp Stats.cpp Coverage.cpp utility.cpp -lbamtools -o BaseStatsSimply2 $(TABIX) 
