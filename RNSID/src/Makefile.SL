# makefile for RNS

EXECTBL  = RNS
EXECTBL2 = RNS_readID
ARCHIVE  = cactus_rns.tgz

SRCS = rnsid_util.c equil.c equil_util.c  hdf5_save.c
OBJS = rnsid_util.o equil.o equil_util.o  hdf5_save.o

ARCFILES = Makefile rnsid.c rnsid_util.c equil.c equil_util.c main_seq.c hdf5_save.c \
           consts.h equil.h equil_util.h nrutil.h  rnsid.h  rnsid_util.h \
           eosSLy

#CC	= pgcc -DRNS_SEQ_COMPILATION -I./include -O0 -g
CC	= gcc-mp-6 -O3 -g
CCFLAGS = -I./include -I/storage5/ExternaLib/gcc/include  -DRNS_SEQ_COMPILATION -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -I/opt/local/include/
LLFLAGS = -L/opt/local/lib -lhdf5_hl -lhdf5 -lz  -lm

####all: $(ARCHIVE) $(EXECTBL) 

all:  $(EXECTBL) $(EXECTBL2)

$(EXECTBL): $(OBJS) RNS.o
	@echo ''
	@echo '-------------------------------------------------' 
	@echo ' Building the executable'
	$(CC) -g -o $(EXECTBL) rnsid_util.o equil.o equil_util.o RNS.o hdf5_save.o $(LLFLAGS) 
	@echo ''
	@echo ' All Done !'

$(EXECTBL2): $(OBJS)  RNS_readID.o
	@echo ''
	@echo '-------------------------------------------------' 
	@echo ' Building the executable'
	$(CC) -g -o $(EXECTBL2) rnsid_util.o equil.o equil_util.o  RNS_readID.o hdf5_save.o $(LLFLAGS) 
	@echo ''
	@echo ' All Done !'

%.o: %.c
	@echo '' 
	@echo '--------------------------------------------------'
	@echo ' Building' $*.o	
	$(CC) -c $(CCFLAGS) $(*).c
	@echo ''	
	@echo ' Done'

$(ARCHIVE):
	@echo ''
	@echo '-------------------------------------------------' 
	@echo ' Building the tar archive'
	tar -zcvf $(ARCHIVE) $(ARCFILES)
	@echo ''	
	@echo ' Done'

# DO NOT DELETE
