FC = gfortran -ffixed-line-length-none
FFLAGS = -O
CC = gcc
CFLAGS = -O -Wno-unused-result

CAP  = cap3D mtdcmp radpttn

SUBS = fft.o Complex.o radiats.o grid3d.o futterman.o sacio.o trap.o

all: $(CAP)

cap3D : %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio

mtdcmp: mtdcmp.o
	$(LINK.f) -o $@ $@.o

radpttn: radpttn.o radiats.o
	$(LINK.c) -o $@ $@.o radiats.o -lm

%/mt.all: %/weight.dat
	for dp in 10 15 20; do cap3D.pl -G${HOME}/data/models/Glib -H0.2 -P300000/65 -Q0.02 -T50/100 -S2/5/0 -D2/1/0.5 -C0.05/0.3/0.02/0.1 -W1 -X10 -Mcus_$$dp $*; head -1 $*/cus_$$dp.out >> $@; done

%/mt.best: %/mt.all
	(depth.pl $< $* > $*/depth.ps) >& $@

clean:
	rm -f $(CAP) *.o
