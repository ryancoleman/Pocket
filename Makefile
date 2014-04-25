#if having problems, run on sgehead, where both 32bit and 64bit libraries are
#installed. doesn't need these once compiled since it is built statically
#
NAME=pocket
DIR=.
NAMEFUL=$(DIR)/$(NAME)
FC = gfortran
FFLAGS = -c -u -O -m32 -static
CC = gcc
CFLAGS = -c -O -static -m32 -L../../GMP/gmp-5.0.5/32bit -I../../GMP/gmp-5.0.5/32bit
LDFLAGS = -O -static -m32 -L../../GMP/gmp-5.0.5/32bit -I../../GMP/gmp-5.0.5/32bit
LIBRARIES= -lgmp

.c.o :
	$(CC) $(CFLAGS) $<

.f.o :
	$(FC) $(FFLAGS) $<

OBJECTS = \
$(NAME).o \
protein_label.o readrad.o \
delcx.o truncate_real.o alpha_fixed.o \
sos_minor_gmp.o alf_tetra_gmp.o tetra_flow_gmp.o \
delaunay_flow.o flow.o depth.o \
measure_vol.o spacefill_vol.o \
pocket_info.o mouth_info.o \
rasmol_script.o kin_script.o alpha_out.o

$(NAMEFUL) : $(OBJECTS)
	$(FC) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTS) $(LIBRARIES)

clean:
	rm -f *.o $(NAMEFUL)

$(OBJECTS) : pocket.h constants.h
