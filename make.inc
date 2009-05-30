#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=0

USEICC=1
USEGCC=0
USECCC=0




SRCS = \
bounds.c boundsint.c coord.c diag.c dudp_calc_3vel.c dudp_calc.c dump.c dump_ener.c fixup.c gaussj.c \
image.c initbase.c init.c interp.c lubksb.c ludcmp.c main.c metric.c \
mnewt.c nrutil.c phys.c ranc.c restart.c set_arrays.c set_grid.c \
step_ch.c tensor.c utoprim.c utoprim_2d.c utoprim_1d.c utoprim_1d_opt.c vchar.c transforms.c fail.c utoprim_ldz.c \
init_mpi.c boundmpi.c boundmpiint.c sources.c utoprim_1d_final.c utoprim_2d_final.c utoprim_5d2_final.c
 
OBJS = \
bounds.o boundsint.o coord.o diag.o dudp_calc_3vel.o dudp_calc.o dump.o dump_ener.o fixup.o gaussj.o \
image.o initbase.o init.o interp.o lubksb.o ludcmp.o main.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o utoprim_2d.o utoprim_1d.o utoprim_1d_opt.o vchar.o transforms.o fail.o utoprim_ldz.o \
init_mpi.o boundmpi.o boundmpiint.o sources.o utoprim_1d_final.o utoprim_2d_final.o  utoprim_5d2_final.o

all:	$(PREP) $(CMD)


$(PREP):
	( sh ./makedecs.h.sh )
	( sh ./makempidecs.h.sh )

$(CMD):	$(OBJS) makefile
	$(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)
	$(CC)  $(CFLAGS) -o $(CMD) $(OBJS) freespace.o $(LDFLAGS)

clean:
	rm *.o

cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~

# dependencies
$(OBJS) : global.h defs.h mympi.h mpidefs.h makefile metric.h
initbase.o : step_ch.h
utoprim_2d.o : utoprim_2d.h utoprim_1d2d.h
utoprim_1d.o : utoprim_1d.h utoprim_1d2d.h
utoprim_1d_opt.o : utoprim_1d_opt.h
step_ch.o : step_ch.h
interp.o : interp.h
utoprim_2d_final.o : u2p_defs.h
utoprim_1d_final.o : u2p_defs.h
utoprim_5d2_final.o : u2p_defs.h



