

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
reconstructeno.o : reconstructeno.h
advance.o : step_ch.h

