

ifeq ($(USEMPI),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif # endif usempich==1

ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)
endif # endif usempich==0

include maketail.harm.inc

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



include maketail.harm.deps.inc

