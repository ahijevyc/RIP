OBJS = 	average_ripdp_by_time_of_day.o \
	kwm_date_utilities.o \
	readrip.o

CMD  = 	avgrip

F90  =	pgf90
FFLAGS = -Mfree -byteswapio

all: $(CMD)

$(CMD):$(OBJS)
	$(F90) -o $(CMD) $(OBJS)

.F.o:
	$(F90) -c $(FFLAGS) $(<)

clean:
	rm -f $(OBJS) *~ *.mod $(CMD)

average_ripdp_by_time_of_day.o:	readrip.o
average_ripdp_by_time_of_day.o:	kwm_date_utilities.o
