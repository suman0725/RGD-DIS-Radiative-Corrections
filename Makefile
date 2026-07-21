# RHEL5 (32 bit)
CERN_ROOT = /u/scigroup/cvmfs/scicomp/sw/el9/cernlib/2023
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       := gfortran
########################################

externals_all_objs = externals_all.o sigmodel_calc_2013.o fy.o target_info.o \
	             F1F209.o f_to_sig.o y_calc.o sig_bar_df.o nform_xem.o
externals_all_srcs = externals_all.f sigmodel_calc_2013.f fy.f target_info.f \
                     F1F209.f f_to_sig.f y_calc.f sig_bar_df.f nform_xem.f

########################################

none: externals_all

all:  externals_all

##externals_all:
##	$(F77) $(FFLAGS) -o jtest $(externals_all_srcs) $(OTHERLIBS)

externals_all.o: externals_all.f
		 $(F77) $(FFLAGS) -c $< -o $@

fy.o: fy.f
	$(F77) $(FFLAGS) -c $< -o $@

sigmodel_calc_2013.o: sigmodel_calc_2013.f
	$(F77) $(FFLAGS) -c $< -o $@

bdisnew4he3.o: bdisnew4he3.f
	$(F77) $(FFLAGS) -c $< -o $@

target_info.o: target_info.f
	$(F77) $(FFLAGS) -c $< -o $@

#f1f2in06.o: f1f2in06.f
#	$(F77) $(FFLAGS) -c $< -o $@

F1F209.o: F1F209.f
	$(F77) $(FFLAGS) -c $< -o $@

#f1f209newmecnof1cor.o: f1f209newmecnof1cor.f
#	$(F77) $(FFLAGS) -c $< -o $@

f_to_sig.o: f_to_sig.f
	$(F77) $(FFLAGS) -c $< -o $@

y_calc.o: y_calc.f
	$(F77) $(FFLAGS) -c $< -o $@

sig_bar_df.o: sig_bar_df.f
	$(F77) $(FFLAGS) -c $< -o $@

nform_xem.o: nform_xem.f
	$(F77) $(FFLAGS) -c $< -o $@

externals_all: $(externals_all_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(externals_all_objs) $(OTHERLIBS)


clean:
	rm -f *.o externals_all
