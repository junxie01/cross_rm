FC=gfortran
#FFLAG=-ffixed-line-length-none -fbounds-check
saclib1=/home/junxie/opt/sac/lib/sacio.a
saclib2=/home/junxie/opt/sac/lib/libsac.a
#saclib1=/Users/junxie/opt/sac/lib/sacio.a
#saclib2=/Users/junxie/opt/sac/lib/libsac.a
all:cross_rm_win
objects=cross_rm_win.o cross.o filter.o usage.o gauss.o rmean.o rtr.o taper.o
cross_rm_win: $(objects)
	$(FC) $(objects) $(saclib1) $(saclib2) -o $@ $(FFLAG)
%.o: %.f
	$(FC) $(FFLAG) -c $< 
%.o: %.f90
	$(FC) $(FFLAG) -c $< 
clean:
	-rm $(objects)
install:
	cp cross_rm_win ../bin
