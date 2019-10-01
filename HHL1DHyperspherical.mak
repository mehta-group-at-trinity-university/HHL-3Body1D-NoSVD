HHL1DHyperspherical.x: modules_qd.o besselnew.o Bsplines.o  ../bspllib_22/bspline90_22.o
	gfortran -O -ffixed-line-length-132 Bsplines.o modules_qd.o ../bspllib_22/bspline90_22.o besselnew.o -L/usr/local/lib/ -L/Users/mehtan/Code/ARPACK/ARPACK -larpack_OSX -framework accelerate -lm HHL1DHyperspherical.f -o HHL1DHyperspherical.x

#1Dpot.o:	1Dpot.f
#	gfortran -O  -ffixed-line-length-132 -c 1Dpot.f

#HHL1DHyperspherical.o:	HHL1DHyperspherical.f
#	gfortran -O  -ffixed-line-length-132 -c HHL1DHyperspherical.f

Bsplines.o:	Bsplines.f
	gfortran -O  -ffixed-line-length-132 -c Bsplines.f	

besselnew.o :	besselnew.f
	gfortran -O  -ffixed-line-length-132 -c besselnew.f	

modules_qd.o :	modules_qd.f90
	gfortran -O  -c modules_qd.f90	

#matrix_stuff.o:	matrix_stuff.f
#	gfortran -O -ffixed-line-length-132 -c matrix_stuff.f
