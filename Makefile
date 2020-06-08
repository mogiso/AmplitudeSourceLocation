##Makefile for AmplitudeSourceLocation_Pulsewidth.F90
##Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
##Copyright: (c) Masashi Ogiso 2020
##License  : MIT License (https://opensource.org/licenses/MIT)

FC = ifort
DEFS = -DDOUBLE -DV_MEA1D
INCDIR = -I${NETCDF_FORTRAN_INC}
LIBDIR = 
LIBS = -lnetcdff
OPTS = -assume byterecl -qopenmp -O3 -xHOST

#FC = gfortran
#DEFS = -DDOUBLE -DV_CONST -DWIN
#INCDIR = -I/usr/include
#LIBDIR =
#LIBS = -lnetcdff
#OPTS = -fopenmp -O3

asl_pw: nrtype.F90 constants.F90 greatcircle.f90 calc_bpf_coef.f90 calc_bpf_order.f90 tandem.f90 \
	itoa.F90 linear_interpolation.F90 rayshooting.F90 read_sacfile.F90 grdfile_io.F90 set_velocity_model.F90 \
	m_util.f90 m_win.f90 m_winch.f90 AmplitudeSourceLocation_PulseWidth.F90 
	$(FC) $(INCDIR) $^ $(DEFS) $(LIBDIR) $(LIBS) $(OPTS) -o $@

clean:
	rm -f *.mod asl_pw
