##Makefile for AmplitudeSourceLocation_Pulsewidth.F90
##Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
##Copyright: (c) Masashi Ogiso 2020
##License  : MIT License (https://opensource.org/licenses/MIT)

## -DDOUBLE: calculate in double precision (recommended)
## -DV_MEA1D, -DV_CONST: velocity structure defined in set_velocity_model.F90
## -DWIN: using win format waveform file as an input (otherwise sac binary format files are used)
## -DSTDP_COR: convert the depth of station from meters to kilometers (origin: sea level, downward positive)
## -DTESTDATA: do not apply bandpass filter to input waveforms (to be used for synthetic data)
## -DOUT_AMPLITUDE: output text file of observed amplitude (without site correction) for asl_masterevent

#FC = ifort
#FFLAGS =
#DEFS_PW = -DDOUBLE -DV_MEA1D -DWIN -DSTDP_COR
#DEFS_MASTEREVENT = -DDOUBLE -DV_MEA1D 
#DEFS_SYNTH = -DDOUBLE -DV_MEA1D -DSTDP_COR
#INCDIR = -I${NETCDF_FORTRAN_INC}
#LIBDIR = 
#LIBS = -lnetcdff
#OPTS = -assume byterecl -qopenmp -O3 -xHOST

FC = gfortran
FFLAGS = -g -Wall -fbounds-check -fbacktrace
DEFS_PW = -DDOUBLE -DV_MEA1D -DWIN -DSTDP_COR
DEFS_MASTEREVENT = -DDOUBLE -DV_MEA1D 
DEFS_SYNTH = -DDOUBLE -DV_MEA1D -DSTDP_COR
INCDIR = -I/usr/include -I/usr/local/include
LIBDIR = -llapack95 -llapack -lblas
LIBS = -lnetcdff
OPTS = 

asl_pw: nrtype.F90 constants.F90 greatcircle.f90 calc_bpf_coef.f90 calc_bpf_order.f90 tandem.f90 \
	itoa.F90 linear_interpolation.F90 rayshooting.F90 read_sacfile.F90 grdfile_io.F90 set_velocity_model.F90 \
	m_util.f90 m_win.f90 m_winch.f90 AmplitudeSourceLocation_PulseWidth.F90 
	$(FC) $(FFLAGS) $(INCDIR) $^ $(DEFS_PW) $(LIBDIR) $(LIBS) $(OPTS) -o $@

asl_masterevent: nrtype.F90 constants.F90 greatcircle.f90 linear_interpolation.F90 rayshooting.F90 grdfile_io.F90 \
        set_velocity_model.F90 AmplitudeSourceLocation_masterevent.F90 
	$(FC) $(FFLAGS) $(INCDIR) $^ $(DEFS_MASTEREVENT) $(LIBDIR) $(LIBS) $(OPTS) -o $@

asl_synthwave: nrtype.F90 constants.F90 greatcircle.f90 linear_interpolation.F90 rayshooting.F90 read_sacfile.F90 \
        grdfile_io.F90 wavelet.f90 set_velocity_model.F90 xorshift1024star.f90 AmplitudeSourceLocation_synthwave.F90 
	$(FC) $(FFLAGS) $(INCDIR) $^ $(DEFS_SYNTH) $(LIBDIR) $(LIBS) $(OPTS) -o $@

clean:
	rm -f *.mod asl_pw asl_masterevent asl_synthwave
