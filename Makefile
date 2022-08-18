##Makefile for AmplitudeSourceLocation_Pulsewidth.F90
##Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
##Copyright: (c) Masashi Ogiso 2020
##License  : MIT License (https://opensource.org/licenses/MIT)

## -DDOUBLE: calculate in double precision (recommended)
## -DV_MEA1D, -DV_CONST: velocity structure defined in set_velocity_model.F90

## -DWIN: using win format waveform file as an input for AmplitudeSourceLocation_PulseWidth.F90 and
##        AmplitudeSourceLocation_masterevent.F90
## -DSAC: use sac-formatted waveform file as an input waveform for AmplitudeSourceLocation_PulseWidth.F90 (default)
##        and AmplitudeSourceLocation_masterevent.F90
## -DAMP_TXT: using amplitude and station parameters written in textfiles instead of win- or sac-formatted waveforms for
##            AmplitudeSourceLocation_PulseWidth.F90
## -DAMP_RATIO: using ampliude ratio instead of amplitude itself (Taisne et al., 2011; Ichihara and Matsumoto, 2017) in
##              AmplitudeSourceLocation_PulseWidth.F90
## -DTESTDATA: do not apply bandpass filter and site correction to input waveforms (to be used for synthetic data) 
## -DOUT_AMPLITUDE: output text file of observed amplitude (without site correction) for asl_masterevent
## -DWITHOUT_TTIME: neglect travel time from assumed source location (or reference source location) to each station
##                  when extract the amplitude from waveform array for AmpliutudeSourceLocation_PulseWidth.F90
##                  (or AmplitudeSourceLocation_masterevent.F90)
## -DDAMPED: Damped least squares solution (check damp() array) for AmplitudeSourceLocation_masterevent.F90
## -DMKL: use MKL; otherwise use lapack95 in AmplitudeSourceLocation_masterevent.F90 and TraveltimeSourceLocation_masterevent.F90
## -DEACH_ERROR: estimate location errors separately in AmplitudeSourceLocation_masterevent.F90 and
##               TraveltimeSourceLocation_masterevent.F90

FC = ifort
FFLAGS = -g -traceback -assume byterecl -CB
#FFLAGS = -traceback -assume byterecl
DEFS = -DDOUBLE -DMKL -DVEL_MEA1D -DEACH_ERROR -DOUT_AMPLITUDE -DSAC -DTESTDATA -DDAMPED
INCDIR = -I${NETCDF_FORTRAN_INC} -I${MKLROOT}/include/intel64/lp64 -I.
LIBDIR = -L${MKLROOT}/lib/intel64 -L${NETCDF_FORTRAN_LIB}
LIBS = -lnetcdff -liomp5 -lpthread -lmkl_core -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_intel_thread
OPTS = -O3 -xHOST -ipo -mcmodel=large

#FC = gfortran
#FFLAGS = -g -Wall -fbounds-check -fbacktrace
#FFLAGS = -fbacktrace
#DEFS = -DDOUBLE -DV_MEA1D -DWIN
#INCDIR = -I/usr/include -I/usr/local/include
#LIBDIR = 
#LIBS = -lnetcdff -llapack95 -llapack -lblas
#OPTS = -O3
#OPTS = -O3 -fopenmp

TARGET		= asl_pw asl_masterevent asl_synthwave ttime_masterevent asl_masterevent_shmdump \
                  asl_masterevent_ttimecor asl_dd
SRCS		= AmplitudeSourceLocation_PulseWidth.F90 AmplitudeSourceLocation_masterevent.F90 \
                  AmplitudeSourceLocation_synthwave.F90 TraveltimeSourceLocation_masterevent.F90 \
                  AmplitudeSourceLocation_masterevent_shmdump.F90 AmplitudeSourceLocation_masterevent_ttimecor.F90 \
                  AmplitudeSourceLocation_DoubleDifference.F90 \
                  wavelet.f90 xorshift1024star.f90 nrtype.F90 linear_interpolation.F90 \
		  greatcircle.f90 grdfile_io.F90 set_velocity_model.F90 m_win.f90 \
		  m_winch.f90 calc_bpf_order.f90 calc_bpf_coef.f90 tandem.f90 calc_env_amplitude.F90 \
		  itoa.F90 grdfile_io.F90 m_util.f90 constants.F90 rayshooting.F90 read_sacfile.F90 read_shmdump.F90 \
                  raybending.F90 def_gridpoint.F90 lsqr_kinds.f90 lsqrblas.f90 lsqr.f90
OBJS		= $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(SRCS)))
MODS		= $(patsubst %.F90,%.mod,$(patsubst %.f90,%.mod,$(SRCS)))

.PHONY: all clean
.SUFFIXES: .f90 .F90
%.o : %.mod


all: $(TARGET)

asl_pw: AmplitudeSourceLocation_PulseWidth.o nrtype.o constants.o rayshooting.o read_sacfile.o \
	set_velocity_model.o linear_interpolation.o itoa.o grdfile_io.o m_win.o m_util.o m_winch.o \
        calc_bpf_order.o calc_bpf_coef.o tandem.o greatcircle.o raybending.o def_gridpoint.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_masterevent: AmplitudeSourceLocation_masterevent.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
	linear_interpolation.o greatcircle.o grdfile_io.o m_win.o m_util.o m_winch.o calc_bpf_order.o calc_bpf_coef.o tandem.o \
	read_sacfile.o raybending.o def_gridpoint.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_dd: AmplitudeSourceLocation_DoubleDifference.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
	linear_interpolation.o greatcircle.o grdfile_io.o raybending.o def_gridpoint.o lsqr_kinds.o lsqrblas.o lsqr.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_synthwave: AmplitudeSourceLocation_synthwave.o nrtype.o constants.o rayshooting.o read_sacfile.o set_velocity_model.o \
		linear_interpolation.o greatcircle.o grdfile_io.o xorshift1024star.o wavelet.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

ttime_masterevent: TraveltimeSourceLocation_masterevent.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
		linear_interpolation.o greatcircle.o grdfile_io.o raybending.o def_gridpoint.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

calc_env_amplitude: calc_env_amplitude.o nrtype.o constants.o calc_bpf_order.o calc_bpf_coef.o tandem.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_masterevent_shmdump: AmplitudeSourceLocation_masterevent_shmdump.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
	read_shmdump.o greatcircle.o linear_interpolation.o grdfile_io.o m_win.o m_util.o m_winch.o \
        calc_bpf_order.o calc_bpf_coef.o tandem.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_masterevent_ttimecor: AmplitudeSourceLocation_masterevent_ttimecor.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
	linear_interpolation.o greatcircle.o grdfile_io.o m_win.o m_util.o m_winch.o calc_bpf_order.o calc_bpf_coef.o tandem.o \
	read_sacfile.o raybending.o def_gridpoint.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

calc_env_amplitude: calc_env_amplitude.o nrtype.o constants.o calc_bpf_order.o calc_bpf_coef.o tandem.o

AmplitudeSourceLocation_PulseWidth.o: nrtype.o constants.o rayshooting.o read_sacfile.o greatcircle.o \
					set_velocity_model.o linear_interpolation.o itoa.o grdfile_io.o \
					calc_bpf_coef.o calc_bpf_order.o tandem.o m_win.o m_winch.o m_util.o raybending.o \
					def_gridpoint.o

AmplitudeSourceLocation_masterevent.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o m_win.o m_winch.o m_util.o \
					calc_bpf_order.o calc_bpf_coef.o tandem.o read_sacfile.o raybending.o def_gridpoint.o

AmplitudeSourceLocation_DoubleDifference.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o raybending.o def_gridpoint.o \
                                        lsqr_kinds.o lsqrblas.o lsqr.o

AmplitudeSourceLocation_synthwave.o: nrtype.o constants.o rayshooting.o read_sacfile.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o xorshift1024star.o wavelet.o

TraveltimeSourceLocation_masterevent.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o raybending.o def_gridpoint.o

AmplitudeSourceLocation_masterevent_shmdump.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o m_winch.o m_util.o \
					calc_bpf_order.o calc_bpf_coef.o tandem.o read_shmdump.o

AmplitudeSourceLocation_masterevent_ttimecor.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o m_win.o m_winch.o m_util.o \
					calc_bpf_order.o calc_bpf_coef.o tandem.o read_sacfile.o raybending.o def_gridpoint.o

rayshooting.o: nrtype.o constants.o greatcircle.o
raybending.o: nrtype.o constants.o def_gridpoint.o greatcircle.o 
m_win.o: m_util.o
m_winch.o: m_win.o
calc_env_amplitude.o: nrtype.o constants.o calc_bpf_coef.o calc_bpf_order.o tandem.o
lsqr_kinds.o: nrtype.o
lsqrblas.o: nrtype.o lsqr_kinds.o
lsqr.o: nrtype.o lsqr_kinds.o lsqrblas.o



.f90.o:
	$(FC) -c $< $(DEFS) $(INCDIR) $(OPTS) $(FFLAGS)
.F90.o:
	$(FC) -c $< $(DEFS) $(INCDIR) $(OPTS) $(FFLAGS)



clean:
	rm -f core* $(TARGET) $(OBJS) $(MODS)

