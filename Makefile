##Makefile for AmplitudeSourceLocation_Pulsewidth.F90
##Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
##Copyright: (c) Masashi Ogiso 2020
##License  : MIT License (https://opensource.org/licenses/MIT)

## -DDOUBLE: calculate in double precision (recommended)
## -DV_MEA1D, -DV_CONST: velocity structure defined in set_velocity_model.F90
## -DWIN: using win format waveform file as an input (otherwise sac binary format files are used) 
##        for AmplitudeSourcelocation_Pulsewidth.F90 
## -DSTDP_COR: convert the depth of station from meters to kilometers (origin: sea level, downward positive)
##             when sac files are uses as an input (if -DWIN is set, station location is given in altitude (m)
##             at channel files, so this variable have no meanings)
## -DTESTDATA: do not apply bandpass filter to input waveforms (to be used for synthetic data) 
##             in AmplitudeSourcelocation_PulseWidth.F90
## -DOUT_AMPLITUDE: output text file of observed amplitude (without site correction) for asl_masterevent
##                  in AmplitudeSourceLocation_PulseWidth.F90
## -DMKL: use MKL; otherwise use lapack95 in AmplitudeSourceLocation_masterevent.F90

#FC = ifort
#FFLAGS = -traceback -assume byterecl -qopenmp
#DEFS = -DDOUBLE -DV_MEA1D -DTESTDATA -DMKL
#INCDIR = -I${NETCDF_FORTRAN_INC} -I${MKLROOT}/include/intel64/lp64
#LIBDIR = -L${MKLROOT}/lib/intel64
#LIBS = -lnetcdff -liomp5 -lpthread -lmkl_core -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_intel_thread
#OPTS = -O3 -xHOST

FC = gfortran
FFLAGS = -g -Wall -fbounds-check -fbacktrace
#FFLAGS = -fbacktrace
DEFS = -DDOUBLE -DV_MEA1D -DWIN 
INCDIR = -I/usr/include -I/usr/local/include
LIBDIR = 
LIBS = -lnetcdff -llapack95 -llapack -lblas
#OPTS = -O3 -fopenmp

TARGET		= asl_pw asl_masterevent asl_synthwave
SRCS		= AmplitudeSourceLocation_PulseWidth.F90 AmplitudeSourceLocation_masterevent.F90 \
                  AmplitudeSourceLocation_synthwave.F90 \
                  wavelet.f90 xorshift1024star.f90 nrtype.F90 linear_interpolation.F90 \
		  greatcircle.f90 grdfile_io.F90 set_velocity_model.F90 m_win.f90 \
		  m_winch.f90 calc_bpf_order.f90 calc_bpf_coef.f90 tandem.f90 \
		  itoa.F90 grdfile_io.F90 m_util.f90 constants.F90 rayshooting.F90 read_sacfile.F90
OBJS		= $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(SRCS)))
MODS		= $(patsubst %.F90,%.mod,$(patsubst %.f90,%.mod,$(SRCS)))

.PHONY: all clean
.SUFFIXES: .f90 .F90


all: $(TARGET)

asl_pw: AmplitudeSourceLocation_PulseWidth.o nrtype.o constants.o rayshooting.o read_sacfile.o \
	set_velocity_model.o linear_interpolation.o itoa.o grdfile_io.o m_win.o m_util.o m_winch.o \
        calc_bpf_order.o calc_bpf_coef.o tandem.o greatcircle.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_masterevent: AmplitudeSourceLocation_masterevent.o nrtype.o constants.o rayshooting.o set_velocity_model.o \
		linear_interpolation.o greatcircle.o grdfile_io.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

asl_synthwave: AmplitudeSourceLocation_synthwave.o nrtype.o constants.o rayshooting.o read_sacfile.o set_velocity_model.o \
		linear_interpolation.o greatcircle.o grdfile_io.o xorshift1024star.o wavelet.o
	$(FC) -o $@ $^ $(OPTS) $(LIBDIR) $(LIBS)

AmplitudeSourceLocation_PulseWidth.o: nrtype.o constants.o rayshooting.o read_sacfile.o greatcircle.o \
					set_velocity_model.o linear_interpolation.o itoa.o grdfile_io.o \
					calc_bpf_coef.o calc_bpf_order.o tandem.o m_win.o m_winch.o m_util.o

AmplitudeSourceLocation_masterevent.o: nrtype.o constants.o rayshooting.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o

AmplitudeSourceLocation_synthwave.o: nrtype.o constants.o rayshooting.o read_sacfile.o set_velocity_model.o \
					linear_interpolation.o greatcircle.o grdfile_io.o xorshift1024star.o wavelet.o

rayshooting.o: greatcircle.o
m_win.o: m_util.o


.f90.o:
	$(FC) -c $< $(DEFS) $(INCDIR) $(OPTS) $(FFLAGS)
.F90.o:
	$(FC) -c $< $(DEFS) $(INCDIR) $(OPTS) $(FFLAGS)



clean:
	rm -f core* $(TARGET) $(OBJS) $(MODS)
