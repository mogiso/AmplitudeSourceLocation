# AmplitudeSourceLocation

## Description
A fortran program package that contains the Amplitude Source Location method using depth-dependent
1-D velocity structure and 3-D attenuation structure (AmplitudeSourceLocation_PulseWidth.F90), relative source location method
using seismic amplitudes (AmplitudeSourceLocation_masterevent.F90), and relative source location method using arrival times of
seismic waves (master event method, TraveltimeSourceLocation_masterevent.F90).

Language: Fortran 90

Require: [NetCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran/index.html "NetCDF-Fortran"),
and [LAPACK95](http://www.netlib.org/lapack95/ "LAPACK95")

These packages are used in AmplitudeSourceLocation_masterevent.F90 and TraveltimeSourceLocation_masterevent.F90.

NetCDF-Fortran depends on [NetCDF](https://www.unidata.ucar.edu/software/netcdf/),
and LAPACK95 depends on [LAPACK](http://www.netlib.org/lapack/). [Intel Math Kernel Library](https://software.intel.com/mkl) can
be used instead of LAPACK95.

## AmplitudeSourceLocation_PulseWidth.F90
### Description
This program estimates source locations of seismic events from seimic amplitudes with 1-D velocity and 3-D seismic attenuation structures
(Kumagai et al., 2019). This program calculates seismic amplitude at each station from SAC- or WIN- or WIN32-formatted waveform file(s) or reads pre-calculated
seismic amplitudes from text file. This program makes travel time and pulse width (i.e., attenuation) tables using ray shooting method. Then, this program conducts
a grid search in a given region to find suitable seismic source locations sequentially. A parallelization using OpenMP is available. I recommend to use OpenMP because
making travel time and pulse width tables is time-consuming work.

Reference:

### Compile
    $ make asl_pw

### Compile options
- `-DDOUBLE`: double precision froating point numbers are used in calculations (This option should be set)
- `-DSAC`: use SAC binary files (NVHDR should be 6) as input waveform files (Default)
- `-DWIN`: use [WIN-formatted](http://wwweic.eri.u-tokyo.ac.jp/WIN/man.en/winformat.html) binary file or
             [WIN32-formatted](https://www.hinet.bosai.go.jp/faq/?LANG=en#Q09) binary file as an input waveform files
- `-DAMP_TXT`: use text-formatted file as observed amplitude data, and do not apply bandpass filter
- `-DTESTDATA`: do not apply band pass filter (When -DSAC or -DWIN is set, this program applys band pass filter to each waveform)
- `-DOUT_AMPLITUDE`: output seismic amplitude at each station into a text-formatted file. This file can be used in AmplitudeSourceLocation_masterevent.F90.
- `-DWITHOUT_TTIME`: do not consider travel time from assumed seismic location to each station when calculating seismic amplitudes

### Executable commands
When `-DSAC` is set,

    $ ./asl_pw (name of digital elevetion map file) (name of station parameter file) (prefix of sac files) (component) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (length_of_rms_time_window) (dirname of results) (filename of result)

or when `-DWIN` is set,

    $ ./asl_pw (name of digital elevetion map file) (name of station parameter file) (WIN- or WIN32-formatted binary file) (filename of channel table) (component) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (length_of_rms_time_window) (dirname of results) (filename of result)

or when `-DAMP_TXT` is set,

    $ ./asl_pw (name of digital elevetion map file) (name of station parameter file) (filename of amplitudes) (frequency) (dirname of results) (filename of result)

### Explanation of arguments
- A digital elevation map file is netcdf-formatted 2-D gridded file (i.e., longitude (degree), latitude (degree), and altitude (meters)).
- The format of station parameter file is as follows:
 
<pre>
    143.9775 43.3797  -0.68 V.MEAB .true. 0.0 1.0    #(longitude, latitude, depth, station name, use_flag, traveltime correction term, and site amplification term of 1st station)
    143.9867 43.3955  -0.74 V.MEAA .true. 0.0 0.738  #(longitude, latitude, depth, station name, use_flag, traveltime correction term, and site amplification term of 2nd station)
    144.0017 43.3818  -1.27 V.PMNS .true. 0.0 2.213 
    144.0042 43.3903  -1.28 V.NSYM .true. 0.0 1.487
    144.0160 43.3695  -1.10 V.MNDK .true. 0.0 2.761
</pre>

The number of stations is automatically determined from station parameter file. Use_flag is either ".true." or ".false.". If use_flag equals .false., the station will
not use in the estimation of source locations.
- When `-DSAC` is set, filename of each SAC file is set to be `trim(prefix of sac files) // "." // trim(station name) // "." // trim(component) // ".sac"`.
Prefix of sac files and component are given in arguments while station name is given by the station parameter file.
- When `-DAMP_TXT` is set, the format of amplitude file is as follows:

<pre>
    # V.MEAB V.MEAA V.PMNS V.NSYM V.MNDK
    0.1621667E+00   0.5075075E-01   0.5697313E+00   0.2616340E+00   0.6868658E+00 305.0
    0.3181189E+00   0.9961827E-01   0.1068634E+01   0.4623199E+00   0.1426192E+01 320.0
    0.3297203E+00   0.1030206E+00   0.1135651E+01   0.4856689E+00   0.1504472E+01 335.0
    0.2320133E+00   0.7649654E-01   0.9502479E+00   0.4577495E+00   0.1055511E+01 350.0
    0.2997828E+00   0.8911001E-01   0.1118786E+01   0.5874434E+00   0.9242919E+00 365.0
    0.4260274E+00   0.1247084E+00   0.1472605E+01   0.7220987E+00   0.1189853E+01 380.0
    0.6150443E+00   0.1840327E+00   0.2174562E+01   0.9410463E+00   0.1964966E+01 395.0
    0.7942829E+00   0.2412873E+00   0.2978328E+01   0.1199161E+01   0.2700312E+01 410.0
    0.1046234E+01   0.2917264E+00   0.3683274E+01   0.1501370E+01   0.3492446E+01 425.0
    0.1226578E+01   0.3338429E+00   0.4143188E+01   0.1755798E+01   0.4092168E+01 440.0
    0.1369444E+01   0.3949369E+00   0.5131268E+01   0.1944450E+01   0.4765992E+01 455.0
</pre>

1st row is comment; it is required but the program read nothing. Seismic amplitudes of each subevent should be written in the following rows.
The order of amplitudes should be the same as the order of station parameter file. For the above example, the order of stations is
V.MEAB-V.MEAA-V.PMNS-V.NSYM-V.MNDK in the station parameter file so that amplitudes at V.MEAB are written in 1st column, those at V.MEAA in 2nd column,
those at V.PMNS in 3rd column, ...
- When `-DWIN` is set, the program will find channel ids from station names given by station parameter files and component given by the argument from
channel table file, then read waveforms of each channel id from the binary file.
- `fl` and `fh` (Hz) are the lower and upper passband edge frequency, and `fs` (Hz) is the stopband edge frequency, respectively. Frequency for the intrinsic
attenuation is the arithmetic mean of `fl` and `fh`. When `-DAMP_TXT` is set, the program expects that each amplitude has already been filtered so that only
frequency for the intrinsic attenuation is required in the arguments.
- `ot_begin and `ot_end` (s) determine the time range of estimating source locations. The assumed origin time begins at `ot_begin` and ends at `ot_end`,
shifting every `ot_shift` s. 

### Other settings
- Some parameters are hard-coded. Velocity and attenuation structures are given in set_velocity_model.F90. Please modify it for different structures.
- The variable `wavetype` in AmplitudeSourceLocation_PulseWidth.F90 indicates the type of body waves. 1 for P-waves and 2 for S-waves, respectively.
- The variables in "Search range" section define the region of grid search while those in "structure range" define the region for ray tracing. If ray
approached the edge of the region defined by them (or approarched to surface), ray shooting is terminated and restarted with another takeoff angle.

### Result files
- One is a text-formatted file with the filename given in arguments. The format is:
<pre>
    # OT min_lon min_lat min_dep source_amp residual
    305.0 144.0040 43.3750 .20   0.7376450E+00   0.1110586E-01
    320.0 144.0020 43.3720 .00   0.1485148E+01   0.6166780E-02
    335.0 144.0020 43.3720 .00   0.1608392E+01   0.6876275E-02
    350.0 144.0070 43.3770 .10   0.1187008E+01   0.1333692E-01
    365.0 144.0040 43.3780 .10   0.1182522E+01   0.2161646E-01
    380.0 144.0010 43.3770 .00   0.1532389E+01   0.2092869E-01
    395.0 144.0010 43.3760 -.10   0.2203310E+01   0.1279828E-01
</pre>
From left to right, assumed origin time, latitude, logitude and depth where the residual is minimum, source amplitude, and residual.
- The others are netcdf-formatted 2-D gridded files. In each assumed origin time, horizontal and vertical slices of residual distibution (i.e., three files) are presented.
Each slice contain the estimated source location.

## AmplitudeSourceLocation_masterevent.F90
### Description
This program estimates relative source locations of seismic events as well as relative amplitude of source radiation from seimic amplitudes. The formulation is given by
Ogiso and Yomogida (2021, submitted to Earth, Planets and Space). Although this
program links set_velocity_model.F90, i.e., 1-D velocity and 3-D attenuation strucutures are used, the formulation requires only velocity and attenuation factors at
the location of reference event. 

### Compile
    $ make asl_masterevent

### Compile options
Some options are same as those of AmplitudeSourceLocation_PulseWidth.F90. The default format of amplitude data is text (see explanation of `-DAMP_TXT`).
- `-DDOUBLE`: double precision froating point numbers are used in calculations (This option should be set)
- `-DDAMPED`: Add additional constraints for model parameters (relative source amplitude and locations)
- `-DSAC`: use SAC binary files (NVHDR should be 6) as input waveform files
- `-DWIN`: use [WIN-formatted](http://wwweic.eri.u-tokyo.ac.jp/WIN/man.en/winformat.html) binary file or
             [WIN32-formatted](https://www.hinet.bosai.go.jp/faq/?LANG=en#Q09) binary file as an input waveform files
- `-DTESTDATA`: do not apply band pass filter (When `-DSAC` or `-DWIN` is set, this program applys band pass filter to each waveform)
- `-DEACH_ERROR`: Estimation errors of model parameter is calculated for each subevent (i.e., variance of data is calculated for each subevent). In default, variance of data
is calculated using all subevents so that estimation errors of locations would be the same for all subevents.
- `-DWITHOUT_TTIME`: do not consider travel time from assumed seismic location to each station when calculating seismic amplitudes
- `-DMKL`: use Intel Math Kernel Library instead of original lapack95/lapack subroutine libraries.

### Executable commands
When `-DSAC` is set,

    $ ./asl_masterevent (name of digital elevetion map file) (name of station parameter file) (filename of amplitude data of a reference event) (prefix of sac files)
    (component) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (length_of_rms_time_window) (filename of result)

or when `-DWIN` is set,

    $ ./asl_masterevent (name of digital elevetion map file) (name of station parameter file) (filename of amplitude data of a reference event)
    (WIN- or WIN32-formatted binary file) (filename of channel table) (component) (fl) (fh) (fs) (ot_begin) (ot_end) (ot_shift) (length_of_rms_time_window) (filename of result)

or default is

    $ ./asl_masterevent (name of digital elevetion map file) (name of station parameter file) (filename of amplitude data of a reference event)
    (filename of amplitude data of subevents) (filename of result)

### Explanation of arguments
- The most of the arguments are the same as those for AmplitudeSourceLocation_PulseWidth.F90.
- The format of the amplitude data of the reference event is as follows:

<pre>
    # V.MEAB V.MEAA V.PMNS V.NSYM V.MNDK
    144.0130 43.3810 0.30
    0.1621667E+00   0.5075075E-01   0.5697313E+00   0.2616340E+00   0.6868658E+00
</pre>

1st row is comment; it is required but the program read nothing. The location (longitude, latitude and depth) of the reference event is given in 2nd row.
Amplitudes at each station is given in 3rd row. The order of stations should be the same as that of the station parameter file.

### Other settings
- The same as those in case of AmplitudeSourceLocation_PulseWidth.F90.
- WIN- (WIN32-) or sac-formatted file(s) is used for continuous waves. To apply this program to event files, please make text-formatted file of amplitude data
from each event files, then use it in this program.

### Result files
- Text-formatted file which filename is given in arguments. The format is:
<pre>
   # amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth residual_sum ot_tmp
 0.14033281E+01  0.94988695E-01  0.14401078E+03  0.64774758E-03  0.43377751E+02  0.74969229E-03  0.14679086E+00  0.11586054E+00  0.43900310E-02     0.0
 0.14646312E+01  0.13915351E+00  0.14401010E+03  0.90919905E-03  0.43378310E+02  0.10522919E-02  0.30514866E-01  0.16262552E+00  0.86491601E-02    15.0
 0.12826188E+01  0.13706276E+00  0.14401027E+03  0.10226216E-02  0.43378508E+02  0.11835652E-02 -0.11308208E-01  0.18291305E+00  0.10941728E-01    30.0
 0.12159722E+01  0.61448032E-01  0.14400896E+03  0.48359011E-03  0.43377395E+02  0.55969917E-03 -0.14821482E+00  0.86498216E-01  0.24468696E-02    45.0
 0.13641865E+01  0.69949849E-01  0.14400868E+03  0.49068877E-03  0.43377529E+02  0.56791503E-03 -0.19671440E+00  0.87767929E-01  0.25192323E-02    60.0
</pre>
From left to right, relative source amplitude, latitude, logitude and depth with its errors. The 10th column is a comment; if `-DSAC` or `-DWIN` is set, assumed origin
times are given in this column. Otherwise, this column will be blank.


## TraveltimeSourceLocation_masterevent.F90
### Description
This program estimates relative source locations of seismic events from travel times of seismic phases (e.g., Aoki, 1974). This program can treat travel times of either
P-waves or S-waves.

### Compile
    $ make ttime_masterevent

### Compile options
Some options are same as those of AmplitudeSourceLocation_PulseWidth.F90.
- `-DDOUBLE`: double precision froating point numbers are used in calculations (This option should be set)
- `-DMKL`: use Intel Math Kernel Library instead of original lapack95/lapack subroutine libraries.
- `-DEACH_ERROR`: Estimation errors of model parameter is calculated for each subevent (i.e., variance of data is calculated for each subevent). In default, variance of data
is calculated using all subevents so that estimation errors of locations would be the same for all subevents.

This program accepts only text-formatted files as travel time data of the reference and the other subevents.

### Executable commands

    $ ./ttime_masterevent (name of digital elevetion map file) (name of station parameter file) (filename of travel time data of a reference event)
    (filename of travel time data of subevents) (filename of result)

### Other settings
- The same as those in case of AmplitudeSourceLocation_masterevent.F90.
- In default, the variable "wavetype" is set to be 1, i.e., using velocity of P-waves.

### Result files
- Text-formatted file which filename is given in arguments. The format is:
<pre>
# OTdiff sigma_OTdiff longitude sigma_lon latitude sigma_lat depth sigma_depth
 0.8174785E-02  0.4855984E-01  0.1440098E+03  0.1153894E-02  0.4338137E+02  0.1177685E-02 -0.3155910E+00  0.1955303E+00
 0.7186357E-03  0.2852916E-01  0.1440076E+03  0.6779187E-03  0.4338348E+02  0.6918960E-03 -0.4451443E+00  0.1148751E+00
-0.1033108E-01  0.1572529E-01  0.1440109E+03  0.3736693E-03  0.4338275E+02  0.3813735E-03 -0.2985357E+00  0.6331922E-01
-0.5672537E-02  0.2319235E-01  0.1440123E+03  0.5511038E-03  0.4338324E+02  0.5624663E-03 -0.2943046E+00  0.9338595E-01
</pre>
From left to right, difference of origin time, latitude, logitude and depth with its errors. The 10th column is a comment; if -DSAC or -DWIN is set, assumed origin
times are given in this column. Otherwise, this column will be blank.

## License
MIT License except m_util.f90, m_win.f90, m_winch.f90 written by Takuto Maeda, calc_bpf_order.f90, calc_bpf_coef.f90, tandem1.f90 taken from Saito (1978), and xorshift1024star.f90 written by Shun Sakuraba.

## References
Aoki (1974), Plate Tectonics of Arc-junction at Central Japan, Journal of Physics of the Earth, 22(1), 141-161,
[doi: 10.4294/jpe1952.22.141](https://doi.org/10.4294/jpe1952.22.141).

Kumagai et al. (2019), Amplitude Source Location Method With Depth-Dependent Scattering and Attenuation Structures: Application at Nevado del Ruiz Volcano, Colombia,
J. Geophys. Res., 124, [doi: 10.1029/2019JB018156](https://doi.org/10.1029/2019JB018156).

Saito (1978), An automatic design algorithm for band selective recursive digital filters, Geophysical Prospecting (Butsuri Tanko), 31(4), 240-263. (In Japanese)

## Acknowledgments
I utilize a part of fwin source code written by Takuto Maeda (https://github.com/tktmyd/fwin).

