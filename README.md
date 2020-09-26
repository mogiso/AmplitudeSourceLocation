# AmplitudeSourceLocation

## Description
Amplitude Source Location method considering depth-dependent 1D velocity structure and 3D attenuation structure.
Waveform format: SAC binary file or WIN format file

Language: Fortran 90

Require: netcdf-fortran

## Search parameters
Almost all the parameters are hard coded in AmplitudeSourceLocation_PulseWidth.F90 except velocity and attenuation
structure in set_velocity_model.F90.

All the stations used in this analysis should be within the search range, because the ray tracing will be stopped when the ray
reaches to the border of the search range.

## Compile
    $ make asl_pw

To use WIN-format file as an input waveform, set -DWIN at $DEFS in Makefile.

## Usage
If sac-formatted waveform files are used,

    $ ./asl_pw sacfile_index dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir result_file_name(txt)

or

    $ ./asl_pw winfile win_chfile dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir result_file_name(txt)

when win-formatted file are used. "ot_begin" and "ot_end" are the time in second measured from the beginning of the waveforms.
These two variables define the search range in time dimension.

In case of sac files, filenames are generated as trim(sacfile_index) // trim(station_name) // trim(sacfile_extension).
To use win file as the input waveform file, give channel ids to st_winch in AmplitudeSourceLocation_PulseWidth.F90.

Many netcdf-format grdfiles are generated in resultdir. These files give three slices (horizontal, longitudinal and latitudinal)
of residual distribution in each assumed origin time. The perl script "plot_min_err.pl" plots them.

## Parallelization
I have parallelized where travel time and pulse width table are made by using OpenMP. If homogeneous velocity
and attenuation structures are adopted, calculation time seems to be reasonably small even if no parallelization,
however, I recommend to use OpenMP when 1D velocity structure is adopted. In this program, I have adopted grid search with ray
shooting as ray tracing method, so making travel time table is time-consuming work.

Grid search section to find minimum residual is also parallelized using OpenMP.

# AmplitudeSourceLocation_masterevent
## Description
Relative location estimation using seismic amplitude.

Input files format: txt files except topography file (netcdf format)

Language: Fortran 90

Require: Lapack (dgels, dgetrf, dgetri), netcdf-fortran

## Compile
    $ make asl_masterevent

## Usage
    $ ./asl_masterevent (dem_grdfile) (station_param_file) (masterevent_param_file) (subevent_param_file) (result_file)

### station_param_file
example:

    5                              #(number of stations)
    143.9775 43.3797  -0.68 V.MEAB #(longitude, latitude, and depth of 1st station)
    143.9867 43.3955  -0.74 V.MEAA #(longitude, latitude, and depth of 2nd station)
    144.0017 43.3818  -1.27 V.PMNS
    144.0042 43.3903  -1.28 V.NSYM
    144.0160 43.3695  -1.10 V.MNDK

The program reads 1st column of 1st row, and initial three columns of 2nd, 3rd, ... rows. Station names shown in 4th column of the example are not required.

### masterevent_param_file
example:

    # V.MEAB V.MEAA V.PMNS V.NSYM V.MNDK 
    144.0040 43.3750 0.20
    0.1621667E+00   0.5075075E-01   0.5697313E+00   0.2616340E+00   0.6868658E+00 305.0

1st row is a comment row. The row is required but the program read nothing from it. 2nd row is the location (longitude, latitude, and depth) of the reference event. 3rd row is the observed amplitudes of the reference event at stations. The order of the colomn must be the same as the order in the station_param_file. In this example, I use five stations, and 6th column is a comment. The program does not read 6th (in this case) column.

### subevent_param_file
example:

    # V.MEAB V.MEAA V.PMNS V.NSYM V.MNDK
    73
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

1st row is a comment raw. Same as the masterevent_param_file, this row is required. 2nd row is the number of subevents. 3rd, 4th, 5th, ..., row is the observed amplitudes of each subevent. Same as the masterevent_param_file, the order of the amplitudes musbe the same as the order in the station_param_file. 6th (in this case) row is a comment.

If AMP_TXT is defined when you compile AmplitudeSourceLocation_Pulsewidth.F90, the program outputs observed amplitude at each station as a txt file so that you can make the masterevent_param_file and subevent_param_file by editing it. Note that each observed amplitude made by AmplitudeSourceLocation_Pulsewidth.F90 depends on the location of each event so that if the estimated location by AmplitudeSourceLocation_Pulsewidth.F90 is not appropriate, the text file made by the program would not be appropriate. Consider defining WITHOUT_TTIME when compile the program and taking a long time window for calculating RMS amplitude to neglect the effect of travel time. 

# TraveltimeSourceLocation_masterevent
## Description
Relative location estimation using either P- or S-wave arrival times.

Input files format: txt files except topography file (netcdf format)

Language: Fortran 90

Require: Lapack (dgels, dgetrf, dgetri), netcdf-fortran

## Compile
    $ make ttime_masterevent

## Usage
    $ ./ttime_masterevent (dem_grdfile) (station_param_file) (masterevent_param_file) (subevent_param_file) (result_file)

The formats of input files are the same as AmplitudeSourceLocation_masterevent.F90 except arrival times of P- or S-waves. Please check set_velocity_model.F90 to make observation and velocity model consistent with each other.

## License
MIT License except m_util.f90, m_win.f90, m_winch.f90 written by Takuto Maeda, calc_bpf_order.f90, calc_bpf_coef.f90, tandem1.f90 taken from Saito (1978).

### References
Kumagai et al., Amplitude Source Location Method With Depth-dependent Scattering and Attenuation Structures:
Application at Nevado del Ruiz Volcano, Colombia, JGR, 2019, doi: 10.1029/2019JB018156.

Saito, An automatic design algorithm for band selective recursive digital filters, Geophysical Prospecting (Butsuri Tanko), 31(4), 240-263, 1978. (In Japanese)

### Acknowledgments
I utilize a part of fwin source code written by Takuto Maeda (https://github.com/tktmyd/fwin).
Also, this project takes advantage of netCDF software developed by UCAR/Unidata (http://doi.org/10.5065/D6H70CW6).

I appreciate their efforts.
