# AmplitudeSourceLocation

## Description
Amplitude Source Location method considering depth-dependent 1D velocity structure and 3D attenuation structure.
Waveform format: SAC binary file or WIN format file

Language: Fortran 90

Require: netcdf-fortran

## Search parameters
Almost all the parameters are hard coded in AmplitudeSourceLocation_PulseWidth.F90 except velocity and attenuation structure in set_velocity_model.F90.

## Compile
$ make asl_pw

To use WIN-format file as an input waveform, set -DWIN at $DEFS on Makefile.

## Usage
If the sac-formatted waveform files are used,

$ ./asl_pw sacfile_index dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir result_file_name(txt)

or

$ ./asl_pw winfile win_chfile dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir result_file_name(txt)

for the case of win-formatted file are used.

In the case of sac files, filenames are generated as trim(sacfile_index) // trim(station_name) // trim(sacfile_extension).
To use win file as the input waveform file, give channel ids to st_winch in AmplitudeSourceLocation_PulseWidth.F90.

Many netcdf-format grdfiles are generated in resultdir. These files give three slices (horizontal, longitudinal and latitudinal) of residual distribution in each assumed origin time. The perl script "plot_min_err.pl" plots them.

## License
MIT License except m_util.f90, m_win.f90, m_winch.f90 written by Takuto Maeda.

### References
Kumagai et al., Amplitude Source Location Method With Depth‐Dependent Scattering and Attenuation Structures: Application at Nevado del Ruiz Volcano, Colombia, JGR, 2019, doi: 10.1029/2019JB018156.

### Acknowledgments
I utilize a part of fwin source code written by Takuto Maeda (https://github.com/tktmyd/fwin). Also, this project takes advantage of netCDF software developed by UCAR/Unidata (http://doi.org/10.5065/D6H70CW6). I appreciate their efforts.
