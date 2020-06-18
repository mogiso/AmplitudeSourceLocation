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

    $ ./asl_pw sacfile_index dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir
 result_file_name(txt)

or

    $ ./asl_pw winfile win_chfile dem_grdfile ot_begin(sec) ot_end(sec) ot_shift(sec) rms_time_window(sec) resultdir
 result_file_name(txt)

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

## License
MIT License except m_util.f90, m_win.f90, m_winch.f90 written by Takuto Maeda.

### References
Kumagai et al., Amplitude Source Location Method With Depth-dependent Scattering and Attenuation Structures:
Application at Nevado del Ruiz Volcano, Colombia, JGR, 2019, doi: 10.1029/2019JB018156.

### Acknowledgments
I utilize a part of fwin source code written by Takuto Maeda (https://github.com/tktmyd/fwin).
Also, this project takes advantage of netCDF software developed by UCAR/Unidata (http://doi.org/10.5065/D6H70CW6).

I appreciate their efforts.
