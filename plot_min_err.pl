#!/usr/bin/env perl
#
#plot ASL result

$yr = $ARGV[0];
$mo = $ARGV[1];
$dy = $ARGV[2];
$hh = $ARGV[3];
$mm = $ARGV[4];
$ss = $ARGV[5];
$result = $ARGV[6];
$demfile = $ARGV[7];

$begin_sec = $hh * 3600 + $mm * 60 + $ss;

##read result file
open IN, "<", $result;
<IN>;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @sec_from_begin, $tmp[0];
  push @min_lon, $tmp[1];
  push @min_lat, $tmp[2];
  push @min_dep, $tmp[3];
  push @source_amp, $tmp[4];
  push @residual, $tmp[5];
}
close IN;

##for GMT
$lon_w = 143.975;
$lon_e = 144.03;
$lat_s = 43.37;
$lat_n = 43.41;
$dep_min = -1.5;
$dep_max = 3.0;

$mapsize_x = 10.0;
$mapsize_y = `echo $lon_e $lat_n | gmt mapproject -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
@tmp = split /\s+/, $mapsize_y;
$mapsize_y = $tmp[1];
$mapsize_z = 4.0;

$cpt = "error.cpt";
#for($i = 0; $i <= $#sec_from_begin; $i++){
for($i = 0; $i <= 0; $i++){
  $sec_from_day = $begin_sec + $sec_from_begin;
  $hh_tmp = int($sec_from_day / 3600.0);
  $mm_tmp = int(($sec_from_day - $hh_tmp * 3600.0) / 60.0);
  $ss_tmp = $sec_from_day - $hh_tmp * 3600.0 - $mm_tmp * 60.0;

  $hh_tmp = sprintf "%02d", $hh_tmp;
  $mm_tmp = sprintf "%02d", $mm_tmp;
  $ss_tmp = sprintf "%02.1f", $ss_tmp;
  $date = "$yr/$mo/$dy $hh_tmp:$mm_tmp:$ss_tmp";

  $output_index = sprintf "%04d", $i;
  $min_err_lonlat = "min_err_lon-lat_${output_index}.grd";
  $min_err_londep = "min_err_lon-dep_${output_index}.grd";
  $min_err_deplat = "min_err_dep-lat_${output_index}.grd";
  $out = "min_err_${output_index}.ps";

  print stderr "output ps = $out\n";
  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -K -P -X3c -Y8c > $out";
  close OUT;

  system "gmt grdimage $min_err_lonlat -C$cpt -JM$mapsize -R$lon_w/$lon_e/$lat_s/$lat_n -O -K -P >> $out";
  system "gmt psbasemap -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -Ba0.02/a0.02::WSen -O -K -P >> $out";
  system "gmt pscontour $demfile -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -C50 -bi3f -W0.5p,black -O -K -P >> $out";


  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -P >> $out";
  close OUT;

}

