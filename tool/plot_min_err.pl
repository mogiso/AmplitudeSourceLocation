#!/usr/bin/env perl
use File::Basename;
#
#plot the results of AmplitudeSourceLocation_PulseWidth.F90 using GMT5
#Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
#Copyright: (c) Masashi Ogiso 2020
#License  : MIT License https://opensource.org/licenses/MIT

$stationparam = $ARGV[0];
$yr = $ARGV[1];
$mo = $ARGV[2];
$dy = $ARGV[3];
$hh = $ARGV[4];
$mm = $ARGV[5];
$ss = $ARGV[6];
$result = $ARGV[7];
$dem_lonlat = $ARGV[8];
$dem_londep = $ARGV[9];
$dem_deplat = $ARGV[10];

$plate = "PB2002_boundaries.dig.txt";

$argc = $#ARGV;
if($argc != 10){
  print stderr "usage: perl plot_min_err.pl stationparam yr mo day hh mm ss resultfile dem_grd(lon-lat) dem_txt(lon-dep) dem_txt(dep-lat)\n";
  die;
}

$resultdir = dirname $result;

$begin_sec = $hh * 3600 + $mm * 60 + $ss;

open IN, "<", $stationparam;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @stlon, $tmp[0];
  push @stlat, $tmp[1];
  push @stname, $tmp[3];
  push @use_flag, $tmp[4]
}
close IN;


#@stlon = (143.9775, 143.9867, 144.0017, 144.0042);
#@stlat = (43.3797, 43.3955, 43.3818, 43.3903);
#@stname = ("V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM");

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
  push @hypo_lon_w, $tmp[7];
  push @hypo_lon_e, $tmp[8];
  push @hypo_lat_s, $tmp[9];
  push @hypo_lat_n, $tmp[10];
  push @hypo_dep_min, $tmp[11];
  push @hypo_dep_max, $tmp[12];
}
close IN;

##for GMT
$lon_w = 135.2;
$lon_e = 137.5;
$lat_s = 32.5;
$lat_n = 34.0;
$dep_min = -1.0;
$dep_max = 21.0;

$mapsize_x = 7.0;
$mapsize_y = `echo $lon_e $lat_n | gmt mapproject -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
@tmp = split /\s+/, $mapsize_y;
$mapsize_y = $tmp[1];
$mapsize_z = 3.0;

$scale_x = $mapsize_x - 0.4;
$scale_y = $mapsize_y - 0.2;
$scale_dist = "40k";

$cpt_x = $mapsize_x / 2;
$cpt_y = -1.0;
$cpt_len_x = $mapsize_x;
$cpt_len_y = "0.3h";

$symbolsize = 0.4;
$symbolsize_st = 0.25;
$title_x = -0.2;
$title_y = $mapsize_y + 0.3;

$annot_a_map = "0.5";
$annot_f_map = "0.25";
$annot_a_dep = "5";
$annot_f_dep = "5";

system "gmt set PS_LINE_JOIN round";
system "gmt set FORMAT_GEO_MAP +D";
system "gmt set FONT_LABEL 12p,Helvetica";
system "gmt set MAP_FRAME_WIDTH 3p";
system "gmt set MAP_FRAME_PEN thin";
system "gmt set MAP_LABEL_OFFSET 6p";
system "gmt set MAP_ANNOT_OFFSET_PRIMARY 6p";
system "gmt set FONT_ANNOT_PRIMARY 11p";

$cpt = "error.cpt";
#for($i = 0; $i <= $#sec_from_begin; $i++){
for($i = 50; $i <= 54; $i++){
#for($i = 0; $i <= 60; $i++){
#for($i = 0; $i <= 0; $i++){
  $sec_from_day = $begin_sec + $sec_from_begin[$i];
  $hh_tmp = int($sec_from_day / 3600.0);
  $mm_tmp = int(($sec_from_day - $hh_tmp * 3600.0) / 60.0);
  $ss_tmp = $sec_from_day - $hh_tmp * 3600.0 - $mm_tmp * 60.0;

  $hh_tmp = sprintf "%02d", $hh_tmp;
  $mm_tmp = sprintf "%02d", $mm_tmp;
  $ss_tmp = sprintf "%.1f", $ss_tmp;
  if($ss_tmp < 10.0){
    $ss_tmp = "0$ss_tmp";
  }
  $date = "$yr-$mo-$dy $hh_tmp:$mm_tmp:$ss_tmp";

  $output_index = sprintf "%04d", $i;
  $min_err_lonlat = "${resultdir}/min_err_lon-lat_${output_index}.grd";
  $min_err_londep = "${resultdir}/min_err_lon-dep_${output_index}.grd";
  $min_err_deplat = "${resultdir}/min_err_dep-lat_${output_index}.grd";
  $out = "${resultdir}/min_err_${output_index}.ps";

  print stderr "output ps = $out\n";
  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -K -P -X3c -Y12c > $out";
  close OUT;

  ##lon-lat
  system "gmt grdimage $min_err_lonlat -C$cpt -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -O -K -P >> $out";
  system "gmt grdcontour $dem_lonlat -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -C200 -W0.4p,dimgray -O -K -P >> $out";
  system "gmt pscoast -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W0.8p,black -Gwhite -O -K -P >> $out";
  system "gmt psbasemap -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                        -Bxya${annot_a_map}f${annot_f_map} -BWeSn -Lx$scale_x/$scale_y+c$lat_s+w$scale_dist+jRT \\
                        -F+gwhite -O -K -P >> $out";
  open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                         -Si$symbolsize_st -W0.6p,black -Gwhite -O -K -P >> $out";
  for($j = 0; $j <= $#stlon; $j++){
    if($use_flag[$j] ne ".false."){
      print OUT "$stlon[$j] $stlat[$j]\n";
    }
  }
  close OUT;

  #open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
  #                       -W1.3p,white -O -K -P >> $out";
  #print OUT "$hypo_lon_w[$i] $min_lat[$i]\n";
  #print OUT "$hypo_lon_e[$i] $min_lat[$i]\n";
  #print OUT ">\n";
  #print OUT "$min_lon[$i] $hypo_lat_s[$i]\n";
  #print OUT "$min_lon[$i] $hypo_lat_n[$i]\n";
  #close OUT;
  open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                         -W0.3p,black -t100 -O -K -P >> $out";
  print OUT "$hypo_lon_w[$i] $min_lat[$i]\n";
  print OUT "$hypo_lon_e[$i] $min_lat[$i]\n";
  print OUT ">\n";
  print OUT "$min_lon[$i] $hypo_lat_s[$i]\n";
  print OUT "$min_lon[$i] $hypo_lat_n[$i]\n";
  close OUT;

  open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                         -Sa$symbolsize -W1.4p,white -O -K -P >> $out";
  print OUT "$min_lon[$i] $min_lat[$i]\n";
  close OUT;
  open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                         -Sa$symbolsize -W0.4p,black -t100 -O -K -P >> $out";
  print OUT "$min_lon[$i] $min_lat[$i]\n";
  close OUT;




  open OUT, " | gmt pstext -JX$mapsize_x/$mapsize_y -R0/$mapsize_x/0/$mapsize_y \\
                           -F+jBL+f12p,Helvetica,black -N -O -K -P >> $out";
  print OUT "$title_x $title_y $date\n";
  close OUT;


  system "gmt psxy $plate -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -W0.8p,black,- -O -K -P >> $out";
  system "gmt psscale -Dx$cpt_x/$cpt_y/$cpt_len_x/$cpt_len_y -C$cpt -Q -Bx+l\"Normalized residual\" -O -K -P >> $out";


  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-6c >> $out";
  close OUT;

  ##lon-dep
  system "gmt grdimage $min_err_londep -C$cpt -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max -O -K -P >> $out";
  system "gmt psbasemap -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                        -Bxf${annot_f_map} -Bya${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BWsen -O -K -P >> $out";




  #open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
  #                       -W1.3p,white -O -K -P >> $out";
  #print OUT "$hypo_lon_w[$i] $min_dep[$i]\n";
  #print OUT "$hypo_lon_e[$i] $min_dep[$i]\n";
  #print OUT ">\n";
  #print OUT "$min_lon[$i] $hypo_dep_min[$i]\n";
  #print OUT "$min_lon[$i] $hypo_dep_max[$i]\n";
  #close OUT;
  open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -W0.3p,black -t100 -O -K -P >> $out";
  print OUT "$hypo_lon_w[$i] $min_dep[$i]\n";
  print OUT "$hypo_lon_e[$i] $min_dep[$i]\n";
  print OUT ">\n";
  print OUT "$min_lon[$i] $hypo_dep_min[$i]\n";
  print OUT "$min_lon[$i] $hypo_dep_max[$i]\n";
  close OUT;

  open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -Sa$symbolsize -W1.4p,white  -O -K -P >> $out";
  print OUT "$min_lon[$i] $min_dep[$i]\n";
  close OUT;
  open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -Sa$symbolsize -W0.4p,black -t100 -O -K -P >> $out";
  print OUT "$min_lon[$i] $min_dep[$i]\n";
  close OUT;




  #system "gmt psxy $dem_londep -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max -W0.3p,black -O -K -P >> $out";

  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -X8c -Y6c >> $out";
  close OUT;

  ##lon-dep
  system "gmt grdimage $min_err_deplat -C$cpt -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n -O -K -P >> $out";
  system "gmt psbasemap -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                        -Byf${annot_f_map} -Bxa${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BwSen -O -K -P >> $out";




  #open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
  #                       -W1.3p,white -O -K -P >> $out";
  #print OUT "$hypo_dep_min[$i] $min_lat[$i]\n";
  #print OUT "$hypo_dep_max[$i] $min_lat[$i]\n";
  #print OUT ">\n";
  #print OUT "$min_dep[$i] $hypo_lat_s[$i]\n";
  #print OUT "$min_dep[$i] $hypo_lat_n[$i]\n";
  #close OUT;
  open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -W0.3p,black -t100 -O -K -P >> $out";
  print OUT "$hypo_dep_min[$i] $min_lat[$i]\n";
  print OUT "$hypo_dep_max[$i] $min_lat[$i]\n";
  print OUT ">\n";
  print OUT "$min_dep[$i] $hypo_lat_s[$i]\n";
  print OUT "$min_dep[$i] $hypo_lat_n[$i]\n";
  close OUT;

  open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -Sa$symbolsize -W1.4p,white -O -K -P >> $out";
  print OUT "$min_dep[$i] $min_lat[$i]\n";
  close OUT;
  open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -Sa$symbolsize -W0.4p,black -t100 -O -K -P >> $out";
  print OUT "$min_dep[$i] $min_lat[$i]\n";
  close OUT;





  #system "gmt psxy $dem_deplat -JX$mapsize_z/$mapsize_y -R/$dep_min/$dep_max/$lat_s/$lat_n -W0.3p,black -O -K -P >> $out";


  open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -P >> $out";
  close OUT;

  system "gmt psconvert $out -A -Te";

}

unlink "gmt.conf";
unlink "gmt.history";

