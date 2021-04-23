#!/usr/bin/env perl
use File::Basename;
#
#plot the results of AmplitudeSourceLocation_Pulsewidth.F90 using GMT5
#plot all the locations in the same figure
#Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
#Copyright: (c) Masashi Ogiso 2020
#License  : MIT License https://opensource.org/licenses/MIT

$out = $ARGV[0];
$result = $ARGV[1];
$dem_lonlat = $ARGV[2];
$dem_londep = $ARGV[3];
$dem_deplat = $ARGV[4];
$stationparam = $ARGV[5];

$argc = $#ARGV;
if($argc != 5){
  print stderr "usage: perl plot_result_masterevent.pl (out_ps) (result_txt) dem_grd(lon-lat) dem_txt(lon-dep) dem_txt(dep-lat) stationparam\n";
  die;
}



#@stlon = (143.9775, 143.9867, 144.0017, 144.0042, 144.0160);
#@stlat = (43.3797, 43.3955, 43.3818, 43.3903, 43.3695);
#@stname = ("V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM", "V.MNDK");
#
open IN, "<", $stationparam;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @stlon, $tmp[0];
  push @stlat, $tmp[1];
  push @stname, $tmp[3];
}
close IN;

##read result file
open IN, "<", $result;
<IN>;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  next if (substr($_, 0, 1)) == "#";
  @tmp = split /\s+/, $_;
  push @origintime, $tmp[0];
  push @evlon, $tmp[1];
  push @evlat, $tmp[2];
  push @evdep, $tmp[3];
}
close IN;

##for GMT
$lon_w = 135.45;
$lon_e = 137.7;
$lat_s = 32.45;
$lat_n = 34.0;
$dep_min = -1;
$dep_max = 21;
$ot_min = "2020-12-1T00:00:00";
$ot_max = "2021-2-1T00:00:00";

$size_x = 11.0;
$size_y = 7.0;
$dx = $size_x + 1.0;
$dy = $size_y + 1.0;;

$size_hist = 5.0;
$freq_min = 0;
$freq_max = 2500;
$bin_deg = 0.1;
$bin_depth = 2;


$symbolsize = 0.3;
$symbolwidth = "0.4p";
$symbolsize_st = 0.4;
$title_x = 0.0;
$title_y = $mapsize_y + 0.4;

$annot_a_map = "0.5";
$annot_f_map = "0.1";
$annot_a_dep = "5";
$annot_f_dep = "1";
$annot_a_time = "10d";
$annot_a_time2 = "1O";
$annot_f_time = "1d";
$annot_a_hist = 500;
$annot_f_hist = 500;

$color = "black";

system "gmt set PS_LINE_JOIN round";
system "gmt set FORMAT_GEO_MAP +D";
system "gmt set FONT_LABEL 17p,Helvetica";
system "gmt set MAP_FRAME_PEN thin";
system "gmt set FONT_ANNOT_PRIMARY 11p,Helvetica";
system "gmt set FONT_ANNOT_SECONDARY 16p,Helvetica";
system "gmt set FORMAT_DATE_MAP yyyy-mm";


print stderr "output ps = $out\n";
open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -K -P -X3c -Y20c > $out";
close OUT;

##time-lon

open OUT, " | gmt psxy -JX${size_x}T/$size_y -R$ot_min/$ot_max/$lon_w/$lon_e \\
                       -Sa$symbolsize -W${symbolwidth},$color \\
                       -BWesn -Bpx${annot_a_time}f${annot_f_time} -By${annot_a_map}f${annot_f_map}+l\"Longitude\" \\
                       -Bsx${annot_a_time2} \\
                       -O -K -P >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$origintime[$j] $evlon[$j]\n";
}
close OUT;

#histogram
open OUT, " | gmt pshistogram -JX$size_hist/$size_y -A -R$lon_w/$lon_e/$freq_min/$freq_max -W$bin_deg \\
                              -L0.5p,black -Gdimgray \\
                              -Bya${annot_a_hist}f${annot_f_hist} -Bxa${annot_a_map}f${annot_f_map} -BwSen \\
                              -O -K -P -X$dx >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evlon[$j]\n";
}
close OUT;

open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-$dy -X-$dx >> $out";
close OUT;

##time-lat
open OUT, " | gmt psxy -JX${size_x}T/$size_y -R$ot_min/$ot_max/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},$color \\
                       -BWesn -Bpx${annot_a_time}f${annot_f_time} -By${annot_a_map}f${annot_f_map}+l\"Latitude\" \\
                       -Bsx${annot_a_time2} \\
                       -O -K -P >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$origintime[$j] $evlat[$j]\n";
}
close OUT;

#histogram
open OUT, " | gmt pshistogram -JX$size_hist/$size_y -A -R$lat_s/$lat_n/$freq_min/$freq_max -W$bin_deg \\
                              -L0.5p,black -Gdimgray \\
                              -Bya${annot_a_hist}f${annot_f_hist} -Bxa${annot_a_map}f${annot_f_map} -BwSen \\
                              -O -K -P -X$dx >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evlat[$j]\n";
}
close OUT;


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-$dy -X-$dx >> $out";
close OUT;

##dep-time

open OUT, " | gmt psxy -JX${size_x}T/-$size_y -R$ot_min/$ot_max/$dep_min/$dep_max \\
                       -Sa$symbolsize -W${symbolwidth},$color \\
                       -BWeSn -Bpx${annot_a_time}f${annot_f_time} -By${annot_a_dep}f${annot_f_dep}+l\"Depth\" \\
                       -Bsx${annot_a_time2} \\
                       -O -K -P >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$origintime[$j] $evdep[$j]\n";
}
close OUT;

#histogram
open OUT, " | gmt pshistogram -JX$size_hist/-$size_y -A -R$dep_min/$dep_max/$freq_min/$freq_max -W$bin_depth \\
                              -L0.5p,black -Gdimgray \\
                              -Bya${annot_a_hist}f${annot_f_hist} -Bxa${annot_a_dep}f${annot_f_dep} -BwSen \\
                              -O -K -P -X$dx >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evdep[$j]\n";
}
close OUT;


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -P >> $out";
close OUT;

system "gmt psconvert $out -A -Tg";

unlink "gmt.conf";
unlink "gmt.history";

