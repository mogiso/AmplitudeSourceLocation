#!/usr/bin/env perl
use File::Basename;
#
#plot the results of (Amplitude|Traveltime)SourceLocation_masterevent.F90 using GMT5
#Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
#Copyright: (c) Masashi Ogiso 2020
#License  : MIT License https://opensource.org/licenses/MIT

$out = $ARGV[0];
$result = $ARGV[1];
$stationparam = $ARGV[2];
$dem_lonlat = $ARGV[3];
$dem_londep = $ARGV[4];
$dem_deplat = $ARGV[5];

$argc = $#ARGV;
if($argc != 5){
  print stderr "usage: perl plot_result_masterevent.pl (out_ps) (result_txt) (stationparam_txt) dem_grd(lon-lat) dem_txt(lon-dep) dem_txt(dep-lat)\n";
  die;
}



#@stlon = (143.9775, 143.9867, 144.0017, 144.0042, 144.0160);
#@stlat = (43.3797, 43.3955, 43.3818, 43.3903, 43.3695);
#@stname = ("V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM", "V.MNDK");
open IN, "<", $stationparam;
while(<IN>){
  chomp;
  @tmp = split /\s+/, $_;
  if($tmp[4] eq ".true."){
    push @stlon, $tmp[0];
    push @stlat, $tmp[1];
    push @stname, $tmp[3];
  }
}
close IN;

##read result file
open IN, "<", $result;
<IN>;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @amp_ratio, $tmp[0];
  push @sigma_amp_ratio, $tmp[1];
  push @evlon, $tmp[2];
  push @sigma_evlon, $tmp[3];
  push @evlat, $tmp[4];
  push @sigma_evlat, $tmp[5];
  push @evdep, $tmp[6];
  push @sigma_evdep, $tmp[7];
}
close IN;

##for GMT
$lon_w = 135.2;
$lon_e = 137.7;
$lat_s = 32.5;;
$lat_n = 34.0;
$dep_min = 0.0;
$dep_max = 20.0;

$mapsize_x = 10.0;
$mapsize_y = `echo $lon_e $lat_n | gmt mapproject -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
@tmp = split /\s+/, $mapsize_y;
$mapsize_y = $tmp[1];
$mapsize_z = 4.5;
$dx = 11.5;
$dy = 5.7;

$scale_x = 9.0;
$scale_y = $mapsize_y + 0.9;
$scale_dist = "20k";

$cpt_x = $mapsize_x / 2;
$cpt_y = -1.0;
$cpt_len_x = $mapsize_x;
$cpt_len_y = "0.3h";

$symbolsize = 0.5;
$symbolwidth = "0.5p";
$symbolsize_st = 0.4;
$title_x = 0.0;
$title_y = $mapsize_y + 0.4;

$annot_a_map = "0.5";
$annot_f_map = "0.25";
$annot_a_dep = "5";
$annot_f_dep = "5";

$color = "black";
$color_ref = "red";


system "gmt set PS_LINE_JOIN round";
system "gmt set FORMAT_GEO_MAP +D";
system "gmt set FONT_LABEL 17p,Helvetica";
system "gmt set MAP_FRAME_PEN thin";
system "gmt set FONT_ANNOT_PRIMARY 16p,Helvetica";


print stderr "output ps = $out\n";
open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -K -P -X3c -Y12c > $out";
close OUT;

##lon-lat
system "gmt psbasemap -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                      -Bxya${annot_a_map}f${annot_f_map} -BWeSn -Lx$scale_x/$scale_y+c$lat_s+w$scale_dist  \\
                      -O -K -P >> $out";
system "gmt pscoast -JM${mapsize_x} -R$lon_w/$lon_e/$lat_s/$lat_n -W0.8p,black -Df -O -K -P >> $out";
#system "gmt grdcontour $dem_lonlat -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -C50 -W0.6p,dimgray -O -K -P >> $out";

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Si$symbolsize_st -W0.9p,black -Gwhite -O -K -P >> $out";
for($j = 0; $j <= $#stlon; $j++){
  print OUT "$stlon[$j] $stlat[$j]\n";
}
close OUT;

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
for($i = 0; $i <= $#amp_ratio; $i++){
  if($amp_ratio[$i] == 1.0){
    $ref_lon = $evlon[$i];
    $ref_lat = $evlat[$i];
    $ref_dep = $evdep[$i];
  }else{
    print OUT "$evlon[$i] $evlat[$i]\n";
  }
}
close OUT;

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},$color_ref -O -K -P >> $out";
  print OUT "$ref_lon $ref_lat\n";
close OUT;


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                      -Bx${annot_f_map} -Bya${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BWsen -O -K -P >> $out";

open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
  for($i = 0; $i <= $#amp_ratio; $i++){
    if($amp_ratio[$i] != 1.0){
      print OUT "$evlon[$i] $evdep[$i]\n";
    }
  }
close OUT;

open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -Sa$symbolsize -W${symbolwidth},$color_ref -O -K -P >> $out";
  print OUT "$ref_lon $ref_dep\n";
close OUT;

system "gmt psxy $dem_londep -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max -W0.6p,black -O -K -P >> $out";


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -X$dx -Y$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                      -By${annot_f_map} -Bxa${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BwSen -O -K -P >> $out";

open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
  for($i = 0; $i <= $#amp_ratio; $i++){
    if($amp_ratio[$i] != 1.0){
      print OUT "$evdep[$i] $evlat[$i]\n";
    }
  }
close OUT;
open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -Sa$symbolsize -W${symbolwidth},$color_ref -O -K -P >> $out";
  print OUT "$ref_dep $ref_lat\n";
close OUT;

system "gmt psxy $dem_deplat -JX$mapsize_z/$mapsize_y -R/$dep_min/$dep_max/$lat_s/$lat_n -W0.6p,black -O -K -P >> $out";


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -P >> $out";
close OUT;

system "gmt psconvert $out -A -Tg";

unlink "gmt.conf";
unlink "gmt.history";

