#!/usr/bin/env perl
use File::Basename;
#
#plot the results of (Amplitude|Traveltime)SourceLocation_masterevent.F90 using GMT5
#Author   : Masashi Ogiso (masashi.ogiso@gmail.com)
#Copyright: (c) Masashi Ogiso 2020
#License  : MIT License https://opensource.org/licenses/MIT

$out = $ARGV[0];
$result = $ARGV[1];
$dem_lonlat = $ARGV[2];
$dem_londep = $ARGV[3];
$dem_deplat = $ARGV[4];

$argc = $#ARGV;
if($argc != 4){
  print stderr "usage: perl plot_result_masterevent.pl (out_ps) (result_txt) dem_grd(lon-lat) dem_txt(lon-dep) dem_txt(dep-lat)\n";
  die;
}



@stlon = (143.9775, 143.9867, 144.0017, 144.0042, 144.0160);
@stlat = (43.3797, 43.3955, 43.3818, 43.3903, 43.3695);
@stname = ("V.MEAB", "V.MEAA", "V.PMNS", "V.NSYM", "V.MNDK");

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
$lon_w = 144.0;
$lon_e = 144.02;
$lat_s = 43.375;
$lat_n = 43.39;
$dep_min = -2.0;
$dep_max = 1.5;

$mapsize_x = 10.0;
$mapsize_y = `echo $lon_e $lat_n | gmt mapproject -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n`;
@tmp = split /\s+/, $mapsize_y;
$mapsize_y = $tmp[1];
$mapsize_z = 4.5;
$dx = 11.5;
$dy = 5.9;

$scale_x = 9.0;
$scale_y = $mapsize_y + 0.9;
$scale_dist = "1.0k";

$cpt_x = $mapsize_x / 2;
$cpt_y = -1.0;
$cpt_len_x = $mapsize_x;
$cpt_len_y = "0.3h";

$symbolsize = 0.5;
$symbolwidth = "1.0p";
$symbolsize_st = 0.4;
$title_x = 0.0;
$title_y = $mapsize_y + 0.4;

$color_p1 = "red";
$color_p2 = "green4";
$color_p3 = "blue";

$color_p1 = "black";
$color_p2 = $color_p1;
$color_p3 = $color_p1;


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
                      -Bxya0.02f0.01 -BWeSn -Lx$scale_x/$scale_y+c$lat_s+w$scale_dist  \\
                      -O -K -P >> $out";
system "gmt grdcontour $dem_lonlat -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -C50 -W0.6p,dimgray -O -K -P >> $out";

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Si$symbolsize_st -W0.9p,black -Gwhite -O -K -P >> $out";
for($j = 0; $j <= $#stlon; $j++){
  print OUT "$stlon[$j] $stlat[$j]\n";
}
close OUT;

for($j = 0; $j <= $#amp_ratio; $j++){
  if($j <= 16){
    $color = $color_p1;
  }elsif($j >= 17 && $j <= 33){
    $color = $color_p2;
  }else{
    $color = $color_p3;
  }
  open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                         -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
  print OUT "$evlon[$j] $evlat[$j] $sigma_evlon[$j] $sigma_evlat[$j]\n";
  close OUT;
}

##96-1 crator
#open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -SE -Gblack -O -K -P >> $out";
#  print OUT "144.0097 43.3826 90.0 0.3 0.15\n";
#close OUT;

##error plot
open OUT, " | gmt pslegend -J -R -Dx0.1c/0.1c+w1c/1.2c+jBL -F+p0.6p,black+gwhite -O -K -P >> $out";
print OUT "N 1\n";
close OUT;
open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},black -E -O -K -P >> $out";
print OUT "143.979 43.363 $sigma_evlon[1] $sigma_evlat[1]\n";
close OUT;


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                      -Bxf0.01 -Bya1f0.5+l\"Depth (km)\" -BWsen -O -K -P >> $out";
for($j = 0; $j <= $#amp_ratio; $j++){
  if($j <= 16){
    $color = $color_p1;
  }elsif($j >= 17 && $j <= 33){
    $color = $color_p2;
  }else{
    $color = $color_p3;
  }
  open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                         -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
    print OUT "$evlon[$j] $evdep[$j] $sigma_evlon[$j] $sigma_evdep[$j]\n";
  close OUT;
}

system "gmt psxy $dem_londep -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max -W0.6p,black -O -K -P >> $out";

##error plot
open OUT, " | gmt pslegend -J -R -Dx0.1c/0.1c+w1c/1.2c+jBL -F+p0.6p,black+gwhite -O -K -P >> $out";
print OUT "N 1\n";
close OUT;
open OUT, " | gmt psxy -J -R \\
                       -Sa$symbolsize -W${symbolwidth},black -E -O -K -P >> $out";
print OUT "143.979 1.4 $sigma_evlon[1] $sigma_evdep[1]\n";
close OUT;

open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -X$dx -Y$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                      -Byf0.01 -Bxa1f0.5+l\"Depth (km)\" -BwSen -O -K -P >> $out";
for($j = $#amp_ratio; $j >= 0; $j--){
  if($j <= 16){
    $color = $color_p1;
  }elsif($j >= 17 && $j <= 33){
    $color = $color_p2;
  }else{
    $color = $color_p3;
  }
  open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                         -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
  print OUT "$evdep[$j] $evlat[$j] $sigma_evdep[$j] $sigma_evlat[$j]\n";
  close OUT;
}

system "gmt psxy $dem_deplat -JX$mapsize_z/$mapsize_y -R/$dep_min/$dep_max/$lat_s/$lat_n -W0.6p,black -O -K -P >> $out";

##error plot
open OUT, " | gmt pslegend -J -R -Dx4.4c/0.1c+w1.2c/1.2c+jBR -F+p0.6p,black+gwhite -O -K -P >> $out";
print OUT "N 1\n";
close OUT;
open OUT, " | gmt psxy -J -R \\
                       -Sa$symbolsize -W${symbolwidth},black -E -O -K -P >> $out";
print OUT "1.4 43.363 $sigma_evdep[1] $sigma_evlat[1]\n";
close OUT;

#open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -X$dx -Y$dy >> $out";
#close OUT;

#open OUT, " | gmt pslegend -J -R -Dx0/-3+w$mapsize_x/3c+jTL -O -K >> $out";
#print OUT "N 1\n";
#print OUT "S 0.1c a $symbolsize - ${symbolwidth},$color_p1 0.6c phase 1\n";
#print OUT "S 0.1c a $symbolsize - ${symbolwidth},$color_p2 0.6c phase 2\n";
#print OUT "S 0.1c a $symbolsize - ${symbolwidth},$color_p3 0.6c phase 3\n";
#close OUT;


open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -P >> $out";
close OUT;

system "gmt psconvert $out -A -Tg";

unlink "gmt.conf";
unlink "gmt.history";

