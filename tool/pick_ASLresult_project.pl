#!/usr/bin/env perl
use File::Basename;

if($#ARGV != 8){
  print stderr "usage: perl pick_ASLresult_project.pl (output_ps) (input_asl_result) (stationparam_file) (project_lon1) (project_lat1) (project_lon2) (project_lat2) (project_maxwidth) (project_output_file)\n";
  die;
}

$out = $ARGV[0];
$result = $ARGV[1];
$stationparam = $ARGV[2];
$proj_lon_b = $ARGV[3];
$proj_lat_b = $ARGV[4];
$proj_lon_e = $ARGV[5];
$proj_lat_e = $ARGV[6];
$proj_width = $ARGV[7];
$proj_out = $ARGV[8];

$argc = $#ARGV;
#if($argc != 5){
#  print stderr "usage: perl plot_result_masterevent.pl (out_ps) (result_txt) dem_grd(lon-lat) dem_txt(lon-dep) dem_txt(dep-lat) stationparam\n";
#  die;
#}


open IN, "<", $stationparam;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
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
  next if (substr($_, 0, 1)) == "#";
  @tmp = split /\s+/, $_;
  push @origintime, $tmp[0];
  push @evlon, $tmp[1];
  push @evlat, $tmp[2];
  push @evdep, $tmp[3];
  push @sourceamp, $tmp[4];
  push @residual, $tmp[5];
  push @nsta_use, $tmp[6];
}
close IN;


##project
unlink "$proj_out";
open OUT, " | gmt project -C$proj_lon_b/$proj_lat_b -E$proj_lon_e/$proj_lat_e -Lw -W-${proj_width}/$proj_width -Fxyz -Q \\
            > $proj_out ";
for($i = 0; $i <= $#origintime; $i++){
  print OUT "$evlon[$i] $evlat[$i] $evdep[$i] $origintime[$i] $sourceamp[$i] $residual[$i] $nsta_use[$i]\n";
}
close OUT;

#for($i = 0; $i <= $#projlist; $i++){
#  print stderr "$projlist[$i]\n";
#}

##for GMT
$lon_w = 135.2;
$lon_e = 137.7;
$lat_s = 32.5;
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

$symbolsize = 0.4;
$symbolwidth = "0.4p";
$symbolsize_proj = 0.6;
$symbolwidth_proj = "0.4p";
$symbolsize_st = 0.4;
$title_x = 0.0;
$title_y = $mapsize_y + 0.3;

$annot_a_map = "0.5";
$annot_f_map = "0.25";
$annot_a_dep = "5";
$annot_f_dep = "5";

$color = "gray";
$color_proj = "red";

system "gmt set PS_LINE_JOIN round";
system "gmt set FORMAT_GEO_MAP +D";
system "gmt set FONT_LABEL 17p,Helvetica";
system "gmt set MAP_FRAME_PEN thin";
system "gmt set FONT_ANNOT_PRIMARY 16p,Helvetica";


print stderr "output ps = $out\n";
open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -K -P -X3c -Y12c > $out";
close OUT;

##lon-lat
#system "gmt grdcontour $dem_lonlat -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -C200 -W0.6p,dimgray -O -K -P >> $out";

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
$count = 0;
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evlon[$j] $evlat[$j]\n";
  $count++;
}
close OUT;

#system " gmt psxy $proj_out -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
#                            -Sa$symbolsize_proj -W${symbolwidth_proj},$color_proj -O -K -P >> $out";

system "gmt pscoast -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -Df -W0.8p,black -O -K -P >> $out";
system "gmt psbasemap -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                      -Bxya${annot_a_map}f${annot_f_map} -BWeSn -Lx$scale_x/$scale_y+c$lat_s+w$scale_dist  \\
                      -O -K -P >> $out";

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n -W1p,black -O -K -P >> $out";
  print OUT "$proj_lon_b $proj_lat_b\n$proj_lon_e $proj_lat_e\n";
close OUT;

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Si$symbolsize_st -W0.9p,black -Gwhite -O -K -P >> $out";
for($j = 0; $j <= $#stlon; $j++){
  print OUT "$stlon[$j] $stlat[$j]\n";
}
close OUT;

open OUT, " | gmt psxy -JM$mapsize_x -R$lon_w/$lon_e/$lat_s/$lat_n \\
                       -Sa$symbolsize_proj -W${symbolwidth_proj},$color_proj -O -K -P >> $out";
open IN, "<", $proj_out;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @origintime_proj, $tmp[3];
  push @evlon_proj, $tmp[0];
  push @evlat_proj, $tmp[1];
  push @evdep_proj, $tmp[2];
  push @sourceamp_proj, $tmp[4];
  push @residual_proj, $tmp[5];
  push @nsta_use_proj, $tmp[6];
  print OUT "$tmp[0] $tmp[1]\n";
}
close IN;
close OUT;
$proj_count = $#origintime_proj + 1;

open OUT, ">", $proj_out;
for($i = 0; $i < $proj_count; $i++){
  print OUT "$origintime_proj[$i] $evlon_proj[$i] $evlat_proj[$i] $evdep_proj[$i] $sourceamp_proj[$i] $residual_proj[$i] $nsta_use_proj[$i]\n";
}
close OUT;
  

open OUT, " | gmt pstext -JX$mapsize_x/$mapsize_y -R0/$mapsize_x/0/$mapsize_y -N -F+f+jLB -O -K -P >> $out";
  print OUT "$title_x $title_y 13p,Helvetica,black N=${count} @;red;N=${proj_count}@;;\n";
close OUT;



open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -Y-$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                      -Bxf${annot_f_map} -Bya${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BWsen -O -K -P >> $out";

open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                       -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evlon[$j] $evdep[$j]\n";
}
close OUT;
open OUT, " | gmt psxy -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max \\
                       -Sa$symbolsize_proj -W${symbolwidth_proj},$color_proj -O -K -P >> $out";
for($i = 0; $i < $proj_count; $i++){
  print OUT "$evlon_proj[$i] $evdep_proj[$i]\n";
}
close OUT;


#system "gmt psxy $dem_londep -JX$mapsize_x/-$mapsize_z -R$lon_w/$lon_e/$dep_min/$dep_max -W0.6p,black -O -K -P >> $out";

open OUT, " | gmt psxy -JX1/1 -R0/1/0/1 -Sc0.1 -O -K -P -X$dx -Y$dy >> $out";
close OUT;

##lon-dep
system "gmt psbasemap -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                      -Byf${annot_f_map} -Bxa${annot_a_dep}f${annot_f_dep}+l\"Depth (km)\" -BwSen -O -K -P >> $out";


open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                       -Sa$symbolsize -W${symbolwidth},$color -O -K -P >> $out";
for($j = 0; $j <= $#origintime; $j++){
  print OUT "$evdep[$j] $evlat[$j]\n";
}
close OUT;

open OUT, " | gmt psxy -JX$mapsize_z/$mapsize_y -R$dep_min/$dep_max/$lat_s/$lat_n \\
                       -Sa$symbolsize_proj -W${symbolwidth_proj},$color_proj -O -K -P >> $out";
for($j = 0; $j < $proj_count; $j++){
  print OUT "$evdep_proj[$j] $evlat_proj[$j]\n";
}
close OUT;
#system "gmt psxy $dem_deplat -JX$mapsize_z/$mapsize_y -R/$dep_min/$dep_max/$lat_s/$lat_n -W0.6p,black -O -K -P >> $out";

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

