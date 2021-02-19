#!/usr/bin/env perl
#
#pick up local maxima from the result of AmplitudeSourceLocation_masterevent.F90,
#then write hypocentral distribution


$in = $ARGV[0];
$out_list = $ARGV[1];
$sourceamp_min = $ARGV[2];

$lon_w = 136.5;
$lon_e = 137.5;
$lat_s = 32.7;
$lat_n = 34.7;
$dep_min = 0.0;
$dep_max = 16.0;
$residual_max = 0.05;

##read result file
open IN, "<", $in;
while(<IN>){
  next if substr $_, 0, 1 eq "#";
  chomp $_;
  push @result_all, $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @sourceamp, $tmp[4];
  push @hypo_lon, $tmp[1];
  push @hypo_lat, $tmp[2];
  push @hypo_dep, $tmp[3];
  push @residual, $tmp[5];
}
close IN;

for($i = 1; $i <= $#sourceamp - 1; $i++){
  $diff_back = $sourceamp[$i] - $sourceamp[$i - 1];
  $diff_forward = $sourceamp[$i] - $sourceamp[$i + 1];
  if($sourceamp[$i] > $sourceamp_min && $diff_back > 0.0 && $diff_forward > 0.0 && 
     #$hypo_lon[$i] > $lon_w && $hypo_lon[$i] < $lon_e &&
     #$hypo_lat[$i] > $lat_s && $hypo_lat[$i] < $lat_n &&
     #$hypo_dep[$i] >= $dep_min && $hypo_dep[$i] <= $dep_max){
     $residual[$i] < $residual_max){
    push @result_localmax_picked, $result_all[$i];
  }
}

open OUT, ">", $out_list;
print OUT "# amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth ot_tmp\n";
for($i = 0; $i <= $#result_localmax_picked; $i++){
  print OUT "$result_localmax_picked[$i]\n";
}
close OUT;


