#!/usr/bin/env perl
#
#pick up local maxima from the result of AmplitudeSourceLocation_(PulseWidth|masterevent).F90,


$in = $ARGV[0];
$out_list = $ARGV[1];
$sourceamp_min = $ARGV[2];
$subevent_amplitude = $ARGV[3];
$subevent_amplitude_out = $ARGV[4];

$lon_w = 136.5;
$lon_e = 137.5;
$lat_s = 32.7;
$lat_n = 34.7;
$dep_min = 0.0;
$dep_max = 16.0;
$residual_max = 0.05;

##read result file
$read_comment = 0;
open IN, "<", $in;
while(<IN>){
  chomp;
  if((substr $_, 0, 1) eq "#"){
    if($read_comment == 0){
      $read_comment = 1;
      $comment1 = $_;
    }else{
      next;
    }
  }
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

$read_comment = 0;
if(-f $subevent_amplitude && $subevent_amplitude ne ""){
  open IN, "<", $subevent_amplitude;
  while(<IN>){
    chomp;
    if((substr $_, 0, 1) eq "#"){
      if($read_comment == 0){
        $read_comment = 1;
        $comment2 = $_;
      }else{
        next;
      }
    }
    push @subevent_amp, $_;
  }
}

for($i = 1; $i <= $#sourceamp - 1; $i++){
  $diff_back = $sourceamp[$i] - $sourceamp[$i - 1];
  $diff_forward = $sourceamp[$i] - $sourceamp[$i + 1];
  if($sourceamp[$i] > $sourceamp_min && $diff_back > 0.0 && $diff_forward > 0.0 && 
     #$hypo_lon[$i] > $lon_w && $hypo_lon[$i] < $lon_e &&
     #$hypo_lat[$i] > $lat_s && $hypo_lat[$i] < $lat_n &&
     #$hypo_dep[$i] >= $dep_min && $hypo_dep[$i] <= $dep_max){
     $residual[$i] < $residual_max){
    push @result_localmax_picked, $result_all[$i];
    if(@subevent_amp){
      push @subevent_amp_picked, $subevent_amp[$i];
    }
  }
}

open OUT, ">", $out_list;
print OUT "$comment1\n";
for($i = 0; $i <= $#result_localmax_picked; $i++){
  print OUT "$result_localmax_picked[$i]\n";
}
close OUT;

if(@subevent_amp_picked){
  open OUT, ">", $subevent_amplitude_out;
  print OUT "$comment2\n";
  for($i = 0; $i <= $#subevent_amp_picked; $i++){
    print OUT "$subevent_amp_picked[$i]\n";
  }
  close OUT;
}


