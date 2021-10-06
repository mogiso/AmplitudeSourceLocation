#!/usr/bin/env perl
#
#pick up local maxima from the result of AmplitudeSourceLocation_(PulseWidth|masterevent).F90,

if($#ARGV != 5 && $#ARGV != 7){
  print stderr "usage: perl pick_local_maxima.pl (asl_result_file) (output_picked_asl_result_file) (minimum_sourceamp) (maximum_residual) (minimum_num_of_used_station) (maximum_num_of_used_station) (subevent_amplitude_file optional) (output_picked_subevent_amplitude_file optional)\n";
  die;
}

$in = $ARGV[0];
$out_list = $ARGV[1];
$sourceamp_min = $ARGV[2];
$residual_max = $ARGV[3];
$nsta_use_min = $ARGV[4];
$nsta_use_max = $ARGV[5];
$subevent_amplitude = $ARGV[6];
$subevent_amplitude_out = $ARGV[7];

$lon_w = 134.9;
$lon_e = 137.6;
$lat_s = 32.4;
$lat_n = 33.8;
$dep_min = -1.0;
$dep_max = 21.0;
$dlon = 0.02;
$dlat = 0.02;
$dz = 2.0;
$delta_hypo = 3.0;
$compare_index = 1;

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
  push @nsta_use, $tmp[6];
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

for($i = $compare_index; $i <= $#sourceamp - $compare_index; $i++){
  $diff_back2    = $sourceamp[$i] - $sourceamp[$i - $compare_index];
  $diff_forward2 = $sourceamp[$i] - $sourceamp[$i + $compare_index];
  $diff_back1    = $sourceamp[$i] - $sourceamp[$i - 1];
  $diff_forward1 = $sourceamp[$i] - $sourceamp[$i + 1];
  if($sourceamp[$i] > $sourceamp_min                              && 
     $diff_back1    > 0.0                                         &&
     $diff_forward1 > 0.0                                         && 
#     $diff_back2    > 0.0                                         &&
#     $diff_forward2 > 0.0                                         && 
     $hypo_lon[$i]  > $lon_w                                      &&
     $hypo_lon[$i]  < $lon_e                                      &&
     $hypo_lat[$i]  > $lat_s                                      &&
     $hypo_lat[$i]  < $lat_n                                      &&
     $hypo_dep[$i]  > $dep_min                                    &&
     $hypo_dep[$i]  < $dep_max                                    &&
     $nsta_use[$i] >= $nsta_use_min                               &&
     $nsta_use[$i] <= $nsta_use_max                               &&
     $residual[$i]  < $residual_max                               &&
     abs($hypo_lon[$i] - $hypo_lon[$i - 1]) < $delta_hypo * $dlon &&
     abs($hypo_lon[$i] - $hypo_lon[$i + 1]) < $delta_hypo * $dlon &&
     abs($hypo_lat[$i] - $hypo_lat[$i - 1]) < $delta_hypo * $dlat &&
     abs($hypo_lat[$i] - $hypo_lat[$i + 1]) < $delta_hypo * $dlat 
#     abs($hypo_dep[$i] - $hypo_dep[$i - 1]) <= $delta_hypo * $dz   &&
#     abs($hypo_dep[$i] - $hypo_dep[$i + 1]) <= $delta_hypo * $dz   
  ){
    if($result_localmax_picked[$#result_localmax_picked] eq $result_all[$i - 2]){
      if($residual[$i - 2] > $residual[$i]){
        pop @result_localmax_picked;
        push @result_localmax_picked, $result_all[$i];
        if(@subevent_amp){
          pop @subevent_amp_picked;
          push @subevent_amp_picked, $subevent_amp[$i] ;
        }
      }
    }else{
      push @result_localmax_picked, $result_all[$i];
      push @subevent_amp_picked, $subevent_amp[$i] if (@subevent_amp);
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


