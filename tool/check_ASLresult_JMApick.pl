#!/usr/bin/env perl
use Time::Local qw(timelocal);

##compare result of ASL method and JMA-formatted hypocenter file
##Origintime in the resultfile of ASL method should be converted to absolute time
#
if($#ARGV != 3 && $#ARGV != 5){
  print stderr "usage: perl check_ASLresult_JMApick.pl (asl_result_file) (JMA-formatted_hypocenter) (not_matched_asl_result) (matched_asl_result) (subevent_amplitude_file optional) (output_subevent_amplitude_file optional)\n";
  die;
}

$aslresult = $ARGV[0];
$jmapick = $ARGV[1];
$asl_passed = $ARGV[2];
$asl_matched = $ARGV[3];
$subevent_amplitude = $ARGV[4];
$subevent_amplitude_out = $ARGV[5];

$ot_diff_threshold = 60;

$tt_exe = "/home/mogiso/kii_shallow_LF/tool/calc_traveltime_jma2001";
$tmp_stlon = 136.5;
$tmp_stlat = 33.0;

##read result of asl method
$read_comment_flag = 0;
open IN, "<", $aslresult;
while(<IN>){
  chomp;
  if((substr $_, 0, 1) eq "#" && $read_comment_flag == 0){
    $comment = $_;
    $read_comment_flag = 1;
    next;
  }
  push @asl_result, $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @evlon, $tmp[1];
  push @evlat, $tmp[2];
  push @evdep, $tmp[3];
  @date = split /[-T:]/, $tmp[0];
  $month = $date[1] - 1;
  $year = $date[0] - 1900;
  push @asl_unixtime, (timelocal($date[5], $date[4], $date[3], $date[2], $month, $year));
}
close IN;

if($#ARGV == 5){
  $read_comment_flag = 0;
  open IN, "<", $subevent_amplitude;
  while(<IN>){
    chomp;
    if((substr $_, 0, 1) eq "#" && $read_comment_flag == 0){
      $comment2 = $_;
      $read_comment_flag = 1;
      next;
    }
    push @subevent_amp_list, $_;
  }
  close IN;
}
  


##read jma-formatted hypocenter catalogue
open IN, "<", $jmapick;
while(<IN>){
  chomp;
  next if ((substr $_, 0, 1) ne "J");
  $flag = substr $_, 95, 1;
  if($flag !~ /[as]/){
    push @jma_catalogue, $_;

    $yy_catalogue = substr $_, 1, 4;
    $mo_catalogue = substr $_, 5, 2;
    $dy_catalogue = substr $_, 7, 2;
    $hh_catalogue = substr $_, 9, 2;
    $mm_catalogue = substr $_, 11, 2;
    $ss_catalogue = int((substr $_, 13, 4) / 100.0);

    $mo_catalogue = $mo_catalogue - 1;
    $yy_catalogue = $yy_catalogue - 1900;
    push @jma_unixtime, (timelocal($ss_catalogue, $mm_catalogue, $hh_catalogue, $dy_catalogue, $mo_catalogue, $yy_catalogue));

    $hypo_lat = substr $_, 21, 3;
    $hypo_lat_m = substr $_, 24, 4;
    $hypo_lat = $hypo_lat + $hypo_lat_m / 100.0 / 60.0;
    $hypo_lon = substr $_, 32, 4;
    $hypo_lon_m = substr $_, 36, 4;
    $hypo_lon = $hypo_lon + $hypo_lon_m / 100.0 / 60.0;
    $hypo_dep = substr $_, 44, 5;
    $hypo_dep = $hypo_dep / 100.0;

    $traveltime = `$tt_exe $hypo_lon $hypo_lat $hypo_dep $tmp_stlon $tmp_stlat`;
    $traveltime =~ s/^\s*(.*?)\s*$/$1/;
    @ttime_array = split /\s+/, $traveltime;
    push @ttime_p, $ttime_array[0];
    push @ttime_s, $ttime_array[1];
    #print stdout "$ttime_array[0] $ttime_array[1]\n";

    #print stderr "$yy_catalogue $mo_catalogue $dy_catalogue $hh_catalogue $mm_catalogue $ss_catalogue $jma_unixtime[$#jma_unixtime]\n";
  }
}
close IN;

##check
for($j = 0; $j <= $#asl_result; $j++){
  $match_flag[$j] = 0;
  for($i = 0; $i <= $#jma_unixtime; $i++){
    #if($jma_unixtime[$i] - $asl_unixtime[$j] > -15.0 && $jma_unixtime[$i] - $asl_unixtime[$j] < $ot_diff_threshold){
    if((($jma_unixtime[$i] + $ttime_p[$i] > $asl_unixtime[$j])        && 
       ($jma_unixtime[$i] + $ttime_p[$i] < $asl_unixtime[$j] + 60.0)) ||
       (($jma_unixtime[$i] + $ttime_s[$i] > $asl_unixtime[$j])        &&
       ($jma_unixtime[$i] + $ttime_s[$i] < $asl_unixtime[$j]))){
      print stderr "$jma_catalogue[$i]\n";
      print stderr "$asl_result[$j]\n\n";
      $match_flag[$j] = 1;
      #last;
    }
  }
}
#
#
#
open OUT, ">", $asl_passed;
print OUT "$comment\n";
open OUT2, ">", $asl_matched;
print OUT2 "$comment\n";
if($#ARGV == 5){
  open OUT3, ">", $subevent_amplitude_out;
  print OUT3 "$comment2\n";
}
for($j = 0; $j <= $#asl_result; $j++){
  if($match_flag[$j] == 0){
    print OUT "$asl_result[$j]\n";
    if($#ARGV == 5){
      print OUT3 "$subevent_amp_list[$j]\n";
    }
  }else{
    print OUT2 "$asl_result[$j]\n";
  }
}
close OUT;
close OUT2;
if($#ARGV == 5){
  close OUT3;
}




