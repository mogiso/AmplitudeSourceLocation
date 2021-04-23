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
    #print stderr "$yy_catalogue $mo_catalogue $dy_catalogue $hh_catalogue $mm_catalogue $ss_catalogue $jma_unixtime[$#jma_unixtime]\n";
  }
}
close IN;
#check

for($j = 0; $j <= $#asl_result; $j++){
  $match_flag[$j] = 0;
  for($i = 0; $i <= $#jma_unixtime; $i++){
    if($jma_unixtime[$i] - $asl_unixtime[$j] > 0.0 && $jma_unixtime[$i] - $asl_unixtime[$j] < $ot_diff_threshold){
      print stderr "$jma_catalogue[$i]\n";
      print stderr "$asl_result[$j]\n\n";
      $match_flag[$j] = 1;
      last;
    }
  }
}



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




