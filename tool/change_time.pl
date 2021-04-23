#!/usr/bin/env perl
#
##change format of origintime
#

if($#ARGV != 8){
  print stderr "usage: perl change_time.pl (input_file) (output_file) (index of column of the time) (ref_year) (ref_month) (ref_day) (ref_hour) (ref_minute) (ref_sec)\n";
  die;
}

$in = $ARGV[0];
$out = $ARGV[1];
$num_timeindex = $ARGV[2];
$yy = $ARGV[3];
$mo = $ARGV[4];
$dy = $ARGV[5];
$hh = $ARGV[6];
$mm = $ARGV[7];
$ss = $ARGV[8];

open IN, "<", $in;
open OUT, ">", $out;

print stderr "change time format of $in\n";

while(<IN>){
  chomp $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  $time_index = $num_timeindex;
  if((substr $_, 0, 1) eq "#"){
    print OUT "$_\n";
  }else{
    $sec_tmp = $tmp[$time_index];
    $sec_from_day = $hh * 3600 + $mm * 60 + $ss + $sec_tmp;
    $hour = int($sec_from_day / 3600);
    $hour = sprintf "%02i", $hour;
    $minute = int(($sec_from_day - $hour * 3600) / 60);
    $minute = sprintf "%02i", $minute;
    $sec = $sec_from_day - $hour * 3600 - $minute * 60;
    $sec = sprintf "%.2f", $sec;
    for($i = 0; $i <= $#tmp; $i++){
      if($i == $time_index){
        print OUT "${yy}-${mo}-${dy}T${hour}:${minute}:${sec} ";
      }else{
        print OUT "$tmp[$i] ";
      }
    }
    print OUT "\n";
  }
}
close IN;
close OUT;


