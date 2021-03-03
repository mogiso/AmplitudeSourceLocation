#!/usr/bin/env perl
#
##change format of origintime
#
$in = $ARGV[0];
$out = $ARGV[1];
$yy = $ARGV[2];
$mo = $ARGV[3];
$dy = $ARGV[4];
$hh = $ARGV[5];
$mm = $ARGV[6];
$ss = $ARGV[7];

open IN, "<", $in;
open OUT, ">", $out;

print stderr "change time format of $in\n";

while(<IN>){
  chomp $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  $time_index = 0;
  #$time_index = $#tmp;
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


