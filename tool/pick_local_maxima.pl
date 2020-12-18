#!/usr/bin/env perl
#
#pick up local maxima from the result of AmplitudeSourceLocation_masterevent.F90,
#then write hypocentral distribution


$in = $ARGV[0];
$out_ps = $ARGV[1];
$out_list = $ARGV[2];
$sourceamp_min = $ARGV[3];

##read result file
open IN, "<", $in;
while(<IN>){
  next if substr $_, 0, 1 eq "#";
  chomp $_;
  push @result_all, $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @sourceamp, $tmp[0];
}
close IN;

for($i = 1; $i <= $#sourceamp - 1; $i++){
  $diff_back = $sourceamp[$i] - $sourceamp[$i - 1];
  $diff_forward = $sourceamp[$i] - $sourceamp[$i + 1];
  if($sourceamp[$i] > $sourceamp_min && $diff_back > 0.0 && $diff_forward > 0.0){
    push @result_localmax_picked, $result_all[$i];
  }
}

open OUT, ">", $out_list;
print OUT "# amp_ratio sigma_ampratio longitude sigma_lon latitude sigma_lat depth sigma_depth ot_tmp\n";
for($i = 0; $i <= $#result_localmax_picked; $i++){
  print OUT "$result_localmax_picked[$i]\n";
}
close OUT;


