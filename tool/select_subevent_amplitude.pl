#!/usr/bin/env perl
#
#select amplitude data from subevent_amplitude_time.txt

if($#ARGV != 3){
  print stderr "usage: perl select_subevent_amplitude.pl (stationparam_file) (ref_asl_result) (input_subevent) (output_subevent)\n";
  die;
}

$stationparam = $ARGV[0];
$asl_selected = $ARGV[1];
$subevent_amp_file_in = $ARGV[2];
$subevent_amp_file_out = $ARGV[3];

##read station parameter
open IN, "<", $stationparam;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  if($tmp[4] eq ".true."){
    $stuse{$tmp[3]} = 1;
  }
}
close IN;

#read selected asl result
open IN, "<", $asl_selected;
while(<IN>){
  chomp;
  next if (substr $_, 0, 1) eq "#";
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  for($i = 0; $i <= $#tmp; $i++){
    if($tmp[$i] =~ /[:T]/){
      print stderr "$tmp[$i]\n";
      $date_index = $i;
      last;
    }
  }
  $asl_result_date{$tmp[$date_index]} = 1;
}
close IN;

#pick up amplitude frmm subevent_file_in
open IN, "<", $subevent_amp_file_in;
open OUT, ">", $subevent_amp_file_out;
$_ = <IN>;
chomp $_;
$_ =~ s/^\s*(.*?)\s*$/$1/;
@stname = split /\s+/, $_;
shift @stname;
print OUT "#";
for($i = 0; $i < $#stname; $i++){
  if($stuse{$stname[$i]} == 1){
    print OUT " $stname[$i]";
  }
}
print OUT " time\n";

while(<IN>){
  chomp $_;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @amp = split /\s+/, $_;
  if($asl_result_date{$amp[$#amp]} == 1){
    for($i = 0; $i < $#amp; $i++){
      if($stuse{$stname[$i]} == 1){
        print OUT "$amp[$i] ";  
      }
    }
    print OUT " $amp[$#amp]\n";
  }
}
close IN;
close OUT;

