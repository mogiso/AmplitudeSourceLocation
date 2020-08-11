#!/usr/bin/env perl

##run hypomh

$hypomh_exe = "/usr/local/win/bin/hypomh";
$hypomh_in = "hypomh_in";
$struct = "../struct_105.tbl";
$finaldir = "../final.man";

@pick = @ARGV;
foreach $pickfile (@pick){
  open IN, "<", $pickfile;
  open OUT, ">", $hypomh_in;
  while(<IN>){
    next if (substr $_, 0, 2) ne "#s";
    chomp $_;
    $stname = substr $_, 3, 6;
    if($stname ne "" && substr($stname, 2, 1) ne "/" && $stname ne "V.MEAB" && $stname ne "V.KNGM"){
      substr ($_, 31, 12, "00.000 0.000");
    }
    $new = substr $_, 3, 100;
    #if($new ne ""){
      print "$new\n";
      print OUT "$new\n";
    #}
  }
  close IN;
  close OUT;
  $finalfile = "$finaldir/$pickfile";
  system "$hypomh_exe $struct $hypomh_in $finalfile /dev/null";
  unlink $hypomh_in;
}

