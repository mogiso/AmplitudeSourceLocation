#!/usr/bin/env perl
use File::Basename;

##make traveltimelist from sac

$out = $ARGV[0];


@station = ("V.MEAB", "V.MEAA", "V.NSYM", "V.MNDK", "V.KNGM", "V.PNMM");
#@station = ("V.MEAB", "V.MEAA", "V.NSYM", "V.MNDK", "V.KNGM", "V.PNMM", "V.MEAK");


$out_p = "${out}_P.txt";
$out_s = "${out}_S.txt";

open OUT_P, ">", $out_p;
open OUT_S, ">", $out_s;

print OUT_P "# @station event_index\n";
print OUT_S "# @station event_index\n";

for ($i = 1; $i <= $#ARGV; $i++){
  $subevent_index = basename $ARGV[$i];
  for($j = 0; $j <= $#station; $j++){
    $sacfile = "${ARGV[$i]}.$station[$j].U.sac";
    print stderr "reading $sacfile\n";
    open SAC, "<", $sacfile;
    seek SAC, 32, SEEK_SET;
    read SAC, $buf, 4;
    $ptime[$j] = sprintf "%.2f", (unpack "f", $buf);
    seek SAC, 40, SEEK_SET;
    read SAC, $buf, 4;
    $stime[$j] = sprintf "%.2f", (unpack "f", $buf);
    close SAC;
  }
  next if $ptime[0] == -12345.0;
  print OUT_P "@ptime $subevent_index\n";
  print OUT_S "@stime $subevent_index\n";
}


close OUT_P;
close OUT_S;
