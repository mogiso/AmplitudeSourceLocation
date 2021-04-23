#!/usr/bin/env perl
use Math::Trig qw(tan atan pi deg2rad great_circle_distance);
#pick up locations from the result of asl

if($#ARGV != 4){
  print stderr "usage: perl pick_ASLresult.pl (input_asl_result) (output_ASL_result) (ref_lon) (ref_lat) (max_distance)\n";
  die;
}

$earth_radius = 6370.291;

$asl_in = $ARGV[0];
$asl_out = $ARGV[1];
$ref_lon = $ARGV[2];
$ref_lat = $ARGV[3];
$max_dist = $ARGV[4];

$ref_lat = geograph2geocen($ref_lat);
@ref_lonlat = NESW($ref_lon, $ref_lat);

##read asl result
open IN, "<", $asl_in;
while(<IN>){
  chomp;
  push @asl_result, $_;
  next if($#asl_result == 0);
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @evlon, $tmp[1];
  push @evlat, $tmp[2];
}
close IN;


##calculate distance and output
open OUT, ">", $asl_out;
print OUT "$asl_result[0]\n";
for($i = 0; $i <= $#evlon; $i++){
  $lat_tmp = geograph2geocen($evlat[$i]);
  @evlonlat = NESW($evlon[$i], $lat_tmp);
  $epdist = great_circle_distance(@ref_lonlat, @evlonlat, $earth_radius);
  if($epdist <= $max_dist){
    print OUT "$asl_result[$i + 1]\n";
  }
}
close OUT;

sub NESW{
  deg2rad($_[0]), deg2rad(90 - $_[1])
}
sub geograph2geocen{
  my $r1 = 6378137.0;
  my $f = 1.0 / 298.257222101;
  my $r2 = $r1 * (1.0 - $f);
  return atan(($r2 / $r1) ** 2 * tan($_[0] * pi / 180.0)) * 180.0 / pi;
}


