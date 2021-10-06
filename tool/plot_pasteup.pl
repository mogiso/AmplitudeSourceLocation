#!/usr/bin/env perl
use Math::Trig qw(tan atan pi deg2rad great_circle_distance);
use List::Util qw(max min);

$out = $ARGV[0];
$winchfile = $ARGV[1];
$stparam = $ARGV[2];
$eqdate = $ARGV[3];
$eqlon = $ARGV[4];
$eqlat = $ARGV[5];
$eqdep = $ARGV[6];

$filter = "-fl 2 -fh 8 -fs 10";
#$filter = "";
$stcomp = "-comp Z";

$dumpwin = "/home/mogiso/dumpwin/dumpwin";
$wavedir = "/home/mogiso/kii_shallow_LF";

##read station parameter
open IN, "<", $stparam;
$i = 0;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  if($tmp[4] eq ".true."){
    push @stlon, $tmp[0];
    push @stlat, $tmp[1];
    push @stdep, $tmp[2];
    push @stname, $tmp[3];
  }
  #($stlon[$i], $stlat[$i], $stdep[$i], $stname[$i], $stflag[$i], $cor1[$i], $cor2[$i], $site[$i], $noise[$i]) = split /\s+/, $_;
  #$i++;
}
close IN;
($eqyy, $eqmo, $eqdy, $eqhh, $eqmm, $eqss) = split /[-T: ]/, $eqdate;
$eqyy = sprintf "%04i", $eqyy;
$eqmo = sprintf "%02i", $eqmo;
$eqdy = sprintf "%02i", $eqdy;
$eqhh = sprintf "%02i", $eqhh;
$eqmm = sprintf "%02i", $eqmm;

$waveyy_b = $eqyy;
$wavemo_b = $eqmo;
$wavedy_b = $eqdy;
$wavehh_b = $eqhh;
$wavemm_b = $eqmm;
$wavess_b = $eqss - 20;

if($wavess_b < 0.0){
  $wavess_b = $wavess_b + 60;
  $wavemm_b = $wavemm_b - 1;
}
if($wavemm_b < 0){
  $wavemm_b = $wavemm_b + 60;
  $wavehh_b = $wavehh_b - 1;
}
if($wavehh_b < 0){
  $wavehh_b = $wavehh_b + 24;
  $wavedy_b = $wavedy_b - 1;
}
if($wavedy_b == 0){
  $wavedy_b = 31;
  $wavemo_b = $wavemo_b - 1;
}
if($wavemo_b == 0){
  $wavemo_b = 12;
  $waveyy_b = $waveyy_b - 1; 
}



$waveyy_e = $eqyy;
$wavemo_e = $eqmo;
$wavedy_e = $eqdy;
$wavehh_e = $eqhh;
$wavemm_e = $eqmm;
$wavess_e = $eqss + 90;
if($wavess_e >= 60.0 && $wavess_e < 120.0){
  $wavess_e = $wavess_e - 60.0;
  $wavemm_e = $wavemm_e + 1;
}elsif($wavess_e >= 120.0){
  $wavess_e = $wavess_e - 120.0;
  $wavemm_e = $wavemm_e + 2;
}
if($wavemm_e >= 60){
  $wavemm_e = $wavemm_e - 60;
  $wavehh_e = $wavehh_e + 1;
}
if($wavehh_e >= 24){
  $wavehh_e = $wavehh_e - 24;
  $wavedy_e = $wavedy_e + 1;
}
if($wavedy_e >= 32){
  $wavedy_e = 1;
  $wavemo_e = $wavemo_e + 1;
}
if($wavemo_e > 12){
  $wavemo_e = 1;
  $waveyy_e = $waveyy_e + 1;
}

$wavehh_b = sprintf "%02d", $wavehh_b;
$wavedy_b = sprintf "%02d", $wavedy_b;

$winfile = "$wavedir/${waveyy_b}${wavemo_b}${wavedy_b}/${waveyy_b}${wavemo_b}${wavedy_b}.${wavehh_b}0000.win32";
for($i = 0; $i <= $#stlon; $i++){
  $stlat_c = geograph2geocen($stlat[$i]);
  @point0 = NESW($stlon[$i], $stlat_c);
  $eqlat_c = geograph2geocen($eqlat);
  @point1 = NESW($eqlon, $eqlat_c);

  $distance[$i] = great_circle_distance(@point0, @point1, 6370.291);
  $delta = $distance[$i] / 6370.291;
  $distance[$i] = sqrt((6370.291 - $stdep[$i]) ** 2 + (6370.291 - $eqdep) ** 2 
                        - 2.0 * (6370.291 - $stdep[$i]) * (6370.291 - $eqdep) * cos($delta));
  #print stderr "$stname[$i] $stlon[$i] $stlat[$i] $stdep[$i] $eqlon $eqlat $eqdep $distance[$i]\n";
}
$maxdist = max @distance;
#$maxdist = 50;
$mindist = min @distance;

$mindist = int($mindist / 10.0 - 1.0) * 10.0;
$maxdist = int($maxdist / 10.0 + 1.0) * 10.0;

$range_x = "${waveyy_b}-${wavemo_b}-${wavedy_b}T${wavehh_b}:${wavemm_b}:${wavess_b}/${waveyy_e}-${wavemo_e}-${wavedy_e}T${wavehh_e}:${wavemm_e}:${wavess_e}";
$range_y = "$mindist/$maxdist";

##GMT
$size_x = 23;
$size_y = 15;
if($filter eq ""){
  $amp = 50;
}else{
  #$amp = 0.5;
  $amp = 3;
}
$scale_x = 8.6;
$scale_y = 6.6;
$annot_x = "a10sf5s";
$annot_x2 = "a1M";
$annot_y = "a20f10";
$label_x = "Time";
$label_y = "Hypocentral distance (km)";
$pen = "0.5p,black";

$txt_x ="${waveyy_e}-${wavemo_e}-${wavedy_e}T${wavehh_e}:${wavemm_e}:${wavess_e}";
$txt_pen = "12p,Helvetica";

$title_x = 0;
$title_y = $size_y + 0.6;
$title_pen = "14p,Helvetica";

system "gmt set PS_LINE_JOIN round";
system "gmt set MAP_FRAME_PEN thin";
system "gmt set FORMAT_CLOCK_MAP hh:mm";

open OUT, " | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -K -X3c -Y3c > $out";
close OUT;

open OUT, " | gmt pswiggle -JX${size_x}T/-${size_y} -R$range_x/$range_y -Z$amp -W$pen \\
                           -Bpx${annot_x}+l\"$label_x\" -Bpy${annot_y}+l\"$label_y\" -BWSen \\
                           -Bsx${annot_x2} \\
                           -O -K >> $out";
  $j = -1;
  open IN, "$dumpwin $winfile $winchfile $stcomp $filter @stname | ";
  while(<IN>){
    chomp;
    $_ =~ s/^\s*(.*?)\s*$/$1/;
    if($_ =~ />/){
      print OUT "$_\n";
      #print stdout "$_\n";
      $j++;
      $k = -1;
    }else{
      $k++;
      @tmp = split /\s+/, $_;
      if($k % 2 == 0){ 
        print OUT "${tmp[0]}-${tmp[1]}-${tmp[2]}T${tmp[3]}:${tmp[4]}:${tmp[5]} $distance[$j] $tmp[6]\n";
        #print stdout "${tmp[0]}-${tmp[1]}-${tmp[2]}T${tmp[3]}:${tmp[4]}:${tmp[5]} $distance[$j] $tmp[6]\n";
      }
    }
  }
  close IN;
       
  #for($i = 0; $i <= $#stname; $i++){
  #  next if $stflag[$i] eq ".false.";
  #  open IN, "$dumpwin $winfile $winchfile $stname[$i] $stcomp $filter | ";
  #  while(<IN>){
  #    $j++;
  #    $_ =~ s/^\s*(.*?)\s*$/$1/;
  #    @tmp = split /\s+/, $_;
  #    if($j % 2 == 0){ 
  #      print OUT "${tmp[0]}-${tmp[1]}-${tmp[2]}T${tmp[3]}:${tmp[4]}:${tmp[5]} $distance[$i] $tmp[6]\n";
  #    }
  #  }
  #  close IN;
  #  print OUT ">\n";
  #}
close OUT;

open OUT, " | gmt pstext -JX${size_x}T/-${size_y} -R$range_x/$range_y -F+f+j -N -O -K -P >> $out";
  for($i = 0; $i <= $#stname; $i++){
    next if $stflag[$i] eq ".false.";
    print OUT "$txt_x $distance[$i] $txt_pen LM $stname[$i]\n";
  }
close OUT;

$eqlon = sprintf "%.2f", $eqlon;
$eqlat = sprintf "%.2f", $eqlat;
$eqdep = sprintf "%i", $eqdep;
open OUT, " | gmt pstext -JX${size_x}/${size_y} -R0/$size_x/0/$size_y -F+f+j -N -O -K -P >> $out";
  print OUT "$title_x $title_y $title_pen LB $eqyy-$eqmo-$eqdy $eqhh:$eqmm:$eqss ${eqlat}N ${eqlon}E ${eqdep}km\n";
close OUT;

    
open OUT, " | gmt psxy -JX1c -R0/1/0/1 -Sc0.1 -O >> $out";
close OUT;

system "gmt psconvert $out -Tf";

unlink "gmt.conf";
unlink "gmt.history";


##距離計算のためのサブルーチン
sub NESW{
  deg2rad($_[0]), deg2rad(90 - $_[1])
}
##地理緯度を地心緯度に変換
sub geograph2geocen{

  ##楕円体の定数(GRS80)
  my $r1 = 6378137.0;
  my $f = 1.0 / 298.257222101;
  my $r2 = $r1 * (1.0 - $f);

  return atan(($r2 / $r1) ** 2 * tan($_[0] * pi / 180.0)) * 180.0 / pi;
}


