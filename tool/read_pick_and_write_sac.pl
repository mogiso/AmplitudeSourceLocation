#!/usr/bin/env perl
use File::Basename;

##winで作ったpickファイルを読み込む

if($#ARGV != 2){
  die "usage: perl read_pick_and_write_sac.pl (win_pick_file) (sacfile_dir)\n";
}

$pick = $ARGV[0];
$hypodir = $ARGV[1];
$sacdir = $ARGV[2];

open IN, "<", $pick;
$count = 0;
while(<IN>){
  chomp;
  $flag = substr $_, 0, 2;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  if($flag eq "#p" && (substr $_, 19, 3) eq "win"){
    $file_ymd = substr $_, 3, 8; 
    $file_hms = substr $_, 12, 6; 
  }
  next if($flag eq "#p");           ##オリジナル読み取り行は飛ばす
  if($flag eq "#s"){
    $count++;
    if($count == 1){                ##1行目: 基準時刻
      @tmp = split /\s+/, $_;
      @date_tmp = split "/", $tmp[1];
      @time_tmp = split /\:/, $tmp[2];
      print stderr "$date_tmp[0] $date_tmp[1] $date_tmp[2] $time_tmp[0] $time_tmp[1]\n";
      $yy = $date_tmp[0];
      $mo = $date_tmp[1];
      $dy = $date_tmp[2];
      $hh = $time_tmp[0];
      $mm = $time_tmp[1];
    }else{   ##観測点行
      @tmp = split /\s+/, $_;
      ##最終行ならやめ
      next if ($#tmp == 0);
      $stname = $tmp[1];
      ##PかSが読み取られていないと止め
      next if ($tmp[3] == 0.0 || $tmp[5] == 0.0);
      $ptime{$stname} = $tmp[3];
      $stime{$stname} = $tmp[5];
      #print stderr "$stname $ptime{$stname} $stime{$stname}\n";
    }
  }
  last if($flag eq "#f");
}
close IN;

$hypofile = basename $pick;
$hypofile = "$hypodir/$hypofile";
open IN, "<", $hypofile;
$_ = <IN>;
close IN;
chomp $_;
$_ =~ s/^\s*(.*?)\s*$/$1/;
@tmp = split /\s+/, $_;
$hypo_yr = $tmp[0];
$hypo_mo = $tmp[1];
$hypo_dy = $tmp[2];
$hypo_hh = $tmp[3];
$hypo_mm = $tmp[4];
$hypo_ss = $tmp[5];
$evlat = $tmp[6];
$evlon = $tmp[7];
$evdep = $tmp[8];

$ot_from_day = $hypo_hh * 60 * 60 + $hypo_mm * 60 + $hypo_ss;
$hypo_ss_int = int($hypo_ss);
$hypo_msec = int(($hypo_ss - $hypo_ss_int) * 1000.0 + 0.5);

@keys = keys %ptime;
for($i = 0; $i <= $#keys; $i++){
  #print stderr "$keys[$i] $ptime{$keys[$i]} $stime{$keys[$i]}\n";

  @filename = ();
  push @filename, "${sacdir}/$file_ymd.$file_hms.$keys[$i].N.sac";
  push @filename, "${sacdir}/$file_ymd.$file_hms.$keys[$i].E.sac";
  push @filename, "${sacdir}/$file_ymd.$file_hms.$keys[$i].U.sac";

  ##ptime, stimeそれぞれ一日の最初からの秒を計算
  $ptime_from_day = $hh * 60 * 60 + $mm * 60 + $ptime{$keys[$i]};
  $stime_from_day = $hh * 60 * 60 + $mm * 60 + $stime{$keys[$i]};


  ##SACファイルごとに書き込む
  for($j = 0; $j <= $#filename; $j++){
    print stderr "$filename[$j]\n";
    next if (!-f $filename[$j]);

    open SAC, "+<", $filename[$j];
      seek SAC, 5 * 4, SEEK_SET;
      read SAC, $begin, 4;
      $begin = unpack "f", $begin;
      seek SAC, 28, SEEK_SET;
      read SAC, $origin, 4;
      $origin = unpack "f", $origin;
      seek SAC, 288, SEEK_SET;
      read SAC, $sachour, 4;
      $sachour = unpack "i", $sachour;
      read SAC, $sacmin, 4;
      $sacmin = unpack "i", $sacmin;
      read SAC, $sacsec, 4;
      $sacsec = unpack "i", $sacsec;
      read SAC, $sacmsec, 4;
      $sacmsec = unpack "i", $sacmsec;
      $begin_from_day = $sachour * 60 * 60 + $sacmin * 60  + $sacsec + $sacmsec / 1000.0 + $begin;

      $ot_sac = 0.0;
      $begin_sac = sprintf "%.2f", ($begin_from_day - $ot_from_day);
      $ptime_sac = sprintf "%.2f", ($ptime_from_day - $ot_from_day);
      $stime_sac = sprintf "%.2f", ($stime_from_day - $ot_from_day);

      print stderr "begin_org and new = $begin, $begin_sac\n";
      print stderr "ptime_sac = $ptime_sac, $stime_sac\n";

      print stderr "writing begin of $filename[$j]\n";
      seek SAC, 5 * 4, SEEK_SET;
      $buf = pack "f", $begin_sac;
      print SAC $buf;
      seek SAC, 28, SEEK_SET;
      $buf = pack "f", $ot_sac;
      print SAC $buf;

      print stderr "writing ptime on A field of $filename[$j]...\n";
      seek SAC, 32, SEEK_SET;
      $buf = pack "f", $ptime_sac;
      print SAC $buf;
      seek SAC, 480, SEEK_SET;
      print SAC "Pmanu   ";
      print stderr "writing stime on T0 field of $filename[$j]...\n";
      seek SAC, 40, SEEK_SET;
      $buf = pack "f", $stime_sac;
      print SAC $buf;
      seek SAC, 488, SEEK_SET;
      print SAC "Smanu   ";

      print stderr "writing NZHOUR, NZMIN, NZSEC, NZMSEC\n";
      seek SAC, 288, SEEK_SET;
      $buf = pack "i", $hypo_hh;
      print SAC $buf;
      $buf = pack "i", $hypo_mm;
      print SAC $buf;
      $buf = pack "i", $hypo_ss_int;
      print SAC $buf;
      $buf = pack "i", $hypo_msec;
      print SAC $buf;
    close SAC;
  }
}

