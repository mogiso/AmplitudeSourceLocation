#!/usr/bin/env perl
##Convert JMA-formatted travel time data for TraveltimeSourceLocation_masterevent

if($#ARGV != 3){
  die "usage: perl JMApick2masterevent.pl (station_param) (JMA-formatted pick file) (traveltime file of P-waves) (traveltime file of S-waves\n";
}

$stationparam = $ARGV[0];
$jmapick = $ARGV[1];
$ttdata_p = $ARGV[2];
$ttdata_s = $ARGV[3];

open IN, "<", $stationparam;
<IN>;
while(<IN>){
  chomp;
  $_ =~ s/^\s*(.*?)\s*$/$1/;
  @tmp = split /\s+/, $_;
  push @stname, $tmp[3];
  $stflag{$tmp[3]} = 1;
}
close IN;

open IN, "<", $jmapick;
open OUT1, ">", $ttdata_p;
open OUT2, ">", $ttdata_s;

print OUT1 "# ";
print OUT2 "# ";
for($i = 0; $i <= $#stname; $i++){
  print OUT1 "$stname[$i] ";
  print OUT2 "$stname[$i] ";
  if($i == $#stname){
    print OUT1 "\n";
    print OUT2 "\n";
  }
}

%ttime_p = ();
%ttime_s = ();
$eqflag = 0;
$pflag = 0;
$sflag = 0;
while(<IN>){
  chomp;
  $pick_flag = substr $_, 0, 1;

  if($pick_flag eq "J"){
    $date = substr $_, 1, 8;
    $day = substr $_, 7, 2;
    $time = substr $_, 9, 8;
    $hour = substr $_, 9, 2;
    $minute = substr $_, 11, 2;
    $sec = (substr $_, 13, 4) / 100.0;
    $sec_from_day = $hour * 3600.0 + $minute * 60.0 + $sec;

  }elsif($pick_flag eq "E"){
    if($eqflag == 1){
      $pflag = 1;
      $sflag = 1;
      for($i = 0; $i <= $#stname; $i++){
        $ttime_p{$stname[$i]} = sprintf "%.2f", $ttime_p{$stname[$i]};
        $ttime_s{$stname[$i]} = sprintf "%.2f", $ttime_s{$stname[$i]};
        $pflag = 0 if $ttime_p{$stname[$i]} == 0.0;
        $sflag = 0 if $ttime_s{$stname[$i]} == 0.0;
      }
      if($pflag == 1){
        for($i = 0; $i <= $#stname; $i++){
          print OUT1 "$ttime_p{$stname[$i]} ";
        }
        print OUT1 "${date}.${time}\n";
      }
      if($sflag == 1){
        for($i = 0; $i <= $#stname; $i++){
          print OUT2 "$ttime_s{$stname[$i]} ";
        }
        print OUT2 "${date}.${time}\n";
      }
    }
    $eqflag = 0;
    $pflag = 0;
    $sflag = 0;
    %ttime_p = ();
    %ttime_s = ();

  }elsif($pick_flag eq "_"){
    $station_tmp = substr $_, 1, 6;
    $day_tmp = substr $_, 13, 2;
    $phase1 = substr $_, 15, 4;
    $hour_tmp = substr $_, 19, 2;
    $minute_tmp1 = substr $_, 21, 2;
    $sec_tmp1 = (substr $_, 23, 4) / 100.0;
    $phase2 = substr $_, 27, 4;
    $minute_tmp2 = substr $_, 31, 2;
    $sec_tmp2 = (substr $_, 33, 4) / 100.0;

    $eqflag = 1 if $stflag{$station_tmp} == 1;

    if($phase1 =~ /P/){
      $ttime_p{$station_tmp} = ($hour_tmp * 3600 + $minute_tmp1 * 60 + $sec_tmp1) - $sec_from_day;
      $ttime_s{$station_tmp} = 0.0;
      $ttime_p{$station_tmp} = $ttime_p{$station_tmp} + 86400.0 if $day_tmp != $day;
    }elsif($phase1 =~ /S/){
      $ttime_p{$station_tmp} = 0.0;
      $ttime_s{$station_tmp} = ($hour_tmp * 3600 + $minute_tmp1 * 60 + $sec_tmp1) - $sec_from_day;
      $ttime_s{$station_tmp} = $ttime_s{$station_tmp} + 86400.0 if $day_tmp != $day;
    }
    if($phase2 =~ /S/){
      $ttime_s{$station_tmp} = ($hour_tmp * 3600 + $minute_tmp2 * 60 + $sec_tmp2) - $sec_from_day;
      $ttime_s{$station_tmp} = $ttime_s{$station_tmp} + 86400.0 if$day_tmp != $day;
    }
  }
}

close IN;
close OUT1;
close OUT2;

