#!/usr/bin/perl -w

use Getopt::Long;

if($#ARGV != 3 && $#ARGV != 4){
    print " peakdensity.pl <peakfile> <output> <windowsize> <genome_table>\n";
    print "     -l: sum total peak length instead of peak number\n";
    exit;
}

$peakfile=$ARGV[0];
$outputfile=$ARGV[1];
$width=$ARGV[2];
$gt=$ARGV[3];
$lensum=0;
GetOptions('len' => \$lensum);

open(File, $gt) ||die "error: can't open $gt.\n";
while(<File>){
    next if($_ eq "\n");
    chomp;
    my @clm = split(/\t/, $_);
    $chrlen{$clm[0]} = $clm[1];
}
close (File);

foreach $chr (keys(%chrlen)){
    my $binnum = int($chrlen{$chr}/$width) +1;
    for($i=0;$i<$binnum;$i++){
	$array{$chr}[$i]=0;
    }
}
open(ListFile, $peakfile) ||die "cannot open $peakfile.";
while(<ListFile>){
    next if($_ eq "\n" || $_ =~ "chromosome");
    chomp($_);
    my @clm = split(/\t/, $_);
    my $chr = $clm[0];
    if(!exists($chrlen{$chr})){
	print "invalid chr: $chr\n"; exit;
    }
    my $s = $clm[1];
    my $e = $clm[2];
    my $summit = ($s + $e)/2;
    my $bin = int($summit/$width);
    if($lensum){$array{$chr}[$bin] += $e - $s;}
    else{$array{$chr}[$bin]++;}
}
close (ListFile);

foreach $chr (keys(%chrlen)){
    open(OUT, ">$outputfile.$width.$chr.xls") ||die "error: can't open $outputfile.$width.chr$chr.xls.\n";
    my $binnum = int($chrlen{$chr}/$width) +1;
    for($i=0;$i<$binnum;$i++){
	printf OUT "%d\t%d\n", $i*$width, $array{$chr}[$i];
    }
    close(OUT);
}

