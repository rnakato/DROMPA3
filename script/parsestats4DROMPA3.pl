#!/usr/bin/perl -w

$file="";
if($ARGV[0] =~ /(.+)\.([0-9]+)\.xls/){
    $file=$1;
}
if($file =~ /(.+)\/(.+)/){
    $file=$2;
}

$pcrbiasthre="";
$FRiP="";
$gcsummit=""; 
open(File, $ARGV[0]) ||die "error: can't open $ARGV[0].\n";
while(<File>){
    next if($_ eq "\n");
    chomp;
    if($_ =~ /Redundancy threshold: (.+)/){
	$pcrbiasthre = $1;
    }elsif($_ =~ /Library complexity: (.+) \((.+)\)/){
	$tested_complexity = $1;
	$tested_reads = $2;
    }elsif($_ =~ /FRiP: (.+)/){
	$FRiP = $1;
    }elsif($_ =~ /GC summit: (.+)/){
	$gcsummit = $1;
    }elsif($_ =~ /Whole genome/){
	chomp;
	my @clm= split(/\t/, $_);
	$total_reads = $clm[4];
	$plus = $clm[5];
	$minus = $clm[6];
	$total_remained = $clm[8];
	$total_filtered = $clm[11];
	if($gcsummit eq ""){
	    $depth = $clm[14];
	    $total_gc_base = $clm[17];
	}else{
	    $depth = $clm[17];
	    $total_gc_base = $clm[20];
	}
    }
}
close (File);

print STDERR "Sample\tMapped reads\t + strand\t - strand\tRedundancy threshold\tNonredundant\tRedundant\tComplexity for10M\tRead depth\tFRiP\tGenome coverage\tTested_reads\tGC summit\n";

print "$file\t$total_reads\t$plus\t$minus\t$pcrbiasthre\t$total_remained\t$total_filtered\t$tested_complexity\t$depth\t$FRiP\t$total_gc_base\t$tested_reads\t$gcsummit\n";
