#!/usr/bin/perl -w

if($#ARGV != 1){
    print " GCcount.pl <chromosome.fa> <window size>\n";
    exit;
}

$fasta=$ARGV[0];
$window=$ARGV[1];

open(ListFile, $fasta) ||die "cannot open $fasta";
while(<ListFile>){
    next if($_ =~ />/);
    chomp;
    $seq .= $_;
}
close (ListFile);

$len = length($seq);
$nwin = int($len/$window) + 1;
for($i=0;$i<$nwin;$i++){
    $subseq = substr($seq, $i*$window, $window);
    $gc_count = ($subseq =~ tr/cgCG/cgCG/);
    $acgt_count = ($subseq =~ tr/acgtACGT/acgtACGT/);
    if($acgt_count >= 20){
	$GCcontent = $gc_count/$acgt_count*100;
    }else{
	$GCcontent = 0;
    }
    printf("%d\t%.2f\n", $i*$window, $GCcontent);
}
