#!/usr/bin/perl

die "makeinteraction.pl <bowtie1> <bowtie2>\n" if($#ARGV !=1);

open(IN, $ARGV[0]) || die "cannot open $ARGV[0].";
while(<IN>) {
    next if($_ eq "\n");
    chomp;
    my @clm = split(/\t/, $_);
    $Hash{$clm[0]} = $_;
}
close IN;


open(IN, $ARGV[1]) || die "cannot open $ARGV[1].";
while(<IN>) {
    next if($_ eq "\n");
    chomp;
    my @clm = split(/\t/, $_);
    if(exists($Hash{$clm[0]})){
	my @clm1= split(/\t/, $Hash{$clm[0]});
	my $chr1 = $clm1[2];
	my $strand1 = $clm1[1];
	my $s1 = $clm1[3];
	my $e1 = $s1 + length($clm1[4]);
	my @clm2 = split(/\t/, $_);
	my $chr2 = $clm2[2];
	my $strand2 = $clm2[1];
	my $s2 = $clm2[3];
	my $e2 = $s2 + length($clm2[4]);
	if($strand1 ne $strand2){
	    print "$chr1\t$s1\t$e1\t$chr2\t$s2\t$e2\n";
	}
    }
}
close IN;

