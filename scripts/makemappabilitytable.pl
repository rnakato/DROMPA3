#!/usr/bin/perl -w

if($#ARGV != 1){
    print " makemappabilitytable.pl <genometable> <binaryfile prefix>\n";
    exit;
}

$gtfile = $ARGV[0];
$dir = $ARGV[1];

$num=0;
open(InputFile,$gtfile) ||die "error: can't open $gtfile.\n";
while($line = <InputFile>){
    next if($line eq "\n");
    chomp $line;
    @clm= split(/\t/, $line);
    $name[$num] = $clm[0];
    $len[$num] = $clm[1];
#    print "$num\t$name[$num]\t$len[$num]\n";
    $num++;
}
close (InputFile);

for($i=0; $i<$num; $i++){
    my $filename =  "${dir}_${name[$i]}_binary.txt";
#    print "$filename\n";
    open(InputFile,$filename) ||die "error: can't open $filename.\n";
    while($line = <InputFile>){
#	$count0 = ($line =~ tr/0/0/);
	$count1 = ($line =~ tr/1/1/);
    }
    my $r = $count1 / $len[$i];
#    print "$name[$i]\t$count1\t$len[$i]\t$r\n";
    print "$name[$i]\t$count1\n";
    close (InputFile);
}

