# Author : Shujia Huang
# Date   : 2015-01-06
#!/usr/bin/perl
use strict;
use warnings;

my @vcfInfiles = @ARGV;
die "perl $0 [vcffiles] > outputvcf\n" if @ARGV < 1;

print STDERR "\nCommand Parameter:\nperl $0\n\t", 
             join ("\n\t", @vcfInfiles), "\n\n";

my %recorl;
for my $file (@vcfInfiles) {

    print STDERR "\n[INFO] Loading $file\n";
    my $n = 0;
    open I, ($file =~ /\.gz$/ ? "gzip -dc $file |": $file) or 
        die "Cannot open file $file: $!\n";
    while (<I>) {

        print STDERR "    -- loading $n lines\n" if $n % 100000 == 0;
        ++$n;

        next if /^#/;
        chomp;
        my @tmp = split;
        my $key = join ":", @tmp[0,1];
        $recorl{$key} = [@tmp[5, 7]] if not exists $recorl{$key};
        # Keep the highest VQ variant
        $recorl{$key} = [@tmp[5, 7]] if $recorl{$key}->[0] < $tmp[5]; 
    }
    close I;
    print STDERR "    -- Complete loading $n lines\n";
}

my %score = (0 => 0, 
             1 => 0, 
             2 => 0, 
             3 => 0, 
             4 => 0);
my %culprit;
my $total = 0;
open I, ($vcfInfiles[0] =~ /\.gz$/ ? "gzip -dc $vcfInfiles[0] |": $vcfInfiles[0])
    or die "Cannot open file $vcfInfiles[0]: $!\n";
while (<I>) {

    if (/^#/) { print; next;}
    chomp;

	++$total;
    my @tmp = split;
    my $key = join ":", @tmp[0,1];
    @tmp[5,7] = @{$recorl{$key}};

	++$score{0} if $tmp[5] >= 0;
	++$score{1} if $tmp[5] >= 1;
	++$score{2} if $tmp[5] >= 2;
	++$score{3} if $tmp[5] >= 3;
	++$score{4} if $tmp[5] >= 4;

	my ($culprit) = $tmp[7] =~ /;CU=[^;]+/;
	++$culprit{$culprit};

    print join "\t", @tmp; print "\n";
}

close I;

print STDERR "\n[Summmary] Here is the summary information:\n";
for my $k (sort {$a <=> $b} keys %score) {

	my $r = sprintf "%.2f", 100 * $score{$k} / $total;
	print STDERR "  ** Variant Site score >= $k: $score{$k}\t$r\n";
}

for my $k (sort {$a cmp $b} keys %culprit) {

	my $r = sprintf "%.2f", 100 * $culprit{$k} / $total;
	print STDERR "  ** Culprit by $k: $culprit{$k}\t$r\n";
}

print STDERR "\n*********************** ALL DONE ***********************\n";



