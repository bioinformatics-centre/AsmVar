#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

die qq/\nperl $0 [InVcffile] > OutputVcf\n\n/ if @ARGV == 0;
my ($vcfInfile) = @ARGV;

my ($n, $m) = (0, 0);
open I, ($vcfInfile =~ /\.gz$/ ? "gzip -dc $vcfInfile |": $vcfInfile) or die "$!";
while (<I>) {

	if (/^#/) { print; next; }
	chomp;
	my @tmp = split;
	$tmp[3] =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
	$tmp[4] =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
	++$m;
	next if uc($tmp[3]) eq uc($tmp[4]);
	++$n;
	print join "\t", @tmp; print "\n";
}
close I;

my $f = $m - $n;
print STDERR "[WARNING] You filter $f lines after translating the alleles.\n" if $f;
print STDERR "******** ALL DONE *******\n";



