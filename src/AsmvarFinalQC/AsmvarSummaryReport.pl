# Author : Shujia Huang
# Date   : 2014-11-05
# Modify : 
#	2014-11-05 15:38:57 COPY from ~/Bin/MyPipe/VariantDetect/bin/Genotyping/src/RmSVdupFromVCF.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util;
use File::Basename qw/dirname/;

use lib dirname($0)."/../lib";
use AsmvarVCFtools;
use AsmvarCommon;

my ($vcffile);
my $qualityThd = 2; # No use now
my $filter     = 'ALL';
GetOptions(

    "v=s"   => \$vcffile,
    "q=i"   => \$qualityThd,
    "f=s"   => \$filter,
);
Usage() if (!$vcffile);
print STDERR "\nCommand Parameter: perl $0 -v $vcffile -q $qualityThd -f $filter\n\n";

my (%sample, %col2sample, @info);
SV_SummaryReport($vcffile, $filter, $qualityThd);

print STDERR "\n********************** ALL DONE ********************\n";

#################################################
sub SV_SummaryReport {

    my ($fn, $filter, $qualityThd) = @_;
    my ($total, $pass, $duplic, $false, $lowQ) = (0,0,0,0,0);
    my $passMultiAllelic = 0;

    my (%col2sam, %sizeSpectrum, %summary, %allsvtype);
    my $fh; 
	my $n = 0;
    open($fh, ($fn =~ /\.gz$/) ? "gzip -dc $fn |" : $fn) || die "Cannot open file $fn : $!\n";
    while (<$fh>) {

        chomp;
        my @col = split;
        if (/^#CHROM/) {
            for (my $i = 9; $i < @col; ++$i) {
                $col2sam{$i-9} = $col[$i];
            }
        }
        next if /^#/;
        ++$n;
        next if AsmvarVCFtools::IsNoGenotype(@col[9..$#col]);
        print STDERR "[INFO] Loading $n lines\n" if $n % 100000 == 0;

        my %format; 
        my @f = split /:/, $col[8];
        for (my $i = 0; $i < @f; ++$i) { $format{$f[$i]} = $i; }
        next if not exists $format{QR}; # May be INTERGAP 

        ++$total;
        ++$false if $col[6] eq 'FALSE';
        ++$pass  if $col[6] eq 'PASS';

        my @mult = split /,/, $col[4];
        ++$passMultiAllelic if @mult > 1 and $col[6] eq 'PASS';

        # Record information for summary output
        Summary(\%summary, 
                \%allsvtype, 
                \%sizeSpectrum,
                @col[3,4], 
                \%col2sam,
                $format{VS}, 
                $format{VT}, 
                $format{QR}, 
                @col[9..$#col]) if (uc($filter) eq 'ALL') or 
                                   (uc($col[6]) eq uc($filter));
    }

    my $rf = sprintf "%.3f", $false/$total;
    my $rp = sprintf "%.3f", $pass/$total;
    my $tr = sprintf "%.3f", $total/$n;
    my $mr = sprintf "%.3f", $passMultiAllelic/$total;
    print "\n** Summary **\n\n";
    print "** The whole set of variants in VCF: $n\n";
    print "** The number of useful variants   : $total ($tr)\n";
    print "** PASS variants    : $pass ($rp)\n";
    print "** PASS MultiAllelic: $passMultiAllelic ($mr)\n";
    print "** FALSE variants   : $false ($rf)\n\n";

    print "-- Just for '$filter' variants --\n";
    OutputSummary(\%allsvtype, \%summary);

    print "\n\n-- Size Spectrum for all '$filter' variants --\n";
    OutputSpectrum(\%sizeSpectrum);

    return;
}

sub Summary {
# Calculate the number and length in different SV types for each variant
    my ($summary, 
        $allsvtype, 
        $sizeSpectrum,
        $refseq, 
        $altseq, 
        $col2sam, 
        $vsIndex, 
        $vtIndex, 
        $qrIndex, 
        @samples) = @_;

    my @seq = ($refseq); # First element is REF: [0]=>REF
    push @seq, $_ for (split /,/, $altseq);

    my %svstat;
    my $isempty = 1;
    for (my $i = 0; $i < @samples; ++$i) {

        my $sampleId = $$col2sam{$i};
        my @f = split /:/, $samples[$i];

        # Did not have genotype
        # If the sample could be genotype here 
        # than we'd better treat it get this SV
        next if $f[0] eq './.' or $f[0] eq '0/0';
#next if $f[0] eq './.';
#next if (@f < $qrIndex + 1) or ($f[$qrIndex] eq '.'); # Sample region. Too strict

        #Get the ALT sequence index
        my $ai = AsmvarVCFtools::GetAltIdxByGTforSample($f[0]); # Get the ALT sequence index
        my ($svtype, $svsize) = AsmvarVCFtools::GetSVtypeAndSizeForSample(
                  $seq[0],     # Ref-sequence
                  $seq[$ai],   # Alt-sequence
                  $f[$vsIndex],# Init svsize 
                  (split /#/, $f[$vtIndex])[0]); # Split '#',in case of 'TRANS'
 
        SetValueToSummary(\$$summary{$sampleId}{$svtype}, $svsize);
        if ($svtype !~ /REF_OR_SNP/) { # Don't include such type when calculate total.

            SetValueToSummary(\$$summary{$sampleId}{'0.Total'}, $svsize);

            # Calculate size spectrum
            my $bin = AsmvarCommon::SizeBinSp($svsize);
            $$sizeSpectrum{$sampleId}{$bin} ++; # Just for Variant
        }

        # Use for getting SVforAll(population) in this position
        $svstat{$svtype}->[0] ++;
        $svstat{$svtype}->[1] = [$svtype, $svsize];

        # Record all the svtype using for output
        $$allsvtype{$svtype} = 1;
        $isempty = 0;
    }
    $$allsvtype{'0.Total'} = 1;
    return if $isempty;

    my ($totalsvtype, $totalsvsize) = 
        AsmvarVCFtools::GetSVforAllPerVariantLine(\%svstat);
    
    SetValueToSummary(\$$summary{'~Population'}{$totalsvtype}, $totalsvsize);
    if ($totalsvtype !~ /REF_OR_SNP/) { # Don't include such type when calculate total.

        SetValueToSummary(\$$summary{'~Population'}{'0.Total'}, $totalsvsize);
        # Calculate size spectrum
        my $bin = AsmvarCommon::SizeBinSp($totalsvsize);
        $$sizeSpectrum{'~Population'}{$bin} ++; # Just for Variant
    }

    return;
}

sub SetValueToSummary {
# Input: an array reference => [] and svsize

    my ($record, $svsize) = @_;

    # Check have been inited or not
    if (not defined $$record->[0]) {
    # [num, all_size, min_size, max_size]
        $$record = [0, 0, $svsize, $svsize];
    }

    $$record->[0] ++;         # sv number
    $$record->[1] += $svsize; # add all the svsize up
    $$record->[2]  = $svsize if $svsize < $$record->[2]; # MIN
    $$record->[3]  = $svsize if $svsize > $$record->[3]; # MAX
}
##################

sub OutputSummary {

    my ($allsvtype, $summaryInfo) = @_;

    my $header = "#SampleID";
    for my $svtype (sort {$a cmp $b} keys %$allsvtype) {
        # Do not include the number in front of the $svtype in output header.
        # The format is always be: 'Number.Type'
        $svtype  = (split /\./, $svtype)[-1];
        $header .= "\t$svtype-NUM\t$svtype-LEN\t$svtype-MIN\t$svtype-MAX\t$svtype-MEAN";
    }

    print "$header\n";
    for my $sampleId (sort {$a cmp $b} keys %$summaryInfo) {

        my @outinfo;
        my $mean;
        for my $svtype (sort {$a cmp $b} keys %$allsvtype) {
            if (exists $$summaryInfo{$sampleId}{$svtype}) {
                $mean = $$summaryInfo{$sampleId}{$svtype}->[1] / 
                        $$summaryInfo{$sampleId}{$svtype}->[0];
                push @outinfo, (join "\t", 
                                @{$$summaryInfo{$sampleId}{$svtype}},
                                sprintf("%.2f", $mean));
            } else {
                push @outinfo, (join "\t", (0) x 5);
            }
        }
        print join "\t", $sampleId, @outinfo, "\n";
    }
}

sub OutputSpectrum {
    my ($sizeSpectrum) = @_;

    my %size;
    for my $s (keys %$sizeSpectrum) {
        # Format : 'start-end'
        $size{$_} = 1 for (keys %{$$sizeSpectrum{$s}});
    }

    my @size = sort {my $da = (split /\-/, $a)[0]; 
                     my $db = (split /\-/, $b)[0]; 
                     $da <=> $db;
                    } keys %size;
    print join "\t", "#SampleID", @size, "\n";
    for my $sampleId (sort {$a cmp $b} keys %$sizeSpectrum) {

        my @outinfo;

        for my $s (@size) {

            $$sizeSpectrum{$sampleId}{$s} = 0 
                if not exists $$sizeSpectrum{$sampleId}{$s};

            push @outinfo, $$sizeSpectrum{$sampleId}{$s};
        }
        print join "\t", $sampleId, @outinfo, "\n";
    }
}

#########
sub Usage {

    print STDERR <<U;
Version: 0.0.1 (2014-11-05)
Author : Shujia Huang

        Last Modify: 2014-11-09  Update many things and debug

        Usage: perl $0 [Options] -v [vcfInfile] > output.vcf

        Options:

              -v  [str]  Variants file. [Reguire]
              -q  [int]  Threshold for Variant qsulity score. [$qualityThd]
              -f  [str]  Specific FILTER. e.g: 'PASS'.  [ALL]

U
     exit(0);
}










