# Author : Shujia Huang
# Date   : 2014-11-05
# Modify : 
#	2014-11-05 15:38:57 COPY from ~/Bin/MyPipe/VariantDetect/bin/Genotyping/src/RmSVdupFromVCF.pl
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

die qq/
Usage : perl $0 <command> [<arguments>] \n
Command:
    summary    Summary report for SV. Including SV number 
               and svsize distribution
    duplidist  Distribution of duplication for SV

/ if @ARGV < 1;

my $command = shift @ARGV;
my %func    = ('summary'   => \&SV_SummaryReport,
               'duplidist' => \&SV_DuplicDist);
die "Unknown command $command\n" if not exists $func{$command};
&{$func{$command}};

print STDERR "\n*************** Processing $command DONE ***************\n";
#################################################
sub SV_DuplicDist {

    my $vcffile;
    my $distance = 20; # The distance of variant nearby for duplication
    my $filter   = "ALL";
    GetOptions("v=s" => \$vcffile, "d=i" => \$distance, "f=s" => \$filter);

    print STDERR "\nCommand Parameter: 
          perl $0 duplidist 
               -v $vcffile 
               -d $distance
               -f $filter\n\n";

    my $fh;
    my @position;
    my %pos;   # %pos for checking the sorted status of vcf, 
    my $mi;    # $mi record the median index of positions
    my %dupliDist; # Distribution of duplication variants
    my $dnum;
    my ($n, $total) = (0, 0);
    my %filterStatistic = ("PASS" => 0, 
                           "PASS_MultiAllelic" => 0, 
                           "FALSE" => 0);
    open($fh, ($vcffile =~ /\.gz$/) ? "gzip -dc $vcffile |" : $vcffile) or
        die "Cannot open file $vcffile : $!\n";
    while (<$fh>) {

        next if /^#/;
        chomp;

        my @col = split;
        ++$n;

        next if uc($filter) ne 'ALL' and uc($filter) ne uc($col[6]);
        next if AsmvarVCFtools::IsNoGenotype(\@col[9..$#col]);
        print STDERR "[INFO] Loading $n lines\n" if $n % 100000 == 0;

        ++$total;
        ++$filterStatistic{$col[6]};
        if (AsmvarVCFtools::IsSpInfo('NEGATIVE_TRAIN_SITE', \$col[7])) {
            ++$filterStatistic{$col[6]."->NEGATIVE_TRAIN_SITE"};
        }

        my @mult = split /,/, $col[4];
        ++$filterStatistic{PASS_MultiAllelic} 
            if @mult > 1 and $col[6] eq 'PASS';

        die "[ERROR] VCF file need to be sorted by position\n" 
            if exists $pos{$col[0]} and $pos{$col[0]} > $col[1];
        $pos{$col[0]} = $col[1];

        if (@position > 0) {

            $dnum = @position;
            while (@position and 
                   $mi < @position and 
                   $position[$mi][0] eq $col[0] and 
                   $position[$mi][1] + $distance < $col[1]) {

                while ($mi >= 0 and 
                       $position[$mi][1] > $position[0][1] + $distance) {
                    shift @position;
                    $dnum = @position;
                    --$mi;
                }
#print join "\t", "** $mi dnum=$dnum", $position[$mi][1], (join ",", (map{$_->[1]} @position)), "\n";
                ++$dupliDist{$dnum};

                if (@position > 1) {
                    ++$mi;
                } else {
                    shift @position;
                }
            }

            if (@position > 0 and $position[0][0] ne $col[0]) {
#print join "\t", "** $mi dnum=$dnum", $position[$mi][1], (join ",", (map{$_->[1]} @position)), "\n";

                $dnum              = @position;
                $dupliDist{$dnum} += $dnum;
                @position          = ();
                $mi                = 0;
            }
        }
        $mi = 0 if not defined $mi or $mi < 0;
        push @position, [@col[0,1]]; # REF_ID and Position
    }
    close $fh;

    if (@position > 0) {

        while ($mi >= 0 and 
               $position[$mi][1] > $position[0][1] + $distance) {
            shift @position;
            $dnum = @position;
            --$mi;
        }
        $dnum              = @position;
        $dupliDist{$dnum} += $dnum;
    }

    print "\n** Summary **\n\n";
    _PrintVarStatistciSummary($n, $total, %filterStatistic);
    print "\n-- SV duplication spectrum for '$filter' variants --\n\n";
    print "#Duplication_number\tNumber\tRatio\n";
    for my $n (sort {$a <=> $b} keys %dupliDist) {
        my $r = sprintf "%.3f", $dupliDist{$n} / $total;
        print join "\t", $n, $dupliDist{$n}, $r, "\n";
    }

    return;
}

sub SV_SummaryReport {

# Calculate the SV number and SV size distribution
    use File::Basename qw/dirname/;
    use lib dirname($0)."/../lib";
    use AsmvarVCFtools;
    use AsmvarCommon;

    my ($vcffile, $isgatkvcf);
    my $filter = 'ALL';
    GetOptions(

        "v=s"   => \$vcffile,
        "f=s"   => \$filter,
        "g"     => \$isgatkvcf,
    );

    die qq/
Version: 0.0.1 (2014-11-05)
Author : Shujia Huang

    Last Modify: 2014-11-18  Update many things and debug

    Usage: perl summary $0 [Options] -v [vcfInfile] > output.vcf

    Options:

          -v  [str]  Variants file. [Reguire]
          -f  [str]  Specific FILTER. e.g: 'PASS'.  [ALL]
          -g         Input gatk vcf. [NULL]

\n/ if not defined $vcffile;

    print STDERR "\nCommand Parameter: 
          perl $0 summary 
               -v $vcffile 
               -f $filter\n\n";

    my ($n, $total) = (0, 0);
    my %filterStatistic = ("PASS" => 0, 
                           "PASS_MultiAllelic" => 0, 
                           "FALSE" => 0);

    my (%col2sam, %sizeSpectrum, %numSpectrum);
    my $fh; 
    open($fh, ($vcffile =~ /\.gz$/) ? "gzip -dc $vcffile |" : $vcffile) or 
        die "Cannot open file $vcffile : $!\n";
    while (<$fh>) {

        chomp;
        my @col = split;
        if (/^#CHROM/) {
            # Get Sample ID and build a hash to point to sample
            for (my $i = 9; $i < @col; ++$i) {
                $col2sam{$i-9} = $col[$i];
            }
        }
        next if /^#/;

        ++$n;
        next if AsmvarVCFtools::IsNoGenotype(\@col[9..$#col]);
        print STDERR "[INFO] Loading $n lines\n" if $n % 100000 == 0;

        my %format; 
        my @f = split /:/, $col[8];
        for (my $i = 0; $i < @f; ++$i) { $format{$f[$i]} = $i; }
        next if not $isgatkvcf and not exists $format{QR}; # May be INTERGAP 

        ++$total;
        ++$filterStatistic{$col[6]};
        ++$filterStatistic{$col[6]."_NEGATIVE_TRAIN_SITE"} 
            if AsmvarVCFtools::IsSpInfo('NEGATIVE_TRAIN_SITE', \$col[7]);

        my @mult = split /,/, $col[4];
        ++$filterStatistic{PASS_MultiAllelic} if @mult > 1 and $col[6] eq 'PASS';

        # Record information for summary output
        _SummarySV($isgatkvcf, # GATK VCF or not!
                   \%numSpectrum, 
                   \%sizeSpectrum,
                   @col[3,4], 
                   \%col2sam,
                   $format{VS}, 
                   $format{VT}, 
                   $format{QR}, 
                   @col[9..$#col]) if (uc($filter) eq 'ALL') or 
                                      (uc($col[6]) eq uc($filter));
    }
    close $fh;
#die "";

    print "\n** Summary **\n\n";
    _PrintVarStatistciSummary($n, $total, %filterStatistic);

    print "\n\n-- SV number spectrum for '$filter' variants --\n\n";
    _PrintNumSpectrum(\%numSpectrum);

    print "\n\n-- Size spectrum for '$filter' variants --\n";
    _PrintSizeSpectrum(\%sizeSpectrum);

    return;
}

sub _SummarySV {
    # Calculate the number and length in different SV types for each variant
    my ($isgatkvcf,
        $numSpectrum, 
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
#my ($ii, $t, $s) = AsmvarVCFtools::GetSVforPop($refseq, $altseq, \@samples);
#my () = AsmvarVCFtools::RecalcuSVBreakpoint(1, $refseq, $seq[$ii], $t);

    my %svstat;
    my $isempty = 1;
    for (my $i = 0; $i < @samples; ++$i) {

        my $sampleId = $$col2sam{$i};
        my @f = split /:/, $samples[$i];

        # Did not have genotype
        # If the sample could be genotype here 
        # than we'd better treat it get this SV
        next if $f[0] eq './.' or $f[0] eq '0/0';

        #Get the ALT sequence index
        my $ai = AsmvarVCFtools::GetAltIdxByGTforSample($f[0]);
        my ($svtype, $svsize);
        if ($isgatkvcf) {

            ($svtype, $svsize) = AsmvarVCFtools::GetGATKSVtypeAndSizeForSample(
                    $seq[0],    # Ref-sequence
                    $seq[$ai]); # Alt-sequence
        } else {
            
            ($svtype, $svsize) = AsmvarVCFtools::GetSVtypeAndSizeForSample(
                  $seq[0],     # Ref-sequence
                  $seq[$ai],   # Alt-sequence
                  $f[$vsIndex],# Init svsize 
                  (split /#/, $f[$vtIndex])[0]); # Split '#',in case of 'TRANS'
        }
 
        _SetValueToSummary(\$$numSpectrum{$sampleId}{$svtype}, $svsize);

        # Don't include 'REF_OR_SNP' when calculate total.
        if ($svtype !~ /REF_OR_SNP/) {

            _SetValueToSummary(\$$numSpectrum{$sampleId}{'0.Total(NON_SNP_VAR)'}, 
                               $svsize);

            # Calculate size spectrum
            my $bin = AsmvarCommon::SizeBinSp($svsize);
            #my $bin = AsmvarCommon::SizeBin($svsize, 10);
            $$sizeSpectrum{$sampleId}{$svtype}{$bin} ++; # Just for Variant
            $$sizeSpectrum{$sampleId}{'0.ALLSV'}{$bin} ++; # For all Variant
        }

        # Use for getting SVforAll(population) in this position
        $svstat{$svtype}->[0] ++;
        $svstat{$svtype}->[1] = [$svtype, $svsize];

        $isempty = 0;
    }
    return if $isempty;

    my ($totalsvtype, $totalsvsize) = 
        AsmvarVCFtools::GetSVforAllPerVariantLine(\%svstat);
    
    _SetValueToSummary(\$$numSpectrum{'~Population'}{$totalsvtype}, 
                        $totalsvsize);
    # Don't include 'REF_OR_SNP' when calculate total.
    if ($totalsvtype !~ /REF_OR_SNP/) {

        _SetValueToSummary(\$$numSpectrum{'~Population'}{'0.Total(NON_SNP_VAR)'}, 
                            $totalsvsize);
        # Calculate size spectrum
        my $bin = AsmvarCommon::SizeBinSp($totalsvsize);
        #my $bin = AsmvarCommon::SizeBin($totalsvsize, 10);
        $$sizeSpectrum{'~Population'}{$totalsvtype}{$bin} ++; #Just for Variant
        $$sizeSpectrum{'~Population'}{'0.ALLSV'}{$bin} ++; # For all Variant
    }

    return;
}

sub _SetValueToSummary {
# Input: an array reference => [] and svsize

    my ($record, $svsize) = @_;

    # Check inited or not
    if (not defined $$record->[0]) {
        # [num, all_size, min_size, max_size]
        $$record = [0, 0, $svsize, $svsize];
    }

    $$record->[0] ++;         # sv number
    $$record->[1] += $svsize; # add all the svsize up
    $$record->[2]  = $svsize if $svsize < $$record->[2]; # MIN
    $$record->[3]  = $svsize if $svsize > $$record->[3]; # MAX
}

sub _PrintVarStatistciSummary {

    my ($allvaraint, $usefulvariant, %filterStatistic) = @_;

    my $tr = sprintf "%.3f", $usefulvariant / $allvaraint;
    print "** The whole set of variants in VCF: $allvaraint\n";
    print "** The number of useful variants   : $usefulvariant ($tr)\n";

    for my $k (sort {$a cmp $b} keys %filterStatistic) {

        my $r = sprintf "%.3f", $filterStatistic{$k} / $usefulvariant;
        print "** $k variants: $filterStatistic{$k} ($r) **\n";
    }
}

sub _PrintNumSpectrum {

    my ($svNumSpectrum) = @_;

    my %allsvtype;
    for my $s (keys %$svNumSpectrum) {
        $allsvtype{$_} = 1 for keys %{$$svNumSpectrum{$s}};
    }

    my $type;
    for my $svtype (sort {$a cmp $b} keys %allsvtype) {
        # Do not include the number in front of the $svtype in output header.
        # The format is always be: 'Number.Type'
        $svtype  = (split /\./, $svtype)[-1];
        $type   .= "\t$svtype-NUM\t$svtype-LEN\t$svtype-MIN".
                    "\t$svtype-MAX\t$svtype-MEAN";
    }

    print "#SampleID\t$type\n";
    for my $sampleId (sort {$a cmp $b} keys %$svNumSpectrum) {

        my @outinfo;
        my $mean;
        for my $svtype (sort {$a cmp $b} keys %allsvtype) {

            if (exists $$svNumSpectrum{$sampleId}{$svtype}) {
                $mean = $$svNumSpectrum{$sampleId}{$svtype}->[1] / 
                        $$svNumSpectrum{$sampleId}{$svtype}->[0];
                push @outinfo, (join "\t", 
                                @{$$svNumSpectrum{$sampleId}{$svtype}},
                                sprintf("%.2f", $mean));
            } else {
                push @outinfo, (join "\t", (0) x 5);
            }
        }
        print join "\t", $sampleId, @outinfo, "\n";
    }
}

sub _PrintSizeSpectrum {

    my ($sizeSpectrum) = @_;

    my %size;
    my %type;
    for my $s (keys %$sizeSpectrum) {

        $type{$_} = 1 for keys %{$$sizeSpectrum{$s}}; # all the svtype
        for my $k (keys %type) {
            # Format : 'start-end'
            $size{$_} = 1 for (keys %{$$sizeSpectrum{$s}{$k}});
        }
    }

    my @type = sort {$a cmp $b} keys %type;
    my @size = sort {my $da = (split /\-/, $a)[0]; 
                     my $db = (split /\-/, $b)[0]; 
                     $da <=> $db;
                    } keys %size;

    for my $k (@type) {

        print "\n** Size Spectrum for '", (split /\./, $k)[-1], "' **\n";
        print join "\t", "#SampleID", @size, "\n";

        for my $sampleId (sort {$a cmp $b} keys %$sizeSpectrum) {

            my @outinfo;
            for my $s (@size) {

                $$sizeSpectrum{$sampleId}{$k}{$s} = 0 
                    if not exists $$sizeSpectrum{$sampleId}{$k}{$s};
                push @outinfo, $$sizeSpectrum{$sampleId}{$k}{$s};
            }
            print join "\t", $sampleId, @outinfo, "\n";
        }
    }
}
