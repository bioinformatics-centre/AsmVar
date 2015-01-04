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

my ($vcffile);
my $qualityThd = 2;
my $distance   = 20;
GetOptions(

    "v=s"   => \$vcffile,
    "d=i"	=> \$distance,
    "q=i"   => \$qualityThd,
);
Usage() if (!$vcffile);
print STDERR "\nCommand Parameter: perl $0 -v $vcffile -d $distance -q $qualityThd\n\n";

my (%sample, %col2sample, @info);
LoadVarRegFromVcf($vcffile,    \@info);
RemoveOverlap    ($distance,   \@info); # The element which need to remove, now is recorded in the first element for each line '$$info[$i][0] = -1'!
Output           ($vcffile, $qualityThd, \@info);

print STDERR "\n********************** ALL DONE ********************\n";

#################################################

sub Output {

    my ($fn, $qualityThd, $info) = @_;
    my ($total, $pass, $duplic, $false, $lowQ) = (0,0,0,0,0);

    my (%col2sam, %summary, %allsvtype);
    my $fh; 
    my ($i, $n) = (0, 0);
    open($fh, ($fn =~ /\.gz$/) ? "gzip -dc $fn |" : $fn) || die "Cannot open file $fn : $!\n";
    while (<$fh>) {

        if (/^##fileformat/) {
            print;
            print "##FILTER=<ID=DUPLIC,Description=\"Duplication variants. Don't need to keep them!\">\n";
            print "##FILTER=<ID=FALSE,Description=\"unconfident variant\">\n";
            print "##FILTER=<ID=q$qualityThd,Description=\"Variant quality below $qualityThd\">\n";
            next;
        } elsif (/^#CHROM/) {
            print "##INFO=<ID=SPN,Number=1,Type=Integer,Description=\"The count of assambly which support this variant region\">\n";
            print; 
            chomp;
            my @t = split;
            for (my $i = 9; $i < @t; ++$i) {
                $col2sam{$i-9} = $t[$i];
            }
            next;
        } elsif ($_ =~ /^#/) {
            print; next;
        }
        chomp;

        ++$n;
##[$mark, $nr, $vq, $asmNum, $nummap, $svtype, $svsize, $tId, $tStart, $tEnd, @t]
        my $mark   = $$info[$i][0];
        my $bestNR = $$info[$i][1];
        my $vq     = $$info[$i][2];
        my $asmNum = $$info[$i][3];
        my $mapnum = $$info[$i][4];
        my $svtype = $$info[$i][5];
        my $svsize = $$info[$i][6];
        my $tId    = $$info[$i][7];
        my $tStart = $$info[$i][8];
        my $tEnd   = $$info[$i][9];

        die qq/[BUG] The output order is not right.\n/ if $n > $mapnum;
        next if $n != $mapnum;
        ++$i; # Increan $i if $i == $mapnum

        my @col = split;
        my %format; my @f = split /:/, $col[8];
        for (my $i = 0; $i < @f; ++$i) { $format{$f[$i]} = $i; }
        die qq/[BUG] 'QR' must in FORMAT Field!\n/ if (!exists $format{QR});

        $col[7] =~ s/;SPN=[^;]+//g;
        $col[7] .= ";SPN=$asmNum";

# For the $mark's detail:
# If the vq is all > qualityThd
# then 0: PASS
#      1: PASS But Duplication
#     -2: DUPLIC
#     -1: FALSE. Anyway
        my $dup = 0;
        if ($mark == -1) {
            $col[6] = 'FALSE';
        } elsif ($mark == -2 and $vq >= $qualityThd) {
            $col[6] = ($col[6] eq '.' || $col[6] eq 'PASS') ? 
                'DUPLIC': 'DUPLIC;'.$col[6];
            ++$duplic;
        } elsif ($mark == 1 and $vq >= $qualityThd) {

            $col[6] = 'PASS';
        } elsif ($vq >= $qualityThd) {

            $col[6] = 'PASS';
        } else {
        # Low Quality variant score
            $col[6] = "q$qualityThd";
        }

        ++$total;
        ++$false if $col[6] eq 'FALSE';
        ++$lowQ  if $col[6] eq "q$qualityThd";
        ++$pass  if $col[6] eq 'PASS';
        print join "\t", @col; print "\n";

        # Record information for summary output
#print STDERR join "\t", "[Before Summary]", @col[0..4], "\n";
        Summary(\%summary, \%allsvtype, @col[3,4], \%col2sam,
                $format{VS}, $format{VT}, $format{QR}, @col[9..$#col]) if $col[6] eq 'PASS';
    }

    my $rf = sprintf "%.3f", $false/$total;
    my $rd = sprintf "%.3f", $duplic/$total;
    my $rp = sprintf "%.3f", $pass/$total;
    my $rl = sprintf "%.3f", $lowQ/$total;
    my $tr = sprintf "%.3f", $total/$n;
    print STDERR "\n** Summary **\n\n";
    print STDERR "** The whole set of variants in VCF: $n\n";
    print STDERR "** The number of useful variants   : $total ($tr)\n";
    print STDERR "** PASS variants   : $pass ($rp)\n";
    print STDERR "** q$qualityThd variants     : $lowQ ($rl)\n";
    print STDERR "** FALSE variants  : $false ($rf)\n";
    print STDERR "** DUPLIC variants : $duplic ($rd)\n\n";

    OutputSummary(\%allsvtype, \%summary);

    return;
}

sub Summary {
# Calculate the number and length in different SV types for each variant
    my ($summary, $allsvtype, $refseq, $altseq, $col2sam, 
        $vsIndex, $vtIndex, $qrIndex, @samples) = @_;

#print STDERR join "\n", "[DEbug]", @samples, "\n";
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
        #next if (@f < $qrIndex + 1 or $f[$qrIndex] eq '.'); # Sample region. Too strict

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

    print STDERR "$header\n";
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
        print STDERR join "\t", $sampleId, @outinfo, "\n";
    }
}

sub LoadVarRegFromVcf {

    my ($fn, $info) = @_;

    my $fh; my $nummap = 0; 
    open($fh, ($fn =~ /\.gz$/) ? "gzip -dc $fn |" : $fn) || die "Cannot open file $fn : $!\n";
    while (<$fh>) {

        next if ($_ =~ /^#/);
        chomp;

        my @t = split;
        die "[ERROR] VCF FORMAT ERROR. The 'GT' field MUST BE present and always appear as the first field in $fn.\n" if ($t[8] !~ /^GT:/);
        ++$nummap; # Record the line order. It'll be use in Output
        next if AsmvarVCFtools::IsNoGenotype(\@t);

        my %format;
        my @f = split /:/, $t[8];
        for (my $i = 0; $i < @f; ++$i) { $format{$f[$i]} = $i; }
        next if not exists $format{QR};

        die "[ERROR] VCF FORMAT ERROR. Missing the 'TR' information in the 'FORMAT'fields in your $fn\n" if (!exists $format{TR});
        die "[ERROR] VCF FORMAT ERROR. Missing the 'NR' information in the 'FORMAT'fields in your $fn\n" if (!exists $format{NR});
        die "[ERROR] VCF FORMAT ERROR. Missing the 'VS' information in the 'FORMAT'fields in your $fn\n" if (!exists $format{VS});
        die "[ERROR] VCF FORMAT ERROR. Missing the 'VT' information in the 'FORMAT'fields in your $fn\n" if (!exists $format{VT});

        my ($nr, $svtype, $svsize, $tId, $tStart, $tEnd, $asmNum) = 
            FindBestInSingleVariant($format{TR}, $format{QR}, $format{NR},
                    $format{VT}, $format{VS}, @t[9..$#t]);

        if ($tEnd - $tStart < 10) {
            $tStart -= 5; $tStart = 1 if ($tStart < 0);
            $tEnd   += 5;
        }
        $t[7] =~ s/;SPN=[^;]+//g;
        $t[7] .= ";SPN=$asmNum";

        my $inb= AsmvarVCFtools::GetDataInSpInfo('InbCoeff', \$t[7]); # Get InbCoeff
        # '-1' is Mark for delete if too much 'N', or inb is not in (-.7, .7)
        my $ma = ($nr > 0.5 or $inb <= -0.7 or $inb >= 0.7) ? -1: 0;
        my $vq = AsmvarVCFtools::GetDataInSpInfo('VQ', \$t[7]); # Get variant score
        push(@$info, [$ma, $nr, $vq, $asmNum, $nummap, $svtype, 
                      $svsize, $tId, $tStart, $tEnd, @t[0,1]]);
    }
    close($fh);

    print STDERR "*** Complete loading the vcf file $fn **\n** Start removing duplicate and calculating the best region **\n";
}

sub FindBestInSingleVariant {

    my ($trIndex, $qrIndex, $nrIndex, $svtIndex, $svsIndex, @sample) = @_;

    my ($nr, $tId, $tStart, $tEnd, $svtype, $svsize);
    my %hash;
    my $asmNum = 0;
    for my $sam (@sample) {

        my @f = split /:/, $sam; #./.:52,2:58,2:T:0.02:F:11,16,0,0,0,18,1:11,35:scaffold71302-155-202:0,0,0,0,0,0,0:0,0:11-151472-151472:3:INS
        next if uc $f[$svtIndex] eq 'REFCALL';

        # ignore the FORMAT just have './.'
        ++$asmNum if ((@f >= $qrIndex + 1) and ($f[$qrIndex] ne '.'));
        next if $f[0] eq './.' or $f[$qrIndex] eq '.'; 

        $svtype = uc $f[$svtIndex];  # SV Type
        $svsize = abs $f[$svsIndex]; # SV Size
        $svtype = 'INDEL' if ($svtype =~ m/INS|DEL/);
        $svtype = 'TRANS' if ($svtype =~ m/TRANS/);
        $nr     = $f[$nrIndex];
        next if $nr eq '.';

        ($tId, $tStart, $tEnd) = (split /[=\-]/, $f[$trIndex]);
        push @{$hash{$nr}{$svtype}{$svsize}}, [$f[$nrIndex], $f[$svtIndex], $svsize, $tId, $tStart, $tEnd];
    }

    if ((keys %hash) == 0) {
        my $sI;
        for (my $sam = 0; $sam < @sample; ++$sam) {
            my @f = split /:/, $sample[$sam];
            next if (@f < $qrIndex + 1 or $f[$qrIndex] eq '.');
            $sI = $sam;
        }
        my @f = split /:/, $sample[$sI];
        return ($f[$nrIndex], $f[$svtIndex], $f[$svsIndex], (split /[=\-]/, $f[$trIndex]), $asmNum);
    }

    my $bk = (sort{$a<=>$b} keys %hash)[0]; # Best NR Key
    my $bt; # Best SV-Type
    if (exists $hash{$bk}{INDEL}     ) {
        $bt = 'INDEL';
    } elsif (exists $hash{$bk}{TRANS}) {
        $bt = 'TRANS';
    } elsif (exists $hash{$bk}{INV}  ) {
        $bt = 'INV';
    } else {
        $bt = (keys %{$hash{$bk}})[0];
    }
    my $bs = (sort{$a<=>$b} %{$hash{$bk}{$bt}})[0];

# I can do here after genotyping process, but Donot random select 
# befroe that or we'll catch format error in SVGenotype
    my $acbI = int(rand(scalar(@{$hash{$bk}{$bt}{$bs}}))); # my $acbI = 0;

    return (@{$hash{$bk}{$bt}{$bs}->[$acbI]}, $asmNum);
}

sub RemoveOverlap { # Find the best region from nerby positions(regions) by vcf format

    print STDERR "** Start removing duplicate and calculating the best region **\n";
    my ($distance_delta, $info) = @_;
    my (%prePos, @index, $id, $start, $end, @data);

    my ($rI, $sI, $eI) = (7, 8, 9);
    @$info = sort { 
        my $da = $a->[$rI]; 
        my $db = $b->[$rI]; 
        if($da eq $db) { 
            $a->[$sI] <=> $b->[$sI]; 
        } else { 
            $a->[$rI] cmp $b->[$rI]; 
        }
    } @$info;

    print STDERR "\n** Starting to RemoveOverlap **\n\n";
    for (my $i = 0; $i < @$info; ++$i) {

        next if $$info[$i][0] < 0;

        die "[ERROR]Your region start > end (@{$$info[$i]}). ".
            "This is not allow in RemoveOverlap() function.\n" 
            if ($$info[$i][$eI] < $$info[$i][$sI]);

        if (!exists $prePos{$$info[$i][$rI]}) {

            Rm($info, @index) if (@index > 1);
            @index = ($i);
            ($id, $start, $end) = @{$$info[$i]}[$rI..$eI];
        } else {

            die "[ERROR]Your array hasn't been sorted.\t".
                "$prePos{$$info[$i][$rI]} > ($$info[$i][$sI])\n@{$$info[$i]}\n"
                if ($prePos{$$info[$i][$rI]} > $$info[$i][$sI]);

            if ($end + $distance_delta >= $$info[$i][$sI]) { # [7,8,9]

                if ($$info[$i][1] < $$info[$index[0]][1]) {

                    $$info[$_][0] = -1 for (@index); # Mark for delete!!
                        @index = ($i);
                    ($id, $start, $end) = @{$$info[$i]}[$rI..$eI];
                } elsif ($$info[$i][1] == $$info[$index[0]][1]) {

                    push @index, $i;
                    $end = $$info[$i][$eI] if ($$info[$i][$eI] > $end); # Keep the preview and shift to next
                } else {

                    $$info[$i][0] = -1; # Mark for delete!!
                }
            } else {

                Rm($info, @index) if (@index > 1);
                @index = ($i);
                ($id, $start, $end) = @{$$info[$i]}[$rI..$eI];
            }
        }

        $prePos{$info[$i][$rI]} = $$info[$i][$sI];
    }

    Rm($info, @index) if (@index > 1);

# Set the output order back to the original vcfInfile.
    @$info = sort {$a->[4] <=> $b->[4]} @$info;
    print STDERR "\n** RemoveOverlap Done **\n";
    return;
}

sub Rm {

    my ($info, @index) = @_;
    my %hash;
    my ($rI, $sI, $eI) = (7, 8, 9);
    for my $j (@index) {
        my $nr = $$info[$j][1];
        push @{ $hash{$nr} }, $j;
    }
    my $bk = (sort{ $a<=>$b } keys %hash)[0];      # Best NR Key
        my $kn = scalar(keys %hash);

##[$mark, $nr, $vq, $asmNum, $nummap, $svtype, $svsize, $tId, $tStart, $tEnd, @t]
    my @leftIndex;
    for my $j (sort { $a <=> $b } @{$hash{$bk}}) {

        next if $$info[$j][0] == -1;
        if (@leftIndex > 0) {
# keep the highest quality variants here
            if ($$info[$j][2] > $$info[$leftIndex[0]][2]) {
                $$info[$_][0] = ($$info[$_][1] == 0.0) ? -2 : -1 for (@leftIndex); # Mark for delete!!
                    @leftIndex = ($j);
            } elsif ($$info[$j][2] == $$info[$leftIndex[0]][2]) {
                push @leftIndex, $j;
            } else {
                $$info[$j][0] = ($$info[$j][1] == 0.0) ? -2 : -1; # Mark for delete!!
            }
        } else {
            @leftIndex = ($j);
        }
    }

    my %bI;
    my ($tmpid, $s, $e) = ('.', 0, 0);
    for my $j (@leftIndex) { # Check Overlap and rm overlap
        $bI{$j} = 1;
        if ($tmpid eq '.') {
            ($tmpid, $s, $e) = @{$$info[$j]}[$rI, $sI, $eI] if ($$info[$j][1] == 0);
        } else {

            if ($tmpid ne $$info[$j][$rI]) {
                ($tmpid, $s, $e) = @{$$info[$j]}[$rI, $sI, $eI];
                next;
            }
            if ($$info[$j][$sI] <= $e and $$info[$j][$eI] >= $s) {
                $$info[$j][0] = ($$info[$j][1] == 0.0) ? -2: -1; # Always keep the last good one
                    ($tmpid, $s, $e) = @{$$info[$j]}[$rI, $sI, $eI];
            }
            $s = $$info[$j][$sI] if ($$info[$j][$sI] < $s);
            $e = $$info[$j][$eI] if ($$info[$j][$eI] > $e)
        }
    }
#my $am = (scalar(keys %bI) > 1) ? 1: 0;

    @leftIndex = ();
    for my $j (@index) {
        $$info[$j][0] = -1  if (!exists $bI{$j} and $$info[$j][0] == 0);
        push @leftIndex, $j if $$info[$j][0] == 0; # Not -1 either -2
    }
    $$info[$leftIndex[0]][0] = 1 if @leftIndex > 0; # Duplication but PASS

# Others treat to DUPLIC
    for (my $i = 1; $i < @leftIndex; ++$i) {
        my $j = $leftIndex[$i];
        $$info[$j][0] = ($$info[$j][1] == 0.0) ? -2: -1;
    }
#    StatisticOverlap($info, @index) if (@index > 1);
}

sub StatisticOverlap {

    my ($info, @index) = @_;
    my (@pos, @nr, @type, @size, @spn);

##[$mark, $nr, $vq, $asmNum, $nummap, $svtype, $svsize, $tId, $tStart, $tEnd, @t]
    for my $i (@index) {

        if ($$info[$i][0] < 0) {

            my $m = ($$info[$i][0] == -2) ? '^' : '*'; # -2 => DUPLIC
            push @pos , $m."$$info[$i][10]-$$info[$i][11]";
            push @type, $m."$$info[$i][5]";
            push @size, $m."$$info[$i][6]";
            push @nr  , $m."$$info[$i][1]";
            push @spn , $m."$$info[$i][3]";			
        } else {

            push @pos , "$$info[$i][10]-$$info[$i][11]";
            push @type, "$$info[$i][5]";
            push @size, "$$info[$i][6]";
            push @nr  , "$$info[$i][1]";
            push @spn , "$$info[$i][3]";			
        }
    }
    print STDERR join "\t", "%%",  (join ",", @pos), (join ",", @type), 
          (join ",", @size), (join ",", @nr),  (join ",", @spn), "\n";	
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
              -d  [int]  Variants overlap distance. [$distance]
              -q  [int]  Threshold for Variant qsulity score. [$qualityThd]

U
     exit(0);
}










