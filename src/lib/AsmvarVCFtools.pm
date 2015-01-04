# Author : Shujia Huang
# Date   : 2014-11-08
package AsmvarVCFtools;
# Package for deal with the AsmVar VCF

$VERSION = '0.0.1';
$DATE    = '2014-11-08';
$Modify  = '2014-11-';
$AUTHOR  = 'Shujia Huang';

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw (

    GetAltIdxByGTforSample
    GetSVtypeAndSizeForSample
    GetGATKSVtypeAndSizeForSample
    GetSVforAllPerVariantLine
    GetSVforPop
    RecalcuSVBreakpoint
    IsNoGenotype
    GetDataInSpFormat
    GetDataInSpInfo
    IsSpInfo
);

########################### Functions ################################

sub GetAltIdxByGTforSample {
# Get the ALT sequence according to 'GT' information
# Input : 'GT' field per sample in VCF, 'GT' should not be like './.'
# Output: Index of ALT-Sequence in 'ALT' field.
#     e.g : if  'GT' = '0/1' '0|1' or '1/1' '1|1', will return 1, 
#           and 'GT' = '0/2' '0|2' or '2/2' '2|2', will return 2,
#           index 0 will always means REF-Sequence

    my ($gt) = @_;

    # It could just be [0,0] or [0,1] or [1,1] or [0,2] or [2,2]
    my @gt  = split /[\|\/]/, $gt;
    my $ind = 0;
    for my $i (@gt) {
        if ($i > 0) { $ind = $i; last; }
    }
#print STDERR "[Debug] GetAltIdx(): ind = $ind\t$gt\n";
    return $ind;
}

sub GetSVtypeAndSizeForSample {
# Get the SV type and size for specific sample
# Input : Ref-sequence, 
#         Alt-sequence, 
#         Init svsize get in 'VS'
#         Init svtype get in 'VT'
# Output: ($svtype, $svsize) for the sample

    my ($refseq, $altseq, $initSVsize, $initSVtype) = @_;

    $initSVsize = abs($initSVsize); # Make sure it's positive value

    my $size = length($altseq) - length($refseq);
    my ($svtype, $svsize);
    # The front number of svtype is used for sorting the output order
    if ($initSVtype =~ /INV/) { # Inversion
        ($svtype, $svsize) = ("5.INV", $initSVsize);
    } elsif ($initSVtype =~ /TRANS/) { # Translocation
        ($svtype, $svsize) = ("6.TRANS", $initSVsize);
    } elsif ($initSVtype eq 'COMPLEX' or $initSVtype eq 'REPLACEMENT') {
        ($svtype, $svsize) = ("7.REPLACEMENT", $initSVsize);
    } elsif ($initSVtype eq 'MNP') {
        ($svtype, $svsize) = ("3.MNP", $initSVsize);
    } elsif ($size > 0) { # Insertion
        ($svtype, $svsize) = ("1.INS", abs($size));
    } elsif ($size < 0) { # Deletion
        ($svtype, $svsize) = ("2.DEL", abs($size));
    } elsif ($altseq ne '.') { # $size == 0 and the genotype is not 0/0
        ($svtype, $svsize) = (length($altseq) == 1) ?
                             ("4.SNP", length($altseq)):
                             ("3.MNP", length($altseq)); # May useless here
    } else {
        ($svtype, $svsize) = ("9.REF", length($refseq));
    }

#print STDERR "[Debug] GetSVtypeAndSize(): ($refseq, $altseq, $initSVsize, $initSVtype): size = $size :($svtype, $svsize)\n";
    return ($svtype, $svsize);
}

sub GetGATKSVtypeAndSizeForSample {
# For GATK vcf
# Get the SV type and size for specific sample from VCF creat by GATK
# Input : Ref-sequence, 
#         Alt-sequence
# Output: ($svtype, $svsize) for the sample
    my ($refseq, $altseq) = @_;

    my $size = length($altseq) - length($refseq);
    my ($svtype, $svsize);

    if ($size > 0) { # Insertion
        ($svtype, $svsize) = ("1.INS", abs($size));
    } elsif ($size < 0) { # Deletion
        ($svtype, $svsize) = ("2.DEL", abs($size));
    } elsif ($altseq ne '.') { # $size == 0 and the genotype is not 0/0
        ($svtype, $svsize) = (length($altseq) == 1) ? 
                             ("4.SNP", length($altseq)): 
                             ("3.MNP", length($altseq));
    } else {
        ($svtype, $svsize) = ("9.REF", length($refseq));
    }

    return ($svtype, $svsize);
}

sub GetSVforAllPerVariantLine {
# Get the population SV types and size for per-variant in one VCF line
# Input : hash with svtype, svtype=>[num_support_svtype, [svtype,svsize]]
# Output: The most supported ($svtype, $svsize) for this VCF line
    my ($svstat) = @_;

    my ($svtype, $svsize);
    my $svnum = 0;
    for my $k (keys %$svstat) {

        if ($k eq 'INV' or $k eq 'TRANS') {
            ($svtype, $svsize) = @{$$svstat{$k}->[1]};
            last;
        } else { # INDEL
            if ($svnum < $$svstat{$k}->[0]) {
                $svnum = $$svstat{$k}->[0];
                ($svtype, $svsize) = @{$$svstat{$k}->[1]}
            }
        }
    }

#print STDERR "\n[Debug] GetSVforAll(): $svtype\t$svsize\n";
    return ($svtype, $svsize);
}

sub GetSVforPop {
# Get the population SV types and size for per-variant in one VCF line
# We use GATK style to determine the SV type, it means we don't care 
# about the MNP, INV, REPLACEMENT and TRANSLOCATION
#
# Input : 
#         (1) Sequence of REF field
#         (2) Sequence of ALT field
#         (3) The reference of all samples' fields array
# Output: The most supported ($altIndex, $svtype, $svsize) for this VCF line
# 
    my ($refseq, $altseq, $samples) = @_;

    my @seq = ($refseq); # First element is REF: [0]=>REF
    push @seq, $_ for (split /,/, $altseq);

    my %svstat;
    my $isempty = 1;
	for (my $i = 0; $i < @$samples; ++$i) {

        my @f = split /:/, $$samples[$i];
        # Ignore un-genotype and reference type samples
        next if $f[0] eq './.' or $f[0] eq '0/0' or $f[0] eq '0|0';

        # Get the ALT sequence index
        my $ai = GetAltIdxByGTforSample($f[0]);
        my ($svtype, $svsize);
  
        # We use GATK Style to determine the SV Type.
        # So that the MNP, INV, Translocation will treat to be the 'REF_OR_SNP'
        ($svtype, $svsize) = GetGATKSVtypeAndSizeForSample(
                $seq[0],    # Ref-sequence
                $seq[$ai]); # Alt-sequence

        $svtype = (split /\./, $svtype)[-1];
        # Use for getting SVforAll(population) in this position
        $svstat{"$svtype.$ai"}->[0] ++;
        $svstat{"$svtype.$ai"}->[1] = [$ai, $svtype, $svsize];
        $isempty = 0;
    }
    return (0, '-', 0) if $isempty;

    my ($altIndex, $svtype, $svsize);
    my $svnum = 0;    
    for my $k (keys %svstat) {

        if ($svnum < $svstat{$k}->[0]) {
            $svnum = $svstat{$k}->[0];
            ($altIndex, $svtype, $svsize) = @{$svstat{$k}->[1]};
        }
    }

#print STDERR "\n[Debug] GetSVforPop(): ($altIndex, $svtype, $svsize)\n";
    return ($altIndex, $svtype, $svsize);
}

sub RecalcuSVBreakpoint {
# Recalculation the Break point position according to the REF-Seq and ALT-Seq
# Often, this function should just be called after calling 'GetSVforPop()'
# Input: (1) REF_POS in POS field of VCF
#        (2) REF Sequence in REF field of VCF
#        (3) The identity ALT sequence, get from 'GetSVforPop()'
#        (4) The SV Type get from 'GetSVforPop()'
#
# Output:
#        (1) recalculate breakpoint on REF
#        (2) SVTYPE: No change, still get from 'GetSVforPop()'
#        (3) variant sequence
#
    my ($refpos, $refseq, $oneAltSeq, $svtype) = @_;

    my $varseq; # Record the sequence of variant
    if (length($oneAltSeq) > length($refseq)) { # INSERTION like

        $refpos += length($refseq) - 1;
        $varseq = substr($oneAltSeq, length($refseq));
    } elsif (length($oneAltSeq) < length($refseq)) { # DELETION like

        $refpos += length($oneAltSeq);
        $varseq = substr($refseq, length($oneAltSeq));
    } else { # '=='
        $varseq = $oneAltSeq; # Nothing change
    }

#print STDERR "\n[Debug] RecalcuSVBreakpoint(): ($refpos, $svtype, $varseq)\n";
    return ($refpos, $svtype, $varseq);
} 

sub IsNoGenotype {
# Determining the samples got non-reference genotype or not
# Input : Each vcf line reference, but we just use sample field 
#         with 'FORMAT' format and must contain 'GT'
# Output: boolean

    my ($vcfline) = @_;

    my $isNogt = 1;
    for (my $i = 9; $i < @$vcfline; ++$i) {

        # Just loop all samples
        my $gt = (split /:/, $$vcfline[$i])[0];
        if ($gt ne '0/0' and $gt ne '0|0' and $gt ne './.') {
            $isNogt = 0;
            last;
        }
    }
    return $isNogt;
}

sub GetDataInSpFormat {
# Get the specific data by specific field in FORMAT for each sample
# Input: (1) 'Specific_format_field'. e.g: 'GT', 'TR' ...
#        (2) 'FORMAT' format
#        (3) All samples in VCF
# Output: An array which recording this specific field data
#         and return NULL if it doesn't contain infomation in this field
#
    my ($sf, $format, $samples) = @_;
    my @format = split /:/, $format;
    my %fm;
    for (my $i = 0; $i < @format; ++$i) { 
        $fm{$format[$i]} = $i;
    }

    if (not exists $fm{$sf}) {
        print STDERR "[WARN]  USER ERROR. '$sf' is not in FORMAT($format)! ". 
                     "the program will continue and return empty data here ".
                     "for '$sf' field for all the samples\n";
        return ();
    }

    my @data;
    for (my $i = 0; $i < @$samples; ++$i) {

        my @s = split /:/, $$samples[$i];

        next if @s < $fm{$sf} or $s[$fm{$sf}] eq '.' or $s[$fm{$sf}] eq './.';
        push @data, $s[$fm{$sf}];
    }

    return @data;
}

sub GetDataInSpInfo {
# Get the specific data by specific field in INFO for each sample
# Input: (1) 'Specific_INFO_field'. e.g: 'AC', 'AF', 'NRatio', 'VQ' ...
#        (2) Reference of INFO
# Output: Return the data of this field
#
# e.g. AsmvarVCFtools::GetDataInSpInfo('VQ', \$t[7])
# 
    my ($spInfo, $info) = @_;

    my ($data) = $$info =~ m/;$spInfo=([^;]+)/;
    if (not defined $data) {# Must be the first one
        ($data) = $$info =~ m/^$spInfo=([^;]+)/;
    }

    return $data;
}

sub IsSpInfo {
# Determine whether the specific field is in INFO or not
# just reture a bool value not the value for this field
# Input: (1) 'Specific_INFO_field'. e.g: 'AC', 'AF', 'NRatio', 'VQ' ...
#        (2) Reference of INFO
# Output: Return 1(Yes) or 0(No)
#
# e.g. AsmvarVCFtools::GetDataInSpInfo('VQ', \$t[7])
#
    my ($spInfo, $info) = @_;
    my $isHaveInfo = 0;
    $isHaveInfo = 1 if $$info =~ /;$spInfo/ or $$info =~ /^$spInfo/;

    return $isHaveInfo;
}




