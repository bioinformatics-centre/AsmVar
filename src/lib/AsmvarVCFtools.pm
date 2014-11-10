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
    IsNoGenotype
    GetDataInSpFormat
    GetDataInSpInfo
);

########################### Functions ################################

sub GetAltIdxByGTforSample {
# Get the ALT sequence according to 'GT' information
# Input : 'GT' field per sample in VCF, 'GT' should not be like './.'
# Output: Index of ALT-Sequence in 'ALT' field of VCF.
#     But I put REF-Sequence and ALT-Sequence in one array, 
#     and REF-Seq is the first element of this array. So the 
#     index should count in the REF-Sequence together.
#     e.g : if 'GT' = '1/1', return 1 not 0. and '0/2' will return 2 not 1!! 
#           index 0 will always => REF-Sequence

    my ($gt) = @_;

    # It could just be [0,0] or [0,1] or [1,1] or [0,2] or [2,2]
    my @gt  = split /\//, $gt;
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

    my $size = length($altseq) - length($refseq);
    my ($svtype, $svsize);

    # The front number of svtype is used for sorting the output order
    if ($initSVtype =~ /INV/) { # Inversion
        ($svtype, $svsize) = ("3.INV", $initSVsize);
    } elsif ($initSVtype =~ /TRANS/) { # Translocation
        ($svtype, $svsize) = ("4.TRANS", $initSVsize);
    } elsif ($size > 0) { # Insertion
        ($svtype, $svsize) = ("1.INS", abs($size));
    } elsif ($size < 0) { # Deletion
        ($svtype, $svsize) = ("2.DEL", abs($size));
    } else { # The genotype is 0/0 or it's SNP
        ($svtype, $svsize) = ("9.REF_OR_SNP", length($refseq));
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
    } else { # The genotype is 0/0 or it's SNP
        ($svtype, $svsize) = ("9.REF_OR_SNP", length($refseq));
    }

    return ($svtype, $svsize);
}

sub GetSVforAllPerVariantLine {
# Get the population SV typs and size for per-variant in one VCF line
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

sub IsNoGenotype {
# Determining the samples got non-reference genotype or not
# Input : Samples field with 'FORMAT' format and must contain 'GT'
# Output: boolean

    my (@samples) = @_;

    my $isNogt = 1;
    for (my $i = 0; $i < @samples; ++$i) {

        my $gt = (split /:/, $samples[$i])[0];
        if ($gt ne '0/0' and $gt ne './.') {
            $isNogt = 0;
            last;
        }
    }
    return $isNogt;
}

sub GetDataInSpFormat {
# Get the specific data by specific field in FORMAT of each sample
# Input: (1) 'Specific_format_field'. e.g: 'GT', 'TR' ...
#        (2) 'FORMAT' format
#        (3) All samples in VCF
# Output: An array which recording this specific field data
#         and igonore the one which doesn't contain infomation in this field
#
    my ($sf, $format, @samples) = @_;
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
    for (my $i = 0; $i < @samples; ++$i) {

        my @s = split /:/, $samples[$i];

        next if @s < $fm{$sf} or $s[$fm{$sf}] eq '.' or $s[$fm{$sf}] eq './.';
        push @data, $s[$fm{$sf}];
    }

    return @data;

}

sub GetDataInSpInfo {
# Get the specific data by specific field in INFO of each sample
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





