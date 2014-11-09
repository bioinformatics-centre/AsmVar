# Author : Shujia Huang
# Date   : 2014-11-08
package AsmvarVCFtools;
# Package for deal with the AsmVar VCF

$VERSION = '0.0.1';
$DATE    = '2014-11-08';
$Modify  = '2014-11-08';
$AUTHOR  = 'Shujia Huang';

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw (

    GetAltIdxByGTforSample
    GetSVtypeAndSizeForSample
    GetSVforAllPerVariantLine
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

    # It'll just be [0,0] or [0,1] or [1,1] or [0,2] or [2,2]
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

sub GetSVforAllPerVariantLine {
# Get the population SV typs and size for per-variant in one VCF line
# Input : hash with svtype, svtype=>[num_support_svtype, [svtype,svsize]]
# Output : ($svtype, $svsize) for this variant line(VCF line)
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

