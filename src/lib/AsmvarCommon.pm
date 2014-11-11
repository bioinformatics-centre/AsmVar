# Author : Shujia Huang
# Date   : 2014-11-10 
package AsmvarCommon;
# package

$VERSION = '0.0.1';
$DATE    = '2014-11-10';
$Modify  = '2014-11-';
$AUTHOR  = 'Shujia Huang';

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = ();

our @EXPORT_OK = qw (

    SizeBin
    SizeBinLog
    SizeBinLog2
    SizeBinLog10
    SizeBinSp
);

###################### Function ######################
sub SizeBin {
# Calculate the size bin
# Inout: (1) size
#        (2) The bin size 
# Return: the bin region in string: start-end
#
    my ($size, $bin) = @_;
    
    my $key   = ($size % $bin) ? int($size / $bin) + 1: int($size / $bin);
    my $start = ($key) ? ($key - 1) * $bin + 1 : 0;
    my $end   = $key * $bin;

    return "$start-$end";
}

sub SizeBinLog {
# Calculate the size bin in logbase
# Input: (1) size
#        (2) logan base
# Return: the bin region in string: start-end
    my ($size, $logbase) = @_;

    die "[ERROR] The data is negative when calling ".
        "log() function.\n" if $size < 0;
    my $key   = ($size > 1) ? int(log($size) / log($logbase)): 0;
    my $start = $logbase ** $key + 1;
    my $end   = $logbase ** ($key + 1);
    
    return "$start-$end";
}

sub SizeBinLog2 {

    my ($size) = @_;
    return SizeBinLog($size, 2);
}

sub SizeBinLog10 {

    my ($size) = @_;
    return SizeBinLog($size, 10);
}

sub SizeBinSp {
# The bin is specific
    my ($size) = @_;

    my $logbin = ($size > 1) ? int(log($size) / log(10)): 0;
    # it'll always turn to be '1-$size', if not step back one
    --$logbin if ($logbin and $size == 10 ** $logbin);

    return SizeBin($size, 10 ** $logbin);
}
