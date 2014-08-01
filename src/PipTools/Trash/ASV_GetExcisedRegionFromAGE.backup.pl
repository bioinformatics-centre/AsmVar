# Author : Shujia Huang
# Date   : 2013-05-26
# Date	 : 2013-09-22 Siyang Liu
#Last Modify : 
#				2013-12-20 23:01:43 
#					1. calculate the NRatio together with 100bp flank region.
#					2. Type: oritype ins and AGE ins-> SInS; oritype del and AGE del -> Sdel; Trans -> make a better mark- remember the partner
#				2013/11/07  Rewrite the program update many things and debug
# Input  : age result
# Output : VCF format

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#for b37
my @sort_b37 = ("1_999999_3000000");
#my @sort_b37 = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y", "MT", "GL000207.1", "GL000226.1", "GL000229.1", "GL000231.1", "GL000210.1", "GL000239.1", "GL000235.1", "GL000201.1", "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1", "GL000246.1", "GL000249.1", "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1", "GL000202.1", "GL000234.1", "GL000232.1", "GL000206.1", "GL000240.1", "GL000236.1", "GL000241.1", "GL000243.1", "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1", "GL000204.1", "GL000198.1", "GL000208.1", "GL000191.1", "GL000227.1", "GL000228.1", "GL000214.1", "GL000221.1", "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1", "GL000211.1", "GL000199.1", "GL000217.1", "GL000216.1", "GL000215.1", "GL000205.1", "GL000219.1", "GL000224.1", "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1", "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1", "NC_007605", "hs37d5");
# For b36
my @sort_b36 = ( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "NT_113887", "NT_113947", "NT_113903", "NT_113908", "NT_113940", "NT_113917", "NT_113963", "NT_113876", "NT_113950", "NT_113946", "NT_113920", "NT_113911", "NT_113907", "NT_113937", "NT_113941", "NT_113909", "NT_113921", "NT_113919", "NT_113960", "NT_113945", "NT_113879", "NT_113938", "NT_113928", "NT_113906", "NT_113904", "NT_113873", "NT_113966", "NT_113943", "NT_113914", "NT_113948", "NT_113886", "NT_113932", "NT_113929", "NT_113878", "NT_113927", "NT_113900", "NT_113918", "NT_113875", "NT_113942", "NT_113926", "NT_113934", "NT_113954", "NT_113953", "NT_113874", "NT_113883", "NT_113924", "NT_113933", "NT_113884", "NT_113890", "NT_113870", "NT_113881", "NT_113939", "NT_113956", "NT_113951", "NT_113902", "NT_113913", "NT_113958", "NT_113949", "NT_113889", "NT_113936", "NT_113957", "NT_113961", "NT_113925", "NT_113882", "NT_113916", "NT_113930", "NT_113955", "NT_113944", "NT_113901", "NT_113905", "NT_113872", "NT_113952", "NT_113912", "NT_113935", "NT_113880", "NT_113931", "NT_113923", "NT_113915", "NT_113885", "NT_113888", "NT_113871", "NT_113964", "NT_113877", "NT_113910", "NT_113962", "NT_113899", "NT_113965", "NT_113898", "NC_007605");

my ( @ageFile, $filelist, $ageFile, $tarFaFile, $qryFaFile, $qryId, $tarId, $sampleID, $bV, $outfile, $header );
my $type  = "ALL";
$sampleID = "Sample";
$bV       = 37;
GetOptions(

	"l=s"	=> \$filelist,
	"r=s"	=> \$tarFaFile,
	"q=s"	=> \$qryFaFile,
	"s=s"	=> \$sampleID,
	"o=s"	=> \$outfile,
	"b=i"	=> \$bV,
	"T=s"	=> \$type,
	"h"		=> \$header,
);
@ageFile = @ARGV;
ReadList ( $filelist, \@ageFile ) if ( defined $filelist );
Usage() if ( @ageFile == 0 or !$tarFaFile or !$qryFaFile or !defined $outfile );
die "-b can just be 36 or 37\n" if ( $bV != 36 and $bV != 37 );
print STDERR "# @ageFile\ttarFaFile: $tarFaFile\tqryFaFile: $qryFaFile\toutfile: $outfile\n";
print STDERR "\n*** The main parameter Options: perl $0 -T $type -t $tarFaFile -q $qryFaFile -s $sampleID -b $bV -o $outfile ***\n\n";

my ( %tarFa, %qryFa );
my $localTime = gmtime ( time() );
print STDERR "****** Time Before loading Fa. $localTime ******\n";
ReadFaSeq ( $tarFaFile, \%tarFa    );
ReadFaSeq ( $qryFaFile, \%qryFa    );
$localTime = gmtime ( time() );
print STDERR "****** Time After  loading Fa. $localTime ******\n";

my ( $vcfRefId, $vcfRefPos,$vcfQuePos, $vcfId, $vcfRefAle, $vcfalt, $vcfQuality, $vcfFilter, $vcfInfo, $vcfFormat);
my ( %output, %markTrans );
for my $file ( @ageFile ) {
	my ( $age_id, $AGEInfo ) = ( "N" ) x 2;
	my ( $varTstart, $varTend, $varQstart, $varQend, $varLen, $tlen, $qlen, $nratio, $varType );
	my ( $first, $second, $strand );
	my ( $tLeftStart, $tLeftEnd, $tRightStart, $tRightEnd );
	my ( $qLeftStart, $qLeftEnd, $qRightStart, $qRightEnd );
	my ( $ave_base,$ave_iden,$left_base,$left_iden,$right_base,$right_iden );

	my ( @info, @info2, $ori_t,$ori_s,$ori_e, $qori_s,$qori_e );
	#@info=split(/[\.]/,basename($file));          #Just default if The age file didn't have: # Ins.+.16.53183502.53183502.1-scaffold146175.3234.3234.1.age'
    #@info2=split(/-/,$info[4]);                   #Just default
    #($ori_t,$ori_s,$ori_e)=(@info[0,3],$info2[0]);#Just default
	my $flag = 0; # A marker to monitor the process is accurate or not

	print STDERR ">>>>>>>>>>> Reading $file <<<<<<<<<<<<\n";
	open I, $file or die "Cannot open file $file : $!\n";
	while ( <I> ) {	

		chomp;
		#---
		#1. filteration
		#---
		if ( $_ =~ /^#/ ) {
		
			die "[ERROR]Something unexpected found in you input file $file. maybe cause by the format problem!\n$_\n" if ( $flag );
			@info     = split( /\./, basename((split /\s+/, $_)[1]) ); #  '# Ins.+.16.53183502.53183502.1-scaffold146175.3234.3234.1.age'
        	my $index = 3;
        	if ($info[2] =~ m/^GL/ ) { ++$index; $info[2] .= ".1"; } # Reference ID is GL***.1
			($ori_t, $ori_s, $ori_e, $qori_s,$qori_e) = @info[0, $index, $index+1, -4, -3 ];
			#($ori_t, $ori_s, $ori_e, $qori_s,$qori_e) = @info[0, $index, $index+1, $index+3, $index+4 ];
			die "[ERROR] Original start > Original End ($_) in file $file \n" if ( $ori_s > $ori_e );
		} elsif ($_ =~ m/^MATCH\s*=/ ) { # Start of the AGE mapping information

			die "[ERROR]Something unexpected found in you input file $file .\n$_\n" if ( $flag );
			<I>; # Space line
			chomp( $first  = <I> ); # First  seq [120994137,121012321] =>     18185 nucs 'X dna:chromosome chromosome:GRCh37:X:1:155270560:1'
			chomp( $second = <I> ); # Second seq [     9651,     8650] =>      1002 nucs 'scaffold45598'
			$first =~ m/nucs\s+'(\S+).*'$/; $tarId = $1;
			$second=~ m/nucs\s+'(\S+).*'$/; $qryId = $1;
			die "Doesn't have the target id $tarId in $tarFaFile\n" if ( !exists $tarFa{$tarId} );
			die "Doesn't have the query  id $qryId in $qryFaFile\n" if ( !exists $qryFa{$qryId} );

			$flag = 1; # Now Turn the light on!
			$vcfRefId =$tarId; $vcfId = "."; $vcfQuality = 255; $vcfFilter = "."; $vcfFormat = "GT:AI:CR:EO:NM:RD:SE:SR:UD";
		} elsif($_ =~ m/^Identic:\s*(\d+)\s*\(\s*(\d+)%\)\s*nucs\s*=>\s*(\d+)\s*\(\s*(\d+)%\)\s*(\d+)\s*\(\s*(\d+)%\)/) {
		#Identic:       998 (100%) nucs =>       499 (100%)       499 (100%)		

			die "[ERROR]Light must be on here!\nYou input file $file must missed the pattern 'MATCH'.\n$_\n" if ( !$flag ); # The light is still off!
			($ave_base,$ave_iden,$left_base,$left_iden,$right_base,$right_iden)=($1,$2,$3,$4,$5,$6);
			( $age_id, $vcfFilter ) = ( "T", "ageT" ); # Default value
			( $age_id, $vcfFilter ) = ( "F", "ageF" ) if ($ave_iden<98 || $left_base<30 || $left_iden<95 || $right_base<30 || $right_iden<95);
			$AGEInfo = join( ":", ($age_id,$ave_base,$ave_iden,$left_base,$left_iden,$right_base,$right_iden) );
		} elsif( $_ =~ m/^Alignment:/ ) {
		#---
		#2. Deal with the alignments
		#---
			die "[ERROR]Light still off!\nYou input file $file must missed the pattern 'MATCH'.\n$_\n" if ( !$flag ); # The light is still off!
			chomp( $first  = <I> );	# first  seq =>  [7028678,7029175] EXCISED REGION [7029176,7031069]
			chomp( $second = <I> ); # second seq =>  [    111,    608] EXCISED REGION [    610,   2503]

			if ( ($first =~ m/\s+EXCISED\s+REGION\s+/) && ($second =~ m/\s+EXCISED\s+REGION\s+/) ) {

				( $varTstart, $varTend, $varLen, $nratio ) = ( 0 ) x 4;
				( $varType, $strand ) = (".") x 2;

				$first  =~ m/\s+\[\s*(\d+),\s*(\d+)\]\s+EXCISED\s+REGION\s+\[\s*(\d+),\s*(\d+)\]/;
				( $tLeftStart, $tLeftEnd, $tRightStart, $tRightEnd ) = ( $1, $2, $3, $4 );
				$second =~ m/\s+\[\s*(\d+),\s*(\d+)\]\s+EXCISED\s+REGION\s+\[\s*(\d+),\s*(\d+)\]/;
				( $qLeftStart, $qLeftEnd, $qRightStart, $qRightEnd ) = ( $1, $2, $3, $4 );
				$strand = ($qRightEnd > $qLeftStart) ? "+" : "-";

				$qlen   = abs( $qRightStart - $qLeftEnd );
				$tlen   = $tRightStart - $tLeftEnd;

				( $vcfRefPos, $varTstart, $varTend ) = ( $tLeftEnd, $tLeftEnd + 1, $tRightStart - 1 );
				( $varTstart, $varTend )             = ( $tLeftEnd, $tLeftEnd ) if ( $tlen == 1 ); # The position is continue. The SV type maybe INS

				($varQstart, $varQend) = ( $qLeftEnd, $qRightStart);
				if ( $strand eq '+' ) { $varQstart += 1; $varQend -= 1; $vcfQuePos = $varQstart - 1; }
				if ( $strand eq '-' ) { $varQstart -= 1; $varQend += 1; $vcfQuePos = $varQend;       }

				$varLen = $varTend - $varTstart + 1;
				if ( $ori_t =~ m/^Tran/i && $ori_t !~ m/^Tran.+\-E$/ ) { # Translocation: need to recover to the original region 
					$varType   = "TRANS"; 
					$strand    = $info[1]; #
					($vcfRefPos, $varTstart, $varTend) = ( $ori_s, $ori_s, $ori_e );
					($vcfQuePos, $varQstart, $varQend) = ( $strand eq '+' ) ? ($qori_s, $qori_s, $qori_e) : ($qori_e, $qori_e, $qori_s);
					$qlen      = $qori_e  - $qori_s    + 1;
					$tlen      = $varTend - $varTstart + 1;
					$varLen    = $tlen;

				} elsif ( $ori_t =~ m/^Inv/i && $age_id eq 'T' ) {
					$vcfRefPos = $ori_s; # Changing to original position if it's inversion region
					$varLen    = $ori_e - $ori_s + 1;
					$tlen      = $varLen;
					$varType   = "INV";
				} elsif ( $tlen > $qlen && $qlen == 1     ) {
				# PureDeletion
					$varQstart = ( $strand eq '+' ) ? $qLeftEnd : $qRightStart;
					$varQend   = $varQstart;
					$varType   = ($ori_t =~ m/^Del/i) ? "SDEL": "DEL"; # I should treat the type to be 'SDEL' if they're agree with SVD's original type.
				} elsif ( ($tlen < $qlen) && ($tlen == 1) ) {
				#PureInsertion
					$varLen    = $qlen - 1; # Should substract 1 cause by the boundary
					$varType   = ($ori_t =~ m/^Ins/i) ? "SINS": "INS"; # I should treat the type to be 'SINS' if they're agree with SVD's original type.
				} elsif ( $tlen == $qlen ) {
				# Substitution or SNPs
					$varType   = ( $varLen == 1 ) ? "SNP" : "BSubstitution";
					$vcfRefPos++ if ( $varType eq "SNP" ); # as the same '$tLeftEnd + 1'
				} elsif ( ($tRightStart - $tLeftEnd) != $qlen && ($tRightStart - $tLeftEnd) !=1 && $qlen !=1 ){ 
				#Simul Gap
					$varType   = "SimulGap";
				} else {
					$varType   = "Unknown";
				}

				my ( @CIPos, @CIEnd ); #Confidence interval 
				my $tmpnum = 0;
				while ( $_ !~ /^Identity\s+/ ) { $_ = <I>; ++$tmpnum; die "[ERROR]Your input format error! $_\n" if ($tmpnum>10); }# I don't care these lines 
				#Identity at breakpoints
				chomp( $first  = <I> ); # first  seq =>         3 nucs [14367089,14367091] to [14367159,14367161]
				chomp( $second = <I> ); # second seq =>         2 nucs [29,28] to [28,27]
				if ( $first  =~ m/\s+\[\s*(\d+),\s*(\d+)\]\s+to\s+\[\s*(\d+),\s*(\d+)\]/ ) { push (@CIPos, [$1,$2]); push (@CIEnd, [$3,$4]); }
				chomp ( $_ = <I> ); die "[ERROR]Program ERROR. Cannot find the word 'Identity'.\n$_\n" if ( $_ !~ /^Identity\s+/ );
				#Identity outside breakpoints: 
				chomp( $first  = <I> ); # first  seq =>         3 nucs [14367086,14367088] to [14367159,14367161]
				chomp( $second = <I> ); # second seq =>        11 nucs [38,28] to [26,16]
				if ( $first  =~ m/\s+\[\s*(\d+),\s*(\d+)\]\s+to\s+\[\s*(\d+),\s*(\d+)\]/ ) { push (@CIPos, [$1,$2]); push (@CIEnd, [$3,$4]); }
				chomp ( $_ = <I> ); die "[ERROR]Program ERROR. Cannot find the word 'Identity'.\n$_\n" if ( $_ !~ /^Identity\s+/ );
				#Identity inside breakpoints: 
				chomp( $first  = <I> ); # first  seq =>         1 nucs [14367089,14367089] to [14367158,14367158]
				chomp( $second = <I> ); # second seq =>         0 nucs
				if ( $first  =~ m/\s+\[\s*(\d+),\s*(\d+)\]\s+to\s+\[\s*(\d+),\s*(\d+)\]/ ) { push (@CIPos, [$1,$2]); push (@CIEnd, [$3,$4]); }

				my ( $CIPos, $CIEnd ) = ( "0,0", "0,0" );
				my ( $tmpS , $tmpE  );
				if ( $ori_t !~ m/^Tran/i && $ori_t !~ m/^Inv/i && @CIPos > 0 ) { 
					( $tmpS, $tmpE ) = &Boundary(@CIPos); $tmpS -= $vcfRefPos; $tmpE -= $vcfRefPos; $CIPos = join ",", ( $tmpS, $tmpE ); 
					( $tmpS, $tmpE ) = &Boundary(@CIEnd); $tmpS -= $varTend  ; $tmpE -= $varTend  ; $CIEnd = join ",", ( $tmpS, $tmpE );
				}
				#
				$vcfRefAle = substr( $tarFa{$tarId}, $vcfRefPos - 1, $tlen );
				$vcfalt    = substr( $qryFa{$qryId}, $vcfQuePos - 1, $qlen );
                die "[ERROR] vcfRefAle is empty.\n$tarId\t$vcfRefPos\t$tlen\n" if length($vcfRefAle) == 0;
                die "[ERROR] vcfalt is empty.\n$qryId\t$vcfQuePos\t$qlen\n"    if length($vcfalt)    == 0;

				$vcfRefAle =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
				$vcfalt    =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
				if( $strand eq "-" ) { $vcfalt = uc( scalar reverse($vcfalt) ); $vcfalt =~ tr/ATCGNn/TAGCNn/; }
				$nratio    = NRatio ( substr($qryFa{$qryId}, (($vcfQuePos-100>0) ? $vcfQuePos-100 : 0), $qlen+200) );
				if (uc($vcfRefAle) eq uc($vcfalt)) { #if their the same GATK will catch ERROR when merge with other vcf.
					$vcfalt = ("N") x length($vcfalt);
					$varType= "NULL" if ( $varType eq 'INS' or $varType eq 'DEL' );
				}

				if ( $vcfFilter eq "ageT" and $nratio < 0.2 ) {
					$vcfFilter = "PASS";
				} else {
					$vcfFilter = "FALSE";
				}
				$vcfInfo = "CIPOS=$CIPos;CIEND=$CIEnd;SVTYPE=$varType;SVSIZE=$varLen;TSTART=$varTstart;TEND=$varTend;Q=$qryId;QSTART=$varQstart;QEND=$varQend;STRAND=$strand;NR=$nratio;AGE=$AGEInfo;ORITYPE=$ori_t;ORITSTART=$ori_s;ORITEND=$ori_e";

				# OK. Now For the most difficul part : Change the Format of AGE result and print them to a new file!  2014-04-04 16:15:03
				my (@leftAlign, @rightAlign);
                

			} else { # No excise region! Means The variant we found before was not right!

				$vcfRefPos = $ori_s;
				$vcfRefAle = substr( $tarFa{$tarId}, $ori_s  - 1, $ori_e  - $ori_s  + 1 ); $vcfRefAle =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
				$vcfalt    = substr( $qryFa{$qryId}, $qori_s - 1, $qori_e - $qori_s + 1 ); $vcfalt    =~ tr/WSKMYRVDBHwskmyrvdbh/ACTACAAATAactacaaata/;
                die "[ERROR] vcfRefAle is empty.\n$tarId\t$ori_s  - 1\t$ori_e  - $ori_s  + 1\n" if length($vcfRefAle) == 0;
                die "[ERROR] vcfalt is empty.\n$qryId\t$qori_s - 1\t$qori_e - $qori_s + 1\n"    if length($vcfalt)    == 0;
				$nratio    = NRatio( substr($qryFa{$qryId}, (($qori_s-100>0) ? $qori_s-100 : 0), $qori_e - $qori_s + 200) );
				$vcfalt    = ("N") x length($vcfalt) if ( uc($vcfRefAle) eq uc($vcfalt) );
				$vcfFilter = "FALSE";
				$vcfQuality= 0;
				$vcfInfo   = "CIPOS=0,0;CIEND=0,0;SVTYPE=NULL;SVSIZE=0;TSTART=0;TEND=0;Q=$qryId;QSTART=0;QEND=0;STRAND=NULL;NR=$nratio;AGE=$AGEInfo;ORITYPE=$ori_t;ORITSTART=$ori_s;ORITEND=$ori_e";
			}

			if ( $ori_t =~ m/^Tran.+\-E$/ ) { 
				$markTrans{"$qryId-$qori_s-$qori_e"}->[0] = $vcfFilter; 
				$markTrans{"$qryId-$qori_s-$qori_e"}->[1] = "$varType#$vcfRefId-$varTstart-$varTend";
				$markTrans{"$qryId-$qori_s-$qori_e"}->[2] = "$ori_t#$tarId-$ori_s-$ori_e"; # This is the expected region of translocation
			} elsif( $ori_t =~ m/^Tran/   ) { # Not the expected region
				$markTrans{"$qryId-$qori_s-$qori_e"}->[3] = $vcfFilter;
                $markTrans{"$qryId-$qori_s-$qori_e"}->[4] = "$varType#$vcfRefId-$varTstart-$varTend";
                $markTrans{"$qryId-$qori_s-$qori_e"}->[5] = "$ori_t#$tarId-$ori_s-$ori_e";
			}
			push (@{$output{$vcfRefId}}, [$vcfRefId,$vcfRefPos,$vcfId,$vcfRefAle,$vcfalt,$vcfQuality,$vcfFilter,$vcfInfo,$vcfFormat,"./.:0:0:0:0:0:0:0:0",$ori_t,"$qryId-$qori_s-$qori_e"]);

			($ave_base, $ave_iden, $left_base, $left_iden, $right_base, $right_iden)=("N")x6;
			$flag  = 0; # Turn the light off!
			$tarId = undef; $qryId = undef;
			@info  = ()   ; ($ori_t,$ori_s,$ori_e, $qori_s,$qori_e)=();
			( $ave_base,$ave_iden,$left_base,$left_iden,$right_base,$right_iden ) = ();
		} else { # Other rows in AGE
			next;
		}
	}
	close I;
}

### Output All The variant information
my @sort;
@sort = @sort_b37 if ( $bV == 37 ); 
@sort = @sort_b36 if ( $bV == 36 );

my $H;
my $flag = 0; # Mark
open $H, ">$outfile" or die "Canot write to file $outfile : $!\n";
VcfHeader( $sampleID, $H ) if ( $header );
for my $id ( @sort ) {

	next if ( !exists $output{$id} );
	$flag = 1;
	@{ $output{$id} } = sort { $a->[1] <=> $b->[1] } @{ $output{$id} };
	for ( my $i = 0; $i < @{ $output{$id} }; ++$i ) {

		my $k = pop @{$output{$id}->[$i]}; # "$qryId:$qori_s-$qori_e"
		my $t = pop @{$output{$id}->[$i]}; # "$ori_t"

		#if ( $t !~ m/^Tran.+\-E$/ && exists $markTrans{$k} ) {
		##	$output{$id}->[$i][6] = $markTrans{$k}->[0];
		#	if ( $output{$id}->[$i][6] eq "FALSE" ) {
		#		print O2 "#!";print O2 join "\t",(@{$output{$id}->[$i]}[0,1],(split /:/,$markTrans{$k}->[1]),(split /[-:]/,$markTrans{$k}->[2])),"\n";
		#	} else {
		#		print O2 join "\t", (@{$output{$id}->[$i]}[0,1], (split /:/, $markTrans{$k}->[1]), (split /[-:]/, $markTrans{$k}->[2])),"\n";
		#	}
		#}
		next if ($type ne 'ALL' && $t !~ /^$type/i);

		if ( $t =~ m/^Tran/i && exists $markTrans{$k} ) {
			my $svtypeInfo;
			if ( $t =~ m/^Tran.+\-E$/ ) { # The expected region of translocation
				if ( !$markTrans{$k}->[4] or !$markTrans{$k}->[5] ) {
					print STDERR "[WARNING] You don't have the translocation's expected region. For : \n@{$output{$id}->[$i]}\n";
					$svtypeInfo = "TRANS,$markTrans{$k}->[1],-,-,$k";
				} else {
					$svtypeInfo = "TRANS,$markTrans{$k}->[1],$markTrans{$k}->[4],$markTrans{$k}->[5],$k";
				}
			} else {
				if ( !$markTrans{$k}->[1] or !$markTrans{$k}->[2] ) {
					print STDERR "[WARNING] You don't have the orignal translocation's region. For : \n@{$output{$id}->[$i]}\n";
					$svtypeInfo = "TRANS,$markTrans{$k}->[4],-,-,$k";
				} else {
					$svtypeInfo = "TRANS,$markTrans{$k}->[4],$markTrans{$k}->[1],$markTrans{$k}->[2],$k";
				}
			}
			$output{$id}->[$i][7] =~ s/;SVTYPE=[^;]+;/;SVTYPE=$svtypeInfo;/;
		}
		print $H join "\t", @{$output{$id}->[$i]}; print $H "\n";
	}
}
close $H;
print STDERR "\n[INFO] Output Nothing!! You'd better to check the Target Id to see if they're consist with human b37 or human b36.\n\n" if !$flag;

$localTime = gmtime ( time() );
print STDERR ">>>>>>>>>>> program $0 done <<<<<<<<<<<<\n";
print STDERR "*********** ALL DONE $localTime ********\n";

###############################################################
###############################################################
###############################################################

sub ReadList {
	my ( $filelist, $ageFile ) = @_;
	open  L, $filelist or die "Cannot open $filelist : $!\n";
	while ( <L> ) { chomp; push ( @$ageFile, $_ ); }
	close L;
}

sub ReadFaSeq {

	my ( $file, $fa ) = @_;

	my ( $refId, $seq );
	open I, $file or die "Cannot open file : $file\n";

	$/ = ">"; <I>; $/ = "\n";
	while ( <I> ) {
		chomp;
		$refId = ( split /\s+/ )[0];
		$/ = ">"; chomp( $seq = <I> ); $/ = "\n";
		$seq =~ s/\s+//g;
		$$fa{$refId} = $seq;
	}
	close I;
}

sub NRatio {

	my ( $seq ) = @_;
	my $len = 0;
	++$len while ( $seq =~ m/N/ig );

	return sprintf "%.6f", $len/length($seq);
}

sub Boundary {
	my @reg = @_; # @reg = ( [s1,e1], [s2,e2], [s3,e3] );
	my ( $left, $right ) = @{ $reg[0] };
	for my $r ( @reg ) {

		$left  = $r->[0] if ( $left  > $r->[0] );
		$right = $r->[1] if ( $right < $r->[1] );
	}
	
	return ( $left, $right );
}

sub VcfHeader {
	my ( $sampleID, $H ) = @_;
#GT:AI:CR:EO:NM:RD:SE:SR:UD
	print $H <<V;
##fileformat=VCFv4.1
##FILTER=<ID=FALSE,Description="NONAGE or HIGH N RATIO">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AI,Number=1,Type=String,Description="Abnormal insertsize mapped">
##FORMAT=<ID=CR,Number=1,Type=String,Description="Cross reads. 'I/D' signal">
##FORMAT=<ID=EO,Number=1,Type=String,Description="Erroneous Orientation">
##FORMAT=<ID=NM,Number=1,Type=String,Description="Normal pair-end mapped">
##FORMAT=<ID=RD,Number=1,Type=String,Description="Read depth ratio. The ratio of the depth of SVs region division the depth of flank regions">
##FORMAT=<ID=SE,Number=1,Type=String,Description="Single end reads mapped">
##FORMAT=<ID=SR,Number=1,Type=String,Description="Split reads signal">
##FORMAT=<ID=UD,Number=1,Type=String,Description="Undetermined Reads">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">
##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of SV which was judged from AGE. For tranclocation the format is : TRANS,SelfVarType#RefId-VarTstart-VarTend,MateVarType#VcfRefId-VarTstart-VarTend,MateOrignalVarType#OrignalRefId-OriginalVarTstart-OriginalVarTend,QueryId-QueryStart-QueryEnd">
##INFO=<ID=SVSIZE,Number=1,Type=Integer,Description="The size of SVs">
##INFO=<ID=TSTART,Number=1,Type=Integer,Description="The left breakpoint of structural variant">
##INFO=<ID=TEND,Number=1,Type=Integer,Description="The right breakpoint of structural variant">
##INFO=<ID=Q,Number=1,Type=String,Description="Query sequenes' id">
##INFO=<ID=QSTART,Number=1,Type=Integer,Description="The POS of variant on query">
##INFO=<ID=QEND,Number=1,Type=Integer,Description="The END of variant on query">
##INFO=<ID=STRAND,Number=1,Type=String,Description="The mapped Strand">
##INFO=<ID=NR,Number=1,Type=Float,Description="N ratio of the query sequences">
##INFO=<ID=AGE,Number=1,Type=String,Description="AGE aligment information. (T/F:ave_base:ave_iden:left_base:left_iden:right_base:right_iden)">
##INFO=<ID=ORITYPE,Number=1,Type=String,Description="The original variants' type">
##INFO=<ID=ORITSTART,Number=1,Type=Integer,Description="The original POS of variants">
##INFO=<ID=ORITEND,Number=1,Type=Integer,Description="The original END of variants">
##reference=file:///home/siyang/Database/b36/human_b37_both.fasta
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleID
V

}

sub Usage {
	print STDERR <<U;

Version : 0.0.1 ( 2013-05-26 )
Author  : Shujia Huang
Created : 2013/05/26

Last Modify : 2013/11/07  Update many things and debug

      Usage : perl $0 [Options] -o output [AGE result file]

      Options :

              -r  [str]  Target fa file.
              -q  [str]  Query fa file.
              -l  [str]  AGE file list: Options.
              -s  [str]  Sample ID.     [Sample]
              -o  [str]  Outfile.       [Reguire]
              -b  [int]  Sorted by the chromosome order of b37 or b36. [37]
              -T  [str]  Pick the special SV type to output. [ALL output]
              -h         Output vcf header or not. [No]

U
	 exit(0);
} 


