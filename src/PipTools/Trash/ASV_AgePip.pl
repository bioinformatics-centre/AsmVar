#!perl -w
#liusiyang@genomics.cn
use strict;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;

my ($svdfile,$queryFa,$targetFa,$outdir, $shdir, $nameprefix, $ageAlg_parameter ); #$outdir is the output directory of age
my $split;
my $sampleId = "sample";
GetOptions(

	"svd=s"	=> \$svdfile,
	"sampleId"=> \$sampleId,
	"q=s"	=> \$queryFa,   #query list file
	"t=s"   => \$targetFa,
	"o=s"   => \$outdir,
	"s=s"   => \$shdir,
	"n=s"   => \$nameprefix,
	"p=s"   => \$ageAlg_parameter, #age_align -indel -both -match=1 -mismatch=-1 -go=-10 -ge=-1
	"split" => \$split,
);
die "perl $0 <-split> -svd [svdfile] -q [queryFa] -t [targetFa] -s [shell dir] -o [out_dir]\n" if ( !defined $svdfile or !defined $queryFa or !defined $targetFa or !defined $shdir or !defined $outdir );
system ( "mkdir -p $outdir" ) if ( ! -d $outdir );
system ( "mkdir -p $shdir"  ) if ( ! -d $shdir  );

my $progDir = File::Basename::dirname abs_path($0); $progDir =~ s#/net/home##;
my $ageProg = "$progDir/ASV_AgeAlign";
my $ExcisedRegionFromAGEProg = "$progDir/ASM_GetExcisedRegionFromAGE.pl";
#1. store the absolute path of the query.fa and target.fa
&generate_age_sh($svdfile,$targetFa,$queryFa,$shdir,$outdir,$nameprefix);

sub generate_age_sh{

	my ($svdfile,$target,$query,$shdir,$outdir,$in_name)=@_;

	my %targetId;
	open(I,"$svdfile") or die "Cannot open $svdfile\n";
	open U, ">$outdir/$in_name.AgeUnresolved" or die "Cannot write $outdir/$in_name.AgeUnresolved\n";
	open O, ">$outdir/$in_name.variant.txt"   or die "Cannot write $outdir/$in_name.variant.txt\n";

	while(<I>){
		next if /^#/;
		chomp;
		my ($tchr,$tvs,$tve,$tvsize,$tvn,$tlen,$qchr,$qvs,$qve,$qvsize,$qvn,$qlen,$str,$type)=split(/\s+/,$_);  

		if ($qvn >0.5 || $tvn>0.5){ print U "[TGapOrQGap]\t$_\n"; next; }
		if ($type eq "Clip")      { print U "[Clip]\t$_\n"      ; next; }
		if ($type eq "Nomadic")   { print U "[Nomadic]\t$_\n"   ; next; }

		my ($ts,$te)=( ($tvs-500>0)?($tvs-500):1, ($tve+500<=$tlen)?($tve+500):$tlen );
		my ($qs,$qe)=( ($qvs-500>0)?($qvs-500):1, ($qve+500<=$qlen)?($qve+500):$qlen );
		die "[ERROR] ( $qe-$qs<0 || $te-$ts<0 ) Negative.\t$_\n" if ( $qe-$qs<0 || $te-$ts<0 );

		if (5 * abs($qe-$qs) * abs ($te-$ts) /1000000000 >= 10) { print U "[BIG]\t$_"; next; }

		$targetId{$tchr} = 1;
		print O join "\t", ( $tchr,$ts,$te,$qchr,$qs,$qe,"$type.$str.$tchr.$tvs.$tve.$tvsize-$qchr.$qvs.$qve.$qvsize.age" ), "\n";
	}
	close I; close U; close O;

	if ($split) {
		for my $k (sort {$a cmp $b} keys %targetId) {
			open  S, ">$shdir/$in_name.$k.age.sh" or die "Cannot write to $shdir/$in_name.$k.age.sh\n";
			print S "time $ageProg $ageAlg_parameter --selectId $k --target $target --query $query $outdir/$in_name.variant.txt > $outdir/$in_name.$k.age";
			close S;
			system( "chmod 755 $shdir/$in_name.$k.age.sh" );
		}
	} else {
		open  S, ">$shdir/$in_name.age.sh" or die "Cannot write to $shdir/$in_name.age.sh\n";
		print S "time $ageProg $ageAlg_parameter --target $target --query $query $outdir/$in_name.variant.txt > $outdir/$in_name.age";
		close S;
		system( "chmod 755 $shdir/$in_name.age.sh" );
	}
}

sub fa_store{
	my ($infile,$hash)=@_;

	open(IN1,$infile) or die "can not open $infile";
	while(<IN1>){
		chomp;
		my @tmp  = split;
		my $name = getname_from_fa( $tmp[0] );
		$$hash{$name} = $tmp[0];
	}
	close IN1;
}

sub getname_from_fa {

	my ( $fa ) = @_;

	my $refId;
	my $i = 0;
	open  I, $fa or die "Cannot open file $fa\n";
	$/ = ">"; <I>; $/ = "\n";
	while ( <I> ) {
		chomp;
		++$i;
		die "[ERROR] There more than one fa sequence in $fa.\n" if ( $i > 1 );
		$refId = ( split /\s+/ )[0];
		$/ = ">"; <I>; $/ = "\n";
	}
	close I;

	return $refId;
}

sub close_handle{
	my @H = @_;
	close $_ for ( @H );
}
