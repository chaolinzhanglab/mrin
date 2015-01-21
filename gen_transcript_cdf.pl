#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Data::Dumper;

use Bed;
use Common;


my $verbose = 0;
my $prog = basename ($0);
my $sameStrand = 0;

my @ARGV0 = @ARGV;

GetOptions ("s"=>\$sameStrand, 
		#"step:i"=>\$step,
		"v"=>\$verbose);


if (@ARGV != 3)
{
	print "get tag coverage cdf on transcripts\n";
	print "Usage: $prog [options] <transcript.bed> <in.wig> <out.txt>\n";
	print " <in.wig>                 : use - to use stdin; gz file is accepted\n";
	print " <out.txt>                : use - to use stdout\n";
	print " -s                       : same strand required\n";
	print " -v                       : verbose\n";
	exit (1);
}



print STDERR "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($regionBedFile, $inWiggleFile, $outFile) = @ARGV;

print STDERR "load transcripts from $regionBedFile ...\n" if $verbose;

select (STDERR);
my $transcripts = readBedFile ($regionBedFile, $verbose);
select(STDOUT);

my $n = @$transcripts;

print STDERR "$n transcripts loaded\n" if $verbose;

my %regionHash;
foreach my $ts (@$transcripts)
{
	
	my $chrom = $ts->{"chrom"};
	Carp::croak "no strand info in ", Dumper ($ts), "\n" unless exists $ts->{"strand"};

	my $strand = $sameStrand ? $ts->{"strand"} : 'b';
	
	#get exons for the gene
	my $exons = geneToExon ($ts);
	$ts->{'exons'} = $exons;

	foreach my $e (@$exons)
	{

		push @{$regionHash{$strand}->{$chrom}}, $e;
	}
}

my $fin;

my $strand = 'b';
my $chrom = "";
my $currRegions = [];
my @currScores;

if ( $inWiggleFile eq "-")
{
    $fin = *STDIN;
}
else
{
    if ($inWiggleFile =~/\.gz$/)
    {
        open ($fin, "gunzip -c $inWiggleFile | ") || Carp::croak "cannot open file $inWiggleFile to read\n";
    }
    else
    {
        open ($fin, "<$inWiggleFile") || Carp::croak "cannot open file $inWiggleFile to read\n";
    }
}

my $firstTrack = 1;

while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	next if $line=~/^\#/;

	if ($line=~/^track\s/ || $firstTrack)
	{
		print STDERR "a new track found ($line)...\n" if $verbose;
		
		$firstTrack = 0;
		if (@$currRegions > 0 && @currScores > 0)
		{
			my $nr = @$currRegions;
			my $ns = @currScores;
			print STDERR "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
			extractCDF ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
		}

		$chrom = "";
		$currRegions = [];
		@currScores = ();
		$strand = 'b';

		if ($sameStrand)
	   	{
			if($line=~/\(([\+|\-])\)/)
			{
				$strand = $1;
			}
			else
			{
				Carp::croak "no strand information found in the track\n";
			}
		}

		next if $line =~/^track\s/;
	}
	
	my ($chr, $chromStart, $chromEnd, $score) = split (/\t/, $line);

	$chromEnd -= 1;
	
	if ($chr ne $chrom)
	{
		if (@$currRegions > 0 && @currScores > 0)
		{
			my $nr = @$currRegions;
			my $ns = @currScores;
			print STDERR "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
			extractCDF ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
		}

		print STDERR "processing $chr ...\n" if $verbose;
		$chrom = $chr;
		$currRegions = [];

		$currRegions = $regionHash{$strand}->{$chrom} if exists $regionHash{$strand} && exists $regionHash{$strand}->{$chrom};
		@currScores = ();
	}
	
	push @currScores, {chrom=>$chr, chromStart=>$chromStart, chromEnd=>$chromEnd, score=>$score};
}

if (@$currRegions > 0 && @currScores > 0)
{
	my $nr = @$currRegions;
	my $ns = @currScores;
	print STDERR "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
	extractCDF ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
}


close ($fin) if $inWiggleFile ne '-';


print STDERR "dumping results to $outFile ...\n" if $verbose;
my $fout;

if ($outFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
}

foreach my $ts (@$transcripts)
{
	my $exons = $ts->{'exons'};

	my @CDF;
	my $prev = 0;
	for (my $i = 0; $i < @$exons; $i++)
	{
		my $e = $exons->[$i];
		$prev = $CDF[$#CDF] if @CDF > 0;

		if (exists $e->{'CDF'})
		{
			map {push @CDF, $_+$prev} @{$e->{'CDF'}};
		}
		else
		{
			#no tag on this exon
			for (my $j = 0; $j < geneToSize ($e); $j++)
			{
				push @CDF, $prev;
			}
		}
	}
	
	#verify consistency
	my $name = $ts->{'name'};
	my $tsSize = geneToSize ($ts);
	my $n = @CDF;

	Carp::croak "$name: inconsistency in size (ts length = $tsSize, number of values=$n)\n" if $n != $tsSize;

	my $total = $CDF[$#CDF];	
	
	#make sure the CDF start from the 3' end of the transcript
	if ($ts->{'strand'} eq '+')
	{
		@CDF = map {$total - $_} @CDF;
		@CDF = reverse (@CDF);
	}	

	print $fout join ("\t", $ts->{'name'}, @CDF), "\n";
}

close ($fout) unless $outFile eq '-';


sub extractCDF
{
	my ($regions, $scores) = @_;
	my @regionsSorted = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regions;
	
	my $firstScoreIdx = 0;

	print STDERR "extracting CDF ...\n" if $verbose;

	my $k = 0;
	foreach my $r (@regionsSorted)
	{

		print STDERR "$k ...\n" if $verbose && $k % 5000 == 0;
		$k++;

		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};
		while ($firstScoreIdx < @$scores && $scores->[$firstScoreIdx]->{"chromEnd"} < $chromStart)
		{
			$firstScoreIdx++;
		}

		#print "firstIdx = $firstScoreIdx\n" if $verbose;

		my $i = $firstScoreIdx;
		my $prevPos = $chromStart - 1;
		while ($i < @$scores && $scores->[$i]->{"chromStart"} <= $chromEnd)
		{
			my $overlapStart = max ($chromStart, $scores->[$i]->{"chromStart"});
			my $overlapEnd = min ($chromEnd, $scores->[$i]->{"chromEnd"});
			
			if ($overlapEnd >= $overlapStart)
			{
				#in case there is a gap between this block and the previous one
				for (my $j =$prevPos+1; $j < $overlapStart; $j++)
				{
					if (exists $r->{'CDF'})
					{
						my $n = @{$r->{'CDF'}};
						my $prev = $r->{'CDF'}->[$n-1];
						push @{$r->{'CDF'}}, $prev;
					}
					else
					{
						push @{$r->{'CDF'}}, 0;
					}
				
				}
				
				#get data in this block
				for (my $j = $overlapStart; $j <=$overlapEnd; $j++)
				{
					if (exists $r->{'CDF'})
					{
						my $n = @{$r->{'CDF'}};
						my $prev = $r->{'CDF'}->[$n-1];
						push @{$r->{'CDF'}}, $prev+$scores->[$i]->{"score"};
					}
					else
					{
						push @{$r->{'CDF'}}, $scores->[$i]->{"score"};
					}
				}
				$prevPos = $overlapEnd;
			}
			$i++;
		}

		#in case the last block does not reach the end
		for (my $j = $prevPos+1; $j <=$chromEnd; $j++)
		{
			if (exists $r->{'CDF'})
			{
				my $n = @{$r->{'CDF'}};
				my $prev = $r->{'CDF'}->[$n-1];
				push @{$r->{'CDF'}}, $prev;
				#print join ("\t", $j, $prev), "\n";
			}
			else
			{
				push @{$r->{'CDF'}}, 0;
			}
		} 
		
		my $tsSize = geneToSize ($r);
		my $n = @{$r->{'CDF'}};
		my $name = $r->{'name'};

		Carp::croak "$name: inconsistency in size (region length = $tsSize, number of values=$n)\n" if $n != $tsSize;
	}
}


