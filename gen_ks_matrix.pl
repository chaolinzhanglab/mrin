#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;


my $prog = basename ($0);
my $verbose = 0;

my $minAbundance = 0;

my $base = "";
#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"min-avg-cov:f"=>\$minAbundance,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate ks matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base         [string] : base dir of input data\n";
	print " --min-avg-cov [float]  : min average coverage ($minAbundance)\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;

print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $samples = readConfigFile ($configFile, $base);

print "done.\n" if $verbose;

print "loading data of individual samples ...\n" if $verbose;

my @sampleData;
my $geneId;
my $n = 0;
my $iter = 0;
my $geneInfo;

my @sampleNames = sort {$samples->{$a}->{"id"} <=> $samples->{$b}->{"id"}} keys %$samples;

foreach my $sName (@sampleNames)
{
	my $inputFile = $samples->{$sName}->{"f"};
	$inputFile = "$base/$inputFile" if $base ne '';
		
	my $sdata = readmRINDataFile ($inputFile);
	$geneInfo = $sdata->{"geneInfo"};
	if ($n != 0)
	{
		Carp::croak "data inconsistency detected\n" if @$geneInfo != $n;
	}
	else
	{
		$n = @$geneInfo;
	}
	$sampleData[$iter] = $sdata->{"data"};
	$iter++;
}

print "$iter samples, $n genes loaded.\n" if $verbose;


my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

print $fout join ("\t", "#name", "exon_len", @sampleNames), "\n";

for (my $i = 0; $i < $n; $i++)
{
	my @out;
	for (my $s = 0; $s < @sampleNames; $s++)
	{
		my $sName = $sampleNames[$s];
		my $d = $sampleData[$s][$i];
		$out[$s] = $d->[0] > $minAbundance ? $d->[1] : '';
	}

	#print $fout join ("\t", @{$geneInfo->[$i]}, @out), "\n";
	print $fout join ("\t", $geneInfo->[$i][0], @out), "\n";
}


close ($fout);




sub readConfigFile
{
	my ($configFile, $base) = @_;
	my $fin;
	open ($fin, "<$configFile") || Carp::croak "cannot open file $configFile to read\n";
	my $i = 0;
	my %samples;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		next if $line=~/^\#/;
		my ($fileName, $sampleName) = split (/\t/, $line);
		$samples{$sampleName}->{"id"} = $i++;
		$samples{$sampleName}->{"f"} = $fileName;
	}
	close ($fin);
	return \%samples;
}

sub readmRINDataFile
{
	my ($inputFile) = @_;
	
	my $fin;
	my @data;
	my @geneInfo;
	open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
	
	my $totalCoverage = 0; #reserved
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;
		
		my ($name, $exonLen, $abundance, $maxDev, $ks) = split (/\t/, $line);
		
		push @geneInfo, [$name, $exonLen];
		push @data, [$abundance, $maxDev, $ks];
		
		$totalCoverage += $abundance * $exonLen;
	}
	close ($fin);

	return {geneInfo=>\@geneInfo, data=>\@data};
}


