#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Carp;

use Common;

#ref:
#http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test#Kolmogorov.E2.80.93Smirnov_statistic


my $prog = basename ($0);

my $verbose = 0;

GetOptions ("v"=>\$verbose);

if (@ARGV != 2)
{
	print "calculate ks statistic for devidation from uniform distribution\n";
	print "Usage: $prog [options] <in.txt> <out.txt>\n";
	print " <in.txt> : gz file supported; use - for stdin\n";
	print " <out.txt> : use - for stdout\n";
	print " -v : verbose\n";
	exit (1);
}

my ($inFile, $outFile) = @ARGV;

my ($fin, $fout);

if ( $inFile eq "-")
{
    $fin = *STDIN;
}
elsif ($inFile =~/\.gz$/)
{
	open ($fin, "gunzip -c $inFile | ") || Carp::croak "cannot open file $inFile to read\n";
}
else
{
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
}


if ($outFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
}
print $fout join("\t", "#name", "pos", "abundance", "max_dev", "ks_statistic"), "\n";

my $i = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	print STDERR "$i ...\n" if $i % 1000 == 0 && $verbose;
	$i++;

	my @elems = split ("\t", $line);

	my $name = shift @elems;
	my $total = $elems [$#elems];
	my $n = @elems;

	my $d = 0;
	for (my $i = 0; $i < $n && $total > 0; $i++)
	{
		if (ABS($d) < ABS($elems[$i] / $total - $i/$n))
		{
			$d = $elems[$i] / $total - $i/$n;
		}
	}
	my $ks = sqrt ($total) * $d;
	print $fout join("\t", $name, $n, $total/$n, $d, $ks), "\n";
}



close ($fin) unless $inFile eq '-';
close ($fout) unless $outFile eq '-';
