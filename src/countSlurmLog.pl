#!/bin/perl

# This script counts the number of entries per minute in the slurm log file
use strict;
use warnings;

my %hash;

my $logfile = $ARGV[0];
my $outfile = "/mnt/research2/bms41/Pigs/NuMTs/src/" . $logfile . ".total.txt";

open(FILE, $logfile) or die "$!\n";

while(<FILE>){
	chomp;
	if($_ =~ m/(20[0-9]{2}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2})/){
		my $timestamp = $1;

		unless( exists $hash{$timestamp}){
			$hash{$timestamp} = 1;
		} else{
			$hash{$timestamp}++;
		}
	}
	
}

close FILE;

open(OUT, ">$outfile") or die "$!\n";
foreach my $t (sort keys %hash){
	print OUT $t . " : " . $hash{$t} . "\n";
}

close OUT;
print "Done\n";
