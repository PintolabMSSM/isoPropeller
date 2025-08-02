#!/usr/bin/env perl
=start
Author: ALin
Purpose: To split a .gtf file into different batches.
Change log:
	v1.1	2023-10	Adopted from split_gtf_v1.3.pl. Modified to enforce transcripts from the same gene are always in the same batch.
=cut

use strict;
use Getopt::Long;

my $version = "v1.1";

my $usage = "perl split_gtf_by_gene_${version}.pl
	-i <String> Input gtf file
	-o <String> Output directory
	-b <Integer> Number of batches (Default: 1)
	-s <Boolean> Sort the output
	-h <Boolean> Help
";

my ($gtf, $out, $batch, $sort, $help) = ("", "", 1, 0, 0);

GetOptions(
	'i=s'	=>	\$gtf,
	'o=s'	=>	\$out,
	'b=i'	=>	\$batch,
	's!'	=>	\$sort,
	'h!'	=>	\$help
);

if($help){
        print "$usage";
        exit;
}

unless(-e $gtf){
	print STDERR "<ERROR> .gtf file, $gtf, does not exist!\n$usage";
	exit;
}

unless($out){
	print STDERR "<ERROR> Output directory missing!\n$usage";
	exit;
}

if(-d $out){
	print STDERR "<WARNNING> Ouput directory $out exists! Will overwrite it...\n";
}
else{
	system("mkdir $out");
}

split_gtf($gtf, $out, $batch, $sort);

sub split_gtf{
	my ($temp_gtf, $temp_out, $temp_batch, $temp_sort) = @_;
	my ($temp_fh);
	if($temp_sort){
		open($temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "<ERROR> Cannot open $temp_gtf!\n";
	}
	else{
		open($temp_fh, "$temp_gtf") or die "<ERROR> Cannot open $temp_gtf!\n";
	}
	my %temp_gene = ();
	my @temp_gene = ();
	my %temp_transcript_count = ();
	my $temp_transcript_count = 0;
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[8] =~ /gene_id \"([^\"]+)\"/){
			my $temp_gene_id = $1;
			unless(exists $temp_gene{$temp_gene_id}){
				@{$temp_gene{$temp_gene_id}} = ();
				push(@temp_gene, $temp_gene_id);
				%{$temp_transcript_count{$temp_gene_id}} = ();
			}
			push(@{$temp_gene{$temp_gene_id}}, $temp_line);
			if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
				my $temp_transcript_id = $1;
				unless(exists $temp_transcript_count{$temp_gene_id}{$temp_transcript_id}){
					$temp_transcript_count{$temp_gene_id}{$temp_transcript_id} = 1;
					$temp_transcript_count++;
				}
			}
			elsif($temp_line[2] ne 'gene'){
				print STDERR "<ERROR> Missing transcript_id in $temp_gtf!\n$temp_line\n";
				exit;
			}
		}
	}
	close $temp_fh;
	my $temp_batch_size = int($temp_transcript_count / $temp_batch) + 1;
	if(($temp_transcript_count / $temp_batch) == 1){
		$temp_batch_size = 1;
	}
	my ($temp_count, $temp_batch_count) = (0, 1);
	open($temp_fh, "> ${temp_out}/batch_${temp_batch_count}.gtf") or die "<ERROR> Cannot create ${temp_out}/batch_${temp_batch_count}.gtf!\n";
	foreach my $temp_gene_id (@temp_gene){
		$temp_count += keys %{$temp_transcript_count{$temp_gene_id}};;
		if($temp_count > ($temp_batch_size * $temp_batch_count)){
			close $temp_fh;
			$temp_batch_count++;
			open($temp_fh, "> ${temp_out}/batch_${temp_batch_count}.gtf");
		}
		foreach my $temp_line (@{$temp_gene{$temp_gene_id}}){
			print $temp_fh "$temp_line\n";
		}
	}
	close $temp_fh;

	return 1;
}



