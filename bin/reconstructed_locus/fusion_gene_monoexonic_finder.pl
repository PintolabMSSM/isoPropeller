#!/usr/bin/env perl
=start
Author: ALin
Purpose: To find fusion genes associated with mono-exonic transcripts based on the transcript information table.
Change log:
=cut

use strict;

my $version = "v1.1";

my $usage = "perl fusion_gene_monoexonic_finder.pl <transcript table> <output>\n";

if(@ARGV != 2){
	print "$usage";
	exit;
}

my ($transcript, $out) = @ARGV;

my ($gene_ref, $fusion_ref) = read_transcript($transcript);

parse_transcript($gene_ref, $fusion_ref, $out);

#Read the transcript table
sub read_transcript{
	my ($temp_transcript) = @_;
	open(my $temp_fh, $temp_transcript) or die "Cannot open $temp_transcript!\n";
	my %temp_gene = ();
	my %temp_fusion = ();
	my %temp_header = ();
	my $temp_header = 0;
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if($temp_header == 0){
			for(my $i = 0 ; $i < @temp_line; $i++){
				if(exists $temp_header{$temp_line[$i]}){
					print STDERR "<ERROR> Duplicated column $temp_line[$i] in $temp_transcript!\n";
					exit;
				}
				else{
					$temp_header{$temp_line[$i]} = $i;
				}
			}
			unless(exists $temp_header{'gene_id'}){
				print STDERR "<ERROR> Missing gene_id column in $temp_transcript!\n";
				exit;
			}
			unless(exists $temp_header{'exon'}){
				print STDERR "<ERROR> Missing exon column in $temp_transcript!\n";
				exit;
			}
			unless(exists $temp_header{'status'}){
				print STDERR "<ERROR> Missing status column in $temp_transcript!\n";
				exit;
			}
			$temp_header = 1;
		}
		else{
			my $temp_gene_id = $temp_line[$temp_header{'gene_id'}];
			if($temp_line[$temp_header{'status'}] eq 'fusion'){
				unless(exists $temp_fusion{$temp_gene_id}){
					$temp_fusion{$temp_gene_id} = 1;
				}
			}
			else{
				unless(exists $temp_gene{$temp_gene_id}){
					%{$temp_gene{$temp_gene_id}} = ();
					$temp_gene{$temp_gene_id}{'num_transcript'} = 0;
					$temp_gene{$temp_gene_id}{'num_monoexonic'} = 0;
				}
				$temp_gene{$temp_gene_id}{'num_transcript'}++;
				if($temp_line[$temp_header{'exon'}] == 1){
					$temp_gene{$temp_gene_id}{'num_monoexonic'}++;
				}
			}
		}
	}
	close $temp_fh;

	return (\%temp_gene, \%temp_fusion);
}

#Parse the transcripts and output the result
sub parse_transcript{
	my ($temp_gene_ref, $temp_fusion_ref, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "Cannot create $temp_out!\n";
	foreach my $temp_gene_id (keys %$temp_fusion_ref){
		my @temp_gene_id = split('\|', $temp_gene_id);
		my $temp_num_monoexonic = 0;
		foreach my $temp_component_gene_id (@temp_gene_id){
			if($temp_gene_ref->{$temp_component_gene_id}{'num_transcript'} == $temp_gene_ref->{$temp_component_gene_id}{'num_monoexonic'}){
				$temp_num_monoexonic++;
			}
		}
		if($temp_num_monoexonic > 0){
			print $temp_fh "$temp_gene_id\n";
		}
	}
	close $temp_fh;

	return 1;
}






