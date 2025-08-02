#!/usr/bin/env perl
#Author: ALin
#Purpose: Parse the output from nmd_splice_junction_finder_v*.pl to find the most prioritized splice junctions of NMD based on the following criterion: 5' most in CDS > 3UTR_NMD > 5' most in 3UTR > 5UTR.

use strict;

my $version = "v1.1";

my $usage = "Usage: perl nmd_splice_junction_parser_${version}.pl <input> <output>\n";

if(@ARGV != 2){
	print STDERR "<ERROR> $usage";
	exit;
}

my ($in, $out) = @ARGV;

my $sj_ref = read_sj($in);
parse_sj($sj_ref, $out);

sub read_sj{
	my ($temp_in) = @_;
	my %temp_sj = ();
	my $temp_fh;
	open($temp_fh, $temp_in) or die "Cannot open $temp_in!\n";
	my $temp_header = 0;
	my %temp_header = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if($temp_header == 0){
			for(my $i = 0; $i < @temp_line; $i++){
				if(exists $temp_header{$temp_line[$i]}){
					print STDERR "<ERROR> Duplicated column $temp_line[$i] in $temp_in!\n";
					exit;
				}
				else{
					$temp_header{$temp_line[$i]} = $i;
				}
			}
			$temp_header = 1;
		}
		else{
			my $temp_transcript_id = $temp_line[$temp_header{'transcript_id'}];
			my $temp_sj_id = $temp_line[$temp_header{'sj_id'}];
			my $temp_origin = $temp_line[$temp_header{'origin'}];
			unless(exists $temp_sj{$temp_transcript_id}){
				%{$temp_sj{$temp_transcript_id}} = ();
			}
			unless(exists $temp_sj{$temp_transcript_id}{$temp_origin}){
				@{$temp_sj{$temp_transcript_id}{$temp_origin}} = ();
			}
			push(@{$temp_sj{$temp_transcript_id}{$temp_origin}}, $temp_sj_id);
		}
	}
	close $temp_fh;

	return \%temp_sj;
}

sub parse_sj{
	my ($temp_sj_ref, $temp_out) = @_;
	my $temp_fh;
	open($temp_fh, "> $temp_out") or die "Cannot create $temp_out!\n";
	print $temp_fh "transcript_id\tsj_id\tdruggable_origin\n";
	foreach my $temp_transcript_id (keys %$temp_sj_ref){
		my ($temp_5_most_sj_id, $temp_type) = ("", "");
		if(exists $temp_sj_ref->{$temp_transcript_id}->{'CDS'}){
			$temp_5_most_sj_id = get_5_most_sj($temp_sj_ref->{$temp_transcript_id}->{'CDS'});
			$temp_type = "CDS";
		}
		elsif(exists $temp_sj_ref->{$temp_transcript_id}->{'3UTR_NMD'}){
			$temp_5_most_sj_id = get_5_most_sj($temp_sj_ref->{$temp_transcript_id}->{'3UTR_NMD'});
			$temp_type = "3UTR_NMD";
		}
		elsif(exists $temp_sj_ref->{$temp_transcript_id}->{'3UTR'}){
			$temp_5_most_sj_id = get_5_most_sj($temp_sj_ref->{$temp_transcript_id}->{'3UTR'});
			$temp_type = "3UTR";
		}
		elsif(exists $temp_sj_ref->{$temp_transcript_id}->{'5UTR'}){
			$temp_5_most_sj_id = get_5_most_sj($temp_sj_ref->{$temp_transcript_id}->{'5UTR'});
			$temp_type = "5UTR";
		}
		else{
			print STDERR "<ERROR> Unrecognized origin in $temp_transcript_id!\n";
			exit;
		}
		print $temp_fh "$temp_transcript_id\t$temp_5_most_sj_id\t$temp_type\n";
	}
	close $temp_fh;

	return 1;
}

sub get_5_most_sj{
	my ($temp_sj_id_ref) = @_;
	my $temp_transcript_strand = "";
	my ($temp_5_end, $temp_5_most_sj_id) = (0, "");
	foreach my $temp_sj_id (@$temp_sj_id_ref){
		if($temp_sj_id =~ /([^_]+)_(\d+)_(\d+)_([+-])/){
			my ($temp_chr, $temp_left, $temp_right, $temp_strand) = ($1, $2, $3, $4);
			if($temp_transcript_strand eq ""){
				$temp_transcript_strand = $temp_strand;
			}
			elsif($temp_transcript_strand ne $temp_strand){
				print STDERR "<ERROR> Incorrect strand in $temp_sj_id\n";
				exit;
			}
			if($temp_transcript_strand eq '-'){
				if($temp_5_end == 0){
					$temp_5_end = $temp_right;
					$temp_5_most_sj_id = $temp_sj_id;
				}
				elsif($temp_5_end < $temp_right){
					$temp_5_end = $temp_right;
					$temp_5_most_sj_id = $temp_sj_id;
				}
			}
			else{
				if($temp_5_end == 0){
					$temp_5_end = $temp_left;
					$temp_5_most_sj_id = $temp_sj_id;
				}
				elsif($temp_5_end > $temp_left){
					$temp_5_end = $temp_left;
					$temp_5_most_sj_id = $temp_sj_id;
				}
			}
		}
		else{
			print STDERR "<ERROR> Incorrect strand in $temp_sj_id\n";
			exit;
		}
	}

	return $temp_5_most_sj_id;
}


