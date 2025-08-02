#!/usr/bin/env perl
=start
Author: ALin
Purpose: This is a script for querying last exon against a reference. Both query and reference files are in .gtf format. It is part of the Alternative Splicing Event Finder (ASEF) package.
Change log:
	v1.1	2021-10	This script was adopted from ASEF_FE_v1.2.pl
	v1.2	2022-02	It will now calculate the distance to the nearest gene and transcript start sites.
	v1.3	2022-04	Fixed an issue of the strand for the nearestt TTS and gene TTS.
	v1.4	2022-07	Enabled bedtools to handle the strand. Changed intermediate files to .gtf.
	v1.5	2023-03 Added gene level output. Added gene level information of TTS.
	v1.6	2023-12 Fixed an issue in parsing gtf attribute.
	v1.7	2023-12	Added an option of TTS region to be used.
	v1.8	2023-12	Fixed an issue in using bedtools closest by allowing all for tie.
=cut

use strict;
use Getopt::Long;

my $version = "v1.8";

my $usage = "Usage: perl ASEF_LE_${version}.pl
	-q <String> .gtf of query
	-r <String> .gtf of reference
	-o <String> Output base
	-e <String> TTS region bed
	-b <String> Path for bedtools
	-h <Boolean> Help
";

my ($query, $ref, $out, $tts_region, $btp, $help) = ("", "", "", "", "", 0);

GetOptions(
	'q=s'	=>	\$query,
	'r=s'	=>	\$ref,
	'o=s'	=>	\$out,
	'e=s'	=>	\$tts_region,
	'b=s'	=>	\$btp,
	'h!'	=>	\$help
);

unless($query && $ref && $out){
	print "$usage";
	exit;
}

if($btp eq ""){
	$btp = "bedtools";
}

unless(`$btp`){
	print "Incorrect path for bedtools!\n";
	print "$usage";
	exit;
}

if($help){
	print "$usage";
	exit;
}

my $dir = "${out}__ASEF_LE_temp";
if(-d $dir){
	print STDERR "<WARNING> Temporary directory $dir exists! Will overwrite it...\n";
}
else{
	system("mkdir $dir");
}

my $tts_region_ref;
if($tts_region){
	if(-e $tts_region){
		$tts_region_ref = read_bed($tts_region);
	}
	else{
		print STDERR "<ERROR> $tts_region does not exist!\n";
		exit;
	}
}

my $ref_transcript = read_gtf($ref);
parse_ref($ref_transcript, $dir, $btp);
my %query_transcript_le = ();
my %query_gene_tts = ();
my $query_transcript = read_gtf($query, $tts_region_ref, \%query_transcript_le, \%query_gene_tts);
parse_query($query_transcript, $dir, $btp);
compare_LE(\%query_transcript_le, \%query_gene_tts, $dir, $btp);
print_result(\%query_transcript_le, \%query_gene_tts, $out);
system("rm -r $dir");

sub read_bed{
	my ($temp_bed) = @_;
	my %temp_bed = ();
	open(my $temp_fh, $temp_bed) or die "Cannot open $temp_bed!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if(exists $temp_bed{$temp_line[3]}){
			print STDERR "<ERROR> Duplicated transcript, $temp_line[3], in $temp_bed!\n";
			exit;
		}
		else{
			%{$temp_bed{$temp_line[3]}} = ();
			$temp_bed{$temp_line[3]}{'chr'} = $temp_line[0];
			$temp_bed{$temp_line[3]}{'strand'} = $temp_line[5];
			$temp_bed{$temp_line[3]}{'left'} = $temp_line[1];
			$temp_bed{$temp_line[3]}{'right'} = $temp_line[2];
		}
	}
	close $temp_fh;

	return (\%temp_bed);
}

sub read_gtf{
	my ($temp_gtf, $temp_tts_ref, $temp_le_ref, $temp_gene_tts_ref) = @_;
	my %temp_transcript = ();
	open(my $temp_fh, $temp_gtf) or die "Cannot open $temp_gtf!\n'";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		unless($temp_line[2] eq "exon"){
			next;
		}
		my $temp_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"/){
			$temp_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Missing gene_id in $temp_gtf!\n$temp_line\n";
			exit;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Missing transcript_id in $temp_gtf!\n$temp_line\n";
			exit;
		}
		if(defined $temp_le_ref){
			unless(exists $temp_le_ref->{$temp_transcript_id}){
				$temp_le_ref->{$temp_transcript_id}{'type'} = "novel";
				$temp_le_ref->{$temp_transcript_id}{'nearest_TTS'} = 999999999;
				%{$temp_le_ref->{$temp_transcript_id}{'nearest_TTS_gene_id'}} = ();
				$temp_le_ref->{$temp_transcript_id}{'nearest_gene_TTS'} = 999999999;
				%{$temp_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}} = ();
			}
		}
		if(defined $temp_gene_tts_ref){
			unless(exists $temp_gene_tts_ref->{$temp_gene_id}){
				%{$temp_gene_tts_ref->{$temp_gene_id}} = ();
				$temp_gene_tts_ref->{$temp_gene_id}{'nearest_gene_TTS'} = 999999999;
				%{$temp_gene_tts_ref->{$temp_gene_id}{'nearest_gene_TTS_gene_id'}} = ();
			}
		}
		unless(exists $temp_transcript{$temp_line[0]}){
			%{$temp_transcript{$temp_line[0]}} = ();
		}
		unless(exists $temp_transcript{$temp_line[0]}{$temp_line[6]}){
			%{$temp_transcript{$temp_line[0]}{$temp_line[6]}} = ();
		}
		unless(exists $temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}){
			%{$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}} = ();
		}
		unless(exists $temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}){
			%{$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}} = ();
			@{$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'exon'}} = ();
		}
		my @temp_interval = ($temp_line[3], $temp_line[4]);
		push(@{$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'exon'}}, \@temp_interval);
		if(defined $temp_tts_ref){
			if(exists $temp_tts_ref->{$temp_transcript_id}){
				$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'tts_left'} = $temp_tts_ref->{$temp_transcript_id}{'left'};
				$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'tts_right'} = $temp_tts_ref->{$temp_transcript_id}{'right'};
			}
		}
	}
	close $temp_fh;

	return \%temp_transcript;
}

sub parse_ref{
	my ($temp_ref_transcript_ref, $temp_dir, $temp_btp) = @_;
	my %temp_ref = ();
	open(my $temp_le_fh, "| $temp_btp sort -i > ${temp_dir}/ref_le.gtf") or die "Cannot create ${temp_dir}/ref_le.gtf!\n";
	open(my $temp_tts_fh, "| $temp_btp sort -i > ${temp_dir}/ref_tts.gtf") or die "Cannot create ${temp_dir}/ref_tts.gtf!\n";
	open(my $temp_gene_tts_fh, "| $temp_btp sort -i > ${temp_dir}/ref_gene_tts.gtf") or die "Cannot create ${temp_dir}/ref_gene_tts.gtf!\n";
	foreach my $temp_chr (keys %$temp_ref_transcript_ref){
		unless(exists $temp_ref{$temp_chr}){
			%{$temp_ref{$temp_chr}} = ();
		}
		foreach my $temp_strand (keys %{$temp_ref_transcript_ref->{$temp_chr}}){
			unless(exists $temp_ref{$temp_chr}{$temp_strand}){
				%{$temp_ref{$temp_chr}{$temp_strand}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'le'}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'tts'}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}} = ();
			}
			foreach my $temp_gene_id (keys %{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}}){
				my ($temp_gene_tts_left, $temp_gene_tts_right) = (0, 0);
				foreach my $temp_transcript_id (keys %{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}}){
					my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'exon'}};
					my ($temp_tts_left, $temp_tts_right);
					my $temp_last_exon;
					if($temp_strand eq "-"){
						if(exists $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'}){
							$temp_tts_left = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'};
							$temp_tts_right = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_right'};
						}
						else{
							$temp_tts_left = $temp_sorted_exon[0][0];
							$temp_tts_right = $temp_sorted_exon[0][0];
						}
						$temp_last_exon = "$temp_sorted_exon[0][0]-$temp_sorted_exon[0][1]";
					}
					else{
						if(exists $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'}){
							$temp_tts_left = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'};
							$temp_tts_right = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_right'};
						}
						else{
							$temp_tts_left = $temp_sorted_exon[@temp_sorted_exon-1][1];
							$temp_tts_right = $temp_sorted_exon[@temp_sorted_exon-1][1];
						}
						$temp_last_exon = "$temp_sorted_exon[@temp_sorted_exon-1][0]-$temp_sorted_exon[@temp_sorted_exon-1][1]";
					}
					if(($temp_gene_tts_left == 0) || ($temp_gene_tts_left > $temp_tts_left)){
						$temp_gene_tts_left = $temp_tts_left;
					}
					if(($temp_gene_tts_right == 0) || ($temp_gene_tts_right < $temp_tts_right)){
						$temp_gene_tts_right = $temp_tts_right;
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}{$temp_gene_id}){
						$temp_ref{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}{$temp_gene_id} = 1;
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}{$temp_gene_id}){
						$temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}{$temp_gene_id} = 1;
					}
					delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id};
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}){
					%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}} = ();
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}){
					%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}} = ();
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}{$temp_gene_id}){
					$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}{$temp_gene_id} = 1;
				}
				delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id};
			}
			delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand};
		}
		delete $temp_ref_transcript_ref->{$temp_chr};
	}
	foreach my $temp_chr (keys %temp_ref){
		foreach my $temp_strand (keys %{$temp_ref{$temp_chr}}){
			foreach my $temp_last_exon (keys %{$temp_ref{$temp_chr}{$temp_strand}{'le'}}){
				my @temp_last_exon = split('-', $temp_last_exon);
				my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}});
				print $temp_le_fh "$temp_chr\tASEF\texon\t$temp_last_exon[0]\t$temp_last_exon[1]\t0\t$temp_strand\t.\texon_id \"$temp_chr:$temp_last_exon[0]-$temp_last_exon[1]:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
			}
			foreach my $temp_tts_left (keys %{$temp_ref{$temp_chr}{$temp_strand}{'tts'}}){
				foreach my $temp_tts_right (keys %{$temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}}){
					my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_left}});
					print $temp_tts_fh "$temp_chr\tASEF\tTTS\t$temp_tts_left\t$temp_tts_right\t0\t$temp_strand\t.\ttts_id \"$temp_chr:$temp_tts_left-$temp_tts_right:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
				}
			}
			foreach my $temp_gene_tts_left (keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}}){
				foreach my $temp_gene_tts_right (keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}}){
					my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}});
					print $temp_gene_tts_fh "$temp_chr\tASEF\tTTS\t$temp_gene_tts_left\t$temp_gene_tts_right\t0\t$temp_strand\t.\ttts_id \"$temp_chr:$temp_gene_tts_left-$temp_gene_tts_right:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
					delete $temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right};
				}
				delete $temp_ref{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left};
			}
			delete $temp_ref{$temp_chr}{$temp_strand};
		}
		delete $temp_ref{$temp_chr};
	}
	close $temp_le_fh;
	close $temp_tts_fh;
	close $temp_gene_tts_fh;

	return 1;
}

sub parse_query{
	my ($temp_query_transcript_ref, $temp_dir, $temp_btp) = @_;
	my %temp_query = ();
	open(my $temp_le_fh, "| $temp_btp sort -i > ${temp_dir}/query_le.gtf") or die "Cannot create ${temp_dir}/query_le.gtf!\n";
	open(my $temp_tts_fh, "| $temp_btp sort -i > ${temp_dir}/query_tts.gtf") or die "Cannot create ${temp_dir}/query_tts.gtf!\n";
	open(my $temp_gene_tts_fh, "| $temp_btp sort -i > ${temp_dir}/query_gene_tts.gtf") or die "Cannot create ${temp_dir}/query_gene_tts.gtf!\n";
	foreach my $temp_chr (keys %$temp_query_transcript_ref){
		unless(exists $temp_query{$temp_chr}){
			%{$temp_query{$temp_chr}} = ();
		}
		foreach my $temp_strand (keys %{$temp_query_transcript_ref->{$temp_chr}}){
			unless(exists $temp_query{$temp_chr}{$temp_strand}){
				%{$temp_query{$temp_chr}{$temp_strand}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'le'}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'tts'}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}} = ();
			}
			foreach my $temp_gene_id (keys %{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}}){
				my ($temp_gene_tts_left, $temp_gene_tts_right) = (0, 0);
				foreach my $temp_transcript_id (keys %{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}}){
					my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'exon'}};
					my ($temp_tts_left, $temp_tts_right);
					my $temp_last_exon;
					if($temp_strand eq "-"){
						if(exists $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'}){
							$temp_tts_left = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'};
							$temp_tts_right = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_right'};
						}
						else{
							$temp_tts_left = $temp_sorted_exon[0][0];
							$temp_tts_right = $temp_sorted_exon[0][0];
						}
						$temp_last_exon = "$temp_sorted_exon[0][0]-$temp_sorted_exon[0][1]";
					}
					else{
						if(exists $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'}){
							$temp_tts_left = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_left'};
							$temp_tts_right = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tts_right'};
						}
						else{
							$temp_tts_left = $temp_sorted_exon[@temp_sorted_exon-1][1];
							$temp_tts_right = $temp_sorted_exon[@temp_sorted_exon-1][1];
						}
						$temp_last_exon = "$temp_sorted_exon[@temp_sorted_exon-1][0]-$temp_sorted_exon[@temp_sorted_exon-1][1]";
					}
					if(($temp_gene_tts_left == 0) || ($temp_gene_tts_left > $temp_tts_left)){
						$temp_gene_tts_left = $temp_tts_left;
					}
					if(($temp_gene_tts_right == 0) || ($temp_gene_tts_right < $temp_tts_right)){
						$temp_gene_tts_right = $temp_tts_right;
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}){
						%{$temp_query{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}{$temp_transcript_id}){
						$temp_query{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}{$temp_transcript_id} = 1;
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}){
						%{$temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}){
						%{$temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}{$temp_transcript_id}){
						$temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}{$temp_transcript_id} = 1;
					}
					delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id};
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}){
					%{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}} = ();
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}){
					%{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}} = ();
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}{$temp_gene_id}){
					$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}{$temp_gene_id} = 1;
				}
				delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id};
			}
			delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand};
		}
		delete $temp_query_transcript_ref->{$temp_chr};
	}
	foreach my $temp_chr (keys %temp_query){
		foreach my $temp_strand (keys %{$temp_query{$temp_chr}}){
			foreach my $temp_last_exon (keys %{$temp_query{$temp_chr}{$temp_strand}{'le'}}){
				my @temp_last_exon = split('-', $temp_last_exon);
				my $temp_transcript_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'le'}{$temp_last_exon}});
				print $temp_le_fh "$temp_chr\tASEF\texon\t$temp_last_exon[0]\t$temp_last_exon[1]\t0\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\;\n";
				delete $temp_query{$temp_chr}{$temp_strand}{$temp_last_exon};
			}
			foreach my $temp_tts_left (keys %{$temp_query{$temp_chr}{$temp_strand}{'tts'}}){
				foreach my $temp_tts_right (keys %{$temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}}){
					my $temp_transcript_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right}});
					print $temp_tts_fh "$temp_chr\tASEF\texon\t$temp_tts_left\t$temp_tts_right\t0\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\;\n";
					delete $temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left}{$temp_tts_right};
				}
				delete $temp_query{$temp_chr}{$temp_strand}{'tts'}{$temp_tts_left};
			}
			foreach my $temp_gene_tts_left (keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}}){
				foreach my $temp_gene_tts_right (keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}}){
					my $temp_gene_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right}});
					print $temp_gene_tts_fh "$temp_chr\tASEF\tTTS\t$temp_gene_tts_left\t$temp_gene_tts_right\t0\t$temp_strand\t.\tgene_id \"$temp_gene_id\"\;\n";
					delete $temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left}{$temp_gene_tts_right};
				}
				delete $temp_query{$temp_chr}{$temp_strand}{'gene_tts'}{$temp_gene_tts_left};
			}
			delete $temp_query{$temp_chr}{$temp_strand};
		}
		delete $temp_query{$temp_chr};
	}
	close $temp_le_fh;
	close $temp_tts_fh;
	close $temp_gene_tts_fh;

	return 1;
}

sub compare_LE{
	my ($temp_query_le_ref, $temp_query_gene_tts_ref, $temp_dir, $temp_btp) = @_;
	open(my $temp_fh, "$temp_btp intersect -wa -s -nonamecheck -a ${temp_dir}/query_le.gtf -b ${temp_dir}/ref_le.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		#print "$temp_line\n";
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id = ();
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_transcript_id (@temp_transcript_id){
				$temp_query_le_ref->{$temp_transcript_id}{'type'} = "overlapping";
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_tts.gtf -b ${temp_dir}/ref_tts.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id = ();
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		my @temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"/){
			@temp_ref_gene_id = split(',', $2);
		}
		elsif($temp_line[18] != -1){
			print "Cannot find gene_id at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_transcript_id (@temp_transcript_id){
			$temp_query_le_ref->{$temp_transcript_id}{'nearest_TTS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_tts.gtf -b ${temp_dir}/ref_gene_tts.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id = ();
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		my @temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"/){
			@temp_ref_gene_id = split(',', $2);
		}
		elsif($temp_line[18] != -1){
			print "Cannot find gene_id at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_transcript_id (@temp_transcript_id){
			$temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_gene_tts.gtf -b ${temp_dir}/ref_gene_tts.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"/){
			@temp_query_gene_id = split(',', $2);
		}
		else{
			print "Cannot find gene_id at line:\n$temp_line\n";
			exit;
		}
		my @temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"/){
			@temp_ref_gene_id = split(',', $2);
		}
		elsif($temp_line[18] != -1){
			print "Cannot find gene_id at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_query_gene_id (@temp_query_gene_id){
			$temp_query_gene_tts_ref->{$temp_query_gene_id}{'nearest_gene_TTS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_gene_tts_ref->{$temp_query_gene_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_gene_tts_ref->{$temp_query_gene_id}{'nearest_gene_TTS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;

	return 1;
}

sub print_result{
	my ($temp_query_le_ref, $temp_query_gene_tts_ref, $temp_out) = @_;
	open(my $temp_transcript_fh, ">${temp_out}_transcript.txt") or die "Cannot create ${temp_out}_transcript.txt!\n";
	open(my $temp_gene_fh, ">${temp_out}_gene.txt") or die "Cannot create ${temp_out}_gene.txt!\n";
	print $temp_transcript_fh "transcript_id\tLE\tnearest_TTS\tnearest_TTS_same_gene\tnearest_gene_TTS\tnearest_gene_TTS_same_gene\n";
	foreach my $temp_transcript_id (keys %$temp_query_le_ref){
		my $temp_nearest_TTS_same_gene = 0;
		my $temp_nearest_gene_TTS_same_gene = 0;
		if(exists $temp_query_le_ref->{$temp_transcript_id}{'nearest_TTS_gene_id'}{$temp_query_le_ref->{$temp_transcript_id}{'gene_id'}}){
			$temp_nearest_TTS_same_gene = 1;
		}
		if(exists $temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS_gene_id'}{$temp_query_le_ref->{$temp_transcript_id}{'gene_id'}}){
			$temp_nearest_gene_TTS_same_gene = 1;
		}
		print $temp_transcript_fh "$temp_transcript_id\t$temp_query_le_ref->{$temp_transcript_id}{'type'}\t$temp_query_le_ref->{$temp_transcript_id}{'nearest_TTS'}\t$temp_nearest_TTS_same_gene\t$temp_query_le_ref->{$temp_transcript_id}{'nearest_gene_TTS'}\t$temp_nearest_gene_TTS_same_gene\n";
	}
	print $temp_gene_fh "gene_id\tnearest_gene_TTS\tnearest_gene_TTS_same_gene\n";
	foreach my $temp_gene_id (keys %$temp_query_gene_tts_ref){
		my $temp_nearest_gene_TTS_same_gene = 0;
		if(exists $temp_query_gene_tts_ref->{$temp_gene_id}{'nearest_gene_TTS_gene_id'}{$temp_gene_id}){
			$temp_nearest_gene_TTS_same_gene = 1;
		}
		print $temp_gene_fh "$temp_gene_id\t$temp_query_gene_tts_ref->{$temp_gene_id}{'nearest_gene_TTS'}\t$temp_nearest_gene_TTS_same_gene\n";
	}
	close $temp_transcript_fh;
	close $temp_gene_fh;

	return 1;
}



