#!/usr/bin/env perl
=start
Author: ALin
Purpose: This is a script for querying first exon against a reference. Both query and reference files are in .gtf format. It is part of the Alternative Splicing Event Finder (ASEF) package.
Change log:
	v1.2	2021-10	Unlike the v1.1, which only looked at the first donor site, v1.2 look at the overlap between the first exons. If an query first exon does not overlap with any reference first exon, it will be an potential novel first exon with a novel FE.
	v1.3	2021-10	In v1.3, it allows users to specify the distance between FEs for bi-directional novel FEs, in both up- and down-stream of known FEs. 
	v1.4	2022-02	It will now calculate the distance to the nearest gene and transcript start sites.
	v1.5	2022-04	An issue of the strand for the nearest TTS and gene TTS was fixed.
	v1.6	2022-07	Separated the steps for finding overlapping and bi-directional cases. Enabled bedtools to handle the strand. Changed the intermediate file to .gtf.
	v1.7	2023-03	Added gene level output. Added gene level information of TSS.
	v1.8	2023-03	Fixed a bug in getting the query first exon in the - strand.
	v1.9	2023-12 Fixed an issue in parsing gtf attribute.
	v1.10	2023-12	Added an option of TSS region to be used.
	v1.11	2023-12	Fixed an issue in using bedtools closest by allowing all for tie.
	v1.12	2025-08	Fixed a bug when handling the left position of input TSS bed file.
=cut

use strict;
use Getopt::Long;

my $version = "v1.12";

my $usage = "perl ASEF_FE_${version}.pl
	-q <String> .gtf of query
	-r <String> .gtf of reference
	-d <Integer> Maximal distance defining a bi-directional novel first exon
	-o <String> Output base
	-e <String> TSS region bed
	-b <String> Path for bedtools
	-h <Boolean> Help
";

my ($query, $ref, $bi_dist, $out, $tss_region, $btp, $help) = ("", "", 0, "", "", "", 0);

GetOptions(
	'q=s'	=>	\$query,
	'r=s'	=>	\$ref,
	'd=i'	=>	\$bi_dist,
	'o=s'	=>	\$out,
	'e=s'	=>	\$tss_region,
	'b=s'	=>	\$btp,
	'h!'	=>	\$help,
);

unless($query && $ref && $out){
	print STDERR "$usage";
	exit;
}

if($btp eq ""){
	$btp = "bedtools";
}

unless(`$btp`){
	print STDERR "<ERROR> Incorrect path for bedtools!\n$usage";
	exit;
}

if($help){
	print STDERR "$usage";
	exit;
}

my $dir = "${out}__ASEF_FE_temp";
if(-d $dir){
	print STDERR "<WARNING> Temporary directory $dir exists! Will overwrite it...\n";
}
else{
	system("mkdir $dir");
}

my $tss_region_ref;
if($tss_region){
	if(-e $tss_region){
		$tss_region_ref = read_bed($tss_region);
	}
	else{
		print STDERR "<ERROR> $tss_region does not exist!\n";
		exit;
	}
}

my $ref_transcript_ref = read_gtf($ref);
parse_ref($ref_transcript_ref, $dir, $btp);
my %query_transcript_fe = ();
my %query_gene_tss = ();
my $query_transcript_ref = read_gtf($query, $tss_region_ref, \%query_transcript_fe, \%query_gene_tss);
parse_query($query_transcript_ref, $dir, $btp);
compare_FE(\%query_transcript_fe, \%query_gene_tss, $dir, $btp);
print_result(\%query_transcript_fe, \%query_gene_tss, $out);
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
			$temp_bed{$temp_line[3]}{'left'} = $temp_line[1] + 1;
			$temp_bed{$temp_line[3]}{'right'} = $temp_line[2];
		}
	}
	close $temp_fh;

	return (\%temp_bed);
}

sub read_gtf{
	my ($temp_gtf, $temp_tss_ref, $temp_fe_ref, $temp_gene_tss_ref) = @_;
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
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\;]+)\"/){
			$temp_gene_id = $2;
		}
		else{
			print STDERR "Mising gene_id in $temp_gtf!\n$temp_line\n";
			exit;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\;]+)\"/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "Mising transcript_id in $temp_gtf!\n$temp_line\n";
			exit;
		}
		if(defined $temp_fe_ref){
			unless(exists $temp_fe_ref->{$temp_transcript_id}){
				%{$temp_fe_ref->{$temp_transcript_id}} = ();
				$temp_fe_ref->{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
				$temp_fe_ref->{$temp_transcript_id}{'type'} = "novel_orphan";
				$temp_fe_ref->{$temp_transcript_id}{'nearest_TSS'} = 999999999;
				%{$temp_fe_ref->{$temp_transcript_id}{'nearest_TSS_gene_id'}} = ();
				$temp_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS'} = 999999999;
				%{$temp_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS_gene_id'}} = ();
			}
		}
		if(defined $temp_gene_tss_ref){
			unless(exists $temp_gene_tss_ref->{$temp_gene_id}){
				%{$temp_gene_tss_ref->{$temp_gene_id}} = ();
				$temp_gene_tss_ref->{$temp_gene_id}{'nearest_gene_TSS'} = 999999999;
				%{$temp_gene_tss_ref->{$temp_gene_id}{'nearest_gene_TSS_gene_id'}} = ();
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
		if(defined $temp_tss_ref){
			if(exists $temp_tss_ref->{$temp_transcript_id}){
				$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'tss_left'} = $temp_tss_ref->{$temp_transcript_id}{'left'};
				$temp_transcript{$temp_line[0]}{$temp_line[6]}{$temp_gene_id}{$temp_transcript_id}{'tss_right'} = $temp_tss_ref->{$temp_transcript_id}{'right'};
			}
		}
	}
	close $temp_fh;

	return \%temp_transcript;
}

sub parse_ref{
	my ($temp_ref_transcript_ref, $temp_dir, $temp_btp) = @_;
	my %temp_ref = ();
	open(my $temp_fe_fh, "| $temp_btp sort -i > ${temp_dir}/ref_fe.gtf") or die "Cannot create ${temp_dir}/ref_fe.gtf!\n";
	open(my $temp_fee_fh, "| $temp_btp sort -i > ${temp_dir}/ref_fee.gtf") or die "Cannot create ${temp_dir}/ref_fee.gtf!\n";
	open(my $temp_tss_fh, "| $temp_btp sort -i > ${temp_dir}/ref_tss.gtf") or die "Cannot create ${temp_dir}/ref_tss.gtf!\n";
	open(my $temp_gene_tss_fh, "| $temp_btp sort -i > ${temp_dir}/ref_gene_tss.gtf") or die "Cannot create ${temp_dir}/ref_gene_tss.gtf!\n";
	foreach my $temp_chr (keys %$temp_ref_transcript_ref){
		unless(exists $temp_ref{$temp_chr}){
			%{$temp_ref{$temp_chr}} = ();
		}
		foreach my $temp_strand (keys %{$temp_ref_transcript_ref->{$temp_chr}}){
			unless(exists $temp_ref{$temp_chr}{$temp_strand}){
				%{$temp_ref{$temp_chr}{$temp_strand}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'fe'}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'tss'}} = ();
				%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}} = ();
			}
			foreach my $temp_gene_id (keys %{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}}){
				my ($temp_gene_tss_left, $temp_gene_tss_right) = (0, 0);
				foreach my $temp_transcript_id (keys %{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}}){
					my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'exon'}};
					my ($temp_tss_left, $temp_tss_right);
					my $temp_first_exon;
					if($temp_strand eq "-"){
						if(exists $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'}){
							$temp_tss_left = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'};
							$temp_tss_right = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_right'};
						}
						else{
							$temp_tss_left = $temp_sorted_exon[@temp_sorted_exon-1][1];
							$temp_tss_right = $temp_sorted_exon[@temp_sorted_exon-1][1];
						}
						$temp_first_exon = "$temp_sorted_exon[@temp_sorted_exon-1][0]-$temp_sorted_exon[@temp_sorted_exon-1][1]";
					}
					else{
						if(exists $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'}){
							$temp_tss_left = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'};
							$temp_tss_right = $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_right'};
						}
						else{
							$temp_tss_left = $temp_sorted_exon[0][0];
							$temp_tss_right = $temp_sorted_exon[0][0];
						}
						$temp_first_exon = "$temp_sorted_exon[0][0]-$temp_sorted_exon[0][1]";
					}
					if(($temp_gene_tss_left == 0) || ($temp_gene_tss_left > $temp_tss_left)){
						$temp_gene_tss_left = $temp_tss_left;
					}
					if(($temp_gene_tss_right == 0) || ($temp_gene_tss_right < $temp_tss_right)){
						$temp_gene_tss_right = $temp_tss_right;
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}{$temp_gene_id}){
						$temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}{$temp_gene_id} = 1;
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}){
						%{$temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}} = ();
					}
					unless(exists $temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}{$temp_gene_id}){
							$temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}{$temp_gene_id} = 1;
					}
					delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id};
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}){
					%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}} = ();
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}){
					%{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}} = ();
				}
				unless(exists $temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}{$temp_gene_id}){
					$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}{$temp_gene_id} = 1;
				}
				delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id};
			}
			delete $temp_ref_transcript_ref->{$temp_chr}{$temp_strand};
		}
		delete $temp_ref_transcript_ref->{$temp_chr};
	}
	foreach my $temp_chr (keys %temp_ref){
		foreach my $temp_strand (keys %{$temp_ref{$temp_chr}}){
			foreach my $temp_first_exon (keys %{$temp_ref{$temp_chr}{$temp_strand}{'fe'}}){
				my @temp_first_exon = split('-', $temp_first_exon);
				my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}});
				print $temp_fe_fh "$temp_chr\tASEF\texon\t$temp_first_exon[0]\t$temp_first_exon[1]\t0\t$temp_strand\t.\texon_id \"$temp_chr:$temp_first_exon[0]-$temp_first_exon[1]:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
				if($bi_dist > 0){
					my $temp_exon_length = $temp_first_exon[1] - $temp_first_exon[0] + 1;
					if($temp_strand eq "+"){
						$temp_first_exon[0] -= $bi_dist;
						$temp_first_exon[1] = $temp_first_exon[0] + (2 * $bi_dist) - 1;
					}
					elsif($temp_strand eq "-"){
						$temp_first_exon[1] += $bi_dist;
						$temp_first_exon[0] = $temp_first_exon[1] - (2 * $bi_dist) + 1;
					}
					if($temp_first_exon[0] <= 0){
						$temp_first_exon[0] = 1;
					}
					print $temp_fee_fh "$temp_chr\tASEF\texon\t$temp_first_exon[0]\t$temp_first_exon[1]\t0\t$temp_strand\t.\texon_id \"$temp_chr:$temp_first_exon[0]-$temp_first_exon[1]:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
				}
				delete $temp_ref{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon};
			}
			foreach my $temp_tss_left (keys %{$temp_ref{$temp_chr}{$temp_strand}{'tss'}}){
				foreach my $temp_tss_right (keys %{$temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}}){
					my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}});
					print $temp_tss_fh "$temp_chr\tASEF\tTSS\t$temp_tss_left\t$temp_tss_right\t0\t$temp_strand\t.\ttss_id \"$temp_chr:$temp_tss_left-$temp_tss_right:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
					delete $temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right};
				}
				delete $temp_ref{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left};
			}
			foreach my $temp_gene_tss_left (keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}}){
				foreach my $temp_gene_tss_right (keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}}){
					my $temp_gene_id = join(',', keys %{$temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}});
					print $temp_gene_tss_fh "$temp_chr\tASEF\tTSS\t$temp_gene_tss_left\t$temp_gene_tss_right\t0\t$temp_strand\t.\ttss_id \"$temp_chr:$temp_gene_tss_left-$temp_gene_tss_right:$temp_strand\"\; gene_id \"$temp_gene_id\"\;\n";
					delete $temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right};
				}
				delete $temp_ref{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left};
			}
			delete $temp_ref{$temp_chr}{$temp_strand};
		}
		delete $temp_ref{$temp_chr};
	}
	close $temp_fe_fh;
	close $temp_fee_fh;
	close $temp_gene_tss_fh;
	close $temp_tss_fh;

	return 1;
}

sub parse_query{
	my ($temp_query_transcript_ref, $temp_dir, $temp_btp) = @_;
	my %temp_query = ();
	open(my $temp_fe_fh, "| $temp_btp sort -i > ${temp_dir}/query_fe.gtf") or die "Cannot create ${temp_dir}/query_fe.gtf!\n";
	open(my $temp_tss_fh, "| $temp_btp sort -i > ${temp_dir}/query_tss.gtf") or die "Cannot create ${temp_dir}/query_tss.gtf!\n";
	open(my $temp_gene_tss_fh, "| $temp_btp sort -i > ${temp_dir}/query_gene_tss.gtf") or die "Cannot create ${temp_dir}/query_gene_tss.gtf!\n";
	foreach my $temp_chr (keys %$temp_query_transcript_ref){
		unless(exists $temp_query{$temp_chr}){
			%{$temp_query{$temp_chr}} = ();
		}
		foreach my $temp_strand (keys %{$temp_query_transcript_ref->{$temp_chr}}){
			unless(exists $temp_query{$temp_chr}{$temp_strand}){
				%{$temp_query{$temp_chr}{$temp_strand}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'fe'}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'tss'}} = ();
				%{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}} = ();
			}
			foreach my $temp_gene_id (keys %{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}}){
				my ($temp_gene_tss_left, $temp_gene_tss_right) = (0, 0);
				foreach my $temp_transcript_id (keys %{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}}){
					my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'exon'}};
					my ($temp_tss_left, $temp_tss_right);
					my $temp_first_exon;
					if($temp_strand eq "-"){
						if(exists $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'}){
							$temp_tss_left = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'};
							$temp_tss_right = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_right'};
						}
						else{
							$temp_tss_left = $temp_sorted_exon[@temp_sorted_exon-1][1];
							$temp_tss_right = $temp_sorted_exon[@temp_sorted_exon-1][1];
						}
						$temp_first_exon = "$temp_sorted_exon[@temp_sorted_exon-1][0]-$temp_sorted_exon[@temp_sorted_exon-1][1]";
					}
					else{
						if(exists $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'}){
							$temp_tss_left = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_left'};
							$temp_tss_right = $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id}{'tss_right'};
						}
						else{
							$temp_tss_left = $temp_sorted_exon[0][0];
							$temp_tss_right = $temp_sorted_exon[0][0];
						}
						$temp_first_exon = "$temp_sorted_exon[0][0]-$temp_sorted_exon[0][1]";
					}
					if(($temp_gene_tss_left == 0) || ($temp_gene_tss_left > $temp_tss_left)){
						$temp_gene_tss_left = $temp_tss_left;
					}
					if(($temp_gene_tss_right == 0) || ($temp_gene_tss_right < $temp_tss_right)){
						$temp_gene_tss_right = $temp_tss_right;
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}){
						%{$temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}{$temp_transcript_id}){
						$temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}{$temp_transcript_id} = 1;
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}){
						%{$temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}){
						%{$temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}} = ();
					}
					unless(exists $temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}{$temp_transcript_id}){
						$temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}{$temp_transcript_id} = 1;
					}
					delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id}{$temp_transcript_id};
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}){
					%{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}} = ();
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}){
					%{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}} = ();
				}
				unless(exists $temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}{$temp_gene_id}){
					$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}{$temp_gene_id} = 1;
				}
				delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand}{$temp_gene_id};
			}
			delete $temp_query_transcript_ref->{$temp_chr}{$temp_strand};
		}
		delete $temp_query_transcript_ref->{$temp_chr};
	}
	foreach my $temp_chr (keys %temp_query){
		foreach my $temp_strand (keys %{$temp_query{$temp_chr}}){
			foreach my $temp_first_exon (keys %{$temp_query{$temp_chr}{$temp_strand}{'fe'}}){
				my @temp_first_exon = split('-', $temp_first_exon);
				my $temp_transcript_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon}});
				print $temp_fe_fh "$temp_chr\tASEF\texon\t$temp_first_exon[0]\t$temp_first_exon[1]\t0\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\;\n";
				delete $temp_query{$temp_chr}{$temp_strand}{'fe'}{$temp_first_exon};
			}
			foreach my $temp_tss_left (keys %{$temp_query{$temp_chr}{$temp_strand}{'tss'}}){
				foreach my $temp_tss_right (keys %{$temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}}){
					my $temp_transcript_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right}});
					print $temp_tss_fh "$temp_chr\tASEF\tTSS\t$temp_tss_left\t$temp_tss_right\t0\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\;\n";
					delete $temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left}{$temp_tss_right};
				}
				delete $temp_query{$temp_chr}{$temp_strand}{'tss'}{$temp_tss_left};
			}
			foreach my $temp_gene_tss_left (keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}}){
				foreach my $temp_gene_tss_right (keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}}){
					my $temp_gene_id = join(',', keys %{$temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right}});
					print $temp_gene_tss_fh "$temp_chr\tASEF\tTSS\t$temp_gene_tss_left\t$temp_gene_tss_right\t0\t$temp_strand\t.\tgene_id \"$temp_gene_id\"\;\n";
					delete $temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left}{$temp_gene_tss_right};
				}
				delete $temp_query{$temp_chr}{$temp_strand}{'gene_tss'}{$temp_gene_tss_left};
			}
			delete $temp_query{$temp_chr}{$temp_strand};
		}
		delete $temp_query{$temp_chr};
	}
	close $temp_fe_fh;
	close $temp_tss_fh;
	close $temp_gene_tss_fh;

	return 1;
}

sub compare_FE{
	my ($temp_query_fe_ref, $temp_query_gene_tss_ref, $temp_dir, $temp_btp) = @_;
	open(my $temp_fh, "$temp_btp intersect -wa -s -nonamecheck -a ${temp_dir}/query_fe.gtf -b ${temp_dir}/ref_fe.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		#print "$temp_line\n";
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_transcript_id (@temp_transcript_id){
			$temp_query_fe_ref->{$temp_transcript_id}{'type'} = "overlapping";
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp intersect -wa -S -nonamecheck -a ${temp_dir}/query_tss.gtf -b ${temp_dir}/ref_fee.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		foreach my $temp_transcript_id (@temp_transcript_id){
			if($temp_query_fe_ref->{$temp_transcript_id}{'type'} eq "novel_orphan"){
				$temp_query_fe_ref->{$temp_transcript_id}{'type'} = "novel_bi-directional";
			}
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_tss.gtf -b ${temp_dir}/ref_tss.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript_id at line:\n$temp_line\n";
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
			$temp_query_fe_ref->{$temp_transcript_id}{'nearest_TSS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_fe_ref->{$temp_transcript_id}{'nearest_TSS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_fe_ref->{$temp_transcript_id}{'nearest_TSS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_tss.gtf -b ${temp_dir}/ref_gene_tss.gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my @temp_transcript_id = ();
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			@temp_transcript_id = split(',', $2);
		}
		else{
			print "Cannot find transcript_id at line:\n$temp_line\n";
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
			$temp_query_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;
	open($temp_fh, "$temp_btp closest -s -D b -nonamecheck -a ${temp_dir}/query_gene_tss.gtf -b ${temp_dir}/ref_gene_tss.gtf |");
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
			$temp_query_gene_tss_ref->{$temp_query_gene_id}{'nearest_gene_TSS'} = $temp_line[18];
			foreach my $temp_ref_gene_id (@temp_ref_gene_id){
				unless(exists $temp_query_gene_tss_ref->{$temp_query_gene_id}{'nearest_gene_TSS_gene_id'}{$temp_ref_gene_id}){
					$temp_query_gene_tss_ref->{$temp_query_gene_id}{'nearest_gene_TSS_gene_id'}{$temp_ref_gene_id} = 1;
				}
			}
		}
	}
	close $temp_fh;

	return 1;
}

sub print_result{
	my ($temp_query_fe_ref, $temp_query_gene_tss_ref, $temp_out) = @_;
	open(my $temp_transcript_fh, "> ${temp_out}_transcript.txt") or die "Cannot create ${temp_out}_transcript.txt!\n";
	open(my $temp_gene_fh, "> ${temp_out}_gene.txt") or die "Cannot create ${temp_out}_gene.txt!\n";
	print $temp_transcript_fh "transcript_id\tFE\tnearest_TSS\tnearest_TSS_same_gene\tnearest_gene_TSS\tnearest_gene_TSS_same_gene\n";
	foreach my $temp_transcript_id (keys %$temp_query_fe_ref){
		my $temp_nearest_TSS_same_gene = 0;
		my $temp_nearest_gene_TSS_same_gene = 0;
		if(exists $temp_query_fe_ref->{$temp_transcript_id}{'nearest_TSS_gene_id'}{$temp_query_fe_ref->{$temp_transcript_id}{'gene_id'}}){
			$temp_nearest_TSS_same_gene = 1;
		}
		if(exists $temp_query_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS_gene_id'}{$temp_query_fe_ref->{$temp_transcript_id}{'gene_id'}}){
			$temp_nearest_gene_TSS_same_gene = 1;
		}
		print $temp_transcript_fh "$temp_transcript_id\t$temp_query_fe_ref->{$temp_transcript_id}{'type'}\t$temp_query_fe_ref->{$temp_transcript_id}{'nearest_TSS'}\t$temp_nearest_TSS_same_gene\t$temp_query_fe_ref->{$temp_transcript_id}{'nearest_gene_TSS'}\t$temp_nearest_gene_TSS_same_gene\n";
	}
	print $temp_gene_fh "gene_id\tnearest_gene_TSS\tnearest_gene_TSS_same_gene\n";
	foreach my $temp_gene_id (keys %$temp_query_gene_tss_ref){
		my $temp_nearest_gene_TSS_same_gene = 0;
		if(exists $temp_query_gene_tss_ref->{$temp_gene_id}{'nearest_gene_TSS_gene_id'}{$temp_gene_id}){
			$temp_nearest_gene_TSS_same_gene = 1;
		}
		print $temp_gene_fh "$temp_gene_id\t$temp_query_gene_tss_ref->{$temp_gene_id}{'nearest_gene_TSS'}\t$temp_nearest_gene_TSS_same_gene\n";
	}
	close $temp_transcript_fh;
	close $temp_gene_fh;

	return 1;
}


