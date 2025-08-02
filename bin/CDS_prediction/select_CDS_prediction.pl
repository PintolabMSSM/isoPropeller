#!/usr/bin/env perl
=start
Author: ALin
Purpose: To select CDS preiction from GMST, CPAT and Transdecoder.
Change log:
	v1.1	2023-11	Adopted from compare_CDS_prediction_v1.2.pl.
=cut

use strict;
use Getopt::Long;

my $version = "v1.1";

my ($in, $out, $min_orf, $min_prob, $help) = ("", "", 1, 0, 0);

my $usage = "Usage: perl select_CDS_prediction.pl
	-i <String> Input list of prefix
	-o <String> Output
	-m <Integer> Minimal ORF length (Default: $min_orf)
	-p <Float> Minimal CPAT coding probablity (Default: $min_prob)
	-h <Boolean> Help
";

GetOptions(
	'i=s'	=>	\$in,
	'o=s'	=>	\$out,
	'm=i'	=>	\$min_orf,
	'p=f'	=>	\$min_prob,
	'h!'	=>	\$help
);

unless($in && $out){
	print "$usage";
	exit;
}

if($help){
	print "usage";
	exit;
}

my ($orf_ref) = read_list($in, $min_orf);
my ($gtf_ref) = read_gtf($in);
print_result($orf_ref, $gtf_ref, $min_prob, $out);

sub read_list{
	my ($temp_in, $temp_min_orf) = @_;
	my %temp_orf = ();
	open(my $temp_fh, $temp_in) or die "Cannot open $temp_in!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my ($temp_prefix, $temp_source) = split('\t', $temp_line);
		open(my $temp_orf_fh, "${temp_prefix}.txt") or die "Cannot open ${temp_prefix}.txt!\n";
		my $temp_header = 0;
		my %temp_header = ();
		while(<$temp_orf_fh>){
			chomp;
			$temp_line = $_;
			$temp_line =~ s/\r//;
			my @temp_line = split('\t', $temp_line);
			if($temp_header == 0){
				for(my $i = 0; $i < @temp_line; $i++){
					if(exists $temp_header{$temp_line[$i]}){
						print STDERR "<ERROR> Duplicated column found: $temp_line[$i]!\n";
						exit;
					}
					else{
						$temp_header{$temp_line[$i]} = $i;
					}
				}
				foreach my $temp_col ('transcript_id', 'length', 'cds_length', '5utr_length', '3utr_length', 'canonical_start', 'tis_efficiency', 'canonical_stop'){
					unless(exists $temp_header{$temp_col}){
						print STDERR "<ERROR> Missing column $temp_col in ${temp_prefix}.txt!\n";
						exit;
					}
				}
				$temp_header = 1;
			}
			else{
				my $temp_transcript_id = $temp_line[$temp_header{'transcript_id'}];
				if(($temp_line[$temp_header{'cds_length'}] >= $temp_min_orf) && ($temp_line[$temp_header{'canonical_stop'}] == 1)){
					my ($temp_prob, $temp_score) = (0, 0);
					if(exists $temp_header{'probability'}){
						$temp_prob = $temp_line[$temp_header{'probability'}];
					}
					elsif(exists $temp_header{'score'}){
						$temp_score = $temp_line[$temp_header{'score'}];
					}
					unless(exists $temp_orf{$temp_transcript_id}){
						@{$temp_orf{$temp_transcript_id}} = ();
					}
					my @temp_orf = ($temp_line[$temp_header{'cds_length'}], $temp_line[$temp_header{'5utr_length'}] + 1, $temp_line[$temp_header{'length'}] - $temp_line[$temp_header{'3utr_length'}], $temp_line[$temp_header{'canonical_start'}], $temp_line[$temp_header{'tis_efficiency'}], $temp_line[$temp_header{'canonical_stop'}], $temp_source, $temp_prob, $temp_score);
					push(@{$temp_orf{$temp_transcript_id}}, \@temp_orf);
				}
			}
		}
		close $temp_orf_fh;
	}
	close $temp_fh;

	return (\%temp_orf);
}

sub read_gtf{
	my ($temp_in) = @_;
	my %temp_gtf = ();
	open(my $temp_fh, $temp_in) or die "Cannot open $temp_in!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my ($temp_prefix, $temp_source) = split('\t', $temp_line);
		unless(exists $temp_gtf{$temp_source}){
			%{$temp_gtf{$temp_source}} = ();
		}
		open(my $temp_gtf_fh, "${temp_prefix}.gtf") or die "Cannot open ${temp_prefix}.gtf!\n";
		while(<$temp_gtf_fh>){
			chomp;
			my $temp_line = $_;
			$temp_line =~ s/\r//;
			my $temp_transcript_id;
			if($temp_line =~ /(^|\s)transcript_id \"([^\"]+)\"/){
				$temp_transcript_id = $2;
			}
			else{
				print STDERR "Missing transcript_id in ${temp_prefix}.gtf!\n$temp_line\n";
				exit;
			}
			unless(exists $temp_gtf{$temp_source}{$temp_transcript_id}){
				@{$temp_gtf{$temp_source}{$temp_transcript_id}} = ();
			}
			push(@{$temp_gtf{$temp_source}{$temp_transcript_id}}, $temp_line);
		}
	}

	return (\%temp_gtf);
}

sub print_result{
	my ($temp_orf_ref, $temp_gtf_ref, $temp_min_prob, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "Cannot create $temp_out!\n";
	foreach my $temp_transcript_id (keys %$temp_orf_ref){
		#print "$temp_transcript_id\n";
		my ($temp_num_source, $temp_diversity) = (1, 1);
		#ORF: CDS length, CDS start, CDS end, canonical start, TIS efficiency, canonical stop, source, CPAT coding probability, Transdecoder score
		if(@{$temp_orf_ref->{$temp_transcript_id}} == 0){
			next;
		}
		elsif(@{$temp_orf_ref->{$temp_transcript_id}} == 1){
			#print "Predicted by single tool\n";
			if(($temp_orf_ref->{$temp_transcript_id}[0][7] >= $temp_min_prob) && ($temp_orf_ref->{$temp_transcript_id}[0][8] == 0)){
				#print "Predicted by CPAT\n";
				print_gtf($temp_fh, $temp_gtf_ref->{$temp_orf_ref->{$temp_transcript_id}[0][6]}{$temp_transcript_id}, 'CPAT', $temp_num_source, $temp_diversity);
				next;
			}
		}
		my ($temp_orf_cluster_ref) = orf_cluster($temp_orf_ref->{$temp_transcript_id});
		#The sorting here may not work for more than three tools.
		my @temp_sorted_orf = sort {($temp_orf_cluster_ref->{$b}{'size'} <=> $temp_orf_cluster_ref->{$a}{'size'})} keys %$temp_orf_cluster_ref;
		#print "Sorted ORFs: @temp_sorted_orf\n";
		my $temp_rep_orf_best = find_rep_orf($temp_orf_cluster_ref, $temp_sorted_orf[0]);
		#print "Rep: $temp_rep_orf_best\n";
		if($temp_orf_cluster_ref->{$temp_rep_orf_best}{'size'} > 1){
			my @temp_candidate_orf = ();
			my @temp_source = ();
			foreach my $temp_orf (keys %$temp_orf_cluster_ref){
				my $temp_rep_orf = find_rep_orf($temp_orf_cluster_ref, $temp_orf);
				if($temp_rep_orf == $temp_rep_orf_best){
					#print "$temp_orf: @{$temp_orf_ref->{$temp_transcript_id}[$temp_orf]}\n";
					push(@temp_candidate_orf, \@{$temp_orf_ref->{$temp_transcript_id}[$temp_orf]});
					push(@temp_source, $temp_orf_ref->{$temp_transcript_id}[$temp_orf][6]);
				}
			}
			my @temp_sorted_candidate_orf = sort {
				($b->[3] <=> $a->[3]) || 
				($b->[4] <=> $a->[4]) || 
				($b->[0] <=> $a->[0]) || 
				($a->[2] <=> $b->[2]) || 
				($a->[6] <=> $b->[6])
			} @temp_candidate_orf;
			#print "$temp_sorted_candidate_orf[0][6]\n";
			my $temp_source = join(",", sort {$a cmp $b} @temp_source);
			$temp_num_source = @temp_source;
			$temp_diversity = keys %{$temp_orf_cluster_ref->{$temp_rep_orf_best}{'diversity'}};
			print_gtf($temp_fh, $temp_gtf_ref->{$temp_sorted_candidate_orf[0][6]}{$temp_transcript_id}, $temp_source, $temp_num_source, $temp_diversity);
		}
	}
	close $temp_fh;

	return 1;
}

sub print_gtf{
	my ($temp_fh, $temp_gtf_ref, $temp_source, $temp_num_source, $temp_diversity) = @_;
	foreach my $temp_line (@$temp_gtf_ref){
		print $temp_fh "$temp_line cds_source \"$temp_source\"\; cds_support \"$temp_num_source\"\; cds_diversity \"$temp_diversity\"\;\n";
	}

	return 1;
}

#Cluster ORFs
sub orf_cluster{
	my ($temp_orf_ref) = @_;
	my %temp_orf_cluster = ();
	for(my $i = 0; $i < @$temp_orf_ref; $i++){
		#ORF: CDS length, CDS start, CDS end, canonical start, TIS efficiency, canonical stop, source, CPAT coding probability, Transdecoder score
		#print "@{$temp_orf_ref->[$i]}\n";
		unless(exists $temp_orf_cluster{$i}){
			%{$temp_orf_cluster{$i}} = ();
			$temp_orf_cluster{$i}{'rep'} = $i;
			$temp_orf_cluster{$i}{'size'} = 1;
			my $temp_orf = "$temp_orf_ref->[$i][1]-$temp_orf_ref->[$i][2]";
			%{$temp_orf_cluster{$i}{'diversity'}} = ();
			$temp_orf_cluster{$i}{'diversity'}{$temp_orf} = 1;
			if($temp_orf_ref->[$i][5]){
				$temp_orf_cluster{$i}{'complete'} = 1;
			}
			else{
				$temp_orf_cluster{$i}{'complete'} = 0;
			}
		}
	}
	for(my $i = 0; $i < (@$temp_orf_ref - 1); $i++){
		for(my $j = $i + 1; $j < @$temp_orf_ref; $j++){
			my ($temp_identifical, $temp_contained, $temp_intersected) = compare_orf($temp_orf_ref->[$i], $temp_orf_ref->[$j]);
			#print "$i: @{$temp_orf_ref->[$i]}, $j: @{$temp_orf_ref->[$j]}\n$temp_identifical, $temp_contained, $temp_intersected\n";
			if($temp_identifical || $temp_contained){
				merge_orf(\%temp_orf_cluster, $i, $j);
			}
		}
	}

	return (\%temp_orf_cluster);
}

#Compare two ORFs
sub compare_orf{
	my @temp_orf = @_;
	my @temp_sorted_orf = sort {($a->[1] <=> $b->[1]) || ($a->[2] <=> $b->[2])} @temp_orf;
	my ($temp_identifical, $temp_contained, $temp_intersected) = (0, 0, 0);
	if($temp_sorted_orf[0][2] >= $temp_sorted_orf[1][1]){
		if(($temp_sorted_orf[0][1] == $temp_sorted_orf[1][1]) && ($temp_sorted_orf[0][2] == $temp_sorted_orf[1][2])){
			$temp_identifical = 1;
		}
		elsif($temp_sorted_orf[0][2] == $temp_sorted_orf[1][2]){
			$temp_contained = 1;
		}
		elsif($temp_sorted_orf[0][2] < $temp_sorted_orf[1][2]){
			$temp_intersected = 1;
		}
	}

	return ($temp_identifical, $temp_contained, $temp_intersected);
}

#Find the representative of ORF in the cluster
sub find_rep_orf{
	my ($temp_orf_cluster_ref, $temp_query_orf) = @_;
	if($temp_orf_cluster_ref->{$temp_query_orf}{'rep'} ne $temp_query_orf){
		my $temp_rep_orf = find_rep_orf($temp_orf_cluster_ref, $temp_orf_cluster_ref->{$temp_query_orf}{'rep'});

		return $temp_rep_orf;
	}
	else{

		return $temp_query_orf;
	}
}

#Merge ORFs that are identical to or contained in each 
sub merge_orf{
	my ($temp_orf_cluster_ref, $temp_orf_a, $temp_orf_b) = @_;
	my $temp_rep_orf_a = find_rep_orf($temp_orf_cluster_ref, $temp_orf_a);
	my $temp_rep_orf_b = find_rep_orf($temp_orf_cluster_ref, $temp_orf_b);
	if($temp_rep_orf_a eq $temp_rep_orf_b){

		return 0;
	}
	#print "Merging $temp_orf_a and $temp_orf_b\nRep of $temp_orf_a: $temp_rep_orf_a ($temp_orf_cluster_ref->{$temp_rep_orf_a}{'size'}, $temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'})\nRep of $temp_orf_b: $temp_rep_orf_b ($temp_orf_cluster_ref->{$temp_rep_orf_b}{'size'}, $temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'})\n";
	if($temp_orf_cluster_ref->{$temp_rep_orf_a}{'size'} > $temp_orf_cluster_ref->{$temp_rep_orf_b}{'size'}){
		$temp_orf_cluster_ref->{$temp_rep_orf_b}{'rep'} = $temp_rep_orf_a;
		$temp_orf_cluster_ref->{$temp_rep_orf_a}{'size'} += $temp_orf_cluster_ref->{$temp_rep_orf_b}{'size'};
		foreach my $temp_orf (keys %{$temp_orf_cluster_ref->{$temp_rep_orf_b}{'diversity'}}){
			unless(exists $temp_orf_cluster_ref->{$temp_rep_orf_a}{'diversity'}{$temp_orf}){
				$temp_orf_cluster_ref->{$temp_rep_orf_a}{'diversity'}{$temp_orf} = 1;
			}
		}
		if($temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'} < $temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'}){
			$temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'} = $temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'};
		}
	}
	else{
		$temp_orf_cluster_ref->{$temp_rep_orf_a}{'rep'} = $temp_rep_orf_b;
		$temp_orf_cluster_ref->{$temp_rep_orf_b}{'size'} += $temp_orf_cluster_ref->{$temp_rep_orf_a}{'size'};
		foreach my $temp_orf (keys %{$temp_orf_cluster_ref->{$temp_rep_orf_a}{'diversity'}}){
			unless(exists $temp_orf_cluster_ref->{$temp_rep_orf_b}{'diversity'}{$temp_orf}){
				$temp_orf_cluster_ref->{$temp_rep_orf_b}{'diversity'}{$temp_orf} = 1;
			}
		}
		if($temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'} < $temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'}){
			$temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'} = $temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'};
		}
	}
	$temp_rep_orf_a = find_rep_orf($temp_orf_cluster_ref, $temp_orf_a);
	$temp_rep_orf_b = find_rep_orf($temp_orf_cluster_ref, $temp_orf_b);
	#print "Updated\nRep of $temp_orf_a: $temp_rep_orf_a ($temp_orf_cluster_ref->{$temp_rep_orf_a}{'size'}, $temp_orf_cluster_ref->{$temp_rep_orf_a}{'complete'})\nRep of $temp_orf_b: $temp_rep_orf_b ($temp_orf_cluster_ref->{$temp_rep_orf_b}{'size'}, $temp_orf_cluster_ref->{$temp_rep_orf_b}{'complete'})\n";

	return 1;
}



