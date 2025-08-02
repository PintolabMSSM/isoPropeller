#!/usr/bin/env perl
=start
Author: ALin
Purpose: To calculate the ratio of expression level between the fusion gene annotated by NIAP_annotate_*.pl and the fusion gene plus its parent genes for specified columns of samples.
Change log:
	v1.1	2022-09	Adopted from fusion_gene_exp_diff_v1.2.pl. Computed the ratio of expression level between a fusion gene and the sum of the fusion gene plus its parent genes. A "-1" ratio indicates the expression level of the sum of the fusion gene plus its parent genes is zero.
	v1.2	2022-10	Cluster fusion genes sharing common parent genes into fusion loci based on the NIAP annotated .gtf file. Added columns showing the fusion loci and the number of fusion genes in each fusion locus in the output.
	v1.3	2023-04	Added the count for protein-coding genes.
	v1.4	2023-10	Modified to be compatible with NIAP.
	v1.5	2023-10	Added output of ratio between the fusion and the most upstream parent gene.
=cut

use strict;

my $version = "v1.5";

my $usage = "perl fusion_gene_exp_ratio_${version}.pl <NIAP annotated .gtf> <expression table> <sample list> <output base>\n";

if(@ARGV != 4){
	print "$usage";
	exit;
}

my ($gtf, $exp, $sample, $out) = @ARGV;

my ($sample_ref) = read_sample($sample);

my ($gene_exp_ref) = read_exp($exp, $sample_ref);

my ($gene_ref, $rep_fusion_gene_ref, $fusion_gene_cluster_size_ref) = read_gtf($gtf);

my ($fusion_gene_cluster_id_ref) = fusion_gene_cluster($gene_ref, $rep_fusion_gene_ref, $fusion_gene_cluster_size_ref);

parse_fusion_gene($sample_ref, $gene_exp_ref, $gene_ref, $rep_fusion_gene_ref, $fusion_gene_cluster_size_ref, $fusion_gene_cluster_id_ref, $out);

#Read the sample list
sub read_sample{
	my ($temp_sample) = @_;
	my $temp_fh;
	open($temp_fh, $temp_sample) or die "Cannot open $temp_sample!\n";
	my %temp_sample = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		unless($temp_line){
			next;
		}
		if(exists $temp_sample{$temp_line}){
			print STDERR "<WARNING> Duplicated sample $temp_line specified!\n";
		}
		else{
			$temp_sample{$temp_line} = 1;
		}
	}
	close $temp_fh;
	my @temp_sample = sort {$a cmp $b} keys %temp_sample;

	return (\@temp_sample);
}

#Read the expression table
sub read_exp{
	my ($temp_exp, $temp_sample_ref) = @_;
	my $temp_fh;
	open($temp_fh, $temp_exp) or die "Cannot open $temp_exp!\n";
	my %temp_gene_exp = ();
	my %temp_transcript_id = ();
	my %temp_header = ();
	my $temp_header = 0;
	while(<$temp_fh>){
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_header == 0){
			my $temp_num_col = @temp_line;
			for(my $i = 0 ; $i < $temp_num_col; $i++){
				if(exists $temp_header{$temp_line[$i]}){
					print STDERR "<ERROR> Duplicated column $temp_line[$i] in $temp_exp!\n";
					exit;
				}
				else{
					$temp_header{$temp_line[$i]} = $i;
				}
			}
			unless(exists $temp_header{'transcript_id'}){
				print STDERR "<ERROR> Please rename the column of the transcript ID to transcript_id in $temp_exp!\n";
				exit;
			}
			unless(exists $temp_header{'gene_id'}){
				print STDERR "<ERROR> Please rename the column of the gene ID to gene_id in $temp_exp!\n";
				exit;
			}
			$temp_header = 1;
		}
		else{
			my $temp_transcript_id = $temp_line[$temp_header{'transcript_id'}];
			if(exists $temp_transcript_id{$temp_transcript_id}){
				print STDERR "<ERROR> Duplicated $temp_transcript_id!\n";
				exit;
			}
			else{
				$temp_transcript_id{$temp_transcript_id} = 1;
			}
			my $temp_gene_id = $temp_line[$temp_header{'gene_id'}];
			unless(exists $temp_gene_exp{$temp_gene_id}){
				%{$temp_gene_exp{$temp_gene_id}} = ();
			}
			foreach my $temp_sample (@$temp_sample_ref){
				unless(exists $temp_gene_exp{$temp_gene_id}{$temp_sample}){
					$temp_gene_exp{$temp_gene_id}{$temp_sample} = 0;
				}
				$temp_gene_exp{$temp_gene_id}{$temp_sample} += $temp_line[$temp_header{$temp_sample}];
			}
		}
	}
	close $temp_fh;

	return \%temp_gene_exp;
}

#Read the .gtf file and extract information for gene
sub read_gtf{
	my ($temp_gtf) = @_;
	my $temp_fh;
	open($temp_fh, "$temp_gtf") or die "Cannot open $temp_gtf!\n";
	my %temp_gene = ();
	my %temp_rep_fusion_gene = ();
	my %temp_fusion_gene_cluster_size = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[2] ne "gene"){
			next;
		}
		my $temp_gene_id;
		if($temp_line[8] =~ /gene_id \"([^\"]+)\"\;/){
			$temp_gene_id = $1;
		}
		else{
			print STDERR "<ERROR> Cannot find gene ID at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_gene{$temp_gene_id}){
			%{$temp_gene{$temp_gene_id}} = ();
			$temp_gene{$temp_gene_id}{'left'} = $temp_line[3];
			$temp_gene{$temp_gene_id}{'right'} = $temp_line[4];
			$temp_gene{$temp_gene_id}{'strand'} = $temp_line[6];
			$temp_gene{$temp_gene_id}{'gene_type'} = "";
		}
		if($temp_line[8] =~ /gene_type \"([^\"]+)\"\;/){
			$temp_gene{$temp_gene_id}{'gene_type'} = $1;
		}
		if($temp_gene_id =~ /\|/){
			unless(exists $temp_rep_fusion_gene{$temp_gene_id}){
				$temp_rep_fusion_gene{$temp_gene_id} = $temp_gene_id;
				$temp_fusion_gene_cluster_size{$temp_gene_id} = 1;
				#print "$temp_gene_id\t$temp_rep_fusion_gene{$temp_gene_id}\t$temp_fusion_gene_cluster_size{$temp_gene_id}\n";
			}
		}
	}
	close $temp_fh;

	return (\%temp_gene, \%temp_rep_fusion_gene, \%temp_fusion_gene_cluster_size);
}

#Compare whether two fusion genes share any parent genes
sub compare_fusion_gene{
	my ($temp_fusion_gene_a, $temp_fusion_gene_b) = @_;
	my @temp_fusion_gene_a = split('\|', $temp_fusion_gene_a);
	my $temp_num_gene_a = @temp_fusion_gene_a;
	my @temp_fusion_gene_b = split('\|', $temp_fusion_gene_b);
	my $temp_num_gene_b = @temp_fusion_gene_b;
	for(my $i = 0; $i < $temp_num_gene_a; $i++){
		for(my $j = 0; $j < $temp_num_gene_b; $j++){
			if($temp_fusion_gene_a[$i] eq $temp_fusion_gene_b[$j]){

				return 1;
			}
		}
	}

	return 0;
}

#Find the representative of fusion gene in the cluster
sub find_rep_fusion_gene{
	my ($temp_rep_fusion_gene_ref, $temp_query_fusion_gene) = @_;
	if($temp_rep_fusion_gene_ref->{$temp_query_fusion_gene} ne $temp_query_fusion_gene){
		my $temp_rep_fusion_gene = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_rep_fusion_gene_ref->{$temp_query_fusion_gene});

		return $temp_rep_fusion_gene;
	}
	else{

		return $temp_query_fusion_gene;
	}
}

#Merge fusion genes share any parent genes
sub merge_fusion_gene{
	my ($temp_rep_fusion_gene_ref, $temp_fusion_gene_cluster_size_ref, $temp_fusion_gene_a, $temp_fusion_gene_b) = @_;
	my $temp_rep_fusion_gene_a = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene_a);
	my $temp_rep_fusion_gene_b = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene_b);
	if($temp_rep_fusion_gene_a eq $temp_rep_fusion_gene_b){

		return 0;
	}
	#print "Merging $temp_fusion_gene_a and $temp_fusion_gene_b\nRep of $temp_fusion_gene_a: $temp_rep_fusion_gene_a ($temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_a})\nRep of $temp_fusion_gene_b: $temp_rep_fusion_gene_b ($temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_b})\n";
	if($temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_a} > $temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_b}){
		$temp_rep_fusion_gene_ref->{$temp_rep_fusion_gene_b} = $temp_rep_fusion_gene_a;
		$temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_a} += $temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_b};
	}
	else{
		$temp_rep_fusion_gene_ref->{$temp_rep_fusion_gene_a} = $temp_rep_fusion_gene_b;
		$temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_b} += $temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_a};
	}
	$temp_rep_fusion_gene_a = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene_a);
	$temp_rep_fusion_gene_b = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene_b);
	#print "Updated\nRep of $temp_fusion_gene_a: $temp_rep_fusion_gene_a ($temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_a})\nRep of $temp_fusion_gene_b: $temp_rep_fusion_gene_b ($temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene_b})\n";

	return 1;
}

#Identify fusion gene cluster
sub fusion_gene_cluster{
	my ($temp_gene_ref, $temp_rep_fusion_gene_ref, $temp_fusion_gene_cluster_size_ref) = @_;
	my @temp_fusion_gene = keys %$temp_rep_fusion_gene_ref;
	my $temp_num_fusion_gene = @temp_fusion_gene;
	for(my $i = 0; $i < ($temp_num_fusion_gene - 1); $i++){
		for(my $j = ($i + 1); $j < $temp_num_fusion_gene; $j++){
			if(compare_fusion_gene($temp_fusion_gene[$i], $temp_fusion_gene[$j])){
				merge_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene_cluster_size_ref, $temp_fusion_gene[$i], $temp_fusion_gene[$j]);
			}
		}
	}
	my %temp_fusion_gene_cluster = ();
	foreach my $temp_fusion_gene (keys %$temp_rep_fusion_gene_ref){
		my $temp_rep_fusion_gene = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_fusion_gene);
		unless(exists $temp_fusion_gene_cluster{$temp_rep_fusion_gene}){
			%{$temp_fusion_gene_cluster{$temp_rep_fusion_gene}} = ();
		}
		my @temp_gene = split('\|', $temp_fusion_gene);
		foreach my $temp_gene (@temp_gene){
			unless(exists $temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$temp_gene}){
				if($temp_gene_ref->{$temp_gene}{'strand'} eq "-"){
					$temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$temp_gene} = $temp_gene_ref->{$temp_gene}{'right'};
				}
				else{
					$temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$temp_gene} = $temp_gene_ref->{$temp_gene}{'left'};
				}
			}
			#print "$temp_rep_fusion_gene\t$temp_fusion_gene\t$temp_gene\t$temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$temp_gene}\n";
		}
	}
	my %temp_fusion_gene_cluster_id = ();
	foreach my $temp_rep_fusion_gene (keys %temp_fusion_gene_cluster){
		my @temp_sorted_gene_id = ();
		if($temp_gene_ref->{$temp_rep_fusion_gene}{'strand'} eq "-"){
			@temp_sorted_gene_id = sort {($temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$b} <=> $temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$a}) || ($a cmp $b)} keys %{$temp_fusion_gene_cluster{$temp_rep_fusion_gene}};
		}
		else{
			@temp_sorted_gene_id = sort {($temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$a} <=> $temp_fusion_gene_cluster{$temp_rep_fusion_gene}{$b}) || ($a cmp $b)} keys %{$temp_fusion_gene_cluster{$temp_rep_fusion_gene}};
		}
		unless(exists $temp_fusion_gene_cluster_id{$temp_rep_fusion_gene}){
			$temp_fusion_gene_cluster_id{$temp_rep_fusion_gene} = join("\|", @temp_sorted_gene_id);
		}
	}

	return (\%temp_fusion_gene_cluster_id);
}

#Parse fusion gene and output the result
sub parse_fusion_gene{
	my ($temp_sample_ref, $temp_gene_exp_ref, $temp_gene_ref, $temp_rep_fusion_gene_ref, $temp_fusion_gene_cluster_size_ref, $temp_fusion_gene_cluster_id_ref, $temp_out) = @_;
	my %temp_num_passed_sample = ();
	#my $temp_anno_fh;
	#my $temp_anno_out_fh;
	#open($temp_anno_fh, $temp_annotation) or die "Cannot open $temp_annotation!\n";
	open(my $temp_total_fh, "> ${temp_out}_fusion_gene_ratio.txt") or die "Cannot create ${temp_out}_fusion_gene_ratio.txt!\n";
	open(my $temp_first_fh, "> ${temp_out}_fusion_gene_first_ratio.txt") or die "Cannot create ${temp_out}_fusion_gene_first_ratio.txt!\n";
	#open($temp_anno_out_fh, "> ${temp_out}_reclocus.txt") or die "Cannot create ${temp_out}_reclocus.txt!\n";
	print $temp_total_fh "gene_id\tcluster_id\tcluster_size\tnum_pc";
	print $temp_first_fh "gene_id\tcluster_id\tcluster_size\tnum_pc";
	foreach my $temp_sample (@$temp_sample_ref){
		print $temp_total_fh "\t$temp_sample";
		print $temp_first_fh "\t$temp_sample";
	}
	print $temp_total_fh "\tnum_passed_sample\n";
	print $temp_first_fh "\n";
	foreach my $temp_gene_id (keys %$temp_gene_exp_ref){
		if($temp_gene_id =~ /\|/){
			my $temp_rep_fusion_gene = find_rep_fusion_gene($temp_rep_fusion_gene_ref, $temp_gene_id);
			my @temp_component_gene = split('\|', $temp_gene_ref->{$temp_gene_id}{'gene_type'});
			my $temp_num_pc = 0;
			foreach my $temp_component_gene (@temp_component_gene){
				if($temp_component_gene eq 'protein_coding'){
					$temp_num_pc++;
				}
			}
			print $temp_total_fh "$temp_gene_id\t$temp_fusion_gene_cluster_id_ref->{$temp_rep_fusion_gene}\t$temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene}\t$temp_num_pc";
			print $temp_first_fh "$temp_gene_id\t$temp_fusion_gene_cluster_id_ref->{$temp_rep_fusion_gene}\t$temp_fusion_gene_cluster_size_ref->{$temp_rep_fusion_gene}\t$temp_num_pc";
			if(exists $temp_num_passed_sample{$temp_gene_id}){
				print STDERR "<ERROR> Duplicated gene $temp_gene_id!\n";
				exit;
			}
			else{
				$temp_num_passed_sample{$temp_gene_id} = 0;
			}
			my @temp_gene_id = split('\|', $temp_gene_id);
			foreach my $temp_sample (@$temp_sample_ref){
				my $temp_total = $temp_gene_exp_ref->{$temp_gene_id}{$temp_sample};
				my $temp_first_total = $temp_gene_exp_ref->{$temp_gene_id[0]}{$temp_sample} + $temp_gene_exp_ref->{$temp_gene_id}{$temp_sample};
				my $temp_fusion = $temp_gene_exp_ref->{$temp_gene_id}{$temp_sample};
				foreach my $temp_parent_gene_id (@temp_gene_id){
					$temp_total += $temp_gene_exp_ref->{$temp_parent_gene_id}{$temp_sample};
				}
				my $temp_ratio = -1;
				if($temp_total > 0){
					$temp_ratio = $temp_fusion / $temp_total;
				}
				print $temp_total_fh "\t$temp_ratio";
				if($temp_ratio > 0.5){
					$temp_num_passed_sample{$temp_gene_id}++;
				}
				$temp_ratio = -1;
				if($temp_first_total > 0){
					$temp_ratio = $temp_fusion / $temp_first_total;
				}
				print $temp_first_fh "\t$temp_ratio";
			}
			print $temp_total_fh "\t$temp_num_passed_sample{$temp_gene_id}\n";
			print $temp_first_fh "\n";
		}
	}
	close $temp_total_fh;
	close $temp_first_fh;
	#my $temp_header = 0;
	#my %temp_header = ();
	#while(<$temp_anno_fh>){
		#chomp;
		#my $temp_line = $_;
		#$temp_line =~ s/\r//;
		#print $temp_anno_out_fh "$temp_line\n";
		#my @temp_line = split('\t', $temp_line);
		#if($temp_header == 0){
			#my $temp_num_col = @temp_line;
			#for(my $i = 0; $i < $temp_num_col; $i++){
				#if(exists $temp_header{$temp_line[$i]}){
					#print "Duplicated column $temp_line[$i]\n";	
					#exit;
				#}
				#else{
					#$temp_header{$temp_line[$i]} = $i;
				#}
			#}
			#unless(exists $temp_header{'gene_id'}){
				#print STDERR "<ERROR> Please rename the column of the gene ID to gene_id in $temp_annotation!\n";
				#exit;
			#}
			#unless(exists $temp_header{'status'}){
				#print STDERR "<ERROR> Please rename the column of the transcript annotation type to status in $temp_annotation!\n";
			#}
			#$temp_header = 1;
		#}
		#else{
			#my $temp_gene_id = $temp_line[$temp_header{'gene_id'}];
			#if((exists $temp_num_passed_sample{$temp_gene_id}) && ($temp_num_passed_sample{$temp_gene_id} > 0)){
				#if($temp_line[$temp_header{'status'}] eq 'fusion'){
					#$temp_line[$temp_header{'status'}] = 'reclocus';
				#}
				#else{
					#print STDERR "<ERROR> There seems to be a conflict between the gene_id and status columns at\n$temp_line\n";
					#exit;
				#}
			#}
			#my $temp_line_out = join("\t", @temp_line);
			#print $temp_anno_out_fh "\t$temp_line_out";
		#}
		#print $temp_anno_out_fh "\n";
	#}
	#close $temp_anno_fh;
	#close $temp_anno_out_fh;

	return 1;
}






