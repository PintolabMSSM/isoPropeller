#!/usr/bin/env perl
=start
Author: ALin
Purpose: To revise the gene_id, gene_name and gene_type of fusion gene and the corresponding parent genes annotated by NIAP_annotate_*.pl for the specified list of fusion gene_id.
Change log:
	v1.1	2023-04	Limited to at most one protein-coding gene and simple fusion gene.
=cut

use strict;

my $version = "v1.1";

my $usage = "perl gtf_reclocus.pl <NIAP annotated .gtf> <reference .gtf> <fusion gene list> <output>\n";

if(@ARGV != 4){
	print "$usage";
	exit;
}

my ($gtf, $ref_gtf, $list, $out) = @ARGV;

my ($mapping_ref) = read_list($list);
update_mapping($mapping_ref, $gtf);
update_mapping($mapping_ref, $ref_gtf);
revise_mapping($mapping_ref);
update_gtf($mapping_ref, $gtf, $out);

#Read the fusion gene list
sub read_list{
	my ($temp_list) = @_;
	open(my $temp_fh, $temp_list) or die "Cannot open $temp_list!\n";
	my %temp_mapping = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		if($temp_line =~ /\|/){
			unless(exists $temp_mapping{$temp_line}){
				%{$temp_mapping{$temp_line}} = ();
				$temp_mapping{$temp_line}{'gene_id'} = "";
				$temp_mapping{$temp_line}{'gene_name'} = "";
				$temp_mapping{$temp_line}{'gene_type'} = "";
				my @temp_line = split('\|', $temp_line);
				foreach my $temp_gene_id (@temp_line){
					unless(exists $temp_mapping{$temp_gene_id}){
						%{$temp_mapping{$temp_gene_id}} = ();
						$temp_mapping{$temp_gene_id}{'gene_id'} = "";
						$temp_mapping{$temp_gene_id}{'gene_name'} = "";
						$temp_mapping{$temp_gene_id}{'gene_type'} = "";
					}
				}
			}
		}
		else{
			if($temp_line =~ /^$/){
				next;
			}
			else{
				print STDERR "<ERROR> $temp_line may be not a fusion gene!\n";
				exit;
			}
		}

	}
	close $temp_fh;

	return (\%temp_mapping);
}

#Read the .gtf file and extract information for the listed fusion genes and the corresponding parent genes
sub update_mapping{
	my ($temp_mapping_ref, $temp_gtf) = @_;
	open(my $temp_fh, $temp_gtf) or die "Cannot open $temp_gtf!\n";
	my %temp_mapping = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] ne 'gene') && ($temp_line[2] ne 'transcript') && ($temp_line[2] ne 'exon')){
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
		if(exists $temp_mapping_ref->{$temp_gene_id}){
			if($temp_mapping_ref->{$temp_gene_id}{'gene_id'} eq ''){
				$temp_mapping_ref->{$temp_gene_id}{'gene_id'} = $temp_gene_id;
				if($temp_line[8] =~ /gene_name \"([^\"]+)\"/){
					$temp_mapping_ref->{$temp_gene_id}{'gene_name'} = $1;
				}
				if($temp_line[8] =~ /gene_type \"([^\"]+)\"/){
					$temp_mapping_ref->{$temp_gene_id}{'gene_type'} = $1;
				}
			}
		}
	}
	close $temp_fh;

	return 1;
}

#Revise the gene information based on the fusion gene type
sub revise_mapping{
	my ($temp_mapping_ref) = @_;
	foreach my $temp_gene_id (keys %$temp_mapping_ref){
		if($temp_gene_id =~ /\|/){
			my @temp_gene_id = split('\|', $temp_gene_id);
			my ($temp_num_pc, $temp_pc_gene_id) = (0, "");
			foreach my $temp_parent_gene_id (@temp_gene_id){
				if($temp_mapping_ref->{$temp_parent_gene_id}{'gene_type'} eq 'protein_coding'){
					$temp_num_pc++;
					$temp_pc_gene_id = $temp_parent_gene_id;
				}
			}
			if($temp_num_pc > 1){
				print STDERR "<ERROR> Not supporting $temp_gene_id containing more than one protein-coding gene!\n";
				exit;
			}
			elsif($temp_num_pc == 1){
				$temp_mapping_ref->{$temp_gene_id}{'gene_id'} = $temp_pc_gene_id;
				$temp_mapping_ref->{$temp_gene_id}{'gene_name'} = $temp_mapping_ref->{$temp_pc_gene_id}{'gene_name'};
				$temp_mapping_ref->{$temp_gene_id}{'gene_type'} = $temp_mapping_ref->{$temp_pc_gene_id}{'gene_type'};
				foreach my $temp_parent_gene_id (@temp_gene_id){
					if($temp_parent_gene_id ne $temp_pc_gene_id){
						$temp_mapping_ref->{$temp_parent_gene_id}{'gene_id'} = $temp_pc_gene_id;
						$temp_mapping_ref->{$temp_parent_gene_id}{'gene_name'} = $temp_mapping_ref->{$temp_pc_gene_id}{'gene_name'};
						$temp_mapping_ref->{$temp_parent_gene_id}{'gene_type'} = $temp_mapping_ref->{$temp_pc_gene_id}{'gene_type'};
					}
				}
			}
			else{
				foreach my $temp_parent_gene_id (@temp_gene_id){
					$temp_mapping_ref->{$temp_parent_gene_id}{'gene_id'} = $temp_mapping_ref->{$temp_gene_id}{'gene_id'};
					$temp_mapping_ref->{$temp_parent_gene_id}{'gene_name'} = $temp_mapping_ref->{$temp_gene_id}{'gene_name'};
					$temp_mapping_ref->{$temp_parent_gene_id}{'gene_type'} = $temp_mapping_ref->{$temp_gene_id}{'gene_type'};
				}
			}
		}
	}

	return 1;
}

#Update the gtf file
sub update_gtf{
	my ($temp_mapping_ref, $temp_gtf, $temp_out) = @_;
	open(my $temp_in_fh, $temp_gtf) or die "Cannot open $temp_gtf!\n";
	open(my $temp_out_fh, "| sort -k 1,1 -k 4,4n > $temp_out") or die "Cannot create $temp_out!\n";
	my %temp_gene = ();
	while(<$temp_in_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		if($temp_line =~ /^#/){
			next;
		}
		my @temp_line = split('\t', $temp_line);
		my $temp_gene_id;
		if($temp_line[8] =~ /gene_id \"([^\"]+)\"\;/){
			$temp_gene_id = $1;
		}
		else{
			print STDERR "<ERROR> Cannot find gene ID at line:\n$temp_line\n";
			exit;
		}
		if(exists $temp_mapping_ref->{$temp_gene_id}){
			my ($temp_reclocus_gene_id, $temp_reclocus_gene_name, $temp_reclocus_gene_type) = ($temp_mapping_ref->{$temp_gene_id}{'gene_id'}, $temp_mapping_ref->{$temp_gene_id}{'gene_name'}, $temp_mapping_ref->{$temp_gene_id}{'gene_type'});
			unless(exists $temp_gene{$temp_reclocus_gene_id}){
				%{$temp_gene{$temp_reclocus_gene_id}} = ();
				$temp_gene{$temp_reclocus_gene_id}{'chr'} = $temp_line[0];
				$temp_gene{$temp_reclocus_gene_id}{'source'} = $temp_line[1];
				$temp_gene{$temp_reclocus_gene_id}{'left'} = $temp_line[3];
				$temp_gene{$temp_reclocus_gene_id}{'right'} = $temp_line[4];
				$temp_gene{$temp_reclocus_gene_id}{'strand'} = $temp_line[6];
				$temp_gene{$temp_reclocus_gene_id}{'gene_name'} = $temp_reclocus_gene_name;
				$temp_gene{$temp_reclocus_gene_id}{'gene_type'} = $temp_reclocus_gene_type;
			}
			if($temp_gene{$temp_reclocus_gene_id}{'left'} > $temp_line[3]){
				$temp_gene{$temp_reclocus_gene_id}{'left'} = $temp_line[3];
			}
			if($temp_gene{$temp_reclocus_gene_id}{'right'} < $temp_line[4]){
				$temp_gene{$temp_reclocus_gene_id}{'right'} = $temp_line[4];
			}
			if($temp_line[2] eq 'gene'){
				next;
			}
			else{
				my @temp_attribute = split(';', $temp_line[8]);
				my %temp_attribute = ();
				for(my $i = 0; $i < @temp_attribute; $i++){
					$temp_attribute[$i] =~ s/^ //;
					if($temp_attribute[$i] =~ /^gene_id \"/){
						$temp_attribute[$i] = "gene_id \"$temp_reclocus_gene_id\"\;";
						if(exists $temp_attribute{'gene_id'}){
							print STDERR "<ERROR> Duplicated gene_id in $temp_gtf!\n$temp_line\n";
							exit;
						}
						else{
							$temp_attribute{'gene_id'} = 1;
						}
					}
					elsif($temp_attribute[$i] =~ /^gene_name \"/){
						$temp_attribute[$i] = "gene_name \"$temp_reclocus_gene_name\"\;";
						if(exists $temp_attribute{'gene_name'}){
							print STDERR "<ERROR> Duplicated gene_name in $temp_gtf!\n$temp_line\n";
							exit;
						}
						else{
							$temp_attribute{'gene_name'} = 1;
						}
					}
					elsif($temp_attribute[$i] =~ /^gene_type \"/){
						$temp_attribute[$i] = "gene_type \"$temp_reclocus_gene_type\"\;";
						if(exists $temp_attribute{'gene_type'}){
							print STDERR "<ERROR> Duplicated gene_type in $temp_gtf!\n$temp_line\n";
							exit;
						}
						else{
							$temp_attribute{'gene_type'} = 1;
						}
					}
					else{
						$temp_attribute[$i] .= "\;";
					}
				}
				unless(exists $temp_attribute{'gene_id'}){
					print STDERR "<ERROR> Missing gene_id in $temp_gtf!\n$temp_line\n";
					exit;
				}
				unless(exists $temp_attribute{'gene_name'}){
					push(@temp_attribute, "gene_name \"$temp_reclocus_gene_name\"\;");
				}
				unless(exists $temp_attribute{'gene_type'}){
					push(@temp_attribute, "gene_type \"$temp_reclocus_gene_type\"\;");
				}
				push(@temp_attribute, "reclocus \"$temp_gene_id\"\;");
				$temp_line[8] = join(" ", @temp_attribute);
				$temp_line = join("\t", @temp_line);
				print $temp_out_fh "$temp_line\n";
			}
		}
		else{
			print $temp_out_fh "$temp_line\n";
		}
	}
	close $temp_in_fh;
	foreach my $temp_gene_id (keys %temp_gene){
		print $temp_out_fh "$temp_gene{$temp_gene_id}{'chr'}\t$temp_gene{$temp_gene_id}{'source'}\tgene\t$temp_gene{$temp_gene_id}{'left'}\t$temp_gene{$temp_gene_id}{'right'}\t.\t$temp_gene{$temp_gene_id}{'strand'}\t.\tgene_id \"$temp_gene_id\"\; gene_name \"$temp_gene{$temp_gene_id}{'gene_name'}\"\; gene_type \"$temp_gene{$temp_gene_id}{'gene_type'}\"\; tag \"reconstructed_locus\"\;\n";
	}
	close $temp_out_fh;

	return 1;
}



