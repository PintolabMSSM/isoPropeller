#!/usr/bin/env perl
=start
Author: ALin
Purpose: To convert transdecoder result to gtf format.
Change log:
	v1.1	2023-11	Adopted from transdecoder_gff2gtf_v1.1.pl. Fix an issue in parssing attribute when reading a gtf. Added input of the transcript .fa and TIS efficiency file for getting TIS efficiency. Added minimal ORF length and NMD identification. Allowed multiple entries per transcript and used various metrics, such that score > stop codon > start codon > TIS efficiency > ORF length > 3' UTR length. Added attribute tis_efficiency in the gtf output.
=cut

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $version = "v1.1";

my ($gtf, $fa, $transdecoder, $out, $min_orf, $nmd, $tis, $help) = ("", "", "", "", 1, -1, "", 0);

my $usage = "Usage: perl transdecoder2gtf.pl
	-g <String> Transcript gtf
	-f <String> Transcript .fa
	-r <String> Transdecoder .cds
	-o <String> Output base name
	-m <Integer> Minimal ORF length (Default: $min_orf)
	-n <Integer> NMD cutoff (Default: $nmd)
	-e <String> TIS efficiency
	-h <Boolean> Help
";

GetOptions(
	'g=s'	=>	\$gtf,
	'f=s'	=>	\$fa,
	'r=s'	=>	\$transdecoder,
	'o=s'	=>	\$out,
	'm=i'	=>	\$min_orf,
	'n=i'	=>	\$nmd,
	'e=s'	=>	\$tis,
	'h!'	=>	\$help
);

unless($gtf && $fa && $transdecoder && $out){
	print "$usage";
	exit;
}

if($help){
	print "$usage";
	exit;
}

my $transcript_ref = read_gtf($gtf);
my $fa_ref = read_fa($fa);
my $tis_ref = read_tis($tis);
my $orf_ref = parse_transdecoder($transcript_ref, $fa_ref, $tis_ref, $transdecoder, $min_orf, $nmd);
print_result($transcript_ref, $orf_ref, $out);

sub read_tis{
	my ($temp_tis) = @_;
	my %temp_tis = ();
	if(-e $temp_tis){
		open(my $temp_fh, $temp_tis) or die "Cannot open $temp_tis!\n";
		while(<$temp_fh>){
			chomp;
			my $temp_line = $_;
			$temp_line =~ s/\r//;
			my @temp_line = split('\t', $temp_line);
			$temp_line[0] =~ tr/Uuatcg/TTATCG/;
			if(exists $temp_tis{$temp_line[0]}){
				print STDERR "<ERROR> Duplicated $temp_line[0] in $temp_tis!\n";
				exit;
			}
			else{
				$temp_tis{$temp_line[0]} = $temp_line[1];
			}
		}
		close $temp_fh;
	}

	return \%temp_tis;
}

sub read_fa{
	my ($temp_fa) = @_;
	my %temp_fa = ();
	my $temp_seq_in = Bio::SeqIO->new(-file => $temp_fa, -format => 'fasta',);
	while(my $temp_seq_obj = $temp_seq_in->next_seq){
		my $temp_transcript_id = $temp_seq_obj->id;
		if(exists $temp_fa{$temp_transcript_id}){
			print STDERR "<ERROR> Duplicated $temp_transcript_id in $temp_fa!\n";
			exit;
		}
		else{
			$temp_fa{$temp_transcript_id} = $temp_seq_obj;
		}
	}

	return \%temp_fa;
}

sub read_gtf{
	my ($temp_gtf) = @_;
	my %temp_transcript = ();
	open(my $temp_fh, $temp_gtf) or die "Cannot open $temp_gtf!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		unless($temp_line[2] eq "exon"){
			next;
		}
		my $temp_transcript_id;
		if($temp_line =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "Missing transcript_id in $temp_gtf for line\n$temp_line\n";
			exit;
		}
		my $temp_gene_id;
		if($temp_line =~ /(^|\s)gene_id \"([^\"]+)\"/){
			$temp_gene_id = $2;
		}
		else{
			print STDERR "Missing gene_id in $temp_gtf for line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			$temp_transcript{$temp_transcript_id}{'length'} = 0;
		}
		my @temp_exon = ($temp_line[3], $temp_line[4]);
		$temp_transcript{$temp_transcript_id}{'length'} += ($temp_line[4] - $temp_line[3] + 1);
		push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
	}
	close $temp_fh;

	return \%temp_transcript;
}	

sub parse_transdecoder{
	my ($temp_transcript_ref, $temp_fa_ref, $temp_tis_ref, $temp_transdecoder, $temp_min_orf, $temp_nmd) = @_;
	my %temp_orf = ();
	my %temp_start_codon = ('ATG', 1);
	my %temp_stop_codon = ('TAA', 1, 'TAG', 1, 'TGA', 1);
	my $temp_seq_in = Bio::SeqIO->new(-file => $temp_transdecoder, -format => 'fasta',);
	my %temp_orf = ();
	while(my $temp_seq_obj = $temp_seq_in->next_seq){
                my $temp_id = $temp_seq_obj->id;
                $temp_id =~ /^(.+)\.p\d+$/;
                my $temp_transcript_id = $1;
		my $temp_desc = $temp_seq_obj->desc;
		my $temp_score;
		if($temp_desc =~ /score=([\d\.-]+)/){
			$temp_score = $1;
		}
		else{
			print STDERR "<ERROR> Missing score in $temp_id!\n";
			exit;
		}
		my ($temp_orf_length, $temp_orf_start, $temp_orf_end, $temp_orf_strand);
		if($temp_desc =~ /len:([^\s]+) $temp_transcript_id:(\d+)-(\d+)\(([+-])\)/){
			($temp_orf_length, $temp_orf_start, $temp_orf_end, $temp_orf_strand) = ($1, $2, $3, $4);
		}
		else{
			print STDERR "<ERROR> Missing ORF information in $temp_id!\n";
			exit;
		}
		if($temp_orf_length < $temp_min_orf){
			next;
		}
		unless(exists $temp_orf{$temp_transcript_id}){
			@{$temp_orf{$temp_transcript_id}} = ();
		}
		#ORF: ORF strand, ORF length, 5' UTR length, 3' UTR length, canonical start codon, TIS efficiency, canonical stop codon, NMD, score
		my @temp_orf = ($temp_orf_strand, $temp_orf_length, $temp_orf_start - 1, $temp_transcript_ref->{$temp_transcript_id}{'length'} - $temp_orf_end, 0, 0, 0, 0, $temp_score);
		my $temp_start_codon = $temp_seq_obj->subseq(1, 3);
		$temp_start_codon =~ tr/Uuatcg/TTATCG/;
		my @temp_tis = ();
		if(exists $temp_start_codon{$temp_start_codon}){
			$temp_orf[4] = 1;
			@temp_tis = ($temp_orf_start - 6, $temp_orf_start + 4);
		}
		else{
			@temp_tis = ($temp_orf_start - 4, $temp_orf_start + 3);
		}
		if(exists $temp_fa_ref->{$temp_transcript_id}){
			if($temp_tis[0] > 0){
				my $temp_tis = $temp_fa_ref->{$temp_transcript_id}->subseq(@temp_tis);
				if(exists $temp_tis_ref->{$temp_tis}){
					$temp_orf[5] = $temp_tis_ref->{$temp_tis};
				}
				else{
					print STDERR "<WARNNING> Missing $temp_tis in the TIS file!\n";
				}
			}
		}
		else{
			print STDERR "Missing $temp_transcript_id in the transcript .fa file!\n";
			exit;
		}
		my $temp_stop_codon = $temp_seq_obj->subseq($temp_orf_length - 2, $temp_orf_length);
		$temp_stop_codon =~ tr/Uuatcg/TTATCG/;
		if(exists $temp_stop_codon{$temp_stop_codon}){
			$temp_orf[6] = 1;
		}
		my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_transcript_ref->{$temp_transcript_id}{'exon'}};
		if(@temp_sorted_exon > 1){
			my $temp_d3s;
			if($temp_transcript_ref->{$temp_transcript_id}{'strand'} eq '-'){
				$temp_d3s = $temp_orf[3] - ($temp_sorted_exon[0][1] - $temp_sorted_exon[0][0] + 1);
			}
			else{
				$temp_d3s = $temp_orf[3] - ($temp_sorted_exon[@temp_sorted_exon-1][1] - $temp_sorted_exon[@temp_sorted_exon-1][0] + 1);
			}
			if($temp_d3s > $temp_nmd){
				$temp_orf[7] = 1;
			}
		}
		push(@{$temp_orf{$temp_transcript_id}}, \@temp_orf);
	}

	return \%temp_orf;
}

sub print_result{
	my ($temp_transcript_ref, $temp_orf_ref, $temp_out) = @_;
	my %temp_cds_type = ("0", "Non-NMD", "1", "NMD");
	open(my $temp_gtf_fh, "> ${temp_out}.gtf") or die "Cannot create ${temp_out}.gtf!\n";
	open(my $temp_table_fh, "> ${temp_out}.txt") or die "Cannot create ${temp_out}.txt!\n";
	print $temp_table_fh "transcript_id\tlength\tcds_type\tcds_length\t5utr_length\t3utr_length\tcanonical_start\ttis_efficiency\tcanonical_stop\tscore\n";
	foreach my $temp_transcript_id (keys %$temp_orf_ref){
		#ORF: ORF strand, ORF length, 5' UTR length, 3' UTR length, canonical start codon, TIS efficiency, canonical stop codon, NMD, score
		my @temp_sorted_orf = sort {
			($b->[8] <=> $a->[8]) || 
			($b->[6] <=> $a->[6]) || 
			($b->[4] <=> $a->[4]) || 
			($b->[5] <=> $a->[5]) || 
			($b->[7] <=> $a->[7]) ||
			($b->[1] <=> $a->[1]) ||
			($b->[3] <=> $a->[3])
		} @{$temp_orf_ref->{$temp_transcript_id}};
		my $temp_cds_type = "Non-NMD";
		print $temp_table_fh "$temp_transcript_id\t$temp_transcript_ref->{$temp_transcript_id}{'length'}\t$temp_cds_type{$temp_sorted_orf[0][7]}\t$temp_sorted_orf[0][1]\t$temp_sorted_orf[0][2]\t$temp_sorted_orf[0][3]\t$temp_sorted_orf[0][4]\t$temp_sorted_orf[0][5]\t$temp_sorted_orf[0][6]\t$temp_sorted_orf[0][8]\n";
		my ($temp_strand, $temp_orf_left, $temp_orf_right) = orf_boundary($temp_transcript_ref->{$temp_transcript_id}{'strand'}, $temp_sorted_orf[0][0], $temp_sorted_orf[0][2] + 1, $temp_transcript_ref->{$temp_transcript_id}{'length'} - $temp_sorted_orf[0][3], $temp_transcript_ref->{$temp_transcript_id}{'length'});
		my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_transcript_ref->{$temp_transcript_id}{'exon'}};
		my $temp_moving_left = 0;
		my $temp_moving_right = 0;
		my $temp_moving_orf = 0;
		EXON:foreach my $temp_exon_ref (@temp_sorted_exon){
			my $temp_exon_length = $temp_exon_ref->[1] - $temp_exon_ref->[0] + 1;
			if($temp_moving_left == 0){
				$temp_moving_left = 1;
			}
			else{
				$temp_moving_left = $temp_moving_right + 1;
			}
			$temp_moving_right += $temp_exon_length;
			if(($temp_orf_left >= $temp_moving_left) && ($temp_orf_right <= $temp_moving_right)){
				my $temp_left = $temp_orf_left - $temp_moving_left + $temp_exon_ref->[0];
				my $temp_right = $temp_orf_right - $temp_moving_left + $temp_exon_ref->[0];
				my $temp_frame = check_frame($temp_moving_orf, $temp_strand, $temp_left, $temp_right);
				print $temp_gtf_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\ttransdecoder\tCDS\t$temp_left\t$temp_right\t.\t$temp_strand\t$temp_frame\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; tis_efficiency \"$temp_sorted_orf[0][5]\"\;\n";
				last EXON;
			}
			elsif($temp_orf_left > $temp_moving_right){
				next EXON;
			}
			elsif(($temp_orf_left >= $temp_moving_left) && ($temp_orf_left <= $temp_moving_right)){
				my $temp_left = $temp_orf_left - $temp_moving_left + $temp_exon_ref->[0];
				my $temp_frame = check_frame($temp_moving_orf, $temp_strand, $temp_left, $temp_exon_ref->[1]);
				print $temp_gtf_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\ttransdecoder\tCDS\t$temp_left\t$temp_exon_ref->[1]\t.\t$temp_strand\t$temp_frame\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; tis_efficiency \"$temp_sorted_orf[0][5]\"\;\n";
				$temp_moving_orf += ($temp_exon_ref->[1] - $temp_left + 1);
			}
			elsif(($temp_orf_right >= $temp_moving_left) && ($temp_orf_right <= $temp_moving_right)){
				my $temp_right = $temp_orf_right - $temp_moving_left + $temp_exon_ref->[0];
				my $temp_frame = check_frame($temp_moving_orf, $temp_strand, $temp_exon_ref->[0], $temp_right);
				print $temp_gtf_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\ttransdecoder\tCDS\t$temp_exon_ref->[0]\t$temp_right\t.\t$temp_strand\t$temp_frame\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; tis_efficiency \"$temp_sorted_orf[0][5]\"\;\n";
				last EXON;
			}
			else{
				my $temp_frame = check_frame($temp_moving_orf, $temp_strand, $temp_exon_ref->[0], $temp_exon_ref->[1]);
				print $temp_gtf_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\ttransdecoder\tCDS\t$temp_exon_ref->[0]\t$temp_exon_ref->[1]\t.\t$temp_strand\t$temp_frame\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; tis_efficiency \"$temp_sorted_orf[0][5]\"\;\n";
				$temp_moving_orf += ($temp_exon_ref->[1] - $temp_exon_ref->[0] + 1);
			}
		}
	}
	close $temp_gtf_fh;

	return 1;
}

sub orf_boundary{
	my ($temp_transcript_strand, $temp_orf_strand, $temp_orf_start, $temp_orf_end, $temp_transcript_length) = @_;
	my ($temp_strand, $temp_orf_left, $temp_orf_right) = ("+", $temp_orf_start, $temp_orf_end);
	if($temp_transcript_strand ne $temp_orf_strand){
		$temp_strand = "-";
		$temp_orf_left = $temp_transcript_length - $temp_orf_end + 1;
		$temp_orf_right = $temp_transcript_length - $temp_orf_start  + 1;
	}

	return ($temp_strand, $temp_orf_left, $temp_orf_right);
}

sub check_frame{
	my ($temp_moving_orf, $temp_strand, $temp_left, $temp_right) = @_;
	my $temp_frame;
	if($temp_strand eq "+"){
		$temp_frame = (3 - ($temp_moving_orf % 3)) % 3;
	}
	elsif($temp_strand eq "-"){
		$temp_frame = ($temp_moving_orf + ($temp_right - $temp_left + 1)) % 3;
	}

	return $temp_frame;
}



