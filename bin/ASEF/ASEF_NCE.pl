#!/usr/bin/env perl
=start
Author: ALin
Purpose: This is a script for querying novel coding exon against a reference. Both query and reference files are in .gtf format. It is part of the Alternative Splicing Event Finder (ASEF) package.
Change log:
	v1.1	2022-04	It was adoptted from ASEF_NE_v1.1.pl.
	v1.2	2022-11 It doesn't output gene level information but generate a .gtf file for the novel coding exons. It will consider whether the stop codons is included in the the query and reference CDS features. It enables multithreading.
	v1.3	2023-03	Break down NCE into N5NCE, N3NCE, KENCE and NENCE.
	v1.4	2023-03	Allow existing tempoary and overwrite it.
	v1.5	2023-12 Fixed an issue in parsing gtf attribute.
=cut

use strict;
use Getopt::Long;
use threads;
use Benchmark qw(:hireswallclock);

my $total_start_time = Benchmark->new;

my $version = "v1.5";

my $usage = "Usage: perl ASEF_NCE_${version}.pl
	-q <String> .gtf of query
	-p <Boolean> The stop codoon is included in the CDS feature of the query.
	-r <String> .gtf of reference
	-s <Boolean> The stop codoon is included in the CDS feature of the reference.
	-o <String> Output base
	-t <Integer> Number of threads
	-b <String> Path for bedtools
	-h <Boolean> Help
";

my ($query, $query_stop_codon, $ref, $ref_stop_codon, $out, $num_thread, $btp, $help) = ("", 0, "", 0, "", 1, "", 0);

GetOptions(
	'q=s'	=>	\$query,
	'p!'	=>	\$query_stop_codon,
	'r=s'	=>	\$ref,
	's!'	=>	\$ref_stop_codon,
	'o=s'	=>	\$out,
	't=i'	=>	\$num_thread,
	'b=s'	=>	\$btp,
	'h!'    =>      \$help,
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

open(my $log_fh, "> ${out}.log") or die "Cannot create ${out}.log!\n";
my $date_format = "+%F-%T.%3N";
my $date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tperl ASEF_NCE_${version}.pl -q $query -r $ref -o $out -b $btp -t $num_thread";
if($query_stop_codon){
	print $log_fh " -p";
}
if($ref_stop_codon){
	print $log_fh " -s";
}
print $log_fh "\n";

my $dir = "${out}__ASEF_NCE_temp";
if(-d $dir){
	print STDERR "<WARNING> Temporary directory $dir exists! Will overwrite it..\n";
}
else{
	system("mkdir $dir");
}

my $total_thread = `grep -c \"^processor\" /proc/cpuinfo`;
chomp $total_thread;
$date = `date $date_format`;
chomp $date;
print  $log_fh "[$date]\tA total of $total_thread threads detected\n";

if($num_thread > $total_thread){
	print STDERR "<WARNING> Exceeding the maximal total number of threads! Setting the number of thread to $total_thread...";
	$num_thread = $total_thread;
}

my $start_time = Benchmark->new;
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tStarting preprocessing...\n";
my $chr_ref = split_gtf($ref, $query, $dir);
my $finish_time = Benchmark->new;
my $time_spent = timediff($finish_time, $start_time);
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tUsed ". timestr($time_spent) . " for preprocessing\n";

my $chr_pos = 0;
my @thread = ();
my @running = ();

while(@thread < @$chr_ref){
	@running = threads->list(threads::running);
	if(@running < $num_thread){
		my $thread = threads->create(\&worker, $chr_ref->[$chr_pos], $query_stop_codon, $ref_stop_codon, $log_fh, $dir, $date_format, $btp);
		push(@thread, $thread);
		$chr_pos++;
	}
}

foreach my $thread (@thread){
	my @result = $thread->join();
	$date = `date $date_format`;
	chomp $date;
	print $log_fh "[$date]\tJoined for $result[0]!\n";
}

open(my $out_fh, "> ${out}_transcript.txt") or die "<ERROR> Cannot create ${out}_transcript.txt!\n";
print $out_fh "transcript_id\tN5NCE\tN3NCE\tKENCE\tNENCE\tFSCE\n";
close $out_fh;
system("cat ${dir}/*_out.txt >> ${out}_transcript.txt");
system("cat ${dir}/*_novel_cds.gtf | sort -k 1,1 -k 4,4n > ${out}_cds.gtf");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tGenerated ${out}_transcript.txt and ${out}_cds.gtf!\n";
system("rm -r ${dir}");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tCleaned temporary files!\n";

my $total_finish_time = Benchmark->new;
my $total_time_spent = timediff($total_finish_time, $total_start_time);
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tTotal time used: ". timestr($total_time_spent) . "\n";
close $log_fh;

sub worker{
	my ($temp_chr, $temp_query_stop_codon, $temp_ref_stop_codon, $temp_log_fh, $temp_dir, $temp_date_format, $temp_btp) = @_;
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $temp_log_fh "[$temp_date]\tStarting to process $temp_chr...\n";
	my $temp_start_time = Benchmark->new;
	my $temp_ref_transcript_ref = read_gtf("${temp_dir}/ref_${temp_chr}.gtf");
	unless($temp_ref_stop_codon){
		add_stop_codon($temp_ref_transcript_ref);
	}
	my ($temp_ref_cds_ref, $temp_ref_cds_end_ref) = parse_ref($temp_ref_transcript_ref);
	my %temp_novel = ();
	my $temp_query_transcript_ref = read_gtf("${temp_dir}/query_${temp_chr}.gtf", \%temp_novel);
	unless($temp_query_stop_codon){
		add_stop_codon($temp_query_transcript_ref);
	}
	parse_query($temp_query_transcript_ref, \%temp_novel, $temp_ref_cds_ref, $temp_ref_cds_end_ref, $temp_chr, "${temp_dir}/${temp_chr}_novel_cds.gtf");
	print_result(\%temp_novel, "${temp_dir}/${temp_chr}_out.txt");
	my $temp_finish_time = Benchmark->new;
	my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $temp_log_fh "[$temp_date]\tDone for $temp_chr in " . timestr($temp_time_spent) , "!\n";
	
	return $temp_chr;
}

sub split_gtf{
	my ($temp_ref, $temp_query, $temp_dir) = @_;
	system("cut -f1 $temp_query | sort | uniq > ${temp_dir}/chr_list.txt");
	my @temp_chr = ();
	open(my $temp_fh, "${temp_dir}/chr_list.txt") or die "<ERROR> Cannot open ${temp_dir}/chr_list.txt!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		push(@temp_chr, $temp_line);
	}
	close $temp_fh;
	system("for i in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$i	\" $temp_ref > ${temp_dir}/ref_\$\{i\}.gtf; done");
	system("for i in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$i	\" $temp_query > ${temp_dir}/query_\$\{i\}.gtf; done");

	return \@temp_chr;
}

sub read_gtf{
	my ($temp_gtf, $temp_query_ref) = @_;
	my %temp_transcript = ();
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "<ERROR> Cannot open $temp_gtf!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] ne "exon") && ($temp_line[2] ne "CDS")){
			next;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> transcript_id missing at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			my $temp_gene_id;
			if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"/){
				$temp_gene_id = $2;
			}
			else{
				print STDERR "<ERROR> gene_id missing at line\n$temp_line\n";
				exit;
			}
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			@{$temp_transcript{$temp_transcript_id}{'cds'}} = ();
			@{$temp_transcript{$temp_transcript_id}{'frame'}} = ();
		}
		if(defined $temp_query_ref){
			unless(exists $temp_query_ref->{$temp_transcript_id}){
				%{$temp_query_ref->{$temp_transcript_id}} = ();
				$temp_query_ref->{$temp_transcript_id}{'N5NCE'} = 0;
				$temp_query_ref->{$temp_transcript_id}{'N3NCE'} = 0;
				$temp_query_ref->{$temp_transcript_id}{'KENCE'} = 0;
				$temp_query_ref->{$temp_transcript_id}{'NENCE'} = 0;
				$temp_query_ref->{$temp_transcript_id}{'FSCE'} = 0;
			}
		}
		if($temp_line[2] eq "exon"){
			my @temp_exon = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
		}
		elsif($temp_line[2] eq "CDS"){
			my @temp_cds = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'cds'}}, \@temp_cds);
			push(@{$temp_transcript{$temp_transcript_id}{'frame'}}, $temp_line[7]);
		}

	}
	close $temp_fh;

	return \%temp_transcript;
}

sub add_stop_codon{
	my ($temp_transcript_ref) = @_;
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		my $temp_num_stop_codon_base = 0;
		if(@{$temp_transcript_ref->{$temp_transcript_id}{'cds'}} > 0){
			for(my $i = 0; $i < 3; $i++){
				$temp_num_stop_codon_base += extend_cds_by_one_base($i, \@{$temp_transcript_ref->{$temp_transcript_id}{'exon'}}, \@{$temp_transcript_ref->{$temp_transcript_id}{'cds'}}, \@{$temp_transcript_ref->{$temp_transcript_id}{'frame'}}, $temp_transcript_ref->{$temp_transcript_id}{'strand'});
			}
			if($temp_num_stop_codon_base == 0){
				print STDERR "<WARNING> Incomplete stop codon found in $temp_transcript_id!\n";
			}
		}
	}

	return 1;
}

sub extend_cds_by_one_base{
	my ($temp_base_pos, $temp_exon_ref, $temp_cds_ref, $temp_frame_ref, $temp_strand) = @_;
	if($temp_strand eq "-"){
		for(my $i = 0; $i < @$temp_exon_ref; $i++){
			if(($temp_cds_ref->[0][0] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[0][0] <= $temp_exon_ref->[$i][1])){
				if(($temp_cds_ref->[0][0] - 1)  >= $temp_exon_ref->[$i][0]){
					$temp_cds_ref->[0][0]--;

					return 1;
				}
				else{
					if($i == 0){

						return 0;
					}
					else{
						my @temp_exon = ($temp_exon_ref->[$i-1][1], $temp_exon_ref->[$i-1][1]);
						my $temp_frame = (3 - $temp_base_pos) % 3;
						unshift(@$temp_cds_ref, \@temp_exon);
						unshift(@$temp_frame_ref, $temp_frame);

						return 1;
					}
				}
			}
		}
	}
	else{
		for(my $i = (@$temp_exon_ref - 1); $i >= 0; $i--){
			if(($temp_cds_ref->[@$temp_cds_ref-1][1] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[@$temp_cds_ref-1][1] <= $temp_exon_ref->[$i][1])){
				if(($temp_cds_ref->[@$temp_cds_ref-1][1] + 1) <= $temp_exon_ref->[$i][1]){
					$temp_cds_ref->[@$temp_cds_ref-1][1]++;

					return 1;
				}
				else{
					if($i == @$temp_exon_ref - 1){

						return 0;
					}
					else{
						my @temp_exon = ($temp_exon_ref->[$i+1][0], $temp_exon_ref->[$i+1][0]);
						my $temp_frame = (3 - $temp_base_pos) % 3;
						push(@$temp_cds_ref, \@temp_exon);
						push(@$temp_frame_ref, $temp_frame);

						return 1;
					}
				}
			}
		}
	}

	return 0;
}

sub parse_ref{
	my ($temp_ref_transcript_ref) = @_;
	my %temp_ref_cds = ();
	my %temp_ref_cds_end = ();
	foreach my $temp_transcript_id (keys %$temp_ref_transcript_ref){
		my $temp_strand = $temp_ref_transcript_ref->{$temp_transcript_id}{'strand'};
		my @temp_cds = @{$temp_ref_transcript_ref->{$temp_transcript_id}{'cds'}};
		my @temp_frame = @{$temp_ref_transcript_ref->{$temp_transcript_id}{'frame'}};
		unless(exists $temp_ref_cds{$temp_strand}){
			%{$temp_ref_cds{$temp_strand}} = ();
			%{$temp_ref_cds_end{$temp_strand}} = ();
			%{$temp_ref_cds_end{$temp_strand}{'5p'}} = ();
			%{$temp_ref_cds_end{$temp_strand}{'3p'}} = ();
		}
		for(my $i = 0; $i < @temp_cds; $i++){
			unless(exists $temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}){
				%{$temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}} = ();
			}
			unless(exists $temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}){
				%{$temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}} = ();
			}
			unless(exists $temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}{$temp_frame[$i]}){
				$temp_ref_cds{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}{$temp_frame[$i]} = 1;
			}
			my ($temp_5p, $temp_3p) = ($temp_cds[$i][0], $temp_cds[$i][1]);
			if($temp_strand eq '-'){
				($temp_5p, $temp_3p) = ($temp_cds[$i][1], $temp_cds[$i][0]);
			}
			unless(exists $temp_ref_cds_end{$temp_strand}{'5p'}{$temp_5p}){
				%{$temp_ref_cds_end{$temp_strand}{'5p'}{$temp_5p}} = ();
			}
			unless(exists $temp_ref_cds_end{$temp_strand}{'5p'}{$temp_5p}{$temp_frame[$i]}){
				$temp_ref_cds_end{$temp_strand}{'5p'}{$temp_5p}{$temp_frame[$i]} = 1;
			}
			unless(exists $temp_ref_cds_end{$temp_strand}{'3p'}{$temp_3p}){
				%{$temp_ref_cds_end{$temp_strand}{'3p'}{$temp_3p}} = ();
			}
			unless(exists $temp_ref_cds_end{$temp_strand}{'3p'}{$temp_3p}{$temp_frame[$i]}){
				$temp_ref_cds_end{$temp_strand}{'3p'}{$temp_3p}{$temp_frame[$i]} = 1;
			}
		}
	}

	return (\%temp_ref_cds, \%temp_ref_cds_end);
}

sub parse_query{
	my ($temp_query_transcript_ref, $temp_novel_ref, $temp_ref_cds_ref, $temp_ref_cds_end_ref, $temp_chr, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "<ERROR> Cannot create $temp_out!\n";
	foreach my $temp_transcript_id (keys %$temp_query_transcript_ref){
		my $temp_gene_id = $temp_query_transcript_ref->{$temp_transcript_id}{'gene_id'};
		my $temp_strand = $temp_query_transcript_ref->{$temp_transcript_id}{'strand'};
		my @temp_cds = @{$temp_query_transcript_ref->{$temp_transcript_id}{'cds'}};
		my @temp_frame = @{$temp_query_transcript_ref->{$temp_transcript_id}{'frame'}};
		for(my $i = 0; $i < @temp_cds; $i++){
			my %temp_novel_tag = ();
			if(exists $temp_ref_cds_ref->{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}){
				unless(exists $temp_ref_cds_ref->{$temp_strand}{$temp_cds[$i][0]}{$temp_cds[$i][1]}{$temp_frame[$i]}){
					$temp_novel_ref->{$temp_transcript_id}{'FSCE'}++;
					unless(exists $temp_novel_tag{'FSCE'}){
						$temp_novel_tag{'FSCE'} = 1;
					}
				}
			}
			else{
				my ($temp_5p, $temp_3p) = ($temp_cds[$i][0], $temp_cds[$i][1]);
				if($temp_strand eq '-'){
					($temp_5p, $temp_3p) = ($temp_cds[$i][1], $temp_cds[$i][0]);
				}
				if((exists $temp_ref_cds_end_ref->{$temp_strand}{'5p'}{$temp_5p}) && (exists $temp_ref_cds_end_ref->{$temp_strand}{'3p'}{$temp_3p})){
					$temp_novel_ref->{$temp_transcript_id}{'KENCE'}++;
					unless(exists $temp_novel_tag{'KENCE'}){
						$temp_novel_tag{'KENCE'} = 1;
					}
				}
				elsif(exists $temp_ref_cds_end_ref->{$temp_strand}{'5p'}{$temp_5p}){
					$temp_novel_ref->{$temp_transcript_id}{'N3NCE'}++;
					unless(exists $temp_novel_tag{'N3NCE'}){
						$temp_novel_tag{'N3NCE'} = 1;
					}
				}
				elsif(exists $temp_ref_cds_end_ref->{$temp_strand}{'3p'}{$temp_3p}){
					$temp_novel_ref->{$temp_transcript_id}{'N5NCE'}++;
					unless(exists $temp_novel_tag{'N5NCE'}){
						$temp_novel_tag{'N5NCE'} = 1;
					}
				}
				else{
					$temp_novel_ref->{$temp_transcript_id}{'NENCE'}++;
					unless(exists $temp_novel_tag{'NENCE'}){
						$temp_novel_tag{'NENCE'} = 1;
					}
				}
			}
			my $temp_novel_tag = keys %temp_novel_tag;
			if($temp_novel_tag > 0){
				print $temp_fh "$temp_chr\tASEF\tCDS\t$temp_cds[$i][0]\t$temp_cds[$i][1]\t.\t$temp_strand\t$temp_frame[$i]\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_gene_id\"\;";
				foreach my $temp_tag ("N5NCE", "N3NCE", "KENCE", "NENCE", "FSCE"){
					if(exists $temp_novel_tag{$temp_tag}){
						print $temp_fh " tag \"$temp_tag\"\;";
					}
				}
				print $temp_fh "\n";
			}
		}
	}

	return 1;
}

sub print_result{
	my ($temp_novel_ref, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "<ERROR> Cannot create $temp_out!\n";
	foreach my $temp_transcript_id (keys %$temp_novel_ref){
		print $temp_fh "$temp_transcript_id";
		for my $temp_col ("N5NCE", "N3NCE", "KENCE", "NENCE", "FSCE"){
			print $temp_fh "\t$temp_novel_ref->{$temp_transcript_id}{$temp_col}";
		}
		print $temp_fh "\n";
	}
	close $temp_fh;

	return 1;
}



