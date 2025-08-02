#!/usr/bin/env perl
=start
Author: ALin
Purpose: To get all CDSs from exons based on a .gtf files.
Change log:
	v1.2	2023-05	Adopted form gtf2fasta_exon_v1.4.pl.
	v1.3	2023-05	Added the output of amino acid sequences.
	v1.4	2023-12	Removed Bio::Perl.
=cut

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use threads;
use Benchmark qw(:hireswallclock);

my $total_start_time = Benchmark->new;

my $version = "v1.4";

my $usage = "perl gtf2fasta_CDS_${version}.pl
	-i <String> Input gtf file
	-o <String> Output base
	-g <String> The genome .fa file
	-t <Integer> Number of threads (Default: 1)
	-h <Boolean> Help
";

my ($gtf, $out, $genome, $num_thread, $help) = ("", "", "", 1, 0);

GetOptions(
	'i=s'	=>	\$gtf,
	'o=s'	=>	\$out,
	'g=s'	=>	\$genome,
	't=i'	=>	\$num_thread,
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

unless(-e $genome){
	print STDERR "<ERROR> Genome file, $genome, does ont exist!\n$usage\n";
	exit;
}

unless($out){
	print STDERR "<ERROR> Output file base missing!\n$usage";
	exit;
}

open(my $log_fh, "> ${out}.log") or die "Cannot create ${out}.log!\n";
my $date_format = "+%F-%T.%3N";
my $date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tperl gtf2fasta_exon_${version}.pl -i $gtf -o $out -g $genome -t $num_thread\n";

my $dir = "${out}__gtf2fasta_CDS_temp";
if(-d $dir){
	print STDERR "<WARNNING> Temporary directory $dir exists! Will overwrite it...\n";
}
else{
	system("mkdir $dir");
}

my $total_thread = `grep -c \"^processor\" /proc/cpuinfo`;
chomp $total_thread;
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tA total of $total_thread threads detected\n";

if($num_thread > $total_thread){
	print STDERR "<WARNING> Exceeding the maximal total number of threads! Setting the number of thread to $total_thread...";
	$num_thread = $total_thread;
}

my $start_time = Benchmark->new;
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tStarting preprocessing...\n";
my ($batch_ref) = split_gtf($gtf, $num_thread, $dir);
if($genome){
	split_genome($batch_ref, $gtf, $genome, $dir);
}

my $finish_time = Benchmark->new;
my $time_spent = timediff($finish_time, $start_time);
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tUsed ". timestr($time_spent) . " for preprocessing\n";

my $batch = 1;
my @thread = ();
my @running = ();

while(@thread < $num_thread){
	@running = threads->list(threads::running);
	if(@running < $num_thread){
		my $thread = threads->create(\&worker, $batch, $dir, $date_format, $log_fh);
		push(@thread, $thread);
		$batch++;
	}
}

foreach my $thread (@thread){
	my @result = $thread->join();
	$date = `date $date_format`;
	chomp $date;
	print $log_fh "[$date]\tJoined for batch $result[0]!\n";
}

system("cat ${dir}/*_cds.fa > ${out}_cds.fa");
system("cat ${dir}/*_aa.fa > ${out}_aa.fa");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tGenerated ${out}_cds.fa and ${out}_aa.fa!\n";
system("rm -r $dir");
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
	my ($temp_batch, $temp_dir, $temp_date_format, $temp_log_fh) = @_;
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $temp_log_fh "[$temp_date]\tStarting to process batch $temp_batch...\n";
	my $temp_start_time = Benchmark->new;
	my $temp_seq_in = Bio::SeqIO->new(-file => "${temp_dir}/batch_${temp_batch}.fa", -format => 'fasta');
	my %temp_seq_obj = ();
	while(my $temp_seq_obj = $temp_seq_in->next_seq){
		my $temp_id = $temp_seq_obj->primary_id;
		$temp_seq_obj{$temp_id} = $temp_seq_obj;
	}
	my $temp_transcipt_ref = read_gtf("${temp_dir}/batch_${temp_batch}.gtf");
	get_sequence(\%temp_seq_obj, $temp_transcipt_ref, "${temp_dir}/batch_${temp_batch}_cds.fa", "${temp_dir}/batch_${temp_batch}_aa.fa");
	my $temp_finish_time = Benchmark->new;
	my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $temp_log_fh "[$temp_date]\tDone for batch $temp_batch in " . timestr($temp_time_spent) , "!\n";

	return $temp_batch;
}

sub split_gtf{
	my ($temp_gtf, $temp_num_thread, $temp_dir) = @_;
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "<ERROR> Cannot open $temp_gtf!\n";
	my %temp_transcript = ();
	my @temp_transcript = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] eq 'CDS') || ($temp_line[2] eq 'stop_codon')){
			 if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
				my $temp_transcript_id = $1;
				unless(exists $temp_transcript{$temp_transcript_id}){
					@{$temp_transcript{$temp_transcript_id}} = ();
					push(@temp_transcript, $temp_transcript_id);
				}
				push(@{$temp_transcript{$temp_transcript_id}}, $temp_line);
			}
			else{
				print STDERR "<ERROR> Missing transcript_id in $temp_gtf!\n$temp_line\n";
				exit;
			}
		}
	}
	close $temp_fh;
	my $temp_batch_size = int((keys %temp_transcript) / $temp_num_thread) + 1;
	my ($temp_count, $temp_batch) = (0, 1);
	my %temp_batch = ();
	open($temp_fh, "> ${temp_dir}/batch_${temp_batch}.gtf") or die "<ERROR> Cannot create ${temp_dir}/batch_${temp_batch}.gtf!\n";
	foreach my $temp_transcript_id (@temp_transcript){
		$temp_count++;
		if($temp_count > ($temp_batch_size * $temp_batch)){
			close $temp_fh;
			$temp_batch++;
			unless(exists $temp_batch{$temp_batch}){
				%{$temp_batch{$temp_batch}} = ();
			}
			open($temp_fh, "> ${temp_dir}/batch_${temp_batch}.gtf");
		}
		foreach my $temp_line (@{$temp_transcript{$temp_transcript_id}}){
			print $temp_fh "$temp_line\n";
			my @temp_line = split('\t', $temp_line);
			unless(exists $temp_batch{$temp_batch}{$temp_line[0]}){
				$temp_batch{$temp_batch}{$temp_line[0]} = 1;
			}
		}
	}
	close $temp_fh;

	return (\%temp_batch);
}

sub split_genome{
	my ($temp_batch_ref, $temp_gtf, $temp_genome, $temp_dir) = @_;
	my $temp_seq_in = Bio::SeqIO->new(-file => $temp_genome, -format => 'fasta');
	my %temp_seq_out = ();
	my %temp_chr_gtf = ();
	my %temp_chr_count = ();
	foreach my $temp_batch (keys %$temp_batch_ref){
		$temp_seq_out{$temp_batch} = Bio::SeqIO->new(-file => "> ${temp_dir}/batch_${temp_batch}.fa", -format => 'fasta');
		foreach my $temp_chr (keys %{$temp_batch_ref->{$temp_batch}}){
			unless(exists $temp_chr_gtf{$temp_chr}){
				$temp_chr_gtf{$temp_chr} = 1;
			}
		}
	}
	while(my $temp_seq_obj = $temp_seq_in->next_seq){
		my $temp_chr = $temp_seq_obj->primary_id;
		foreach my $temp_batch (keys %$temp_batch_ref){
			if(exists $temp_batch_ref->{$temp_batch}{$temp_chr}){
				$temp_seq_out{$temp_batch}->write_seq($temp_seq_obj);
				unless(exists $temp_chr_count{$temp_chr}){
					$temp_chr_count{$temp_chr} = 1;
				}
			}
		}
	}
	my $temp_chr_gtf = keys %temp_chr_gtf;
	my $temp_chr_count = keys %temp_chr_count;
	if($temp_chr_gtf != $temp_chr_count){
		print STDERR "<ERROR> Sequences in $temp_genome not matching $temp_gtf!\n";
		exit;
	}

	return 1;
}

sub read_gtf{
	my ($temp_gtf)  = @_;
	open(my $temp_fh, $temp_gtf) or die "<ERROR> Cannot open $temp_gtf!\n";
	my %temp_transcript = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my ($temp_gene_id, $temp_transcript_id);
		if($temp_line[8] =~ /gene_id \"([^\"]+)\"/){
			$temp_gene_id = $1;
		}
		else{
			print STDERR "<ERROR> gene_id missing at line\n$temp_line\n";
			exit;
		}
		if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $1;
		}
		else{
			print STDERR "<ERROR> transcript_id missing at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
		}
		my @temp_exon = ($temp_line[3], $temp_line[4]);
		push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
	}
	close $temp_fh;

	return (\%temp_transcript);
}

sub get_sequence{
	my ($temp_seq_obj_ref, $temp_transcipt_ref, $temp_cds, $temp_aa) = @_;
	my $temp_cds_out = Bio::SeqIO->new(-file => ">$temp_cds", -format => 'fasta',);
	my $temp_aa_out = Bio::SeqIO->new(-file => ">$temp_aa", -format => 'fasta',);
	foreach my $temp_transcript_id (keys %$temp_transcipt_ref){
		my @temp_sorted_exon = sort {$a->[0] <=> $b->[0]} @{$temp_transcipt_ref->{$temp_transcript_id}{'exon'}};
		my $temp_desc = "GeneID=$temp_transcipt_ref->{$temp_transcript_id}{'gene_id'};Strand=$temp_transcipt_ref->{$temp_transcript_id}{'strand'};Pos=$temp_transcipt_ref->{$temp_transcript_id}{'chr'}:";
		my @temp_pos = ();
		my $temp_seq = "";
		for(my $i = 0; $i < @temp_sorted_exon; $i++){
			$temp_seq .= $temp_seq_obj_ref->{$temp_transcipt_ref->{$temp_transcript_id}{'chr'}}->subseq($temp_sorted_exon[$i][0], $temp_sorted_exon[$i][1]);
			push(@temp_pos, "$temp_sorted_exon[$i][0]-$temp_sorted_exon[$i][1]");
		}
		my $temp_pos = join(",", @temp_pos);
		$temp_desc .= $temp_pos;
		if($temp_transcipt_ref->{$temp_transcript_id}{'strand'} eq '-'){
			$temp_seq = reverse_complement($temp_seq);
		}
		my $temp_cds_seq_obj = Bio::Seq->new();
		$temp_cds_seq_obj->id($temp_transcript_id);
		$temp_cds_seq_obj->seq($temp_seq);
		$temp_cds_seq_obj->desc($temp_desc);
		my $temp_aa_seq_obj = $temp_cds_seq_obj->translate;
		$temp_cds_out->write_seq($temp_cds_seq_obj);
		$temp_aa_out->write_seq($temp_aa_seq_obj);
	}

	return 1;
}

sub reverse_complement{
	my ($temp_seq) = @_;
	my $temp_revcomp = reverse($temp_seq);
	$temp_revcomp =~ tr/ACGTacgt/TGCAtgca/;

	return $temp_revcomp;
}




