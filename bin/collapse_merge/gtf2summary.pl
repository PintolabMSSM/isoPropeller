#!/usr/bin/env perl
=start
Author: ALin
Purpose: To summarize information in transcript and gene level.
Change log:
	v1.1	2023-02	This script was adopted from gtf2transcript_summary_v1.7.pl.
	v1.2	2023-02	Added the function to look for external splice junction support. Removed the random seed. Ouput the log to a .log file.
	v1.3	2023-03 Included the gene feature during parsing. Added gene_name, gene_type, gene and gene_biotype as gene columns if they are present in the attribute list.
	v1.4	2023-03	Added the CDS exon number column for coding transcripts.
	v1.5	2023-03	Fixed a bug for having coordinate longer than the sequence length.
	v1.6	2023-04	Added the genomic poisitions of CDS start and end. Fixed a bug in getting the correct CDS length for stop codon not in the CDS feature.
	v1.7	2023-04	Enabled to retrieve tag in gene feature.
	v1.8	2023-07	Fixed a bug for getting attributes only present in transcript but not gene.
	v1.9	2023-09	Restructured the code. Added the number of longest continuous "A"s in the output when -g is specified.
	v1.10	2023-09	Added the GC content in the output  when -g is specified.
	v1.11	2023-10	Parallelize jobs by number of transcripts instead of chromosomes. Fixed an issue of exon order.
	v1.12	2023-11	Fixed an issue in parsing gtf attribute.
	v1.13	2024-03	Allowed attribute to omit the ".
	v1.14	2024-05	Allowed the numbers of chromosomes being different between the gtf and .fa files by only giving a warning. Prevented batch with missing sequences being lost.
=cut

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use threads;
use Benchmark qw(:hireswallclock);

my $total_start_time = Benchmark->new;

my $version = 'v1.14';

my ($gtf, $out, $attribute, $tag, $stop_codon_flag, $nmd, $genome, $window, $sj, $sj_intron, $sj_gtf, $num_thread, $help) = ('', '', '', '', 0, -1, '', 20, '', 0, 0, 1, 0);

my $usage = "Usage: perl gtf2summary.pl
	-i <String> Input gtf file
	-o <String> Output base
	-a <String> The list of attribute to be retrieved
	-d <String> The list of transcript tag to be retrieved
	-s <Boolean> The stop codoon is included in the CDS feature.
	-n <Integer> The distance between the stop codon (exclusive) and the downstream 3' most splice site for defining NMD
	-g <String> The genome .fa file for finding fraction of \"A\"s downstream of transcripts
	-w <Integer> The window size for finding fraction of \"A\"s downstream of transcripts (Default: $window)
	-j <String> The splice junction support file (Default: STAR SJ.out.tab format)
	-r <Boolean> Splice junction as intron
	-f <Boolean> The splice junction support file in gtf format
	-t <Integer> Number of threads (Default: $num_thread)
	-h <Boolean> Help
";

GetOptions(
	'i=s'	=>	\$gtf,
	'o=s'	=>	\$out,
	'a=s'	=>	\$attribute,
	'd=s'	=>	\$tag,
	's!'	=>	\$stop_codon_flag,
	'n=i'	=>	\$nmd,
	'g=s'	=>	\$genome,
	'w=i'	=>	\$window,
	'j=s'	=>	\$sj,
	'r!'	=>	\$sj_intron,
	'f!'	=>	\$sj_gtf,
	't=i'	=>	\$num_thread,
	'h!'	=>	\$help,
);

if($help){
        print "$usage";
        exit;
}

unless(-e $gtf){
	print STDERR "<ERROR> .gtf file, $gtf, does not exist!\n$usage";
	exit;
}

unless($out){
	print STDERR "<ERROR> Output file missing!\n$usage";
	exit;
}

if($sj ne ''){
	unless(-e $sj){
		print STDERR "<ERROR> The splice junction support file, $sj, does not exist!\n$usage";
		exit;
	}
}

open(my $log_fh, "> ${out}.log") or die "Cannot create ${out}.log!\n";
my $date_format = "+%F-%T.%3N";
my $date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tperl gtf2transcript_summary_${version}.pl -i $gtf -o $out";

if($attribute){
	unless(-e $attribute){
		print STDERR "<ERROR> Attribute list, $attribute, does not exist!\n$usage";
		exit;
	}
	print $log_fh " -a $attribute";
}

if($tag){
	unless(-e $tag){
		print STDERR "<ERROR> Tag list, $tag, does not exist!\n$usage";
		exit;
	}
	print $log_fh " -d $tag";
}

if($stop_codon_flag){
	print $log_fh " -s";
}

if($nmd >= 0){
	print $log_fh " -n $nmd";
}

if($genome){
	print $log_fh " -g $genome -w $window";
}

if($sj ne ''){
	print $log_fh " -j $sj";
	if($sj_intron){
		print $log_fh " -r";
	}
	if($sj_gtf){
		print $log_fh " -f";
	}	
}

print $log_fh " -t $num_thread\n";

my $dir = "${out}__gtf2summary_temp";
if(-d $dir){
	print STDERR "<WARNNING> Temporary directory $dir exists! Will overwrite it...\n";
	system("rm ${dir}/*");
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
my $batch_ref = split_gtf($gtf, $num_thread, $dir);
if($genome){
	split_genome($batch_ref, $gtf, $genome, $dir);
}
if($sj){
	split_sj($batch_ref, $sj, $dir);
}
my $finish_time = Benchmark->new;
my $time_spent = timediff($finish_time, $start_time);
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tUsed ". timestr($time_spent) . " for preprocessing\n";
my $batch = 1;
my @thread = ();
my @running = ();

while(@thread < keys %$batch_ref){
	@running = threads->list(threads::running);
	if(@running < $num_thread){
		my $thread = threads->create(\&worker, $batch, $dir, $date_format, $attribute, $tag, $stop_codon_flag, $nmd, $genome, $window, $sj_intron, $sj_gtf);
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

open(my $out_transcript_fh, "> ${out}_transcript.txt") or die "<ERROR> Cannot create ${out}_transcript.txt!\n";
open(my $out_gene_fh, "> ${out}_gene.txt") or die "<ERROR> Cannot create ${out}_gene.txt!\n";
print $out_transcript_fh "transcript_id\tgene_id\tchr\tstrand\tleft\tright\texon\tsj_support\tlength\tcds_type\tcds_exon\tcds_length\t5utr_length\t3utr_length\tstart_codon\tstop_codon\tcds_left\tcds_right\tfraction_downstream_A\tlength_downstream_A\tgc_content";
print $out_gene_fh "gene_id\tchr\tstrand\tleft\tright\ttranscript\tnon_nmd\tnmd\tnon_coding";
my $attribute_ref = read_list($attribute);
my $tag_ref = read_list($tag);
foreach my $attr (sort {$a cmp $b} keys %$attribute_ref){
        print $out_transcript_fh "\t$attr";
	if(($attr eq 'gene_name') || ($attr eq 'gene_type') || ($attr eq 'gene') || ($attr eq 'gene_biotype') || ($attr eq 'level')){
		print $out_gene_fh "\t$attr";
	}
}
foreach my $t (sort {$a cmp $b} keys %$tag_ref){
	print $out_transcript_fh "\t$t";
	print $out_gene_fh "\t$t";
}
print $out_transcript_fh "\n";
print $out_gene_fh "\n";
close $out_transcript_fh;
close $out_gene_fh;
system("cat ${dir}/*_transcript.txt | sort -k 3,3 -k 5,5n >> ${out}_transcript.txt");
system("cat ${dir}/*_gene.txt | sort -k 2,2 -k 4,4n >> ${out}_gene.txt");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tGenerated $out!\n";
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
	my ($temp_batch, $temp_dir, $temp_date_format, $temp_attribute, $temp_tag, $temp_stop_codon_flag, $temp_nmd, $temp_genome, $temp_window, $temp_sj_intron, $temp_sj_gtf) = @_;
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $log_fh "[$temp_date]\tStarting to process batch $temp_batch...\n";
	my $temp_start_time = Benchmark->new;
	my $temp_attribute_ref = read_list($temp_attribute);
	my $temp_tag_ref = read_list($temp_tag);
	my ($temp_transcript_ref, $temp_gene_ref) = read_gtf("${temp_dir}/batch_${temp_batch}.gtf", $temp_attribute_ref, $temp_tag_ref);
	print_summary($temp_transcript_ref, $temp_gene_ref, "${temp_dir}/batch_${temp_batch}.fa", "${temp_dir}/batch_${temp_batch}_transcript.txt", "${temp_dir}/batch_${temp_batch}_gene.txt", $temp_attribute_ref, $temp_tag_ref, $temp_stop_codon_flag, $temp_nmd, $temp_window, "${temp_dir}/batch_${temp_batch}.sj", $temp_sj_intron, $temp_sj_gtf);
	my $temp_finish_time = Benchmark->new;
	my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $log_fh "[$temp_date]\tDone for batch $temp_batch in " . timestr($temp_time_spent) , "!\n";

	return $temp_batch;

}

sub split_gtf{
	my ($temp_gtf, $temp_num_thread, $temp_dir) = @_;
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "<ERROR> Cannot open $temp_gtf!\n";
	my %temp_gene = ();
	my @temp_gene = ();
	my %temp_transcript_count = ();
	my $temp_transcript_count = 0;
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] eq 'gene') || ($temp_line[2] eq 'transcript') || ($temp_line[2] eq 'exon') || ($temp_line[2] eq 'CDS')){
			 if($temp_line[8] =~ /(^|\s)gene_id \"?([^";]+)\"?/){
				my $temp_gene_id = $2;
				unless(exists $temp_gene{$temp_gene_id}){
					@{$temp_gene{$temp_gene_id}} = ();
					push(@temp_gene, $temp_gene_id);
					%{$temp_transcript_count{$temp_gene_id}} = ();
				}
				push(@{$temp_gene{$temp_gene_id}}, $temp_line);
				if($temp_line[8] =~ /(^|\s)transcript_id \"?([^";]+)\"?/){
					my $temp_transcript_id = $2;
					unless(exists $temp_transcript_count{$temp_gene_id}{$temp_transcript_id}){
						$temp_transcript_count{$temp_gene_id}{$temp_transcript_id} = 1;
						$temp_transcript_count++;
					}
				}
			}
			else{
				print STDERR "<ERROR> Missing gene_id in $temp_gtf!\n$temp_line\n";
				exit;
			}
		}
	}
	close $temp_fh;
	my $temp_batch_size = int($temp_transcript_count / $temp_num_thread) + 1;
	my ($temp_count, $temp_batch) = (0, 1);
	my %temp_batch = ();
	open($temp_fh, "> ${temp_dir}/batch_${temp_batch}.gtf") or die "<ERROR> Cannot create ${temp_dir}/batch_${temp_batch}.gtf!\n";
	foreach my $temp_gene_id (@temp_gene){
		$temp_count += keys %{$temp_transcript_count{$temp_gene_id}};
		if($temp_count > ($temp_batch_size * $temp_batch)){
			close $temp_fh;
			$temp_batch++;
			unless(exists $temp_batch{$temp_batch}){
				%{$temp_batch{$temp_batch}} = ();
			}
			open($temp_fh, "> ${temp_dir}/batch_${temp_batch}.gtf");
		}
		foreach my $temp_line (@{$temp_gene{$temp_gene_id}}){
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
		print STDERR "<WARNING> Sequences in $temp_genome not matching $temp_gtf! Transcripts in inconsistent sequences will be lost.\n";
	}

	return 1;
}

sub split_sj{
	my ($temp_batch_ref, $temp_sj, $temp_dir) = @_;
	open(my $temp_in_fh, $temp_sj) or die "Cannot open $temp_sj!\n";
	my %temp_out_fh = ();
	foreach my $temp_batch (keys %$temp_batch_ref){
		open($temp_out_fh{$temp_batch}, "> ${temp_dir}/batch_${temp_batch}.sj") or die "Cannot create ${temp_dir}/batch_${temp_batch}.sj!\n";
		foreach my $temp_chr (keys %{$temp_batch_ref->{$temp_batch}}){
			system("grep \"^${temp_chr}\"\$\'\\t\' $temp_sj >> ${temp_dir}\/batch_${temp_batch}.sj");
		}
		close $temp_out_fh{$temp_batch};
	}

	return 1;
}

sub read_list{
	my ($temp_list) = @_;
	my %temp_list = ();
	if($temp_list){
		open(my $temp_fh, $temp_list) or die "<ERROR> Cannot open $temp_list!\n";
		while(<$temp_fh>){
			chomp;
			my $temp_line = $_;
			$temp_line =~ s/\r//;
			unless(exists $temp_list{$temp_line}){
				$temp_list{$temp_line} = 1;
			}
		}
		close $temp_fh;
	}

	return \%temp_list;
}

sub read_gtf{
	my ($temp_gtf, $temp_attribute_ref, $temp_tag_ref) = @_;
	open(my $temp_fh, $temp_gtf) or die "<ERROR> Cannot open $temp_gtf!\n";
	my %temp_transcript = ();
	my %temp_gene = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] ne 'gene') && ($temp_line[2] ne 'transcript') && ($temp_line[2] ne 'exon') && ($temp_line[2] ne 'CDS')){
			next;
		}
		my $temp_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"?([^";]+)\"?/){
			$temp_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> gene_id missing at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_gene{$temp_gene_id}){
			%{$temp_gene{$temp_gene_id}} = ();
			$temp_gene{$temp_gene_id}{'chr'} = $temp_line[0];
			$temp_gene{$temp_gene_id}{'strand'} = $temp_line[6];
			$temp_gene{$temp_gene_id}{'left'} = $temp_line[3];
			$temp_gene{$temp_gene_id}{'right'} = $temp_line[4];
			$temp_gene{$temp_gene_id}{'transcript'} = 0;
			$temp_gene{$temp_gene_id}{'Non-NMD'} = 0;
			$temp_gene{$temp_gene_id}{'NMD'} = 0;
			$temp_gene{$temp_gene_id}{'Non-coding'} = 0;
			foreach my $temp_attr (keys %$temp_attribute_ref){
				$temp_gene{$temp_gene_id}{$temp_attr} = '';
			}
			foreach my $temp_tag (keys %$temp_tag_ref){
				$temp_gene{$temp_gene_id}{$temp_tag} = 0;
			}
		}
		foreach my $temp_attr (keys %$temp_attribute_ref){
			if($temp_line[8] =~ /(^|\s)$temp_attr \"?([^";]+)\"?/){
				$temp_gene{$temp_gene_id}{$temp_attr} = $2;
			}
		}
		if($temp_gene{$temp_gene_id}{'left'} > $temp_line[3]){
			$temp_gene{$temp_gene_id}{'left'} = $temp_line[3];
		}
		if($temp_gene{$temp_gene_id}{'right'} < $temp_line[4]){
			$temp_gene{$temp_gene_id}{'right'} = $temp_line[4];
		}
		if($temp_line[2] eq 'gene'){
			foreach my $temp_tag (keys %$temp_tag_ref){
				if($temp_line[8] =~ /(^|\s)tag \"?$temp_tag\"?/){
					$temp_gene{$temp_gene_id}{$temp_tag}++;
				}
			}
			next;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"?([^";]+)\"?/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> transcript_id missing at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			$temp_gene{$temp_gene_id}{'transcript'}++;
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
			$temp_transcript{$temp_transcript_id}{'length'} = 0;
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			@{$temp_transcript{$temp_transcript_id}{'cds'}} = ();
			foreach my $temp_attr (keys %$temp_attribute_ref){
				$temp_transcript{$temp_transcript_id}{$temp_attr} = '';
			}
			foreach my $temp_tag (keys %$temp_tag_ref){
				$temp_transcript{$temp_transcript_id}{$temp_tag} = 0;
			}
		}
		foreach my $temp_attr (keys %$temp_attribute_ref){
			if($temp_line[8] =~ /(^|\s)$temp_attr \"?([^";]+)\"?/){
				$temp_transcript{$temp_transcript_id}{$temp_attr} = $2;
				if($temp_gene{$temp_gene_id}{$temp_attr} eq ''){
					$temp_gene{$temp_gene_id}{$temp_attr} = $2;
				}
			}
			elsif($temp_gene{$temp_gene_id}{$temp_attr} ne ''){
				if(($temp_attr eq 'gene_name') || ($temp_attr eq 'gene_type') || ($temp_attr eq 'gene') || ($temp_attr eq 'gene_biotype') || ($temp_attr eq 'level')){
					$temp_transcript{$temp_transcript_id}{$temp_attr} = $temp_gene{$temp_gene_id}{$temp_attr};
				}
			}
		}
		foreach my $temp_tag (keys %$temp_tag_ref){
			if($temp_line[8] =~ /(^|\s)tag \"?$temp_tag\"?/){
				if($temp_transcript{$temp_transcript_id}{$temp_tag} == 0){
					$temp_transcript{$temp_transcript_id}{$temp_tag} = 1;
					$temp_gene{$temp_gene_id}{$temp_tag}++;
				}
			}
		}
		if($temp_transcript{$temp_transcript_id}{'left'} > $temp_line[3]){
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
		}
		if($temp_transcript{$temp_transcript_id}{'right'} < $temp_line[4]){
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
		}
		if($temp_line[2] eq 'exon'){
			$temp_transcript{$temp_transcript_id}{'length'} += ($temp_line[4] - $temp_line[3] + 1);
			my @temp_exon = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
		}
		elsif($temp_line[2] eq 'CDS'){
			my @temp_cds = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'cds'}}, \@temp_cds);
		}
	}
	close $temp_fh;

	return (\%temp_transcript, \%temp_gene);
}

sub read_sj{
	my ($temp_sj, $temp_sj_intron, $temp_sj_gtf) = @_;
	my %temp_sj = ();
	my %temp_strand = ('0' => '.', '1' => '+', '2' => '-', '+' => '+', '-' => '-', '.' => '.');
	open(my $temp_fh, $temp_sj) or die "Cannot open $temp_sj!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		my ($temp_chr, $temp_strand, $temp_left, $temp_right) = ($temp_line[0], $temp_strand{$temp_line[3]}, $temp_line[1], $temp_line[2]);
		if($temp_sj_gtf){
			($temp_strand, $temp_left, $temp_right) = ($temp_strand{$temp_line[6]}, $temp_line[3], $temp_line[4]);
		}
		if($temp_sj_intron){
			$temp_left--;
			$temp_right++;
		}
		unless(exists $temp_sj{$temp_chr}){
			%{$temp_sj{$temp_chr}} = ();
		}
		unless(exists $temp_sj{$temp_chr}{$temp_strand}){
			%{$temp_sj{$temp_chr}{$temp_strand}} = ();
		}
		unless(exists $temp_sj{$temp_chr}{$temp_strand}{$temp_left}){
			%{$temp_sj{$temp_chr}{$temp_strand}{$temp_left}} = ();
		}
		unless(exists $temp_sj{$temp_chr}{$temp_strand}{$temp_left}{$temp_right}){
			$temp_sj{$temp_chr}{$temp_strand}{$temp_left}{$temp_right} = 1;
		}
	}
	close $temp_fh;

	return \%temp_sj;
}

sub print_summary{
	my ($temp_transcript_ref, $temp_gene_ref, $temp_genome, $temp_transcript_out, $temp_gene_out, $temp_attribute_ref, $temp_tag_ref, $temp_stop_codon_flag, $temp_nmd, $temp_window, $temp_sj, $temp_sj_intron, $temp_sj_gtf) = @_;
	my $temp_seq_in;
	my %temp_seq_obj = ();
	if(-e $temp_genome){
		$temp_seq_in = Bio::SeqIO->new(-file => $temp_genome, -format => 'fasta');
		while(my $temp_seq_obj = $temp_seq_in->next_seq){
			my $temp_id = $temp_seq_obj->primary_id;
			if(exists $temp_seq_obj{$temp_id}){
				print STDERR "Duplicated sequence $temp_id in $temp_genome!\n";
				exit;
			}
			else{
				$temp_seq_obj{$temp_id} = $temp_seq_obj;
			}
		}
	}
	my $temp_sj_ref;
	if(-e ${temp_sj}){
		$temp_sj_ref = read_sj($temp_sj, $temp_sj_intron, $temp_sj_gtf);
	}
	open(my $temp_transcript_fh, "> $temp_transcript_out") or die "<ERROR> Cannot open $temp_transcript_out!\n";
	open(my $temp_gene_fh, "> $temp_gene_out") or die "<ERROR> Cannot open $temp_gene_out!\n";
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		my ($temp_gene_id, $temp_chr, $temp_strand, $temp_left, $temp_right, $temp_length) = ($temp_transcript_ref->{$temp_transcript_id}{'gene_id'}, $temp_transcript_ref->{$temp_transcript_id}{'chr'}, $temp_transcript_ref->{$temp_transcript_id}{'strand'}, $temp_transcript_ref->{$temp_transcript_id}{'left'}, $temp_transcript_ref->{$temp_transcript_id}{'right'}, $temp_transcript_ref->{$temp_transcript_id}{'length'});
		my $temp_exon_num = @{$temp_transcript_ref->{$temp_transcript_id}{'exon'}};
		my $temp_sj_support = sj_support($temp_sj_ref->{$temp_chr}{$temp_strand}, $temp_transcript_ref->{$temp_transcript_id}{'exon'});
		my $temp_cds_num = @{$temp_transcript_ref->{$temp_transcript_id}{'cds'}};
		my ($temp_cds_type, $temp_cds_length, $temp_5utr_length, $temp_3utr_length, $temp_start_codon, $temp_stop_codon, $temp_cds_left, $temp_cds_right) = cds_classification($temp_genome, $temp_seq_obj{$temp_chr}, $temp_strand, $temp_transcript_ref->{$temp_transcript_id}{'exon'}, $temp_transcript_ref->{$temp_transcript_id}{'cds'}, $temp_stop_codon_flag, $temp_nmd);
		my ($temp_fraction_downstream_a, $temp_length_downstream_a) = downstream_a($temp_genome, $temp_seq_obj{$temp_chr}, $temp_strand, $temp_left, $temp_right, $temp_window);
		#my $temp_start_time = Benchmark->new;
		my ($temp_gc) = count_gc($temp_genome, $temp_seq_obj{$temp_chr}, $temp_transcript_ref->{$temp_transcript_id}{'exon'});
		$temp_gc /= $temp_length;
		#my $temp_finish_time = Benchmark->new;
		#my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
		#my $temp_date_format = "+%F-%T.%3N";
		#my $temp_date = `date $temp_date_format`;
		#chomp $temp_date;
		#print $log_fh "[$temp_date]\tgc_content done for $temp_transcript_id in " . timestr($temp_time_spent) , "!\n";
		print $temp_transcript_fh "$temp_transcript_id\t$temp_gene_id\t$temp_chr\t$temp_strand\t$temp_left\t$temp_right\t$temp_exon_num\t$temp_sj_support\t$temp_length\t$temp_cds_type\t$temp_cds_num\t$temp_cds_length\t$temp_5utr_length\t$temp_3utr_length\t$temp_start_codon\t$temp_stop_codon\t$temp_cds_left\t$temp_cds_right\t$temp_fraction_downstream_a\t$temp_length_downstream_a\t$temp_gc";
		$temp_gene_ref->{$temp_gene_id}{$temp_cds_type}++;
		foreach my $temp_attr (sort {$a cmp $b} keys %$temp_attribute_ref){
			print $temp_transcript_fh "\t$temp_transcript_ref->{$temp_transcript_id}{$temp_attr}";
		}
		foreach my $temp_tag (sort {$a cmp $b} keys %$temp_tag_ref){
			print $temp_transcript_fh "\t$temp_transcript_ref->{$temp_transcript_id}{$temp_tag}";
		}
		print $temp_transcript_fh "\n";
	}
	close $temp_transcript_fh;
	foreach my $temp_gene_id (keys %$temp_gene_ref){
		print $temp_gene_fh "$temp_gene_id\t$temp_gene_ref->{$temp_gene_id}{'chr'}\t$temp_gene_ref->{$temp_gene_id}{'strand'}\t$temp_gene_ref->{$temp_gene_id}{'left'}\t$temp_gene_ref->{$temp_gene_id}{'right'}\t$temp_gene_ref->{$temp_gene_id}{'transcript'}\t$temp_gene_ref->{$temp_gene_id}{'Non-NMD'}\t$temp_gene_ref->{$temp_gene_id}{'NMD'}\t$temp_gene_ref->{$temp_gene_id}{'Non-coding'}";
		foreach my $temp_attr (sort {$a cmp $b} keys %$temp_attribute_ref){
			if(($temp_attr eq 'gene_name') || ($temp_attr eq 'gene_type') || ($temp_attr eq 'gene') || ($temp_attr eq 'gene_biotype') || ($temp_attr eq 'level')){
				print $temp_gene_fh "\t$temp_gene_ref->{$temp_gene_id}{$temp_attr}";
			}
		}
		foreach my $temp_tag (sort {$a cmp $b} keys %$temp_tag_ref){
			print $temp_gene_fh "\t$temp_gene_ref->{$temp_gene_id}{$temp_tag}";
		}
		print $temp_gene_fh "\n";
	}

	return 1;
}

sub cds_classification{
	my ($temp_genome, $temp_seq_obj_ref, $temp_strand, $temp_exon_ref, $temp_cds_ref, $temp_stop_codon_flag, $temp_nmd) = @_;
	if(($temp_nmd == -1) || (!@$temp_cds_ref)){

		return ("Non-coding", 0, 0, 0, "Non-coding", "Non-coding", 0, 0);
	}
	my $temp_num_stop_codon_base = 0;
	if($temp_stop_codon_flag){
		$temp_num_stop_codon_base = 3;
	}
	else{
		for(my $i = 0; $i < 3; $i++){
			$temp_num_stop_codon_base += extend_cds_by_one_base($temp_exon_ref, $temp_cds_ref, $temp_strand);
		}
	}
	my $temp_cds_length = 0;
	foreach my $temp_exon_ref (@$temp_cds_ref){
		$temp_cds_length += ($temp_exon_ref->[1] - $temp_exon_ref->[0] + 1);
	}
	my ($temp_d3s, $temp_terminal_exon_flag, $temp_5utr_length, $temp_3utr_length, $temp_start_codon, $temp_stop_codon, $temp_cds_left, $temp_cds_right) = (0, 1, 0, 0, "Non-coding", "Non-coding", 0, 0);
	EXON:for(my $i = 0; $i < @$temp_exon_ref; $i++){
		if($temp_strand eq '-'){
			if(($temp_cds_ref->[0][0] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[0][0] <= $temp_exon_ref->[$i][1])){
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_cds_ref->[0][0] - $temp_exon_ref->[$i][0]);
				}
				$temp_3utr_length += ($temp_cds_ref->[0][0] - $temp_exon_ref->[$i][0]);
				last EXON;
			}
			else{
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
				}
				$temp_3utr_length += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
			}
		}
		else{
			if(($temp_cds_ref->[0][0] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[0][0] <= $temp_exon_ref->[$i][1])){
				$temp_5utr_length += ($temp_cds_ref->[0][0] - $temp_exon_ref->[$i][0]);
				last EXON;
			}
			else{
				$temp_5utr_length += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
			}
		}
		if($temp_terminal_exon_flag == 1){
			$temp_terminal_exon_flag = 0;
		}
	}
	$temp_terminal_exon_flag = 1;
	EXON:for(my $i = (@$temp_exon_ref - 1); $i >= 0; $i--){
		if($temp_strand eq '-'){
			if(($temp_cds_ref->[@$temp_cds_ref-1][1] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[@$temp_cds_ref-1][1] <= $temp_exon_ref->[$i][1])){
				$temp_5utr_length += ($temp_exon_ref->[$i][1] - $temp_cds_ref->[@$temp_cds_ref-1][1]);
				last EXON;
			}
			else{
				$temp_5utr_length += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
			}
		}
		else{
			if(($temp_cds_ref->[@$temp_cds_ref-1][1] >= $temp_exon_ref->[$i][0]) && ($temp_cds_ref->[@$temp_cds_ref-1][1] <= $temp_exon_ref->[$i][1])){
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_exon_ref->[$i][1] - $temp_cds_ref->[@$temp_cds_ref-1][1]);
				}
				$temp_3utr_length += ($temp_exon_ref->[$i][1] - $temp_cds_ref->[@$temp_cds_ref-1][1]);
				last EXON;
			}
			else{
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
				}
				$temp_3utr_length += ($temp_exon_ref->[$i][1] - $temp_exon_ref->[$i][0] + 1);
			}
		}
		if($temp_terminal_exon_flag == 1){
			$temp_terminal_exon_flag = 0;
		}
	}
	if((-e $temp_genome) && (defined $temp_seq_obj_ref)){
		($temp_start_codon, $temp_stop_codon) = start_stop_codon($temp_seq_obj_ref, $temp_cds_ref, $temp_strand, $temp_num_stop_codon_base);
		if($temp_num_stop_codon_base == 0){
			$temp_stop_codon = 'Missing';
		}
	}
	if($temp_d3s > $temp_nmd){

		return ('NMD', $temp_cds_length, $temp_5utr_length, $temp_3utr_length, $temp_start_codon, $temp_stop_codon, $temp_cds_ref->[0][0], $temp_cds_ref->[@$temp_cds_ref-1][1]);
	}
	else{

		return ('Non-NMD', $temp_cds_length, $temp_5utr_length, $temp_3utr_length, $temp_start_codon, $temp_stop_codon, $temp_cds_ref->[0][0], $temp_cds_ref->[@$temp_cds_ref-1][1]);
	}
}

sub extend_cds_by_one_base{
	my ($temp_exon_ref, $temp_cds_ref, $temp_strand) = @_;
	if($temp_strand eq '-'){
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
						unshift(@$temp_cds_ref, \@temp_exon);

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
						push(@$temp_cds_ref, \@temp_exon);

						return 1;
					}
				}
			}
		}
	}

	return 0;
}

sub start_stop_codon{
	my ($temp_seq_obj_ref, $temp_cds_ref, $temp_strand, $temp_num_stop_codon_base) = @_;
	my ($temp_codon, $temp_left_codon, $temp_right_codon) = ('', '', '');
	EXON:for(my $i = 0; $i < @$temp_cds_ref; $i++){
		for(my $j = $temp_cds_ref->[$i][0]; $j <= $temp_cds_ref->[$i][1]; $j++){
			$temp_left_codon .= $temp_seq_obj_ref->subseq($j, $j);
			if(length($temp_left_codon) == 3){
				last EXON;
			}
		}
	}
	$temp_codon = '';
	EXON:for(my $i = (@$temp_cds_ref - 1); $i >= 0; $i--){
		for(my $j = $temp_cds_ref->[$i][1]; $j >= $temp_cds_ref->[$i][0]; $j--){
			$temp_codon .= $temp_seq_obj_ref->subseq($j, $j);
			if(length($temp_codon) == 3){
				last EXON;
			}
		}
	}
	$temp_right_codon = reverse($temp_codon);
	if($temp_strand eq '-'){

		return (reverse_complement($temp_right_codon), substr(reverse_complement($temp_left_codon), (3 - $temp_num_stop_codon_base)));
	}
	else{

		return ($temp_left_codon, substr($temp_right_codon, (3 - $temp_num_stop_codon_base)));
	}	
}

sub reverse_complement{
        my ($temp_seq) = @_;
        my $temp_revcomp = reverse($temp_seq);
        $temp_revcomp =~ tr/ACGTacgt/TGCAtgca/;

	return $temp_revcomp;
}

sub downstream_a{
	my ($temp_genome, $temp_seq_obj_ref, $temp_strand, $temp_left, $temp_right, $temp_window) = @_;
	my ($temp_fraction_a, $temp_length_a) = (0, 0);
	if((-e $temp_genome) && (defined $temp_seq_obj_ref)){
		if($temp_strand eq '-'){
			$temp_right = $temp_left - 1;
			$temp_left -= $temp_window;
			if(($temp_right >= 1) && ($temp_left >= 1)){
				my $temp_seq = $temp_seq_obj_ref->subseq($temp_left, $temp_right);
				($temp_fraction_a, $temp_length_a) = count_a_t($temp_seq, 'T');
			}
		}
		else{
			$temp_left = $temp_right + 1;
			$temp_right += $temp_window;
			if(($temp_left <= $temp_seq_obj_ref->length()) && ($temp_right <= $temp_seq_obj_ref->length())){
				my $temp_seq = $temp_seq_obj_ref->subseq($temp_left, $temp_right);
				($temp_fraction_a, $temp_length_a) = count_a_t($temp_seq, 'A');
			}
		}
		$temp_fraction_a /= $temp_window;
	}

	return ($temp_fraction_a, $temp_length_a);
}

sub count_a_t{
	my ($temp_seq, $temp_a_t) = @_;
	my ($temp_count, $temp_cur_count, $temp_longest_count) = (0, 0, 0);
	my @temp_seq = split('', $temp_seq);
	foreach my $temp_nt (@temp_seq){
		if($temp_nt eq $temp_a_t){
			$temp_count++;
			$temp_cur_count++;
		}
		else{
			if($temp_longest_count < $temp_cur_count){
				$temp_longest_count = $temp_cur_count;
			}
			$temp_cur_count = 0;
		}
	}
	if($temp_longest_count < $temp_cur_count){
		$temp_longest_count = $temp_cur_count;
	}

	return ($temp_count, $temp_longest_count);
}

sub count_gc{
	my ($temp_genome, $temp_seq_obj_ref, $temp_exon_ref) = @_;
	my $temp_gc = 0;
	if((-e $temp_genome) && (defined $temp_seq_obj_ref)){
		for(my $i = 0; $i < @$temp_exon_ref; $i++){
			my $temp_seq = $temp_seq_obj_ref->subseq($temp_exon_ref->[$i][0], $temp_exon_ref->[$i][1]);
			my @temp_seq = split('', $temp_seq);
			foreach my $temp_nt (@temp_seq){
				if(($temp_nt eq 'G') || ($temp_nt eq 'C')){
					$temp_gc++;
				}
			}
		}
	}

	return ($temp_gc);
}

sub sj_support{
	my ($temp_sj_ref, $temp_exon_ref) = @_;
	my $temp_sj_support = 0;
	for(my $i = 0; $i < (@$temp_exon_ref - 1); $i++){
		my $temp_left = $temp_exon_ref->[$i][1];
		my $temp_right = $temp_exon_ref->[$i+1][0];
		if(exists $temp_sj_ref->{$temp_left}{$temp_right}){
			$temp_sj_support++;
		}
	}

	return $temp_sj_support;
}



