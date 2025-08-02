#!/usr/bin/env perl
#Author: ALin
#Purpose: To find splice junctions present in NMD, non-NMD and non-coding.
#Change log:
#	v1.1	2022-12	This script was adopted from gtf2transcritp_summary_v1.6.pl. A bug in cds_classification was fixed.

use strict;
use Getopt::Long;
use threads;
use Benchmark qw(:hireswallclock);

my $total_start_time = Benchmark->new;

my $version = "v1.1";
my $random = `tr -dc A-Za-z0-9 </dev/urandom | head -c 10`;

my $usage = "Usage: perl nmd_spilce_junction_finder_${version}.pl
	-i <String> Input gtf file
	-o <String> Output base
	-s <Boolean> The stop codoon is included in the CDS feature.
	-m <Integer> The minimal CDS length (Including the stop codon, default: 24)
	-n <Integer> The distance between the stop codon (exclusive) and the downstream 3' most splice site
	-t <Integer> Number of threads (Default: 1)
	-h <Boolean> Help
";


my ($gtf, $out, $stop_codon, $min, $nmd, $num_thread, $help) = ("", "", 0, 24, -1, 1, 0);

GetOptions(
	'i=s'	=>	\$gtf,
	'o=s'	=>	\$out,
	's!'	=>	\$stop_codon,
   'm=i'	=>	\$min,
	'n=i'	=>	\$nmd,
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

my $date_format = "+%F-%T.%3N";
my $date = `date $date_format`;
chomp $date;
print "[$date]\tperl nmd_splice_junction_finder_${version}.pl -i $gtf -o $out";

if($stop_codon){
	print " -s";
}

if($nmd >= 0){
	print " -n $nmd";
}

print " -t $num_thread\n";

my $dir = "temp_${random}";
if(-d $dir){
	print STDERR "<ERROR> Temporary directory $dir exists! Please remove this directory and re-run!\n";
	exit;
}
elsif(-e $dir){
	print STDERR "<ERROR> Exists file $dir having the same name as the temporary directory! Please remove this file and re-run!\n";
	exit;
}
else{
	system("mkdir $dir");
}

my $total_thread = `grep -c \"^processor\" /proc/cpuinfo`;
chomp $total_thread;
$date = `date $date_format`;
chomp $date;
print "[$date]\tA total of $total_thread threads detected\n";

if($num_thread > $total_thread){
	print "<WARNING> Exceeding the maximal total number of threads! Setting the number of thread to $total_thread...";
	$num_thread = $total_thread;
}

my $start_time = Benchmark->new;
$date = `date $date_format`;
chomp $date;
print "[$date]\tStarting preprocessing...\n";
my $chr_ref = split_gtf($gtf, $random, $dir);
my $finish_time = Benchmark->new;
my $time_spent = timediff($finish_time, $start_time);
$date = `date $date_format`;
chomp $date;
print "[$date]\tUsed ". timestr($time_spent) . " for preprocessing\n";

my $chr_pos = 0;
my @thread = ();
my @running = ();

while(@thread < @$chr_ref){
	@running = threads->list(threads::running);
	if(@running < $num_thread){
		my $thread = threads->create(\&worker, $chr_ref->[$chr_pos], $random, $dir, $date_format, $stop_codon, $min, $nmd);
		push(@thread, $thread);
		$chr_pos++;
	}
}

foreach my $thread (@thread){
	my @result = $thread->join();
	$date = `date $date_format`;
	chomp $date;
	print "[$date]\tJoined for $result[0]!\n";
}

my $out_fh;
open($out_fh, "> ${out}_nmd_sj.txt") or die "<ERROR> Cannot create ${out}_nmd_sj.txt!\n";
print $out_fh "sj_id\tNMD\tNon-NMD\tNon-coding\n";
close $out_fh;
system("cat ${dir}/temp_${random}_*_nmd_sj.txt | sort -k 1 >> ${out}_nmd_sj.txt");
open($out_fh, "> ${out}_nmd_sj_per_transcript.txt") or die "<ERROR> Cannot create ${out}_nmd_sj_per_transcript.txt!\n";
print $out_fh "transcript_id\tsj_id\torigin\n";
close $out_fh;
system("cat ${dir}/temp_${random}_*_nmd_sj_per_transcript.txt | sort -k 1 >> ${out}_nmd_sj_per_transcript.txt");
$date = `date $date_format`;
chomp $date;
print "[$date]\tGenerated ${out}_nmd_sj.txt and ${out}_nmd_sj_per_transcript.txt!\n";
system("rm ${dir}/temp_${random}*");
system("rmdir ${dir}");
$date = `date $date_format`;
chomp $date;
print "[$date]\tCleaned temporary files!\n";

my $total_finish_time = Benchmark->new;
my $total_time_spent = timediff($total_finish_time, $total_start_time);
$date = `date $date_format`;
chomp $date;
print "[$date]\tTotal time used: ". timestr($total_time_spent) . "\n";

sub worker{
	my ($temp_chr, $temp_random, $temp_dir, $temp_date_format, $temp_stop_codon, $temp_min, $temp_nmd) = @_;
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[$temp_date]\tStarting to process $temp_chr...\n";
	my $temp_start_time = Benchmark->new;
	my $temp_transcript_ref = read_gtf("${temp_dir}/temp_${temp_random}_${temp_chr}.gtf");
	my ($temp_sj_ref, $temp_sj_per_transcript_ref) = analyze_sj($temp_chr, $temp_transcript_ref, $temp_stop_codon, $temp_min, $temp_nmd);
	print_result($temp_sj_ref, $temp_sj_per_transcript_ref, "${temp_dir}/temp_${temp_random}_${temp_chr}_nmd_sj.txt", "${temp_dir}/temp_${temp_random}_${temp_chr}_nmd_sj_per_transcript.txt");
	my $temp_finish_time = Benchmark->new;
	my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[$temp_date]\tDone for $temp_chr in " . timestr($temp_time_spent) , "!\n";

	return $temp_chr;

}

sub split_gtf{
	my ($temp_gtf, $temp_random, $temp_dir) = @_;
	system("grep -v \"^#\" $temp_gtf | cut -f1 | sort | uniq > ${temp_dir}/temp_${temp_random}_chr_list.txt");
	my @temp_chr = ();
	my $temp_fh;
	open($temp_fh, "${temp_dir}/temp_${temp_random}_chr_list.txt") or die "<ERROR> Cannot open ${temp_dir}/temp_${temp_random}_chr_list.txt!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		push(@temp_chr, $temp_line);
	}
	close $temp_fh;
	system("for chr in `cat ${temp_dir}/temp_${temp_random}_chr_list.txt`; do grep \"^\$chr	\" $temp_gtf > ${temp_dir}/temp_${temp_random}_\$\{chr\}.gtf; done");

	return \@temp_chr;
}

sub read_gtf{
	my ($temp_gtf) = @_;
	my $temp_fh;
	open($temp_fh, $temp_gtf) or die "<ERROR> Cannot open $temp_gtf!\n";
	my %temp_transcript = ();

	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[2] ne "exon") && ($temp_line[2] ne "transcript") && ($temp_line[2] ne "CDS")){
			next;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $1;
		}
		else{
			print STDERR "<ERROR> transcript_id missing at line\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			my $temp_gene_id;
			if($temp_line[8] =~ /gene_id \"([^\"]+)\"/){
				$temp_gene_id = $1;
			}
			else{
				print STDERR "<ERROR> gene_id missing at line\n$temp_line\n";
				exit;
			}
			$temp_transcript{$temp_transcript_id}{'gene'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
			$temp_transcript{$temp_transcript_id}{'length'} = 0;
			$temp_transcript{$temp_transcript_id}{'exon_num'} = 0;
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			@{$temp_transcript{$temp_transcript_id}{'cds'}} = ();
		}
		if($temp_transcript{$temp_transcript_id}{'left'} > $temp_line[3]){
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
		}
		if($temp_transcript{$temp_transcript_id}{'right'} < $temp_line[4]){
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
		}
		if($temp_line[2] eq "exon"){
			$temp_transcript{$temp_transcript_id}{'exon_num'}++;
			$temp_transcript{$temp_transcript_id}{'length'} += ($temp_line[4] - $temp_line[3] + 1);
			my @temp_exon = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
		}
		elsif(($temp_line[2] eq "CDS") || ($temp_line[2] eq "stop_codon")){
			my @temp_cds = ($temp_line[3], $temp_line[4]);
			push(@{$temp_transcript{$temp_transcript_id}{'cds'}}, \@temp_cds);
		}
	}
	close $temp_fh;

	return \%temp_transcript;
}

sub analyze_sj{
	my ($temp_chr, $temp_transcript_ref, $temp_stop_codon, $temp_min, $temp_nmd) = @_;
	my %temp_sj = ();
	my %temp_sj_per_transcript = ();
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		if(@{$temp_transcript_ref->{$temp_transcript_id}->{'exon'}} > 1){
			my ($temp_cds_type, $temp_nmd_sj_left, $temp_nmd_sj_right) = cds_classification($temp_transcript_ref->{$temp_transcript_id}->{'strand'}, $temp_transcript_ref->{$temp_transcript_id}->{'exon'}, $temp_transcript_ref->{$temp_transcript_id}->{'cds'}, $temp_stop_codon, $temp_min, $temp_nmd);
			my @temp_sorted_exon = sort { $a->[0] <=> $b->[0] } @{$temp_transcript_ref->{$temp_transcript_id}->{'exon'}};
			for(my $i = 0; $i < (@temp_sorted_exon - 1); $i++){
				my $temp_sj_id = "${temp_chr}_$temp_sorted_exon[$i][1]_$temp_sorted_exon[$i+1][0]_$temp_transcript_ref->{$temp_transcript_id}->{'strand'}";
				my $temp_sj_pos = "NC";
				if($temp_cds_type ne "Non-coding"){
					$temp_sj_pos = check_sj_pos($temp_transcript_ref->{$temp_transcript_id}->{'strand'}, $temp_transcript_ref->{$temp_transcript_id}->{'cds'}, $temp_sorted_exon[$i][1], $temp_sorted_exon[$i+1][0]);
				}
				unless(exists $temp_sj{$temp_sj_id}){
					%{$temp_sj{$temp_sj_id}} = ();
					foreach my $temp_type ("NMD", "Non-NMD", "Non-coding"){
						%{$temp_sj{$temp_sj_id}{$temp_type}} = ();
					}
				}
				unless(exists $temp_sj{$temp_sj_id}{$temp_cds_type}{$temp_sj_pos}){
					$temp_sj{$temp_sj_id}{$temp_cds_type}{$temp_sj_pos} = 1;
				}
				unless(exists $temp_sj_per_transcript{$temp_transcript_id}){
					%{$temp_sj_per_transcript{$temp_transcript_id}} = ();
				}
				if(($temp_sorted_exon[$i][1] == $temp_nmd_sj_left) && ($temp_sorted_exon[$i+1][0] == $temp_nmd_sj_right)){
					$temp_sj_pos .= "_${temp_cds_type}";
				}
				unless(exists $temp_sj_per_transcript{$temp_transcript_id}{$temp_sj_pos}){
					%{$temp_sj_per_transcript{$temp_transcript_id}{$temp_sj_pos}} = ();
				}
				unless(exists $temp_sj_per_transcript{$temp_transcript_id}{$temp_sj_pos}{$temp_sj_id}){
					$temp_sj_per_transcript{$temp_transcript_id}{$temp_sj_pos}{$temp_sj_id} = 1;
				}
			}
		}
	}

	return (\%temp_sj, \%temp_sj_per_transcript);
}

sub cds_classification{
	my ($temp_strand, $temp_exon_ref, $temp_cds_ref, $temp_stop_codon, $temp_min, $temp_nmd) = @_;
	if(($temp_nmd == -1) || (!@$temp_cds_ref)){

		return ("Non-coding");
	}
	my @temp_sorted_exon = sort { $a->[0] <=> $b->[0] } @$temp_exon_ref;
	my @temp_sorted_cds = sort { $a->[0] <=> $b->[0] } @$temp_cds_ref;
	my $temp_num_stop_codon_base = 0;
	unless($temp_stop_codon){
		for(my $i = 0; $i < 3; $i++){
			$temp_num_stop_codon_base += extend_cds_by_one_base(\@temp_sorted_exon, \@temp_sorted_cds, $temp_strand);
		}
	}
	my $temp_cds_length = 0;
	foreach my $temp_exon_ref (@$temp_cds_ref){
		$temp_cds_length += ($temp_exon_ref->[1] - $temp_exon_ref->[0] + 1);
	}
	if($temp_cds_length < $temp_min){
		return ("Non-coding", 0, 0);
	}
	my ($temp_d3s, $temp_terminal_exon_flag, $temp_cds_terminal_exon) = (0, 1, 0, 0);
	EXON:for(my $i = 0; $i < @temp_sorted_exon; $i++){
		if($temp_strand eq "-"){
			if(($temp_sorted_cds[0][0] >= $temp_sorted_exon[$i][0]) && ($temp_sorted_cds[0][0] <= $temp_sorted_exon[$i][1])){
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_sorted_cds[0][0] - $temp_sorted_exon[$i][0]);
					$temp_cds_terminal_exon = $i;
				}
				last EXON;
			}
			else{
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_sorted_exon[$i][1] - $temp_sorted_exon[$i][0] + 1);
				}
			}
		}
		if($temp_terminal_exon_flag == 1){
			$temp_terminal_exon_flag = 0;
		}
	}
	$temp_terminal_exon_flag = 1;
	EXON:for(my $i = (@temp_sorted_exon - 1); $i >= 0; $i--){
		if($temp_strand ne "-"){
			if(($temp_sorted_cds[@temp_sorted_cds-1][1] >= $temp_sorted_exon[$i][0]) && ($temp_sorted_cds[@temp_sorted_cds-1][1] <= $temp_sorted_exon[$i][1])){
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_sorted_exon[$i][1] - $temp_sorted_cds[@temp_sorted_cds-1][1]);
					$temp_cds_terminal_exon = $i;
				}
				last EXON;
			}
			else{
				if($temp_terminal_exon_flag == 0){
					$temp_d3s += ($temp_sorted_exon[$i][1] - $temp_sorted_exon[$i][0] + 1);
				}
			}
		}
		if($temp_terminal_exon_flag == 1){
			$temp_terminal_exon_flag = 0;
		}
	}
	if($temp_d3s > $temp_nmd){
		my ($temp_d3s_nmd, $temp_nmd_sj_left, $temp_nmd_sj_right) = (0, 0, 0);
		if($temp_strand eq "-"){
			EXON:for(my $i = $temp_cds_terminal_exon; $i >= 0; $i--){
				if($i == $temp_cds_terminal_exon){
					$temp_d3s_nmd += ($temp_sorted_cds[0][0] - $temp_sorted_exon[$i][0]);
				}
				else{
					$temp_d3s_nmd += ($temp_sorted_exon[$i][1] - $temp_sorted_exon[$i][0] + 1);
				}
				if($temp_d3s_nmd > $temp_nmd){
					($temp_nmd_sj_left, $temp_nmd_sj_right) = ($temp_sorted_exon[$i-1][1], $temp_sorted_exon[$i][0]);
					last EXON;
				}	
			}
		}
		else{
			EXON: for(my $i = $temp_cds_terminal_exon; $i < @temp_sorted_exon; $i++){
				if($i == $temp_cds_terminal_exon){
					$temp_d3s_nmd += ($temp_sorted_exon[$i][1] - $temp_sorted_cds[@temp_sorted_cds-1][1]);
				}
				else{
					$temp_d3s_nmd += ($temp_sorted_exon[$i][1] - $temp_sorted_exon[$i][0] + 1);
				}
				if($temp_d3s_nmd > $temp_nmd){
					($temp_nmd_sj_left, $temp_nmd_sj_right) = ($temp_sorted_exon[$i][1], $temp_sorted_exon[$i+1][0]);
					last EXON;
				}
			}
		}

		return ("NMD", $temp_nmd_sj_left, $temp_nmd_sj_right);
	}
	else{

		return ("Non-NMD", 0, 0);
	}
}

sub extend_cds_by_one_base{
	my ($temp_exon_ref, $temp_cds_ref, $temp_strand) = @_;
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

sub check_sj_pos{
	my ($temp_strand, $temp_cds_ref, $temp_left, $temp_right) = @_;
	my @temp_cds = sort { $a->[0] <=> $b->[0] } @$temp_cds_ref;
	my $temp_sj_pos;
	if($temp_strand eq "-"){
		if($temp_right < $temp_cds[0][0]){
			$temp_sj_pos = "3UTR";
		}
		elsif($temp_left > $temp_cds[@temp_cds-1][1]){
			$temp_sj_pos = "5UTR";
		}
		else{
			$temp_sj_pos = "CDS";
		}
	}
	else{
		if($temp_right < $temp_cds[0][0]){
			$temp_sj_pos = "5UTR";
		}
		elsif($temp_left > $temp_cds[@temp_cds-1][1]){
			$temp_sj_pos = "3UTR";
		}
		else{
			$temp_sj_pos = "CDS";
		}
	}

	return $temp_sj_pos;
}

sub print_result{
	my ($temp_sj_ref, $temp_sj_per_transcript_ref, $temp_sj_out, $temp_per_transcript_out) = @_;
	my ($temp_sj_fh, $temp_per_transcript_fh);
	open($temp_sj_fh, "> $temp_sj_out") or die "Cannot create $temp_sj_out!\n";
	open($temp_per_transcript_fh, "> $temp_per_transcript_out") or die "Cannot create $temp_per_transcript_out!\n";
	foreach my $temp_sj_id (keys %$temp_sj_ref){
		print $temp_sj_fh "$temp_sj_id";
		my $temp_nmd_uniq = 0;
		foreach my $temp_type ("NMD", "Non-NMD", "Non-coding"){
			my @temp_sj_pos = sort {$a cmp $b} keys %{$temp_sj_ref->{$temp_sj_id}->{$temp_type}};
			my $temp_sj_pos = join(",", @temp_sj_pos);
			print $temp_sj_fh "\t$temp_sj_pos";
		}
		print $temp_sj_fh "\n";
	}
	close $temp_sj_fh;
	foreach my $temp_transcript_id (keys %$temp_sj_per_transcript_ref){
		foreach my $temp_sj_pos (keys %{$temp_sj_per_transcript_ref->{$temp_transcript_id}}){
			foreach my $temp_sj_id (keys %{$temp_sj_per_transcript_ref->{$temp_transcript_id}->{$temp_sj_pos}}){
				print $temp_per_transcript_fh "$temp_transcript_id\t$temp_sj_id\t$temp_sj_pos\n";
			}
		}
	}
	close $temp_per_transcript_fh;

	return 1;
}



