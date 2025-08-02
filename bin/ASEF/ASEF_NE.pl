#!/usr/bin/env perl
=start
Author: ALin
Purpose: This is a script for querying novel exon against a reference. Both query and reference files are in .gtf format. It is part of the Alternative Splicing Event Finder (ASEF) package.
Change log:
	v1.1	2022-04	It was adoptted from ASEF_NIE_v1.2.pl.
	v1.2	2022-05	It is now able to find different types of novel exon, including novel donor/acceptor site, intron retention and otherwise novel. Fixed a bug when parsing isoforms with two exons in the query.
	v1.3	2022-11	It now includes comparison of splicing junctions and novel exons not overlapping with each any reference. The codes have been restructured and simply the classification of novel exon types. It doesn't output gene level information but generate a .gtf file for the novel exons. It enables multithreading.
	v1.4	2023-02	Fixed a bug in parsing the referenc donor/acceptor sites. Removed ranndom seed.
	v1.5	2023-02	Added the column showing the number of supported novel splice junctions by external data. Ouput the log to a .log file.
	v1.6	2023-12 Fixed an issue in parsing gtf attribute.
=cut

use strict;
use Getopt::Long;
use threads;
use Benchmark qw(:hireswallclock);

my $total_start_time = Benchmark->new;

my $version = "v1.6";

my $usage="Usage: perl ASEF_NE_${version}.pl
	-q <String> .gtf of query
	-r <String> .gtf of reference
	-o <String> Output base
	-j <String> The splice junction support file (Default: STAR SJ.out.tab format)
	-i <Boolean> Splice junction as intron
	-f <Boolean> The splice junction support file in gtf format
	-t <Integer> Number of threads
	-b <String> Path for bedtools
	-h <Boolean> Help
";

my ($query, $ref, $out, $sj, $sj_intron, $sj_gtf, $num_thread, $btp, $help) = ("", "", "", "", 0, 0, 1, "", 0);

GetOptions(
	'q=s' =>	\$query,
	'r=s'	=>	\$ref,
	'o=s'	=>	\$out,
	'j=s'	=>	\$sj,
	'i!'	=>	\$sj_intron,
	'f!'	=>	\$sj_gtf,
	't=i'	=>	\$num_thread,
	'b=s'	=>	\$btp,
	'h!'	=>	\$help,
);

unless($query && $ref && $out){
	print "$usage";
	exit;
}

if($sj ne ""){
	unless(-e $sj){
		print STDERR "<ERROR> The splice junction support file, $sj, does not exist!\n$usage";
		exit;
	}
}

if($btp eq ""){
	$btp = "bedtools";
}

unless(`$btp`){
	print "Incorrect path for bedtools!\n$usage";
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
print $log_fh "[$date]\tperl ASEF_NE_${version}.pl -q $query -r $ref -o $out -b $btp -t $num_thread\n";

my $dir = "${out}__ASEF_NE_temp";
if(-d $dir){
	print STDERR "<WARNING> Temporary directory $dir exists! Will overwrite it...\n";
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
system("grep \"	exon	\" $ref > ${dir}/ref.gtf");
system("grep \"	exon	\" $query > ${dir}/query.gtf");
my $chr_ref = split_gtf("${dir}/ref.gtf", "${dir}/query.gtf", $dir);
if($sj ne ""){
	split_sj($sj, $dir);
}
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
		my $thread = threads->create(\&worker, $chr_ref->[$chr_pos], $dir, $log_fh, $date_format, $sj_intron, $sj_gtf, $btp);
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
print $out_fh "transcript_id\tNSJ\tSNSJ\tNOE\tNFD\tNFE\tNLA\tNLE\tNDIE\tNAIE\tKEIE\tNEIE\n";
close $out_fh;
system("cat ${dir}/*_transcript.txt >> ${out}_transcript.txt");
system("cat ${dir}/*_novel_exon.gtf | sort -k 1,1 -k 4,4n > ${out}_exon.gtf");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tGenerated ${out}_transcript.txt and ${out}_exon.gtf!\n";
system("rm -r $dir");
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tCleaned temporary files!\n";

my $total_finish_time = Benchmark->new;
my $total_time_spent = timediff($total_finish_time, $total_start_time);
$date = `date $date_format`;
chomp $date;
print $log_fh "[$date]\tTotal time used: ". timestr($total_time_spent) . "\n";

sub worker{
	my ($temp_chr, $temp_dir, $temp_log_fh, $temp_date_format, $temp_sj_intron, $temp_sj_gtf, $temp_btp) = @_;
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print $temp_log_fh "[$temp_date]\tStarting to process $temp_chr...\n";
	my $temp_start_time = Benchmark->new;
	my $temp_ref_transcript_ref = read_gtf("${temp_dir}/ref_${temp_chr}.gtf");
	my $temp_sj_ref;
	if(-e "${temp_dir}/${temp_chr}.sj"){
		$temp_sj_ref = read_sj("${temp_dir}/${temp_chr}.sj", $temp_sj_intron, $temp_sj_gtf);
	}
	my ($temp_ref_exon_ref, $temp_ref_junction_ref, $temp_ref_donor_ref, $temp_ref_acceptor_ref) = parse_ref($temp_ref_transcript_ref);
	print_ref_exon($temp_ref_exon_ref, $temp_chr, "${temp_dir}/ref_${temp_chr}_exon.gtf");
	my %temp_novel = ();
	my $temp_query_transcript_ref = read_gtf("${temp_dir}/query_${temp_chr}.gtf", \%temp_novel);
	my $temp_novel_exon_ref = find_novel_exon(\%temp_novel, "${temp_dir}/query_${temp_chr}.gtf", "${temp_dir}/ref_${temp_chr}_exon.gtf", $temp_btp);
	parse_query($temp_query_transcript_ref, \%temp_novel, $temp_novel_exon_ref, $temp_ref_exon_ref, $temp_ref_junction_ref, $temp_ref_donor_ref, $temp_ref_acceptor_ref, $temp_sj_ref, $temp_chr, "${temp_dir}/${temp_chr}_novel_exon.gtf");
	print_result(\%temp_novel, "${temp_dir}/${temp_chr}_transcript.txt");
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

sub split_sj{
	my ($temp_sj, $temp_dir) = @_;
	open(my $temp_fh, $temp_sj) or die "Cannot open $temp_sj!\n";
	system("for chr in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$chr	\" $temp_sj > ${temp_dir}\/\$\{chr\}.sj; done");

	return 1;
}

sub read_gtf{
	my ($temp_gtf, $temp_query_ref) = @_;
	my %temp_transcript = ();
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "<ERROR> Cannot open $temp_gtf!\n'";
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
			print STDERR "<ERROR> Cannot find gene ID at line:\n$temp_line\n";
			exit;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\;]+)\"/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			if(defined $temp_query_ref){
				unless(exists $temp_query_ref->{$temp_transcript_id}){
					%{$temp_query_ref->{$temp_transcript_id}} = ();
					$temp_query_ref->{$temp_transcript_id}->{'NSJ'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'SNSJ'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NOE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NFD'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NLA'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NFE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NLE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NDIE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NAIE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'KEIE'} = 0;
					$temp_query_ref->{$temp_transcript_id}->{'NEIE'} = 0;
				}
			}
		}
		my @temp_interval = ($temp_line[3], $temp_line[4]);
		push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_interval);
	}
	close $temp_fh;

	return \%temp_transcript;
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
		my ($temp_strand, $temp_left, $temp_right) = ($temp_strand{$temp_line[3]}, $temp_line[1], $temp_line[2]);
		if($temp_sj_gtf){
			($temp_strand, $temp_left, $temp_right) = ($temp_strand{$temp_line[6]}, $temp_line[3], $temp_line[4]);
		}
		if($temp_sj_intron){
			$temp_left--;
			$temp_right++;
		}
		unless(exists $temp_sj{$temp_strand}){
			%{$temp_sj{$temp_strand}} = ();
		}
		unless(exists $temp_sj{$temp_strand}{$temp_left}){
			%{$temp_sj{$temp_strand}{$temp_left}} = ();
		}
		unless(exists $temp_sj{$temp_strand}{$temp_left}{$temp_right}){
			$temp_sj{$temp_strand}{$temp_left}{$temp_right} = 1;
		}
	}
	close $temp_fh;

	return \%temp_sj;
}

sub parse_ref{
	my ($temp_transcript_ref) = @_;
	my %temp_exon = ();
	my %temp_junction = ();
	my %temp_donor = ();
	my %temp_acceptor = ();
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		my $temp_strand = $temp_transcript_ref->{$temp_transcript_id}->{'strand'};
		unless(exists $temp_exon{$temp_strand}){
			%{$temp_exon{$temp_strand}} = ();
			%{$temp_exon{$temp_strand}{'mono'}} = ();
			%{$temp_exon{$temp_strand}{'internal'}} = ();
			%{$temp_exon{$temp_strand}{'external'}} = ();
		}
		unless(exists $temp_junction{$temp_strand}){
			%{$temp_junction{$temp_strand}} = ();
		}
		unless(exists $temp_donor{$temp_strand}){
			%{$temp_donor{$temp_strand}} = ();
			%{$temp_donor{$temp_strand}{'all'}} = ();
			%{$temp_donor{$temp_strand}{'first'}} = ();
		}
		unless(exists $temp_acceptor{$temp_strand}){
			%{$temp_acceptor{$temp_strand}} = ();
			%{$temp_acceptor{$temp_strand}{'all'}} = ();
			%{$temp_acceptor{$temp_strand}{'last'}} = ();
		}
		my @temp_exon = @{$temp_transcript_ref->{$temp_transcript_id}->{'exon'}};
		if(@temp_exon > 1){
			for(my $i = 0; $i < (@temp_exon - 1); $i++){
				if($i > 0){
					unless(exists $temp_exon{$temp_strand}{'internal'}{$temp_exon[$i][0]}){
						%{$temp_exon{$temp_strand}{'internal'}{$temp_exon[$i][0]}} = ();
					}
					unless(exists $temp_exon{$temp_strand}{'internal'}{$temp_exon[$i][0]}{$temp_exon[$i][1]}){
						$temp_exon{$temp_strand}{'internal'}{$temp_exon[$i][0]}{$temp_exon[$i][1]} = 1;
					}
				}
				else{
					unless(exists $temp_exon{$temp_strand}{'external'}{$temp_exon[$i][0]}){
						%{$temp_exon{$temp_strand}{'external'}{$temp_exon[$i][0]}} = ();
					}
					unless(exists $temp_exon{$temp_strand}{'external'}{$temp_exon[$i][0]}{$temp_exon[$i][1]}){
						$temp_exon{$temp_strand}{'external'}{$temp_exon[$i][0]}{$temp_exon[$i][1]} = 1;
					}
				}
				unless(exists $temp_junction{$temp_strand}{$temp_exon[$i][1]}){
					%{$temp_junction{$temp_strand}{$temp_exon[$i][1]}} = ();
				}
				unless(exists $temp_junction{$temp_strand}{$temp_exon[$i][1]}{$temp_exon[$i+1][0]}){
					$temp_junction{$temp_strand}{$temp_exon[$i][1]}{$temp_exon[$i+1][0]} = 1;
				}
				if($temp_strand eq "-"){
					unless(exists $temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i][1]}){
						$temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i][1]} = 1;
					}
					if($i > 0){
						unless(exists $temp_donor{$temp_strand}{'all'}{$temp_exon[$i][0]}){
							$temp_donor{$temp_strand}{'all'}{$temp_exon[$i][0]} = 0;
						}
					}
					if($i == 0){
						unless(exists $temp_acceptor{$temp_strand}{'last'}{$temp_exon[$i][1]}){
							$temp_acceptor{$temp_strand}{'last'}{$temp_exon[$i][1]} = 1;
						}
					}
					if($i == (@temp_exon - 2)){
						unless(exists $temp_donor{$temp_strand}{'all'}{$temp_exon[$i+1][0]}){
							$temp_donor{$temp_strand}{'all'}{$temp_exon[$i+1][0]} = 1;
						}
						unless(exists $temp_donor{$temp_strand}{'first'}{$temp_exon[$i+1][0]}){
							$temp_donor{$temp_strand}{'first'}{$temp_exon[$i+1][0]} = 1;
						}
					}
				}
				else{
					unless(exists $temp_donor{$temp_strand}{'all'}{$temp_exon[$i][1]}){
						$temp_donor{$temp_strand}{'all'}{$temp_exon[$i][1]} = 1;
					}
					if($i > 0){
						unless(exists $temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i][0]}){
							$temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i][0]} = 1;
						}
					}
					if($i == 0){
						unless(exists $temp_donor{$temp_strand}{'first'}{$temp_exon[$i][1]}){
							$temp_donor{$temp_strand}{'first'}{$temp_exon[$i][1]} = 1;
						}
					}
					if($i == (@temp_exon - 2)){
						unless(exists $temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i+1][0]}){
							$temp_acceptor{$temp_strand}{'all'}{$temp_exon[$i+1][0]} = 1;
						}
						unless(exists $temp_acceptor{$temp_strand}{'last'}{$temp_exon[$i+1][0]}){
							$temp_acceptor{$temp_strand}{'last'}{$temp_exon[$i+1][0]} = 1;
						}
					}
				}
			}
			unless(exists $temp_exon{$temp_strand}{'external'}{$temp_exon[@temp_exon-1][0]}){
				%{$temp_exon{$temp_strand}{'external'}{$temp_exon[@temp_exon-1][0]}} = ();
			}
			unless(exists $temp_exon{$temp_strand}{'external'}{$temp_exon[@temp_exon-1][0]}{$temp_exon[@temp_exon-1][1]}){
				$temp_exon{$temp_strand}{'external'}{$temp_exon[@temp_exon-1][0]}{$temp_exon[@temp_exon-1][1]} = 1;
			}
		}
		else{
			unless(exists $temp_exon{$temp_strand}{'mono'}{$temp_exon[@temp_exon-1][0]}){
				%{$temp_exon{$temp_strand}{'mono'}{$temp_exon[@temp_exon-1][0]}} = ();
			}
			unless(exists $temp_exon{$temp_strand}{'mono'}{$temp_exon[@temp_exon-1][0]}{$temp_exon[@temp_exon-1][1]}){
				$temp_exon{$temp_strand}{'mono'}{$temp_exon[@temp_exon-1][0]}{$temp_exon[@temp_exon-1][1]} = 1;
			}
		}
	}

	return (\%temp_exon, \%temp_junction, \%temp_donor, \%temp_acceptor);
}

sub print_ref_exon{
	my ($temp_exon_ref, $temp_chr, $temp_out) = @_;
	open(my $temp_fh, "| sort -k 1,1 -k 4,4n > $temp_out") or die "<ERROR> Cannot create $temp_out!\n";
	my $temp_count = 1;
	foreach my $temp_strand (keys %$temp_exon_ref){
		foreach my $temp_type (keys %{$temp_exon_ref->{$temp_strand}}){
			foreach my $temp_left (keys %{$temp_exon_ref->{$temp_strand}->{$temp_type}}){
				foreach my $temp_right (keys %{$temp_exon_ref->{$temp_strand}->{$temp_type}->{$temp_left}}){
					print $temp_fh "$temp_chr\tASEF\texon\t$temp_left\t$temp_right\t.\t$temp_strand\t.\texon_id \"exon_${temp_count}\"\;\n";
					$temp_count++;
				}
			}
		}
	}
	close $temp_fh;
}

sub find_novel_exon{
	my ($temp_novel_ref, $temp_query, $temp_ref, $temp_btp) = @_;
	my %temp_novel_exon = ();
	open(my $temp_fh, "$temp_btp intersect -wa -v -s -nonamecheck -a $temp_query -b $temp_ref |") or die "<ERROR> Cannot intersect $temp_query and $temp_ref!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		my @temp_line = split('\t', $temp_line);
		my $temp_transcript_id;
		if($temp_line[8] =~ /transcript_id \"([^\"]+)\"/){
			$temp_transcript_id = $1;
		}
		else{
			print STDERR "<ERROR> Cannot find transcript ID at line:\n$temp_line\n";
			exit;
		}
		if(exists $temp_novel_ref->{$temp_transcript_id}){
			$temp_novel_ref->{$temp_transcript_id}->{'NOE'}++;
		}
		else{
			print STEDRR "<ERROR> Transcript $temp_transcript_id was not properly initiated!\n";
			exit;
		}
		unless(exists $temp_novel_exon{$temp_transcript_id}){
			%{$temp_novel_exon{$temp_transcript_id}} = ();
		}
		unless(exists $temp_novel_exon{$temp_transcript_id}{$temp_line[3]}){
			%{$temp_novel_exon{$temp_transcript_id}{$temp_line[3]}} = ();
		}
		unless(exists $temp_novel_exon{$temp_transcript_id}{$temp_line[3]}{$temp_line[4]}){
			$temp_novel_exon{$temp_transcript_id}{$temp_line[3]}{$temp_line[4]} = $temp_line[6];
		}
	}
	close $temp_fh;

	return \%temp_novel_exon;
}

sub parse_query{
	my ($temp_query_transcript_ref, $temp_novel_ref, $temp_novel_exon_ref, $temp_ref_exon_ref, $temp_ref_junction_ref, $temp_ref_donor_ref, $temp_ref_acceptor_ref, $temp_sj_ref, $temp_chr, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "<ERROR> Cannot create $temp_out!\n";
	foreach my $temp_transcript_id (keys %$temp_query_transcript_ref){
		my $temp_gene_id = $temp_query_transcript_ref->{$temp_transcript_id}->{'gene_id'};
		my $temp_strand = $temp_query_transcript_ref->{$temp_transcript_id}->{'strand'};
		my @temp_exon = @{$temp_query_transcript_ref->{$temp_transcript_id}->{'exon'}};
		if(@temp_exon > 1){
			for(my $i = 0; $i < (@temp_exon - 1); $i++){
				my %temp_novel_tag = ();
				unless(exists $temp_ref_junction_ref->{$temp_strand}->{$temp_exon[$i][1]}{$temp_exon[$i+1][0]}){
					$temp_novel_ref->{$temp_transcript_id}->{'NSJ'}++;
					if(exists $temp_sj_ref->{$temp_strand}->{$temp_exon[$i][1]}{$temp_exon[$i+1][0]}){
						$temp_novel_ref->{$temp_transcript_id}->{'SNSJ'}++;
					}
				}
				if($i == 0){
					if($temp_strand eq "-"){
						unless(exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]}){
							$temp_novel_ref->{$temp_transcript_id}->{'NLA'}++;
							unless(exists $temp_novel_tag{'NLA'}){
								$temp_novel_tag{'NLA'} = 1;
							}
						}
						unless(exists $temp_ref_acceptor_ref->{$temp_strand}->{'last'}->{$temp_exon[$i][1]}){
							$temp_novel_ref->{$temp_transcript_id}->{'NLE'}++;
							unless(exists $temp_novel_tag{'NLE'}){
								$temp_novel_tag{'NLE'} = 1;
							}
						}
					}
					else{
						unless(exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]}){
							$temp_novel_ref->{$temp_transcript_id}->{'NFD'}++;
							unless(exists $temp_novel_tag{'NFD'}){
								$temp_novel_tag{'NFD'} = 1;
							}
						}
						unless(exists $temp_ref_donor_ref->{$temp_strand}->{'first'}->{$temp_exon[$i][1]}){
							$temp_novel_ref->{$temp_transcript_id}->{'NFE'}++;
							unless(exists $temp_novel_tag{'NFE'}){
								$temp_novel_tag{'NFE'} = 1;
							}
						}
					}
				}
				else{
					unless(exists $temp_ref_exon_ref->{$temp_strand}->{'internal'}->{$temp_exon[$i][0]}->{$temp_exon[$i][1]}){
						if($temp_strand eq "-"){
							if((exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][0]}) && (exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]})){
								$temp_novel_ref->{$temp_transcript_id}->{'KEIE'}++;
								unless(exists $temp_novel_tag{'KEIE'}){
									$temp_novel_tag{'KEIE'} = 1;
								}
							}
							elsif(exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][0]}){
								$temp_novel_ref->{$temp_transcript_id}->{'NAIE'}++;
								unless(exists $temp_novel_tag{'NAIE'}){
									$temp_novel_tag{'NAIE'} = 1;
								}
							}
							elsif(exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]}){
								$temp_novel_ref->{$temp_transcript_id}->{'NDIE'}++;
								unless(exists $temp_novel_tag{'NDIE'}){
									$temp_novel_tag{'NDIE'} = 1;
								}
							}
							else{
								$temp_novel_ref->{$temp_transcript_id}->{'NEIE'}++;
								unless(exists $temp_novel_tag{'NEIE'}){
									$temp_novel_tag{'NEIE'} = 1;
								}
							}
						}
						else{
							if((exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][0]}) && (exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]})){
								$temp_novel_ref->{$temp_transcript_id}->{'KEIE'}++;
								unless(exists $temp_novel_tag{'KEIE'}){
									$temp_novel_tag{'KEIE'} = 1;
								}
							}
							elsif(exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][0]}){
								$temp_novel_ref->{$temp_transcript_id}->{'NDIE'}++;
								unless(exists $temp_novel_tag{'NDIE'}){
									$temp_novel_tag{'NDIE'} = 1;
								}
							}
							elsif(exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[$i][1]}){
								$temp_novel_ref->{$temp_transcript_id}->{'NAIE'}++;
								unless(exists $temp_novel_tag{'NAIE'}){
									$temp_novel_tag{'NAIE'} = 1;
								}
							}
							else{
								$temp_novel_ref->{$temp_transcript_id}->{'NEIE'}++;
								unless(exists $temp_novel_tag{'NEIE'}){
									$temp_novel_tag{'NEIE'} = 1;
								}
							}
						}
					}
				}
				if(exists $temp_novel_exon_ref->{$temp_transcript_id}->{$temp_exon[$i][0]}->{$temp_exon[$i][1]}){
					unless(exists $temp_novel_tag{'NOE'}){
						$temp_novel_tag{'NOE'} = 1;
					}
				}
				my $temp_novel_tag = keys %temp_novel_tag;
				if($temp_novel_tag > 0){
					print $temp_fh "$temp_chr\tASEF\texon\t$temp_exon[$i][0]\t$temp_exon[$i][1]\t.\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_gene_id\"\;";
					foreach my $temp_tag ("NOE", "NFD", "NFE", "NLA", "NLE", "NDIE", "NAIE", "KEIE", "NEIE"){
						if(exists $temp_novel_tag{$temp_tag}){
							print $temp_fh " tag \"$temp_tag\"\;";
						}
					}
					print $temp_fh "\n";
				}
			}
			my %temp_novel_tag = ();
			if($temp_strand eq "-"){
				unless(exists $temp_ref_donor_ref->{$temp_strand}->{'all'}->{$temp_exon[@temp_exon-1][0]}){
					$temp_novel_ref->{$temp_transcript_id}->{'NFD'}++;
					unless(exists $temp_novel_tag{'NFD'}){
						$temp_novel_tag{'NFD'} = 1;
					}
				}
				unless(exists $temp_ref_donor_ref->{$temp_strand}->{'first'}->{$temp_exon[@temp_exon-1][0]}){
					$temp_novel_ref->{$temp_transcript_id}->{'NFE'}++;
					unless(exists $temp_novel_tag{'NFE'}){
						$temp_novel_tag{'NFE'} = 1;
					}
				}
			}
			else{
				unless(exists $temp_ref_acceptor_ref->{$temp_strand}->{'all'}->{$temp_exon[@temp_exon-1][0]}){
					$temp_novel_ref->{$temp_transcript_id}->{'NLA'}++;
					unless(exists $temp_novel_tag{'NLA'}){
						$temp_novel_tag{'NLA'} = 1;
					}
				}
				unless(exists $temp_ref_acceptor_ref->{$temp_strand}->{'last'}->{$temp_exon[@temp_exon-1][0]}){
					$temp_novel_ref->{$temp_transcript_id}->{'NLE'}++;
					unless(exists $temp_novel_tag{'NLE'}){
						$temp_novel_tag{'NLE'} = 1;
					}
				}
			}
			if(exists $temp_novel_exon_ref->{$temp_transcript_id}->{$temp_exon[@temp_exon-1][0]}->{$temp_exon[@temp_exon-1][1]}){
				unless(exists $temp_novel_tag{'NOE'}){
					$temp_novel_tag{'NOE'} = 1;
				}
			}
			my $temp_novel_tag = keys %temp_novel_tag;
			if($temp_novel_tag > 0){
				print $temp_fh "$temp_chr\tASEF\texon\t$temp_exon[@temp_exon-1][0]\t$temp_exon[@temp_exon-1][1]\t.\t$temp_strand\t.\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_query_transcript_ref->{$temp_transcript_id}->{'gene_id'}\"\;";
				foreach my $temp_tag ("NOE", "NFD", "NFE", "NLA", "NLE", "NDIE", "NAIE", "KEIE", "NEIE"){
					if(exists $temp_novel_tag{$temp_tag}){
						print $temp_fh " tag \"$temp_tag\"\;";
					}
				}
				print $temp_fh "\n"
			};
		}
	}

	return 1;
}

sub print_result{
	my ($temp_novel_ref, $temp_out) = @_;
	open(my $temp_fh, "> $temp_out") or die "<ERROR> Cannot create $temp_out!\n";
	foreach my $temp_transcript_id (keys %$temp_novel_ref){
		print $temp_fh "$temp_transcript_id";
		for my $temp_col ("NSJ", "SNSJ", "NOE", "NFD", "NFE", "NLA", "NLE", "NDIE", "NAIE", "KEIE", "NEIE"){
			print $temp_fh "\t$temp_novel_ref->{$temp_transcript_id}->{$temp_col}";
		}
		print $temp_fh "\n";
	}
	close $temp_fh;

	return 1;
}



