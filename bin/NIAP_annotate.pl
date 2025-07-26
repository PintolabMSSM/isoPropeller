#!/usr/bin/env perl
=start
Author: ALin
Purpose: To compare a query gtf with a refence gtf and annotate the transcripts in the query gtf.
Change logs:
	v1.1	2021-02	This script is a succession of compare_Cufflinks_gtf_with_reference_gtf_v1.6.pl. It can now bettern handle cases of suspected isoform.
	v1.2	2021-03	This script is now able to make detailed annotation of the unkown transcripts, as antisense, intergenic, intronic and overlapping, partial transcripts of the annotated ones, which may be truncated, novel transcripts, including those with alternative ends and novel structures.
	v1.3	2021-03	It now considers the direction of the cases with novel ends.
	v1.4	2021-04	A bug was fixed when the reference .gtf contained entries other than "exon".
	v1.5	2022-02	A bug was fixed for replacing the gene_id and transcript_id in a hit. Removed the operation for "contained_in".
	v1.6	2022-02	Changed the "suspected" category as a supplementary category. Added the random string in the temporary output file names.
	v1.7	2022-03	Fixed the partial match bug to consider query with not mismatch exon(s). Simplified the ranking for transcripts hitting multiple genes by taking the first ranked genes as the annotation. Added the suspected tag to annotated transcripts if they overlap with multiple genes.
	v1.8	2022-06	Restructured the codes. Provided new steps for annotation. First, identify transcripts from annotated genes and novel genes. Second, recluster novel genes. Third, add novel genes to the reference for refining annotation of transcripts from annotated genes. Fourth, classify unannotated transcripts depnding on whether they are parts of the annotated genes or novel. Remove the stranded parameter and leave the unstranded features as they are.
	v1.9	2022-06	Keep the transcript IDs and added ref_transcript_id, instead of using asm_transcript_id. Enables multithreading by splitting the chromosomes. Fiexed a regex bug in matching the codes for the partial case.
	v1.9.6	2022-12 Fixed a bug in comparing terminal exons of multiexonic transcripts with internal exons of multiexonic ones. Employed a five-pass scheme for annotating transcripts in reference gene loci. Added the gene types of monoexonic ncRNA gene for consideration when classifying novel genes, which requires the gene_type attribute in the reference and mainly work with GENCODE. Fixed a bug in sorting gene/transcript id so that the annotation is stable.
	(1) Identify isoforms not overlapping any reference genes and cluster them into reference novel gene loci.
	(2) Group reference and reference novel gene isoforms overlapping with each other (requiring at least exon boundary for multiexonic or overlapping with at least monoexonic exon) into loci.
	(3) Compare each query monoexonic isoform against each reference monoexonic isoform and find matched exons. Require overlapping and TSSs within the  distance (100 bp by default).
	(4) Compare each query multiexonic isoform against each reference multiexonic isoform and find matched exons. An internal exon vs an internal exon: Require matches in both exon boundaries. An internal, or first, exon vs a first exon: Require a match in the donor site. An internal, or last, exon vs a last exon: Require a match in the acceptor site.
	(5) Annotate query multiexonic isoforms that have at least one matched exon against the reference multiexonic isoforms. known: Require a query isoform to match exons of a reference isoform. partial: Require all exons of a query isoform to match the same amount of consecutive exons, in a reference isoform. If the 5’, or 3’, end of the query isoform is longer than the matched reference exon, it will be tagged “partial_5p_extended”, or “partial_3p_extended”. Concurrence of “partial_5p_extended” and “partial_3p_extended” will be tagged “partial_extended”. novel_structure: Require a query isoform to match at least one exon in a reference isoform. fusion: If the query isoform overlaps exons from different reference genes and different loci from step (2). novel_gene: Require a query isoform to match at least one exon in an isoform of a reference novel gene locus.
	(6) Repeat step (3) - (5) with all remaining query isoforms and the isoforms (annotated, partial, novel_structure, novel_gene) from step (5) until no novel isoforms can be annotated.
	(7) Repeat step (3) - (5) with all remaining query isoforms and the fusion isoforms from step (5) until no novel isoforms can be added to the fusion loci . Isoforms will be annotated as “fusion” at this step.
	(8) Compare the TSS of each remaining query isoforms , including monoexonic and multiexonic, against that of each expanded reference isoforms (annotated, partia, lnovel_structure, novel_gene). novel_structure: Require a query TSS is within 100 bp, by default, to an expanded reference TSS, excluding monoexonic noncoding genes. Annotate use the isoforms having the nearest TSS. novel_gene: Require a query TSS is within 100 bp, by default, to an expanded reference novel gene TSS. Annotate use the isoforms having the nearest TSS.
	(9) Compare all remaining query monoexonic isoforms to reference monoexonic isoforms, excluding monoexonic non-coding gene. known: Require overlapping between the two monoexonic isoforms.
	(10) Group all remaining query isoforms, including monoexonic and multiexonic, overlapping with each other into loci. 
	(11) Classify novel genes based on https://pubmed.ncbi.nlm.nih.gov/32747822/, with reference to protein-coding genes. ncRNA_host_gene > divergent (100 bp by default) > intronic > overlapping > antisense > intergenic (sub-type overlapping)
	v1.10	2023-10	Restructured the code. Added ref_gene_name and used gene_name corresponding to the gene_id. Modified to use NIAP_merge_1.6.pl.
	v1.11	2023-11	Added an option to use TSS region instead of TSS position. Modified the functions to output TSS without adding the promoter width. Revised the step for re-arranging fusion gene ID. Switched the order of the 4th and 5th passes. Added query_match and query_overlap during evaluation, and updated the ranking, query_match > query_overlap > ref_match > ref_overlap > overlap. Fixed an issue to keep gene_name updated as gene_id for novel genes. Fixed an issue of getting gtf attributes.
	v1.12	2024-05	Fixed a bug in tss_dist. Fixed a bug in the 1st-pass when the best reference transcript is from a novel gene but there is no match of exon. Fixed the issue when the reference is empty. Minor modifications. Fixed a bug when a partial isoform is from a novel gene.
	v1.13	2024-09	Fixed the bed file left position issue.
=cut

use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use threads;
use Benchmark qw(:hireswallclock);

my $niap_merge = `dirname $0`;
chomp $niap_merge;
$niap_merge .= "/NIAP_merge.pl";
unless(-e $niap_merge){
	print STDERR "$niap_merge not found. Please install the script.\n";
	exit;
}

my $total_start_time = Benchmark->new;

my $version = "v1.13";

my ($ref, $query, $out, $max_p, $max_bd, $missing, $tss_region, $btp, $num_thread, $help) = ("", "", "", 100, 100, 0, "", "", 1, 0);

my $usage = "Usage: perl NIAP_annotate.pl
	-q <String> Input query gtf file
	-r <String> Input reference gtf file
	-o <String> Output file
	-p <Integer> The maximal length (bp) of a promoter. (Default: $max_p)
	-d <Integer> The maximal distance (bp) defining a bi-directional promoter of a novel gene. (Default: $max_bd)
	-m <Boolean> Add those missing reference exons in the query gtf as in the refence gtf file
	-e <String> Bed file of TSS region
	-b <String> Path for bedtools
	-t <Integer> Number of threads (Default: $num_thread)
	-h <Boolean> Help
";



GetOptions(
	'q=s'	=>	\$query,
	'r=s'	=>	\$ref,
	'o=s'	=>	\$out,
	'p=i'	=>	\$max_p,
	'd=i'	=>	\$max_bd,
	'm!'	=>	\$missing,
	'e=s'	=> 	\$tss_region,
	'b!'	=>	\$btp,
	't=i'	=>	\$num_thread,
	'h!'	=>	\$help,
);

unless($btp){
	$btp = "bedtools";
}

unless(`$btp`){
	print STDERR "<ERROR> Incorrect path for bedtools!\n$usage";
	exit;
}

if($help){
	print "$usage";
	exit;
}

unless(-e $query){
	print STDERR "<ERROR> Query file, $query, does not exist!\n$usage";
	exit;
}

unless(-e $ref){
	print STDERR "<ERROR> Reference file, $ref, does not exist!\n$usage";
	exit;
}

unless($out){
	print STDERR "<ERROR> Output file missing!\n$usage";
	exit;
}

if($max_p < 0){
	print STDERR "<ERROR> The maximal length (bp) of a promoter must be greater than 0!\n$usage";
	exit;
}

if($max_bd < 0){
	print STDERR "<ERROR> The maximal distance (bp) defining a bi-directional promoter of a novel gene must be greater than 0!\n$usage";
	exit;
}

open(my $out_fh, "> $out") or die "Cannot create $out!\n";
my $date_format = "+%F-%T.%3N";
my $date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tperl NIAP_annotate_${version}.pl -q $query -r $ref -o $out -b $btp -t $num_thread -d $max_bd -p $max_p";
print $out_fh "#perl NIAP_annotate_${version}.pl -q $query -r $ref -o $out -b $btp -t $num_thread -d $max_bd -p $max_p";
if($missing){
	print " -m";
	print $out_fh " -m";
}
if($tss_region){
	if(-e $tss_region){
		print " -e $tss_region";
		print $out_fh " -e $tss_region";
	}
	else{
		print STDERR "<ERROR> $tss_region does not exist!\n";
		exit;
	}
}
print "\n";
print $out_fh "\n";
close $out_fh;

my $dir = "${out}__NIAP_annotate_temp";
if(-d $dir){
	print STDERR "<WARNING> Temporary directory $dir exists! Will overwrite it...\n";
	system("rm ${dir}/*");
}
elsif(-e $dir){
	print STDERR "<ERROR> File $dir has the same name as the temporary directory! Please remove this file and re-run!\n";
	exit;
}
else{
	system("mkdir $dir");
}

my $total_thread = `grep -c \"^processor\" /proc/cpuinfo`;
chomp $total_thread;
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tA total of $total_thread threads detected\n";

if($num_thread > $total_thread){
	print "<WARNING> Exceeding the maximal total number of threads! Setting the number of thread to $total_thread...";
	$num_thread = $total_thread;
}

my $start_time = Benchmark->new;
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tStarting preprocessing...\n";
system("grep \$\'\\t\'\"exon\"\$\'\\t\' $ref > ${dir}/ref.gtf");
system("grep \$\'\\t\'\"exon\"\$\'\\t\' $query > ${dir}/query.gtf");
my $chr_ref = split_gtf("${dir}/ref.gtf", "${dir}/query.gtf", $dir);
if($tss_region){
	split_bed($tss_region, $dir);
}
my $finish_time = Benchmark->new;
my $time_spent = timediff($finish_time, $start_time);
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tUsed ". timestr($time_spent) . " for preprocessing\n";

my $chr_pos = 0;
my @thread = ();
my @running = ();

while(@thread < @$chr_ref){
	@running = threads->list(threads::running);
	if(@running < $num_thread){
		my $thread = threads->create(\&worker, $chr_ref->[$chr_pos], $niap_merge, $dir, $date_format, $max_p, $max_bd, $missing, $tss_region, $btp);
		push(@thread, $thread);
		$chr_pos++;
	}
}

foreach my $thread (@thread){
	my @result = $thread->join();
	$date = `date $date_format`;
	chomp $date;
	print "[NIAP_annotate $date]\tJoined for $result[0]!\n";
}

system("cat ${dir}/*_out.gtf | $btp sort >> $out");
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tGenerated $out!\n";
system("rm -r ${dir}");
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tCleaned temporary files!\n";

my $total_finish_time = Benchmark->new;
my $total_time_spent = timediff($total_finish_time, $total_start_time);
$date = `date $date_format`;
chomp $date;
print "[NIAP_annotate $date]\tTotal time used: ". timestr($total_time_spent) . "\n";

sub worker{
	my ($temp_chr, $temp_niap_merge, $temp_dir, $temp_date_format, $temp_max_p, $temp_max_bd, $temp_missing, $temp_tss_region, $temp_btp) = @_;
	my %temp_monoexonic_gene_type = ("miRNA", 1, "misc_RNA", 1, "Mt_rRNA", 1, "Mt_tRNA", 1, "ribozyme", 1, "rRNA", 1, "rRNA_pseudogene", 1, "scaRNA", 1, "scRNA", 1, "snoRNA", 1, "snRNA", 1, "sRNA", 1, "vaultRNA", 1);
	my %temp_pc_gene_type = ("protein_coding", 1);
	my %temp_immune_gene_type = ("IG_C_gene", 1, "IG_C_pseudogene", 1, "IG_D_gene", 1, "IG_J_gene", 1, "IG_J_pseudogene", 1, "IG_pseudogene", 1, "IG_V_gene", 1, "IG_V_pseudogene", 1, "TR_C_gene", 1, "TR_D_gene", 1, "TR_J_gene", 1, "TR_J_pseudogene", 1, "TR_V_gene", 1, "TR_V_pseudogene", 1);
	my %temp_novel_cluster_count = ();
	my $temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tStarting to process $temp_chr...\n";
	my $temp_start_time = Benchmark->new;

	my $temp_tss_ref;
	if($temp_tss_region){
		($temp_tss_ref) = read_bed("${temp_dir}/query_${temp_chr}_tss.bed");
	}
	#Initial pass annotation: Identify novel genes not overlapping any reference genes and create the loci information by adding them to the reference
	my ($temp_ref_transcript_ref, $temp_ref_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}.gtf", $temp_tss_ref);
	my ($temp_query_transcript_ref, $temp_query_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}.gtf", $temp_tss_ref);
	classify_novel_initial(
		$temp_ref_gene_ref, 
		$temp_query_transcript_ref, 
		"${temp_dir}/query_${temp_chr}_transcript.gtf", 
		"${temp_dir}/query_${temp_chr}_transcript_tss.gtf", 
		"${temp_dir}/ref_${temp_chr}_gene.gtf", 
		"${temp_dir}/ref_${temp_chr}_gene_tss.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_ref.gtf", 
		"${temp_dir}/query_${temp_chr}_novel.gtf", 
		$temp_btp
	);
	recluster_novel(
		\%temp_novel_cluster_count, 
		"${temp_dir}/query_${temp_chr}_novel_ref.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf"
	);
	system("cat ${temp_dir}/ref_${temp_chr}.gtf ${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
	system("echo ${temp_dir}/ref_${temp_chr}_expanded.gtf > ${temp_dir}/ref_${temp_chr}_gtf_list.txt");
	system("perl $temp_niap_merge -i ${temp_dir}/ref_${temp_chr}_gtf_list.txt -o ${temp_dir}/ref_${temp_chr}_expanded_merged -b $temp_btp");
	my $temp_loci_ref = read_loci_id("${temp_dir}/ref_${temp_chr}_expanded_merged_id.txt");

	#1st-pass annotation: Identify annotated, fusion and novel transcripts of known genes and those of novel genes
	($temp_ref_transcript_ref, $temp_ref_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
	($temp_query_transcript_ref, $temp_query_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel.gtf", $temp_tss_ref);
	my $temp_map_ref = compare_transcript(
		\%temp_monoexonic_gene_type, 
		$temp_max_p, 
		"${temp_dir}/ref_${temp_chr}_expanded.gtf", 
		"${temp_dir}/query_${temp_chr}_novel.gtf", 
		$temp_ref_transcript_ref, 
		$temp_query_transcript_ref, 
		$temp_ref_gene_ref, 
		$temp_btp
	);
	evaluate_transcript(
		1, 
		\%temp_monoexonic_gene_type, 
		$temp_map_ref, $temp_loci_ref, 
		$temp_ref_transcript_ref, 
		$temp_query_transcript_ref, 
		$temp_ref_gene_ref, 
		"${temp_dir}/query_${temp_chr}_annotated.gtf", 
		"${temp_dir}/query_${temp_chr}_fusion.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf", 
		"${temp_dir}/query_${temp_chr}_novel.gtf"
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\t1st-pass annotation done for $temp_chr!\n";

	#2nd-pass annotation: Annotate novel transcripts using the annotated transcripts in addition to the reference
	my $temp_novel_ref_last_size = -1;
	my $temp_novel_last_size = -1;
	my $temp_novel_ref_cur_size = 0;
	my $temp_novel_cur_size = 0;
	my $temp_it = 1;
	while(($temp_novel_ref_last_size != $temp_novel_ref_cur_size) || ($temp_novel_last_size != $temp_novel_cur_size)){
		$temp_novel_ref_last_size = $temp_novel_ref_cur_size;
		$temp_novel_last_size = $temp_novel_cur_size;
		system("cat ${temp_dir}/query_${temp_chr}_annotated.gtf ${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
		my ($temp_ref_expanded_transcript_ref, $temp_ref_expanded_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
		my ($temp_query_novel_transcript_ref, $temp_query_novel_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel.gtf", $temp_tss_ref);
		$temp_map_ref = compare_transcript(
			\%temp_monoexonic_gene_type, 
			$temp_max_p, 
			"${temp_dir}/ref_${temp_chr}_expanded.gtf", 
			"${temp_dir}/query_${temp_chr}_novel.gtf", 
			$temp_ref_expanded_transcript_ref, 
			$temp_query_novel_transcript_ref, 
			$temp_ref_gene_ref, 
			$temp_btp
		);
		evaluate_transcript(
			2, 
			\%temp_monoexonic_gene_type, 
			$temp_map_ref, 
			$temp_loci_ref, 
			$temp_ref_expanded_transcript_ref, 
			$temp_query_novel_transcript_ref, 
			$temp_ref_expanded_gene_ref, 
			"${temp_dir}/query_${temp_chr}_annotated.gtf", 
			"${temp_dir}/query_${temp_chr}_fusion.gtf", 
			"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf", 
			"${temp_dir}/query_${temp_chr}_novel.gtf"
		);
		$temp_novel_ref_cur_size = -s "${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf";
		$temp_novel_cur_size = -s "${temp_dir}/query_${temp_chr}_novel.gtf";
		$temp_date = `date $temp_date_format`;
		chomp $temp_date;
		print "[NIAP_annotate $temp_date]\t2nd-pass annotation iteration $temp_it done for $temp_chr!\n";
		$temp_it++;
	}

	#3rd-pass annotation: Annotate novel transcripts using the fusion transcript
	$temp_novel_ref_last_size = -1;
	$temp_novel_last_size = -1;
	$temp_novel_ref_cur_size = 0;
	$temp_novel_cur_size = 0;
	$temp_it = 1;
	while(($temp_novel_ref_last_size != $temp_novel_ref_cur_size) || ($temp_novel_last_size != $temp_novel_cur_size)){
		$temp_novel_ref_last_size = $temp_novel_ref_cur_size;
		$temp_novel_last_size = $temp_novel_cur_size;
		system("cat ${temp_dir}/query_${temp_chr}_fusion.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
		my ($temp_ref_expanded_transcript_ref, $temp_ref_expanded_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
		my ($temp_query_novel_transcript_ref, $temp_query_novel_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel.gtf", $temp_tss_ref);
		$temp_map_ref = compare_transcript(
			\%temp_monoexonic_gene_type, 
			$temp_max_p, 
			"${temp_dir}/ref_${temp_chr}_expanded.gtf", 
			"${temp_dir}/query_${temp_chr}_novel.gtf", 
			$temp_ref_expanded_transcript_ref, 
			$temp_query_novel_transcript_ref, 
			$temp_ref_gene_ref, 
			$temp_btp
		);
		evaluate_transcript(
			3, 
			\%temp_monoexonic_gene_type, 
			$temp_map_ref, $temp_loci_ref, 
			$temp_ref_expanded_transcript_ref, 
			$temp_query_novel_transcript_ref, 
			$temp_ref_expanded_gene_ref, 
			"${temp_dir}/query_${temp_chr}_annotated.gtf", 
			"${temp_dir}/query_${temp_chr}_fusion.gtf", 
			"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf", 
			"${temp_dir}/query_${temp_chr}_novel.gtf"
		);
		$temp_novel_ref_cur_size = -s "${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf";
		$temp_novel_cur_size = -s "${temp_dir}/query_${temp_chr}_novel.gtf";
		$temp_date = `date $temp_date_format`;
		chomp $temp_date;
		print "[NIAP_annotate $temp_date]\t3nd-pass annotation iteration $temp_it done for $temp_chr!\n";
		$temp_it++;
	}

	#4th-pass annotation: Annotate novel monoexonic transcripts the overlapping monoexonic transcripts
	my ($temp_query_novel_transcript_ref, $temp_query_novel_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel.gtf", $temp_tss_ref);
	$temp_map_ref = compare_transcript(
		\%temp_monoexonic_gene_type, 
		$temp_max_p, 
		"${temp_dir}/ref_${temp_chr}.gtf", 
		"${temp_dir}/query_${temp_chr}_novel.gtf", 
		$temp_ref_transcript_ref, 
		$temp_query_novel_transcript_ref, 
		$temp_ref_gene_ref, 
		$temp_btp
	);
	evaluate_transcript(
		4, 
		\%temp_monoexonic_gene_type, 
		$temp_map_ref, 
		$temp_loci_ref, 
		$temp_ref_transcript_ref, 
		$temp_query_novel_transcript_ref, 
		$temp_ref_gene_ref, 
		"${temp_dir}/query_${temp_chr}_annotated.gtf", 
		"${temp_dir}/query_${temp_chr}_fusion.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf", 
		"${temp_dir}/query_${temp_chr}_novel.gtf"
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\t4th-pass annotation done for $temp_chr!\n";

	#5th-pass annotation: Annotate novel transcripts using TSS of reference, annotated and novel gene transcripts
	$temp_novel_ref_last_size = -1;
	$temp_novel_last_size = -1;
	$temp_novel_ref_cur_size = 0;
	$temp_novel_cur_size = 0;
	$temp_it = 1;
	while(($temp_novel_ref_last_size != $temp_novel_ref_cur_size) || ($temp_novel_last_size != $temp_novel_cur_size)){
		$temp_novel_ref_last_size = $temp_novel_ref_cur_size;
		$temp_novel_last_size = $temp_novel_cur_size;
		system("cat ${temp_dir}/ref_${temp_chr}.gtf ${temp_dir}/query_${temp_chr}_annotated.gtf ${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
		my ($temp_ref_expanded_transcript_ref, $temp_ref_expanded_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
		my ($temp_query_novel_transcript_ref, $temp_query_novel_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel.gtf", $temp_tss_ref);
		compare_tss(
			\%temp_monoexonic_gene_type, 
			$temp_max_p, 
			$temp_loci_ref, 
			$temp_ref_expanded_transcript_ref, 
			$temp_query_novel_transcript_ref, 
			$temp_ref_gene_ref, 
			"${temp_dir}/ref_${temp_chr}_expanded_transcript.gtf", 
			"${temp_dir}/query_${temp_chr}_novel_transcript.gtf", 
			"${temp_dir}/ref_${temp_chr}_expanded_tss.gtf", 
			"${temp_dir}/query_${temp_chr}_novel_tss.gtf", 
			"${temp_dir}/query_${temp_chr}_annotated.gtf", 
			"${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf", 
			"${temp_dir}/query_${temp_chr}_novel.gtf", 
			$temp_btp
		);
		$temp_novel_ref_cur_size = -s "${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf";
		$temp_novel_cur_size = -s "${temp_dir}/query_${temp_chr}_novel.gtf";
		$temp_date = `date $temp_date_format`;
		chomp $temp_date;
		print "[NIAP_annotate $temp_date]\t5th-pass annotation iteration $temp_it done for $temp_chr!\n";
		$temp_it++;
	}

	#Re-cluster novel transcripts
	recluster_novel(
		\%temp_novel_cluster_count, 
		"${temp_dir}/query_${temp_chr}_novel.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_reclustered.gtf"
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tNovel genes re-clustering done for $temp_chr!\n";

	#Classify novel gene loci
	system("cat ${temp_dir}/ref_${temp_chr}.gtf ${temp_dir}/query_${temp_chr}_annotated.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
	system("cat ${temp_dir}/query_${temp_chr}_novel_ref_reclustered.gtf >> ${temp_dir}/query_${temp_chr}_novel_reclustered.gtf");
	my ($temp_ref_expanded_transcript_ref, $temp_ref_expanded_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
	my ($temp_query_novel_transcript_reclustered_ref, $temp_query_novel_gene_reclustered_ref) = read_gtf("${temp_dir}/query_${temp_chr}_novel_reclustered.gtf", $temp_tss_ref);
	classify_novel(
		\%temp_monoexonic_gene_type, 
		\%temp_pc_gene_type, 
		$temp_max_bd, 
		$temp_ref_expanded_transcript_ref, 
		$temp_query_novel_transcript_reclustered_ref, 
		$temp_ref_expanded_gene_ref, 
		$temp_query_novel_gene_reclustered_ref, 
		"${temp_dir}/ref_${temp_chr}_expanded.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_reclustered.gtf", 
		"${temp_dir}/ref_${temp_chr}_expanded_gene.gtf", 
		"${temp_dir}/ref_${temp_chr}_expanded_gene_tss.gtf", 
		"${temp_dir}/ref_${temp_chr}_expanded_transcript.gtf", 
		"${temp_dir}/ref_${temp_chr}_expanded_transcript_tss.gtf", 
		"${temp_dir}/ref_${temp_chr}_expanded_intron.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_gene_reclustered.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_gene_reclustered_tss.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_transcript_reclustered.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_transcript_reclustered_tss.gtf", 
		"${temp_dir}/query_${temp_chr}_novel_classified.gtf", 
		$temp_btp
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tClassifying novel genes done for $temp_chr!\n";

	#Re-arrange the refernece gene IDs for fusion genes
	system("cat ${temp_dir}/ref_${temp_chr}.gtf ${temp_dir}/query_${temp_chr}_annotated.gtf ${temp_dir}/query_${temp_chr}_novel_classified.gtf > ${temp_dir}/ref_${temp_chr}_expanded.gtf");
	($temp_ref_expanded_transcript_ref, $temp_ref_expanded_gene_ref) = read_gtf("${temp_dir}/ref_${temp_chr}_expanded.gtf", $temp_tss_ref);
	my ($temp_query_fusion_transcript_ref, $temp_query_fusion_gene_ref) = read_gtf("${temp_dir}/query_${temp_chr}_fusion.gtf", $temp_tss_ref);
	$temp_map_ref = compare_transcript(
		\%temp_monoexonic_gene_type, 
		$temp_max_p, 
		"${temp_dir}/ref_${temp_chr}_expanded.gtf", 
		"${temp_dir}/query_${temp_chr}_fusion.gtf", 
		$temp_ref_expanded_transcript_ref, 
		$temp_query_fusion_transcript_ref, 
		$temp_ref_expanded_gene_ref, 
		$temp_btp
	);
	rearrange_fusion(
		\%temp_monoexonic_gene_type,
		$temp_map_ref, 
		$temp_ref_expanded_transcript_ref, 
		$temp_ref_expanded_gene_ref, 
		$temp_query_fusion_transcript_ref, 
		"${temp_dir}/query_${temp_chr}_fusion.gtf"
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tRe-arranging gene IDs of fusion genes done for $temp_chr!\n";

	#Print exons
	if($temp_missing){
		print_missing($temp_ref_transcript_ref, "${temp_dir}/ref_${temp_chr}_missing.gtf");
		system("cat ${temp_dir}/query_${temp_chr}_annotated.gtf ${temp_dir}/query_${temp_chr}_fusion.gtf ${temp_dir}/query_${temp_chr}_novel_classified.gtf ${temp_dir}/ref_${temp_chr}_missing.gtf > ${temp_dir}/${temp_chr}_out.gtf");
	}
	else{
		system("cat ${temp_dir}/query_${temp_chr}_annotated.gtf ${temp_dir}/query_${temp_chr}_fusion.gtf ${temp_dir}/query_${temp_chr}_novel_classified.gtf > ${temp_dir}/${temp_chr}_out.gtf");
	}
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tMerged output .gtf for $temp_chr!\n";

	#Print genes and transcripts
	my ($temp_out_transcript_ref, $temp_out_gene_ref) = read_gtf("${temp_dir}/${temp_chr}_out.gtf", $temp_tss_ref);
	print_gene_transcript(
		$temp_out_transcript_ref, 
		$temp_out_gene_ref, 
		"${temp_dir}/${temp_chr}_out.gtf"
	);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tAdded genes and transcripts for $temp_chr!\n";

	my $temp_finish_time = Benchmark->new;
	my $temp_time_spent = timediff($temp_finish_time, $temp_start_time);
	$temp_date = `date $temp_date_format`;
	chomp $temp_date;
	print "[NIAP_annotate $temp_date]\tDone for $temp_chr in " . timestr($temp_time_spent) , "!\n";

	return $temp_chr;
}

#Split annotation
sub split_gtf{
	my ($temp_ref, $temp_query, $temp_dir) = @_;
	system("grep -v \"^#\" $temp_query | cut -f1  | sort | uniq > ${temp_dir}/chr_list.txt");
	my @temp_chr = ();
	open(my $temp_fh, "${temp_dir}/chr_list.txt") or die "Cannot open ${temp_dir}/chr_list.txt!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		push(@temp_chr, $temp_line);
	}
	close $temp_fh;
	system("for i in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$i\"\$\'\\t\' $temp_ref > ${temp_dir}/ref_\$\{i\}.gtf; done");
	system("for i in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$i\"\$\'\\t\' $temp_query > ${temp_dir}/query_\$\{i\}.gtf; done");

	return \@temp_chr;
}

#Split bed file
sub split_bed{
	my ($temp_bed, $temp_dir) = @_;
	system("for i in `cat ${temp_dir}/chr_list.txt`; do grep \"^\$i\"\$\'\\t\' $temp_bed > ${temp_dir}/query_\$\{i\}_tss.bed; done");

	return 1;
}

#Load gtf
sub read_gtf{
	my ($temp_gtf, $temp_tss_ref) = @_;
	open(my $temp_fh, "sort -k 1,1 -k 4,4n $temp_gtf |") or die "Cannot open $temp_gtf!\n";
	my %temp_transcript = ();
	my %temp_gene = ();
	my $temp_gene_name_flag = 0;
	my $temp_gene_type_flag = 0;
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[2] ne 'exon'){
			next;
		}
		my $temp_num_col = @temp_line;
		if($temp_num_col != 9){
			print STDERR "<ERROR> $temp_gtf\n$temp_line\nThe number of columns does not match the gtf format.\n";
			exit;
		}
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find transcript ID in $temp_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find gene ID in $temp_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_gene_name;
		if($temp_line[8] =~ /(^|\s)gene_name \"([^\"]+)\"\;/){
			$temp_gene_name = $2;
			$temp_gene_name_flag = 1;
		}
		else{
			$temp_gene_name = $temp_gene_id;
		}
		my $temp_gene_type;
		if($temp_line[8] =~ /(^|\s)gene_type \"([^\"]+)\"\;/){
			$temp_gene_type = $2;
			$temp_gene_type_flag = 1;
		}
		else{
			$temp_gene_type = 'na';
		}
		my $temp_status;
		if($temp_line[8] =~ /(^|\s)status \"([^\"]+)\"\;/){
			$temp_status = $2;
		}
		my $temp_ref_gene_id;
		if($temp_line[8] =~ /(^|\s)ref_gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		my $temp_ref_transcript_id;
		if($temp_line[8] =~ /(^|\s)ref_transcript_id \"([^\"]+)\"\;/){
			$temp_ref_transcript_id = $2;
		}
		my $temp_ref_gene_name;
		if($temp_line[8] =~ /(^|\s)ref_gene_name \"([^\"]+)\"\;/){
			$temp_ref_gene_name = $2;
		}
		my $temp_ref_gene_type;
		if($temp_line[8] =~ /(^|\s)ref_gene_type \"([^\"]+)\"\;/){
			$temp_ref_gene_type = $2;
		}
		my $temp_asm_gene_id;
		if($temp_line[8] =~ /(^|\s)asm_gene_id \"([^\"]+)\"\;/){
			$temp_asm_gene_id = $2;
		}
		my $temp_novel_gene_tag = 0;
		if($temp_line[8] =~ /(^|\s)tag \"novel_gene\"\;/){
			$temp_novel_gene_tag = 1;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'gene_id'} = $temp_gene_id;
			$temp_transcript{$temp_transcript_id}{'gene_name'} = $temp_gene_name;
			$temp_transcript{$temp_transcript_id}{'gene_type'} = $temp_gene_type;
			$temp_transcript{$temp_transcript_id}{'chr'} = $temp_line[0];
			$temp_transcript{$temp_transcript_id}{'source'} = $temp_line[1];
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
			if(exists $temp_tss_ref->{$temp_transcript_id}){
				$temp_transcript{$temp_transcript_id}{'tss_left'} = $temp_tss_ref->{$temp_transcript_id}{'left'};
				$temp_transcript{$temp_transcript_id}{'tss_right'} = $temp_tss_ref->{$temp_transcript_id}{'right'};
			}
			$temp_transcript{$temp_transcript_id}{'length'} = 0;
			$temp_transcript{$temp_transcript_id}{'strand'} = $temp_line[6];
			$temp_transcript{$temp_transcript_id}{'status'} = $temp_status;
			$temp_transcript{$temp_transcript_id}{'ref_gene_id'} = $temp_ref_gene_id;
			$temp_transcript{$temp_transcript_id}{'ref_transcript_id'} = $temp_ref_transcript_id;
			$temp_transcript{$temp_transcript_id}{'ref_gene_name'} = $temp_ref_gene_name;
			$temp_transcript{$temp_transcript_id}{'ref_gene_type'} = $temp_ref_gene_type;
			$temp_transcript{$temp_transcript_id}{'asm_gene_id'} = $temp_asm_gene_id;
			$temp_transcript{$temp_transcript_id}{'hit'} = 0;
			$temp_transcript{$temp_transcript_id}{'novel_gene'} = $temp_novel_gene_tag;
			@{$temp_transcript{$temp_transcript_id}{'exon'}} = ();
			#print "$temp_transcript_id\t$temp_transcript{$temp_transcript_id}{'gene_id'}\n";
		}
		my $temp_length = $temp_line[4] - $temp_line[3] + 1;
		$temp_transcript{$temp_transcript_id}{'length'} += $temp_length;
		if($temp_transcript{$temp_transcript_id}{'left'} > $temp_line[3]){
			$temp_transcript{$temp_transcript_id}{'left'} = $temp_line[3];
		}
		if($temp_transcript{$temp_transcript_id}{'right'} < $temp_line[4]){
			$temp_transcript{$temp_transcript_id}{'right'} = $temp_line[4];
		}
		my @temp_exon = ($temp_line[3], $temp_line[4]);
		push(@{$temp_transcript{$temp_transcript_id}{'exon'}}, \@temp_exon);
		unless(exists $temp_gene{$temp_gene_id}){
			%{$temp_gene{$temp_gene_id}} = ();
			$temp_gene{$temp_gene_id}{'gene_name'} = $temp_gene_name;
			$temp_gene{$temp_gene_id}{'gene_type'} = $temp_gene_type;
			$temp_gene{$temp_gene_id}{'chr'} = $temp_line[0];
			$temp_gene{$temp_gene_id}{'source'} = $temp_line[1];
			$temp_gene{$temp_gene_id}{'left'} = $temp_line[3];
			$temp_gene{$temp_gene_id}{'right'} = $temp_line[4];
			$temp_gene{$temp_gene_id}{'strand'} = $temp_line[6];
			$temp_gene{$temp_gene_id}{'ref_gene_id'} = $temp_ref_gene_id;
			$temp_gene{$temp_gene_id}{'ref_gene_name'} = $temp_ref_gene_name;
			$temp_gene{$temp_gene_id}{'ref_gene_type'} = $temp_ref_gene_type;
			$temp_gene{$temp_gene_id}{'asm_gene_id'} = $temp_asm_gene_id;
			$temp_gene{$temp_gene_id}{'multiexonic'} = 0;
			%{$temp_gene{$temp_gene_id}{'transcript_id'}} = ();
		}
		if($temp_gene{$temp_gene_id}{'left'} > $temp_line[3]){
			$temp_gene{$temp_gene_id}{'left'} = $temp_line[3];
		}
		if($temp_gene{$temp_gene_id}{'right'} < $temp_line[4]){
			$temp_gene{$temp_gene_id}{'right'} = $temp_line[4];
		}
		unless(exists $temp_gene{$temp_gene_id}{'transcript_id'}{$temp_transcript_id}){
			$temp_gene{$temp_gene_id}{'transcript_id'}{$temp_transcript_id} = 1;
		}
		if(($temp_gene{$temp_gene_id}{'multiexonic'} == 0) && (@{$temp_transcript{$temp_transcript_id}{'exon'}} > 1)){
			$temp_gene{$temp_gene_id}{'multiexonic'} = 1;
		}
	}
	close $temp_fh;
	foreach my $temp_transcript_id (keys %temp_transcript){
		my $temp_gene_id = $temp_transcript{$temp_transcript_id}{'gene_id'};
		unless(exists $temp_transcript{$temp_transcript_id}{'tss_left'}){
			if($temp_transcript{$temp_transcript_id}{'strand'} eq '-'){
				$temp_transcript{$temp_transcript_id}{'tss_left'} = $temp_transcript{$temp_transcript_id}{'right'};
				$temp_transcript{$temp_transcript_id}{'tss_right'} = $temp_transcript{$temp_transcript_id}{'right'};
			}
			else{
				$temp_transcript{$temp_transcript_id}{'tss_left'} = $temp_transcript{$temp_transcript_id}{'left'};
				$temp_transcript{$temp_transcript_id}{'tss_right'} = $temp_transcript{$temp_transcript_id}{'left'};
			}
		}
		if(exists $temp_gene{$temp_gene_id}{'tss_left'}){
			if($temp_gene{$temp_gene_id}{'tss_left'} > $temp_transcript{$temp_transcript_id}{'tss_left'}){
				$temp_gene{$temp_gene_id}{'tss_left'} = $temp_transcript{$temp_transcript_id}{'tss_left'};
			}
			if($temp_gene{$temp_gene_id}{'tss_right'} < $temp_transcript{$temp_transcript_id}{'tss_right'}){
				$temp_gene{$temp_gene_id}{'tss_right'} = $temp_transcript{$temp_transcript_id}{'tss_right'};
			}
		}
		else{
			$temp_gene{$temp_gene_id}{'tss_left'} = $temp_transcript{$temp_transcript_id}{'tss_left'};
			$temp_gene{$temp_gene_id}{'tss_right'} = $temp_transcript{$temp_transcript_id}{'tss_right'};
		}
	}
	if($temp_gene_name_flag == 0){
		print STDERR "<WARNING> gene_name attribute missing in $temp_gtf\n";
	}
	if($temp_gene_type_flag == 0){
		print STDERR "<WARNING> gene_type attribute missing in $temp_gtf\n";
	}
	
	return (\%temp_transcript, \%temp_gene);
}

#Load bed file
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
			print STDERR "<ERROR> Duplicated transcript, $temp_line[3], in the TSS bed file, $temp_bed!\n";
			exit;
		}
		else{
			%{$temp_bed{$temp_line[3]}} = ();
			$temp_bed{$temp_line[3]}{'chr'} = $temp_line[0];
			$temp_bed{$temp_line[3]}{'strand'} = $temp_line[5];
			$temp_bed{$temp_line[3]}{'left'} = $temp_line[1] - 1;
			$temp_bed{$temp_line[3]}{'right'} = $temp_line[2];
		}
	}
	close $temp_fh;

	return (\%temp_bed);
}

#Classify novel loci for the initial pass annotation
sub classify_novel_initial{
	my ($temp_ref_gene_ref, $temp_query_transcript_ref, $temp_query_transcript_gtf, $temp_query_transcript_tss_gtf, $temp_ref_gene_gtf, $temp_ref_gene_tss_gtf, $temp_query_novel_ref_gtf, $temp_query_novel_gtf, $temp_btp) = @_;
	#Create gene annotation for the reference
	print_gene_tss(
		"", 
		"", 
		$temp_ref_gene_ref, 
		$temp_ref_gene_gtf, 
		$temp_ref_gene_tss_gtf
	);
	#Create transcript annotation for the query
	print_transcript_tss(
		"", 
		"", 
		$temp_query_transcript_ref, 
		$temp_query_transcript_gtf, 
		$temp_query_transcript_tss_gtf
	);
	open(my $temp_fh, "$temp_btp intersect -s -v -wa -nonamecheck -a $temp_query_transcript_gtf -b $temp_ref_gene_gtf |");
	open(my $temp_query_novel_ref_fh, "> $temp_query_novel_ref_gtf") or die "Cannot create $temp_query_novel_ref_gtf!\n";
	open(my $temp_query_novel_fh, "> $temp_query_novel_gtf") or die "Cannot create $temp_query_novel_gtf!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		$temp_line =~ s/\r//;
		if($temp_line =~ /(^|\s)transcript_id \"([^\"]+)\"/){
			my $temp_query_transcript_id = $2;
			print_exon(
				$temp_query_novel_ref_fh, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}, 
				$temp_query_transcript_id, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
				"", 
				"", 
				"", 
				"", 
				"", 
				"", 
				"novel_gene"
			);
			delete($temp_query_transcript_ref->{$temp_query_transcript_id});
		}
		else{
			print STDERR "<ERROR> transcript_id missing when inntersecting $temp_query_transcript_gtf and $temp_ref_gene_gtf for\n$temp_line\n";
			exit;
		}
	}
	foreach my $temp_query_transcript_id (keys %$temp_query_transcript_ref){
		print_exon(
			$temp_query_novel_fh, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}, 
			$temp_query_transcript_id, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
			"", 
			"", 
			"", 
			"", 
			"", 
			"", 
			""
		);
		delete $temp_query_transcript_ref->{$temp_query_transcript_id};
	}
	close $temp_query_novel_ref_fh;
	close $temp_query_novel_fh;
	close $temp_fh;

	return 1;
}

#Load the *_loci_count.txt file
sub read_loci_id{
	my ($temp_file) = @_;
	my %temp_loci = ();
	open(my $temp_fh, $temp_file) or die "Cannot open $temp_file!\n";
	my %temp_header = ();
	my $temp_header = 0;
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_header == 0){
			my $temp_num_col = @temp_line;
			for(my $i = 0; $i < $temp_num_col; $i++){
				if(exists $temp_header{$temp_line[$i]}){
					print STDERR "<ERROR> Duplicated column $temp_line[$i] in $temp_file!\n";
					exit;
				}
				else{
					$temp_header{$temp_line[$i]} = $i;
				}
			}
			$temp_header = 1;

			next;
		}
		my @temp_id = split(',', $temp_line[$temp_header{'original_id'}]);
		foreach my $temp_id (@temp_id){
			unless(exists $temp_loci{$temp_id}){
				$temp_loci{$temp_id} = $temp_line[$temp_header{'gene_id'}];
			}
		}
	}
	close $temp_fh;

	return \%temp_loci;
}

#Look for the current exon number
sub cur_exon_num{
	my ($temp_exon_ref, $temp_left, $temp_right, $temp_strand) = @_;
	my $temp_exon_num = @$temp_exon_ref;
	for(my $i = 0; $i < $temp_exon_num; $i++){
		if(($temp_exon_ref->[$i][0] == $temp_left) && ($temp_exon_ref->[$i][1] == $temp_right)){
			my $temp_cur_exon_num;
			if($temp_strand eq "-"){
				$temp_cur_exon_num = $temp_exon_num - $i;
			}
			else{
				$temp_cur_exon_num = $i + 1;
			}

			return $temp_cur_exon_num;
		}
	}

	return 0;
}

#Map query transcirpts to reference
sub compare_transcript{
	my ($temp_monoexonic_gene_type_ref, $temp_max_p, $temp_ref_gtf, $temp_query_gtf, $temp_ref_transcript_ref, $temp_query_transcript_ref, $temp_ref_gene_ref, $temp_btp) = @_;
	my %temp_map = ();
	open(my $temp_fh, "$temp_btp intersect -wao -s -nonamecheck -a $temp_query_gtf -b $temp_ref_gtf |") or die "Cannot interesect $temp_query_gtf and $temp_ref_gtf!\n";
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[18] == 0) || ($temp_line[2] ne "exon") || ($temp_line[11] ne "exon")){
			next;
		}
		my $temp_query_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_query_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query transcript ID when intersecting $temp_query_gtf and $temp_ref_gtf at line:\n$temp_line\n";
			exit;
		}
		if($temp_line[6] eq "\."){
			print STDERR "<WARNING> $temp_query_transcript_id is non-strand-specific!\n";
		}
		#print "$temp_query_transcript_id\n";
		my $temp_query_exon_num = @{$temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}};
		#print "$temp_query_exon_num\n";
		#Get reference transcript name
		#Amendment is needed for alternative reference annotation.
		my $temp_ref_transcript_id;
		if($temp_line[17] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_ref_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference transcript ID when intersecting $temp_query_gtf and $temp_ref_gtf at line:\n$temp_line\n";
			exit;
		}
		#print "$temp_line\n$temp_query_transcript_id\t$temp_ref_transcript_id\n";
		my $temp_ref_exon_num = @{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'exon'}};
		#Mapping query transcript to reference
		unless(exists $temp_map{$temp_query_transcript_id}){
			%{$temp_map{$temp_query_transcript_id}} = ();
		}
		unless(exists $temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}){
			%{$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}} = ();
			@{$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query'}} = (0) x $temp_query_exon_num;
			@{$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref'}} = (0) x $temp_ref_exon_num;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'} = 0;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_match'} = 0;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'overlap'} = 0;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_overlap'} = 0;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'transcript_tss_dist'} = 999999999;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'gene_tss_dist'} = 999999999;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'known'} = 0;
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'partial'} = 0;
		}
		#Get the exon number of the current query exon
		my $temp_cur_query_exon_num = cur_exon_num(
			\@{$temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}}, 
			$temp_line[3], 
			$temp_line[4], 
			$temp_line[6]
		);
		if($temp_cur_query_exon_num == 0){
			print STDERR "<ERROR> Cannot find query exon $temp_line[0]:$temp_line[3]-$temp_line[4]:$temp_line[6] for query transcript $temp_query_transcript_id\n";
			exit;
		}
		#Get the exon number of the current reference exon
		my $temp_cur_ref_exon_num = cur_exon_num(
			\@{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'exon'}}, 
			$temp_line[12], 
			$temp_line[13], 
			$temp_line[15]
		);
		if($temp_cur_ref_exon_num == 0){
			print STDERR "<ERROR> Cannot find reference exon $temp_line[9]:$temp_line[12]:$temp_line[13]:$temp_line[15] for reference transcript $temp_ref_transcript_id\n";
			exit;
		}
		my $temp_match = 0;
		if(($temp_ref_exon_num == 1) && ($temp_query_exon_num == 1)){
			#If both query and reference transcripts are monoexonic, check TSSs within the distance.
			my $temp_dist = tss_dist(
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_left'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_right'}
			);
			#print "$temp_query_transcript_id\t$temp_ref_transcript_id\n$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}, $temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}, $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_left'}, $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_right'}\n$temp_dist\n";
			if($temp_dist <= $temp_max_p){
				$temp_match = 1;
			}
		}
		elsif(($temp_ref_exon_num == 1) && ($temp_cur_query_exon_num == 1) && !(exists $temp_monoexonic_gene_type_ref->{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'}}) && (($temp_ref_gene_ref->{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_id'}}{'multiexonic'} > 0) || ($temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'} eq "protein_coding"))){
			#If a reference transcript is monoexonic, excluding monoexonic noncoding genes and non-protein-coding monoexonic genes, and overlapping the first exon of a multiexonic query transcript, check TSSs within the distance.
			my $temp_dist = tss_dist(
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_left'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_right'}
			);
			if($temp_dist <= $temp_max_p){
				$temp_match = 1;
			}
		}
		elsif(($temp_ref_exon_num > 1) && ($temp_query_exon_num > 1)){
			#Both query and reference transcripts are multiexonic.
			if((($temp_cur_query_exon_num == 1) && ($temp_cur_ref_exon_num != $temp_ref_exon_num)) || (($temp_cur_ref_exon_num == 1) && ($temp_cur_query_exon_num != $temp_query_exon_num))){
				if($temp_line[6] eq "-"){
					#The rightmost case
					if($temp_line[3] == $temp_line[12]){
						$temp_match = 1;
					}
				}
				else{
					#The leftmost case
					if($temp_line[4] == $temp_line[13]){
						$temp_match = 1;
					}
				}
			}
			elsif((($temp_cur_query_exon_num == $temp_query_exon_num) && ($temp_cur_ref_exon_num != 1)) || (($temp_cur_ref_exon_num == $temp_ref_exon_num) && ($temp_cur_query_exon_num != 1))){
				if($temp_line[6] eq '-'){
					#The leftmost case
					if($temp_line[4] == $temp_line[13]){
						$temp_match = 1;
					}
				}
				else{
					#The rightmost case
					if($temp_line[3] == $temp_line[12]){
						$temp_match = 1;
					}
				}
			}
			else{
				#Internal exon
				if(($temp_line[3] == $temp_line[12]) && ($temp_line[4] == $temp_line[13])){
					#Matching both exon ends
					$temp_match = 1;
				}
			}
		}
		if($temp_match == 1){
			if($temp_line[6] eq '-'){
				$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query'}[($temp_query_exon_num - $temp_cur_query_exon_num)] = 1;
				$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref'}[($temp_ref_exon_num - $temp_cur_ref_exon_num)] = 1;
			}
			else{
				$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query'}[($temp_cur_query_exon_num - 1)] = 1;
				$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref'}[($temp_cur_ref_exon_num - 1)] = 1;
			}
			$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'}++;
		}
		$temp_map{$temp_query_transcript_id}{$temp_ref_transcript_id}{'overlap'} += $temp_line[18];
	}
	close $temp_fh;

	return \%temp_map;
}

#Print exons to .gtf
sub print_exon{
	my ($temp_fh, $temp_transcript_ref, $temp_transcript_id, $temp_gene_id, $temp_gene_name, $temp_gene_type, $temp_status, $temp_asm_gene_id, $temp_ref_transcript_id, $temp_ref_gene_id, $temp_ref_gene_name, $temp_ref_gene_type, $temp_tag) = @_;
	my $temp_num_exon = 1;
	if($temp_transcript_ref->{'strand'} eq '-'){
		$temp_num_exon = @{$temp_transcript_ref->{'exon'}};
	}
	foreach my $temp_exon_ref (@{$temp_transcript_ref->{'exon'}}){
		print $temp_fh "$temp_transcript_ref->{'chr'}\tNIAP\texon\t$temp_exon_ref->[0]\t$temp_exon_ref->[1]\t.\t$temp_transcript_ref->{'strand'}\t.\ttranscript_id \"$temp_transcript_id\"\; exon_num \"$temp_num_exon\"\; gene_id \"$temp_gene_id\"\; gene_name \"$temp_gene_name\"\; gene_type \"$temp_gene_type\"\;";
		if($temp_asm_gene_id){
			print $temp_fh " asm_gene_id \"$temp_asm_gene_id\"\;";
		}
		if($temp_status){
			print $temp_fh " status \"$temp_status\"\;";
		}
		if($temp_ref_transcript_id){
			print $temp_fh " ref_transcript_id \"$temp_ref_transcript_id\"\;";
		}
		if($temp_ref_gene_id){
			print $temp_fh " ref_gene_id \"$temp_ref_gene_id\"\;";
		}
		if($temp_ref_gene_name){
			print $temp_fh " ref_gene_name \"$temp_ref_gene_name\"\;";
		}
		elsif($temp_transcript_ref->{'ref_gene_name'}){
			print $temp_fh " ref_gene_name \"$temp_transcript_ref->{'ref_gene_name'}\"\;";
		}
		if($temp_ref_gene_type){
			print $temp_fh " ref_gene_type \"$temp_ref_gene_type\"\;";
		}
		elsif($temp_transcript_ref->{'ref_gene_type'}){
			print $temp_fh " ref_gene_type \"$temp_transcript_ref->{'ref_gene_type'}\"\;";
		}
		if($temp_tag =~ /novel_gene/){
			print $temp_fh " tag \"novel_gene\"\;";
		}
		print $temp_fh "\n";
		if($temp_transcript_ref->{'strand'} eq '-'){
			$temp_num_exon--;
		}
		else{
			$temp_num_exon++;
		}
	}

	return 1;
}

#Add genes and transcripts to .gtf
sub print_gene_transcript{
	my ($temp_transcript_ref, $temp_gene_ref, $temp_gtf) = @_;
	open(my $temp_fh, ">> $temp_gtf") or die "Cannot open $temp_gtf!\n";
	foreach my $temp_gene_id (keys %$temp_gene_ref){
		print $temp_fh "$temp_gene_ref->{$temp_gene_id}{'chr'}\tNIAP\tgene\t$temp_gene_ref->{$temp_gene_id}{'left'}\t$temp_gene_ref->{$temp_gene_id}{'right'}\t.\t$temp_gene_ref->{$temp_gene_id}{'strand'}\t.\tgene_id \"$temp_gene_id\"\;";
		if($temp_gene_ref->{$temp_gene_id}{'ref_gene_id'}){
			print $temp_fh " ref_gene_id \"$temp_gene_ref->{$temp_gene_id}{'ref_gene_id'}\"\;";
		}
		if($temp_gene_ref->{$temp_gene_id}{'asm_gene_id'}){
			print $temp_fh " asm_gene_id \"$temp_gene_ref->{$temp_gene_id}{'asm_gene_id'}\"\;";
		}
		if($temp_gene_ref->{$temp_gene_id}{'gene_name'}){
			print $temp_fh " gene_name \"$temp_gene_ref->{$temp_gene_id}{'gene_name'}\"\;";
		}
		if($temp_gene_ref->{$temp_gene_id}{'gene_type'}){
			print $temp_fh " gene_type \"$temp_gene_ref->{$temp_gene_id}{'gene_type'}\"\;";
		}
		print $temp_fh "\n";
	}
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		print $temp_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\tNIAP\ttranscript\t$temp_transcript_ref->{$temp_transcript_id}{'left'}\t$temp_transcript_ref->{$temp_transcript_id}{'right'}\t.\t$temp_transcript_ref->{$temp_transcript_id}{'strand'}\t.\ttranscript_id \"$temp_transcript_id\"\; gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; gene_name \"$temp_transcript_ref->{$temp_transcript_id}{'gene_name'}\"\; gene_type \"$temp_transcript_ref->{$temp_transcript_id}{'gene_type'}\"\;";
		if($temp_transcript_ref->{$temp_transcript_id}{'asm_gene_id'}){
			print $temp_fh " asm_gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'asm_gene_id'}\"\;";
		}
		if($temp_transcript_ref->{$temp_transcript_id}{'status'}){
			print $temp_fh " status \"$temp_transcript_ref->{$temp_transcript_id}{'status'}\"\;";
		}
		if($temp_transcript_ref->{$temp_transcript_id}{'ref_transcript_id'}){
			print $temp_fh " ref_transcript_id \"$temp_transcript_ref->{$temp_transcript_id}{'ref_transcript_id'}\"\;";
		}
		if($temp_transcript_ref->{$temp_transcript_id}{'ref_gene_id'}){
			print $temp_fh " ref_gene_id \"$temp_transcript_ref->{$temp_transcript_id}{'ref_gene_id'}\"\;";
		}
		if($temp_transcript_ref->{$temp_transcript_id}{'ref_gene_name'}){
			print $temp_fh " ref_gene_name \"$temp_transcript_ref->{$temp_transcript_id}{'ref_gene_name'}\"\;";
		}
		if($temp_transcript_ref->{$temp_transcript_id}{'ref_gene_type'}){
			print $temp_fh " ref_gene_type \"$temp_transcript_ref->{$temp_transcript_id}{'ref_gene_type'}\"\;";
		}
		print $temp_fh "\n";
	}
	close $temp_fh;

	return 1;
}

#Get the 5' position in the query of a hit
sub get_5p{
	my ($temp_ref_transcript_ref, $temp_query_transcript_ref, $temp_ref_code_ref, $temp_query_code_ref) = @_;
	my ($temp_ref_pos, $temp_query_pos);
	if($temp_ref_transcript_ref->{'strand'} eq '-'){
		if((@{$temp_ref_transcript_ref->{'exon'}} == 1) && (@{$temp_query_transcript_ref->{'exon'}} == 1)){
			$temp_ref_pos = $temp_ref_transcript_ref->{'exon'}[0][1];
			$temp_query_pos = $temp_query_transcript_ref->{'exon'}[0][1];

			return min($temp_ref_pos, $temp_query_pos);
		}
		for(my $i = (@{$temp_ref_transcript_ref->{'exon'}} - 1); $i >= 0; $i--){
			if($temp_ref_code_ref->[$i] == 1){
				$temp_ref_pos = $temp_ref_transcript_ref->{'exon'}[$i][1];
				for(my $i = (@{$temp_query_transcript_ref->{'exon'}} - 1); $i >= 0; $i--){
					if($temp_query_code_ref->[$i] == 1){
						$temp_query_pos = $temp_query_transcript_ref->{'exon'}[$i][1];

						return min($temp_ref_pos, $temp_query_pos);
					}
				}
			}
		}
	}
	else{
		if((@{$temp_ref_transcript_ref->{'exon'}} == 1) && (@{$temp_query_transcript_ref->{'exon'}} == 1)){
			$temp_ref_pos = $temp_ref_transcript_ref->{'exon'}[0][0];
			$temp_query_pos = $temp_query_transcript_ref->{'exon'}[0][0];

			return max($temp_ref_pos, $temp_query_pos);
		}
		for(my $i = 0; $i < @{$temp_ref_transcript_ref->{'exon'}}; $i++){
			if($temp_ref_code_ref->[$i] == 1){
				$temp_ref_pos = $temp_ref_transcript_ref->{'exon'}[$i][0];
				for(my $i = 0; $i < @{$temp_query_transcript_ref->{'exon'}}; $i++){
					if($temp_query_code_ref->[$i] == 1){
						$temp_query_pos = $temp_query_transcript_ref->{'exon'}[$i][0];

						return max($temp_ref_pos, $temp_query_pos);
					}
				}
			}
		}
	}

	return 0;
}

#Evaluate transcript mapping
sub evaluate_transcript{
	my ($temp_pass, $temp_monoexonic_gene_type_ref, $temp_map_ref, $temp_loci_ref, $temp_ref_transcript_ref, $temp_query_transcript_ref, $temp_ref_gene_ref, $temp_query_annotated_gtf, $temp_query_fusion_gtf, $temp_query_novel_ref_gtf, $temp_query_novel_gtf) = @_;
	open(my $temp_annotated_fh, ">> $temp_query_annotated_gtf") or die "Cannot open $temp_query_annotated_gtf!\n";
	open(my $temp_fusion_fh, ">> $temp_query_fusion_gtf") or die "Cannot open $temp_query_fusion_gtf!\n";
	open(my $temp_novel_ref_fh, ">> $temp_query_novel_ref_gtf") or die "Cannot open $temp_query_novel_ref_gtf!\n";
	open(my $temp_novel_fh, "> $temp_query_novel_gtf") or die "Cannot create $temp_query_novel_gtf!\n";
	foreach my $temp_query_transcript_id (keys %$temp_map_ref){
		my $temp_known_flag = 0;
		my %temp_gene_id = ();
		my %temp_loci = ();
		my $temp_query_exon_num = @{$temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}};
		#Loop through the map and identify the best annotation
		foreach my $temp_ref_transcript_id (keys %{$temp_map_ref->{$temp_query_transcript_id}}){
			my $temp_ref_gene_id = $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_id'};
			my $temp_ref_gene_type = $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'};
			my $temp_ref_exon_num = @{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'exon'}};
			#print "T:$temp_query_transcript_id\t$temp_ref_transcript_id\t$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}\t$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}\t$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_left'}\t$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_right'}\n";
			$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'transcript_tss_dist'} = tss_dist(
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}, 
					$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_left'}, 
					$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'tss_right'}
			);
			#print "$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'transcript_tss_dist'}\n";
			foreach my $temp_all_ref_transcript_id (keys %{$temp_ref_gene_ref->{$temp_ref_gene_id}{'transcript_id'}}){
				#print "G:$temp_query_transcript_id\t$temp_all_ref_transcript_id\t$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}\t$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}\t$temp_ref_transcript_ref->{$temp_all_ref_transcript_id}{'tss_left'}\t$temp_ref_transcript_ref->{$temp_all_ref_transcript_id}{'tss_right'}\n";
				my $temp_dist = tss_dist(
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_left'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'tss_right'}, 
					$temp_ref_transcript_ref->{$temp_all_ref_transcript_id}{'tss_left'}, 
					$temp_ref_transcript_ref->{$temp_all_ref_transcript_id}{'tss_right'}
				);
				#print "$temp_dist\n";
				if($temp_dist < $temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'gene_tss_dist'}){
					$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'gene_tss_dist'} = $temp_dist;
					#print "Got value\n"
				}
			}
			my $temp_query_code = join('', @{$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query'}});
			my $temp_ref_code = join('', @{$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref'}});
			$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query_match'} = $temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'} / $temp_query_exon_num;
			$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_match'} = $temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'} / $temp_ref_exon_num;
			$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query_overlap'} = $temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'overlap'} / $temp_query_transcript_ref->{$temp_query_transcript_id}{'length'};
			$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_overlap'} = $temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'overlap'} / $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'length'};
			if(($temp_query_code eq $temp_ref_code) && !($temp_query_code =~ /0/) && !($temp_ref_code =~ /0/)){
				$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'known'} = 1;
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'hit'}++;
				$temp_known_flag = 1;
			}
			elsif(!($temp_query_code =~ /0/) && ($temp_ref_code =~ /$temp_query_code/)){
				$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'partial'} = 1;
			}

			#print "$temp_pass pass: $temp_query_transcript_id\t$temp_query_code\t$temp_ref_transcript_id\t$temp_ref_code\t$temp_ref_gene_id\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'known'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'partial'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_match'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'overlap'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref_overlap'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'transcript_tss_dist'}\t$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'gene_tss_dist'}\n";

			#Get the location of 5' for each gene and locus ID for the 1st-pass, excluding monoexonic noncoding genes.
			if($temp_pass == 1){
				if(((($temp_query_exon_num == 1) && ($temp_ref_exon_num == 1)) || ($temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'} > 0)) && !(exists $temp_monoexonic_gene_type_ref->{$temp_ref_gene_type})){
					unless(exists $temp_gene_id{$temp_ref_gene_id}){
						$temp_gene_id{$temp_ref_gene_id} = 1;
					}
					if(exists $temp_loci_ref->{$temp_ref_transcript_id}){
						unless(exists $temp_loci{$temp_loci_ref->{$temp_ref_transcript_id}}){
							
							$temp_loci{$temp_loci_ref->{$temp_ref_transcript_id}} = 1;
						}
					}
					#print "Added $temp_ref_gene_id\t$temp_gene_id{$temp_ref_gene_id}\t$temp_loci_ref->{$temp_ref_transcript_id}\n";
				}
			}
			elsif($temp_pass == 4){
				if(($temp_query_exon_num == 1) && ($temp_ref_exon_num == 1) && !(exists $temp_monoexonic_gene_type_ref->{$temp_ref_gene_type})){
					$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'known'} = 1;
					$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'hit'}++;
					$temp_known_flag = 1;
					#print "Hit in the 4th-pass\n";
				}
			}
		}
		#The sorting of transcripts is arbitary if the best hits are equivalent.
		my @temp_sorted_transcript_id = sort {
			($temp_map_ref->{$temp_query_transcript_id}{$b}{'known'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'known'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'partial'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'partial'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'match'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'match'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'query_match'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'query_match'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'query_overlap'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'query_overlap'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'ref_match'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'ref_match'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'ref_overlap'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'ref_overlap'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$b}{'overlap'} <=> $temp_map_ref->{$temp_query_transcript_id}{$a}{'overlap'}) 
			|| ($temp_map_ref->{$temp_query_transcript_id}{$a}{'gene_tss_dist'} <=> $temp_map_ref->{$temp_query_transcript_id}{$b}{'gene_tss_dist'}) 
			|| ($a cmp $b) 
		} keys %{$temp_map_ref->{$temp_query_transcript_id}};
		#print "Sorted:@temp_sorted_transcript_id\n";
		my $temp_ref_code = join('', @{$temp_map_ref->{$temp_query_transcript_id}{$temp_sorted_transcript_id[0]}{'ref'}});
		my $temp_query_code = join('', @{$temp_map_ref->{$temp_query_transcript_id}{$temp_sorted_transcript_id[0]}{'query'}});
		if($temp_pass == 1){
			#Known, partial and fusion transcripts of known genes and those of novel genes will be annotated in the 1st-pass annotation.
			my $temp_num_gene = keys %temp_gene_id;
			my $temp_num_loci = keys %temp_loci;
			if(($temp_num_gene > 1) && ($temp_num_loci > 1)){
				#Sort the gene IDs alphabetically.
				my @temp_sorted_gene_id = sort {($a cmp $b)} keys %temp_gene_id;
				my $temp_ref_gene_id = join("\|", @temp_sorted_gene_id);
				my @temp_sorted_gene_name = ();
				my @temp_sorted_gene_type = ();
				foreach my $temp_gene_id (@temp_sorted_gene_id){
					if($temp_ref_gene_ref->{$temp_gene_id}{'gene_name'}){
						push(@temp_sorted_gene_name, $temp_ref_gene_ref->{$temp_gene_id}{'gene_name'});
					}
					if($temp_ref_gene_ref->{$temp_gene_id}{'gene_type'}){
						push(@temp_sorted_gene_type, $temp_ref_gene_ref->{$temp_gene_id}{'gene_type'});
					}
				}
				my $temp_ref_gene_name = join("\|", @temp_sorted_gene_name);
				my $temp_ref_gene_type = join("\|", @temp_sorted_gene_type);
				print_exon(
					$temp_fusion_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_ref_gene_id, 
					$temp_ref_gene_name, 
					$temp_ref_gene_type, 
					"fusion", 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					"", 
					"", 
					"", 
					"", 
					""
				);
			}
			elsif($temp_known_flag == 1){
				print_exon(
					$temp_annotated_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
					"known", 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					$temp_sorted_transcript_id[0], 
					"", 
					"", 
					"", 
					""
				);
			}
			elsif($temp_map_ref->{$temp_query_transcript_id}{$temp_sorted_transcript_id[0]}{'partial'} == 1){
				if($temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'novel_gene'} == 1){
					#Isoforms of novel gene identified in the initial pass
					print_exon(
						$temp_novel_ref_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						"", 
						"", 
						"", 
						"", 
						"novel_gene"
					);
				}
				else{
					my $temp_status;
					#Novel isoform being partial of known isoforms
					my $temp_partial5e = 0;
					my $temp_partial3e = 0;
					my $temp_first_seg_start = 0;
					if($temp_ref_code =~ /(^0*)1+0*$/){
						my $temp_first_seg = $1;
						$temp_first_seg_start = length($temp_first_seg);
					}
					if($temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}[0][0] < $temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'exon'}[$temp_first_seg_start][0]){
						if($temp_query_transcript_ref->{$temp_query_transcript_id}{'strand'} eq "-"){
							$temp_partial3e = 1;
						}
						else{
							$temp_partial5e = 1;
						}
					}
					my $temp_last_seg_start = (length($temp_ref_code) - 1);
					if($temp_ref_code =~ /(^0*1+)0*$/){
						my $temp_last_seg = $1;
						$temp_last_seg_start = (length($temp_last_seg) - 1);
					}
					if($temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}[-1][1] > $temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'exon'}[$temp_last_seg_start][1]){
						if($temp_query_transcript_ref->{$temp_query_transcript_id}{'strand'} eq "-"){
							$temp_partial5e = 1;
						}
						else{
							$temp_partial3e = 1;
						}
					}
					if(($temp_partial5e == 1) && ($temp_partial3e == 1)){
						#Novel isoforms with varied transcript ends
						$temp_status = "partial_extended";
					}
					elsif(($temp_partial5e == 1) && ($temp_partial3e == 0)){
						#Novel isoform with varied 5' end
						$temp_status = "partial_5p_extended";
					}
					elsif(($temp_partial5e == 0) && ($temp_partial3e == 1)){
						#Novel isoform with varied 3' end
						$temp_status = "partial_3p_extended";
					}
					else{
						#Partial transcript
						$temp_status = "partial";
					}
					print_exon(
						$temp_annotated_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						$temp_status, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						$temp_sorted_transcript_id[0], 
						"", 
						"", 
						"", 
						""
					);
				}
			}
			elsif($temp_query_code =~ /1/){
				if($temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'novel_gene'} == 1){
					#Isoforms of novel gene identified in the initial pass
					print_exon(
						$temp_novel_ref_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						"", 
						"", 
						"", 
						"", 
						"novel_gene"
					);
				}
				else{
					#Novel isoform with novel transcript structure
					print_exon(
						$temp_annotated_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"novel_structure", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						$temp_sorted_transcript_id[0], 
						"", 
						"", 
						"", 
						""
					);
				}
			}
			else{
				#Isoform of novel gene
				print_exon(
					$temp_novel_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
					"", 
					"", 
					"", 
					"", 
					"", 
					"", 
					""
				);
			}
		}
		elsif(($temp_pass == 2) || ($temp_pass == 3)){
			#Annotate for 2nd and 3rd passes
			if($temp_query_code =~ /1/){
				if($temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'} =~ /\|/){
					#Novel isoform of fusion
					print_exon(
						$temp_fusion_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"fusion", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						$temp_sorted_transcript_id[0], 
						"", 
						"", 
						"", 
						""
					);
				}
				elsif($temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'novel_gene'} == 1){
					#Isoforms of novel gene identified in the initial pass
					print_exon(
						$temp_novel_ref_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						"", 
						"", 
						"", 
						"", 
						"novel_gene"
					);
				}
				else{
					#Novel isoform with novel transcript structure
					print_exon(
						$temp_annotated_fh, 
						$temp_query_transcript_ref->{$temp_query_transcript_id}, 
						$temp_query_transcript_id, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
						$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
						"novel_structure", 
						$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
						$temp_sorted_transcript_id[0], 
						"", 
						"", 
						"", 
						""
					);
				}
			}
			else{
				#Isoform of novel gene
				print_exon(
					$temp_novel_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
					"", 
					"", 
					"", 
					"", 
					"", 
					"", 
					""
				);
			}
		}
		elsif($temp_pass == 4){
			#Annotation for 4th pass
			if(($temp_query_exon_num == 1) && ($temp_known_flag == 1)){
				print_exon(
					$temp_annotated_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
					"known", 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					$temp_sorted_transcript_id[0], 
					"", 
					"", 
					"", 
					""
				);
			}
			else{
				#Isoform of novel gene
				print_exon(
					$temp_novel_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, $temp_query_transcript_id, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
					"", 
					"", 
					"", 
					"", 
					"", 
					"", 
					""
				);
			}
		}
		delete $temp_map_ref->{$temp_query_transcript_id};
		delete $temp_query_transcript_ref->{$temp_query_transcript_id};
	}
	foreach my $temp_query_transcript_id (keys %$temp_query_transcript_ref){
		#Isoform of novel gene
		print_exon(
			$temp_novel_fh, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}, 
			$temp_query_transcript_id, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
			"", 
			"", 
			"", 
			"", 
			"", 
			"", 
			""
		);
		delete $temp_query_transcript_ref->{$temp_query_transcript_id};
	}
	close $temp_annotated_fh;
	close $temp_fusion_fh;
	close $temp_novel_fh;

	return 1;
}

#Print TSS of transcript to .gtf
sub print_transcript_tss{
	my ($temp_monoexonic_gene_type_ref, $temp_ref_gene_type_ref, $temp_transcript_ref, $temp_transcript_gtf, $temp_transcript_tss_gtf, $temp_max_p) = @_;
	open(my $temp_transcript_fh, "| sort -k 1,1 -k 4,4n > $temp_transcript_gtf") or die "Cannot create $temp_transcript_gtf!\n";
	open(my $temp_tss_fh, "| sort -k 1,1 -k 4,4n > $temp_transcript_tss_gtf") or die "Cannot create $temp_transcript_tss_gtf!\n";
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		my $temp_gene_type = $temp_transcript_ref->{$temp_transcript_id}{'gene_type'};
		if($temp_monoexonic_gene_type_ref && $temp_ref_gene_type_ref){
			#Skip if not monoexonic or protein-coding gene
			unless((exists $temp_monoexonic_gene_type_ref->{$temp_gene_type}) || (exists $temp_ref_gene_type_ref->{$temp_gene_type})){
				next;
			}
		}
		print $temp_transcript_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\t$temp_transcript_ref->{$temp_transcript_id}{'source'}\ttranscript\t$temp_transcript_ref->{$temp_transcript_id}{'left'}\t$temp_transcript_ref->{$temp_transcript_id}{'right'}\t.\t$temp_transcript_ref->{$temp_transcript_id}{'strand'}\t.\tgene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; transcript_id \"$temp_transcript_id\"\;\n";
		print $temp_tss_fh "$temp_transcript_ref->{$temp_transcript_id}{'chr'}\t$temp_transcript_ref->{$temp_transcript_id}{'source'}\ttss\t$temp_transcript_ref->{$temp_transcript_id}{'tss_left'}\t$temp_transcript_ref->{$temp_transcript_id}{'tss_right'}\t.\t$temp_transcript_ref->{$temp_transcript_id}{'strand'}\t.\tgene_id \"$temp_transcript_ref->{$temp_transcript_id}{'gene_id'}\"\; transcript_id \"$temp_transcript_id\"\;\n";
	}
	close $temp_transcript_fh;
	close $temp_tss_fh;

	return 1;
}

#Calculate the distance between TSSs
sub tss_dist{
	my ($temp_a_left, $temp_a_right, $temp_b_left, $temp_b_right) = @_;
	my $temp_dist = 0;
	my @temp_a = ($temp_a_left, $temp_a_right);
	my @temp_b = ($temp_b_left, $temp_b_right);
	my @temp_tss = (\@temp_a, \@temp_b);
	my @temp_sorted_tss = sort {$a->[0] <=> $b->[0]} @temp_tss;
	#foreach my $temp_tss_ref (@temp_tss){
		#print "$temp_tss_ref->[0]-$temp_tss_ref->[1]\n";
	#}
	if($temp_sorted_tss[1][0] > $temp_sorted_tss[0][1]){
		#print "$temp_sorted_tss[1][0] > $temp_sorted_tss[0][1]\n";
		$temp_dist = $temp_sorted_tss[1][0] - $temp_sorted_tss[0][1];
	}

	return $temp_dist;
}

#Map query transcirpts to reference
sub compare_tss{
	my ($temp_monoexonic_gene_type_ref, $temp_max_p, $temp_loci_ref, $temp_ref_transcript_ref, $temp_query_transcript_ref, $temp_ref_gene_ref, $temp_ref_transcript_gtf, $temp_query_transcript_gtf, $temp_ref_transcript_tss_gtf, $temp_query_transcript_tss_gtf, $temp_query_annotated_gtf, $temp_query_novel_ref_gtf, $temp_query_novel_gtf, $temp_btp) = @_;
	print_transcript_tss(
		"", 
		"", 
		$temp_ref_transcript_ref, 
		$temp_ref_transcript_gtf, 
		$temp_ref_transcript_tss_gtf
	);
	print_transcript_tss(
		"", 
		"", 
		$temp_query_transcript_ref, 
		$temp_query_transcript_gtf, 
		$temp_query_transcript_tss_gtf
	);
	open(my $temp_fh, "$temp_btp closest -s -nonamecheck -D a -a $temp_query_transcript_tss_gtf -b $temp_ref_transcript_tss_gtf |");
	open(my $temp_annotated_fh, ">> $temp_query_annotated_gtf") or die "Cannot open $temp_query_annotated_gtf!\n";
	open(my $temp_novel_ref_fh, ">> $temp_query_novel_ref_gtf") or die "Cannot create $temp_query_novel_ref_gtf!\n";
	open(my $temp_novel_fh, "> $temp_query_novel_gtf") or die "Cannot create $temp_query_novel_gtf!\n";
	my %temp_dist = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my $temp_query_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_query_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query transcript ID when intersecting $temp_query_transcript_tss_gtf and $temp_ref_transcript_tss_gtf at line:\n$temp_line\n";
			exit;
		}
		if(abs($temp_line[18]) > $temp_max_p){
			next;
		}
		my $temp_ref_transcript_id;
		if($temp_line[17] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_ref_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference transcript ID when intersecting $temp_query_transcript_tss_gtf and $temp_ref_transcript_tss_gtf at line:\n$temp_line\n";
			next;
		}
		if((exists $temp_monoexonic_gene_type_ref->{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'}}) || ($temp_ref_gene_ref->{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_id'}}{'multiexonic'} == 0)){
			next;
		}
		unless(exists $temp_dist{$temp_query_transcript_id}){
			%{$temp_dist{$temp_query_transcript_id}} = ();
		}
		unless(exists $temp_dist{$temp_query_transcript_id}{$temp_ref_transcript_id}){
			$temp_dist{$temp_query_transcript_id}{$temp_ref_transcript_id} = $temp_line[18];
		}
	}
	close $temp_fh;
	foreach my $temp_query_transcript_id (keys %$temp_query_transcript_ref){
		if(exists $temp_dist{$temp_query_transcript_id}){
			my @temp_sorted_transcript_id = sort {($temp_dist{$temp_query_transcript_id}{$a} <=> $temp_dist{$temp_query_transcript_id}{$b}) || ($a cmp $b)} keys %{$temp_dist{$temp_query_transcript_id}};
			if($temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'novel_gene'} == 1){
				#Isoform of novel gene identified in the initial pass
				print_exon(
					$temp_novel_ref_fh, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_query_transcript_id, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, 
					$temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, 
					"", 
					$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
					"", 
					"", 
					"", 
					"", 
					"novel_gene"
				);
			}
			else{
				#Novel isoform with novel structure
				print_exon($temp_annotated_fh, $temp_query_transcript_ref->{$temp_query_transcript_id}, $temp_query_transcript_id, $temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_id'}, $temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_name'}, $temp_ref_transcript_ref->{$temp_sorted_transcript_id[0]}{'gene_type'}, "novel_structure", $temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, $temp_sorted_transcript_id[0], "", "", "", "");
			}
		}
		else{
			#Isoform of novel gene
			print_exon(
				$temp_novel_fh, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}, 
				$temp_query_transcript_id, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_name'}, 
				$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
				"", 
				"", 
				"", 
				"", 
				"", 
				"", 
				""
			);
		}
		delete $temp_query_transcript_ref->{$temp_query_transcript_id};
	}
	close $temp_annotated_fh;
	close $temp_novel_ref_fh;
	close $temp_novel_fh;

	return 1;
}

#Merge gene 
sub cluster_transcript{
	my ($temp_sref, $temp_start) = @_;
	if(!defined $temp_start){
		if(wantarray){
			my @temp_sets = map {[@{$_}]} @{$temp_sref};
			$temp_sref = \@temp_sets;
		}
		@{$temp_sref} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$temp_sref};
		$temp_start = 0;
	}
	my $temp_last = $temp_sref->[$temp_start];
	$temp_start++;

	if(@{$temp_last}){
		#print "@{$temp_last}\n";
		for(my $i = $temp_start; $i < @{$temp_sref}; $i++){
			my $temp_cur = $temp_sref->[$i];
			if(!@{$temp_cur}){
				next;
			}

			if (($temp_cur->[0] >= $temp_last->[0]) && ($temp_cur->[0] <= $temp_last->[1])){
				if ($temp_cur->[1] > $temp_last->[1]){
					$temp_last->[1] = $temp_cur->[1];
				}
				push(@{$temp_last->[2]}, @{$temp_cur->[2]});
				@{$temp_cur} = ();
			}
			else{
				last;
			}
		}
	}
	if($temp_start < @{$temp_sref}){
		cluster_transcript($temp_sref, $temp_start);
	}
	if(wantarray){
		return sort {$a->[0] <=> $b->[0]} map {@{$_} ? $_ : () } @{$temp_sref};
	}
}

#Re-cluster novel transcripts
sub recluster_novel{
	my ($temp_novel_cluster_count_ref, $temp_in_gtf, $temp_out_gtf) = @_;
	#Read gtf of novel transcripts. Cannot use the general function read_gtf for this.
	open(my $temp_fh, $temp_in_gtf) or die "Cannot open $temp_in_gtf!\n";
	my %temp_transcript = ();
	my %temp_gene = ();
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my $temp_transcript_id;
		if($temp_line[8] =~ /(^|\s)transcript_id \"([^\"]+)\"\;/){
			$temp_transcript_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query transcritp ID in $temp_in_gtf at line:\n$temp_line\n";
			exit;
		}
		my ($temp_asm_gene_id, $temp_gene_id);
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_asm_gene_id = $2;
			$temp_gene_id = $2;
			$temp_gene_id =~ s/_\d+$//;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID in $temp_in_gtf at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_transcript{$temp_transcript_id}){
			%{$temp_transcript{$temp_transcript_id}} = ();
			$temp_transcript{$temp_transcript_id}{'asm_gene_id'} = $temp_asm_gene_id;
			@{$temp_transcript{$temp_transcript_id}{'original'}} = ();
		}
		push(@{$temp_transcript{$temp_transcript_id}{'original'}}, $temp_line);
		unless(exists $temp_gene{$temp_gene_id}){
			%{$temp_gene{$temp_gene_id}} = ();
		}
		unless(exists $temp_gene{$temp_gene_id}{$temp_transcript_id}){
			@{$temp_gene{$temp_gene_id}{$temp_transcript_id}} = ($temp_line[3], $temp_line[4]);
		}
		if($temp_gene{$temp_gene_id}{$temp_transcript_id}[0] > $temp_line[3]){
			$temp_gene{$temp_gene_id}{$temp_transcript_id}[0] = $temp_line[3];
		}
		if($temp_gene{$temp_gene_id}{$temp_transcript_id}[1] < $temp_line[4]){
			$temp_gene{$temp_gene_id}{$temp_transcript_id}[1] = $temp_line[4];
		}
	}
	close $temp_fh;
	#Recluster novel transcripts
	open($temp_fh, "> $temp_out_gtf") or die "Cannot create $temp_out_gtf!\n";
	foreach my $temp_gene_id (keys %temp_gene){
		my @temp_transcript = ();
		foreach my $temp_transcript_id (keys %{$temp_gene{$temp_gene_id}}){
			my @temp_transcript_id = ($temp_transcript_id);
			my @temp_interval = ($temp_gene{$temp_gene_id}{$temp_transcript_id}[0], $temp_gene{$temp_gene_id}{$temp_transcript_id}[1], \@temp_transcript_id);
			push(@temp_transcript, \@temp_interval);
		}
		my @temp_transcript_cluster = cluster_transcript(\@temp_transcript);
		unless(exists $temp_novel_cluster_count_ref->{$temp_gene_id}){
			$temp_novel_cluster_count_ref->{$temp_gene_id} = 1;
		}
		foreach my $temp_cluster_ref (@temp_transcript_cluster){
			foreach my $temp_transcript_id (@{$temp_cluster_ref->[2]}){
				my $temp_novel_gene_id = "${temp_gene_id}_$temp_novel_cluster_count_ref->{$temp_gene_id}";
				foreach my $temp_line (@{$temp_transcript{$temp_transcript_id}{'original'}}){
					if($temp_line =~ /^(gene_id \"[^\"]+\")\;/){
						my $temp_gene_id = $1;
						$temp_gene_id =~ s/\|/\\\|/g;
						$temp_line =~ s/$temp_gene_id/gene_id \"$temp_novel_gene_id\"/;
					}
					elsif($temp_line =~ /( gene_id \"[^\"]+\")\;/){
						my $temp_gene_id = $1;
						$temp_gene_id =~ s/\|/\\\|/g;
						$temp_line =~ s/$temp_gene_id/ gene_id \"$temp_novel_gene_id\"/;
					}
					else{
						print STDERR "<ERROR> Cannot find query gene ID at line:\n$temp_line\n";	
						exit;
					}
					print $temp_fh "$temp_line asm_gene_id \"$temp_transcript{$temp_transcript_id}{'asm_gene_id'}\"\;\n";
				}
			}
			$temp_novel_cluster_count_ref->{$temp_gene_id}++;
		}
	}
	close $temp_fh;
	
	return 1;
}

#Print gene and gene TSS to .gtf
sub print_gene_tss{
	my ($temp_monoexonic_gene_type_ref, $temp_ref_gene_type_ref, $temp_gene_ref, $temp_gene_gtf, $temp_tss_gtf) = @_;
	open(my $temp_gene_fh, "| sort -k 1,1 -k 4,4n > $temp_gene_gtf") or die "Cannot create $temp_gene_gtf!\n";
	open(my $temp_tss_fh, "| sort -k 1,1 -k 4,4n > $temp_tss_gtf") or die "Cannot create $temp_tss_gtf!\n";
	foreach my $temp_gene_id (keys %$temp_gene_ref){
		my $temp_gene_type = $temp_gene_ref->{$temp_gene_id}{'gene_type'};
		if($temp_monoexonic_gene_type_ref && $temp_ref_gene_type_ref){
			#Skip if not monoexonic or protein-coding gene
			unless((exists $temp_monoexonic_gene_type_ref->{$temp_gene_type}) || (exists $temp_ref_gene_type_ref->{$temp_gene_type})){
				next;
			}
		}
		print $temp_gene_fh "$temp_gene_ref->{$temp_gene_id}{'chr'}\t$temp_gene_ref->{$temp_gene_id}{'source'}\tgene\t$temp_gene_ref->{$temp_gene_id}{'left'}\t$temp_gene_ref->{$temp_gene_id}{'right'}\t.\t$temp_gene_ref->{$temp_gene_id}{'strand'}\t.\tgene_id \"$temp_gene_id\"\;";
		if($temp_gene_type){
			print $temp_gene_fh " gene_type \"$temp_gene_type\"\;";
		}
		print $temp_gene_fh "\n";
		print $temp_tss_fh "$temp_gene_ref->{$temp_gene_id}{'chr'}\t$temp_gene_ref->{$temp_gene_id}{'source'}\ttss\t$temp_gene_ref->{$temp_gene_id}{'tss_left'}\t$temp_gene_ref->{$temp_gene_id}{'tss_right'}\t.\t$temp_gene_ref->{$temp_gene_id}{'strand'}\t.\tgene_id \"$temp_gene_id\"\;";
		if($temp_gene_type){
			print $temp_tss_fh " gene_type \"$temp_gene_type\"\;";
		}
		print $temp_tss_fh "\n";
	}
	close $temp_gene_fh;
	close $temp_tss_fh;

	return 1;
}

#Merge exons in a gene
sub merge_exon{
	my ($temp_exon_ref, $temp_start) = @_;

	if(!defined $temp_start) {
		my @temp_exon = map {[@{$_}]} @{$temp_exon_ref};
		$temp_exon_ref = \@temp_exon;
		@{$temp_exon_ref} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$temp_exon_ref};
		$temp_start = 0;
	}
	my $temp_last = $temp_exon_ref->[$temp_start];
	$temp_start++;

	if(@{$temp_last}){
		for my $i ($temp_start .. @{$temp_exon_ref}-1){
			my $temp_cur = $temp_exon_ref->[$i];
			if(!@{$temp_cur}){
				next;
			}

			if(($temp_cur->[0] >= $temp_last->[0]) && ($temp_cur->[0] <= ($temp_last->[1] + 1))){
				if($temp_cur->[1] > $temp_last->[1]){
					$temp_last->[1] = $temp_cur->[1];
				}
				@{$temp_cur} = ();
			}
			else{
				last;
			}
		}
	}
	if($temp_start < @{$temp_exon_ref}){
		merge_exon($temp_exon_ref, $temp_start);
	}

	return sort {$a->[0] <=> $b->[0]} map {@{$_} ? $_ : () } @{$temp_exon_ref};
}


#Print reference intron to .gtf
sub print_intron{
	my ($temp_ref_gene_type_ref, $temp_transcript_ref, $temp_intron_gtf) = @_;
	open(my $temp_fh, "> $temp_intron_gtf") or die "Cannot create $temp_intron_gtf!\n";
	my %temp_gene_exon = ();
	foreach my $temp_transcript_id (keys %$temp_transcript_ref){
		unless(exists $temp_ref_gene_type_ref->{$temp_transcript_ref->{$temp_transcript_id}{'gene_type'}}){
			next;
		}
		my $temp_gene_id = $temp_transcript_ref->{$temp_transcript_id}{'gene_id'};
		unless(exists $temp_gene_exon{$temp_gene_id}){
			%{$temp_gene_exon{$temp_gene_id}} = ();
			$temp_gene_exon{$temp_gene_id}{'chr'} = $temp_transcript_ref->{$temp_transcript_id}{'chr'};
			$temp_gene_exon{$temp_gene_id}{'strand'} = $temp_transcript_ref->{$temp_transcript_id}{'strand'};
			$temp_gene_exon{$temp_gene_id}{'gene_type'} = $temp_transcript_ref->{$temp_transcript_id}{'gene_type'};
			@{$temp_gene_exon{$temp_gene_id}{'exon'}} = ();
		}
		foreach my $temp_exon_ref (@{$temp_transcript_ref->{$temp_transcript_id}{'exon'}}){
			push(@{$temp_gene_exon{$temp_gene_id}{'exon'}}, $temp_exon_ref);
		}
	}
	foreach my $temp_gene_id (keys %temp_gene_exon){
		my @temp_merged_exon = merge_exon($temp_gene_exon{$temp_gene_id}{'exon'});
		my $temp_num_exon = @temp_merged_exon;
		for(my $i = 1; $i < $temp_num_exon; $i++){
			my $temp_left = $temp_merged_exon[$i-1][1] + 1;
			my $temp_right = $temp_merged_exon[$i][0] - 1;
			print $temp_fh "$temp_gene_exon{$temp_gene_id}{'chr'}\tNIAP\tintron\t$temp_left\t$temp_right\t.\t$temp_gene_exon{$temp_gene_id}{'strand'}\t.\tgene_id \"$temp_gene_id\"\; gene_type \"$temp_gene_exon{$temp_gene_id}{'gene_type'}\"\;\n";
		}
	}
	close $temp_fh;

	return 1;
}

#Classify novel loci
sub classify_novel{
	my ($temp_monoexonic_gene_type_ref, $temp_ref_gene_type_ref, $temp_max_bd, $temp_ref_transcript_ref, $temp_query_transcript_ref, $temp_ref_gene_ref, $temp_query_gene_ref, $temp_ref_gtf, $temp_query_gtf, $temp_ref_gene_gtf, $temp_ref_gene_tss_gtf, $temp_ref_transcript_gtf, $temp_ref_transcript_tss_gtf, $temp_ref_intron_gtf, $temp_query_gene_gtf, $temp_query_gene_tss_gtf, $temp_query_transcript_gtf, $temp_query_transcript_tss_gtf, $temp_query_classified_gtf, $temp_btp) = @_;
	#Create gene annotation for the reference
	print_gene_tss(
		$temp_monoexonic_gene_type_ref, 
		$temp_ref_gene_type_ref, 
		$temp_ref_gene_ref, 
		$temp_ref_gene_gtf, 
		$temp_ref_gene_tss_gtf
	);
	#Create transcript annotation for the reference
	print_transcript_tss(
		$temp_monoexonic_gene_type_ref, 
		$temp_ref_gene_type_ref, 
		$temp_ref_transcript_ref, 
		$temp_ref_transcript_gtf, 
		$temp_ref_transcript_tss_gtf
	);
	#Create intron annotation for reference
	print_intron(
		$temp_ref_gene_type_ref, 
		$temp_ref_transcript_ref, 
		$temp_ref_intron_gtf
	);
	#Create gene annotation for the query
	print_gene_tss(
		"", 
		"", 
		$temp_query_gene_ref, 
		$temp_query_gene_gtf, 
		$temp_query_gene_tss_gtf
	);
	#Create transcript annotation for reference
	print_transcript_tss(
		"", 
		"", 
		$temp_query_transcript_ref, 
		$temp_query_transcript_gtf, 
		$temp_query_transcript_tss_gtf
	);
	#Classify unannotated genes in the query
	my %temp_class = ();

	#print "Analyzing for nc_host_gene...\n";
	open(my $temp_fh, "$temp_btp intersect -wao -s -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_gene_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		#print "$temp_line\n";
		my @temp_line = split('\t', $temp_line);
		if($temp_line[18] == 0){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_monoexonic_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		#print "$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}\n$temp_query_gene_id\t$temp_ref_gene_id\n";
		#print "$temp_query_gene_id\t$temp_ref_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'ncRNA_host_gene'){
				#print "Already ncRNA_host_gene\n";
				unless(exists $temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id}){
					$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
					#print "Added reference gene\n";
				}
			}
		}
		else{
			#print "Initiated ncRNA_host_gene\nAdded reference gene\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'ncRNA_host_gene';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
	}
	close $temp_fh;

	#print "Analyzing for divergent...\n";
	open($temp_fh, "$temp_btp closest -S -nonamecheck -D a -a $temp_query_transcript_tss_gtf -b $temp_ref_transcript_tss_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[9] eq '.'){
			next;
		}
		elsif(abs($temp_line[18]) > $temp_max_bd){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_transcript_tss_gtf and $temp_ref_transcript_tss_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene when intersecting $temp_query_transcript_tss_gtf and $temp_ref_transcript_tss_gtf ID at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_ref_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		#print "$temp_query_gene_id\t$temp_ref_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'divergent'){
				#print "Already divergent\n";
				unless(exists $temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id}){
					$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
					#print "Added reference gene\n";
				}
			}
		}
		else{
			#print "Initiated divergent\nAdded reference gene\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'divergent';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
		
	}
	close $temp_fh;

	#print "Analyzing for intronic...\n";
	open($temp_fh, "$temp_btp intersect -s -wao -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_intron_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if(($temp_line[3] < $temp_line[12]) || ($temp_line[4] > $temp_line[13])){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_gene_gtf and $temp_ref_intron_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene ID when intersecting $temp_query_gene_gtf and $temp_ref_intron_gtf at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_ref_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		#print "$temp_query_gene_id\t$temp_ref_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existsing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'intronic'){
				#print "Already intronic\n";
				unless(exists $temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id}){
					$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
					#print "Added reference gene\n";
				}
			}
		}
		else{
			#print "Initiated intronic\nAdded reference gene\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'intronic';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
	}

	#print "Analyzing for overlapping...\n";
	open($temp_fh, "$temp_btp intersect -wao -s -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_gene_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[18] == 0){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_ref_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		#print "$temp_query_gene_id\t$temp_ref_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existsing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'overlapping'){
				#print "Already overlapping\n";
				unless(exists $temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id}){
					$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
					#print "Added reference gene\n";
				}
			}
		}
		else{
			#print "Initiated overlapping\nAdded reference gene\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'overlapping';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
	}
	close $temp_fh;

	#print "Analyzing for antisense...\n";
	open($temp_fh, "$temp_btp intersect -wao -S -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_gene_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[18] == 0){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		unless(exists $temp_ref_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		#print "$temp_query_gene_id\t$temp_ref_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'antisense'){
				#print "Already antisense\nAdded reference gene\n";
				$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
			}
		}
		else{
			#print "Initiated antisense\nAdded reference gene\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'antisense';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
	}
	close $temp_fh;

	#print "Analyzing for intergenic...\n";
	open($temp_fh, "$temp_btp intersect -s -v -wa -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_gene_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find query gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		#print "$temp_query_gene_id\n";
		if(exists $temp_class{$temp_query_gene_id}){
			#print "Existing\n";
			if($temp_class{$temp_query_gene_id}{'class'} eq 'intronic'){
				#print "Already intronic\nInitiated intergenic\n";
				$temp_class{$temp_query_gene_id}{'class'} = 'intergenic';
				%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
			}
			#elsif($temp_class{$temp_query_gene_id}{'class'} eq 'intergenic'){
				#print "Already intergenic\n";
			#}
		}
		else{
			#print "Initiated intergenic\n";
			%{$temp_class{$temp_query_gene_id}} = ();
			$temp_class{$temp_query_gene_id}{'class'} = 'intergenic';
			%{$temp_class{$temp_query_gene_id}{'ref_gene_id'}} = ();
		}
	}
	close $temp_fh;

	#print "Analyzing for intergenic sub-type...\n";
	print_gene_tss(
		"", 
		"", 
		$temp_ref_gene_ref, 
		$temp_ref_gene_gtf, 
		$temp_ref_gene_tss_gtf
	);
	open($temp_fh, "$temp_btp intersect -wao -nonamecheck -a $temp_query_gene_gtf -b $temp_ref_gene_gtf |");
	while(<$temp_fh>){
		chomp;
		my $temp_line = $_;
		my @temp_line = split('\t', $temp_line);
		if($temp_line[18] == 0){
			next;
		}
		my $temp_query_gene_id;
		if($temp_line[8] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_query_gene_id = $2;
		}
		unless($temp_class{$temp_query_gene_id}{'class'} eq 'intergenic'){
			next;
		}
		my $temp_ref_gene_id;
		if($temp_line[17] =~ /(^|\s)gene_id \"([^\"]+)\"\;/){
			$temp_ref_gene_id = $2;
		}
		else{
			print STDERR "<ERROR> Cannot find reference gene ID when intersecting $temp_query_gene_gtf and $temp_ref_gene_gtf at line:\n$temp_line\n";
			exit;
		}
		if(exists $temp_ref_gene_type_ref->{$temp_ref_gene_ref->{$temp_ref_gene_id}{'gene_type'}}){
			next;
		}
		unless(exists $temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id}){
			$temp_class{$temp_query_gene_id}{'ref_gene_id'}{$temp_ref_gene_id} = 1;
		}
	}
	close $temp_fh;

	#Print
	open($temp_fh, "> $temp_query_classified_gtf");
	foreach my $temp_query_transcript_id (keys %$temp_query_transcript_ref){
		my %temp_gene_5p = ();
		my @temp_ref_gene_id = keys %{$temp_class{$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}}{'ref_gene_id'}};
		foreach my $temp_gene_id (@temp_ref_gene_id){
			unless(exists $temp_gene_5p{$temp_gene_id}){
				if($temp_ref_gene_ref->{$temp_gene_id}{'strand'} eq "-"){
					$temp_gene_5p{$temp_gene_id} = $temp_ref_gene_ref->{$temp_gene_id}{'right'};
				}
				else{
					$temp_gene_5p{$temp_gene_id} = $temp_ref_gene_ref->{$temp_gene_id}{'left'};
				}
			}
		}
		my @temp_sorted_gene_id = ();
		if($temp_query_transcript_ref->{$temp_query_transcript_id}{'strand'} eq "-"){
			@temp_sorted_gene_id = sort {($temp_gene_5p{$b} <=> $temp_gene_5p{$a}) || ($a cmp $b)} keys %temp_gene_5p;
		}
		else{
			@temp_sorted_gene_id = sort {($temp_gene_5p{$a} <=> $temp_gene_5p{$b}) || ($a cmp $b)} keys %temp_gene_5p;
		}
		my $temp_ref_gene_id = join("\|", @temp_sorted_gene_id);
		my @temp_sorted_gene_name = ();
		my @temp_sorted_gene_type = ();
		foreach my $temp_gene_id (@temp_sorted_gene_id){
			if($temp_ref_gene_ref->{$temp_gene_id}{'gene_name'}){
				push(@temp_sorted_gene_name, $temp_ref_gene_ref->{$temp_gene_id}{'gene_name'});
			}
			if($temp_ref_gene_ref->{$temp_gene_id}{'gene_type'}){
				push(@temp_sorted_gene_type, $temp_ref_gene_ref->{$temp_gene_id}{'gene_type'});
			}
		}
		my $temp_ref_gene_name = join("\|", @temp_sorted_gene_name);
		my $temp_ref_gene_type = join("\|", @temp_sorted_gene_type);
		print_exon(
			$temp_fh, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}, 
			$temp_query_transcript_id, $temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_type'}, 
			$temp_class{$temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'}}{'class'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'asm_gene_id'}, 
			"", 
			$temp_ref_gene_id, 
			$temp_ref_gene_name, 
			$temp_ref_gene_type, 
			""
		);
	}
	close $temp_fh;
	
	return 1;
}

sub rearrange_fusion{
	my ($temp_monoexonic_gene_type_ref, $temp_map_ref,$temp_ref_transcript_ref, $temp_ref_gene_ref, $temp_query_transcript_ref, $temp_query_gtf) = @_;
	open(my $temp_fh, "> $temp_query_gtf") or die "Cannot open $temp_query_gtf!\n";
	my %temp_gene_5p = ();
	foreach my $temp_query_transcript_id (keys %$temp_map_ref){
		my $temp_gene_id = $temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'};
		unless(exists $temp_gene_5p{$temp_gene_id}){
			%{$temp_gene_5p{$temp_gene_id}} = ();
		}
		foreach my $temp_ref_transcript_id (keys %{$temp_map_ref->{$temp_query_transcript_id}}){
			my $temp_ref_gene_id = $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_id'};
			my $temp_ref_gene_type= $temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'};
			if((((@{$temp_query_transcript_ref->{$temp_query_transcript_id}{'exon'}} == 1) && (@{$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'exon'}} == 1)) || ($temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'match'} > 0)) && !(exists $temp_monoexonic_gene_type_ref->{$temp_ref_gene_type})){
				my $temp_5p = get_5p(
					$temp_ref_transcript_ref->{$temp_ref_transcript_id}, 
					$temp_query_transcript_ref->{$temp_query_transcript_id}, 
					$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'ref'}, 
					$temp_map_ref->{$temp_query_transcript_id}{$temp_ref_transcript_id}{'query'}
				);
				if(exists $temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}){
					if($temp_ref_gene_ref->{$temp_ref_gene_id}{'strand'} eq '-'){
						if($temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} < $temp_ref_gene_ref->{$temp_ref_gene_id}{'right'}){
							$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} = $temp_ref_gene_ref->{$temp_ref_gene_id}{'right'};
						}
						if($temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'} < $temp_5p){
							$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'} = $temp_5p;
						}
					}
					else{
						if($temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} > $temp_ref_gene_ref->{$temp_ref_gene_id}{'left'}){
							$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} = $temp_ref_gene_ref->{$temp_ref_gene_id}{'left'};
						}
						if($temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'} > $temp_5p){
							$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'} = $temp_5p;
						}
					}
				}
				else{
					%{$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}} = ();
					if($temp_ref_gene_ref->{$temp_ref_gene_id}{'strand'} eq "-"){
						$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} = $temp_ref_gene_ref->{$temp_ref_gene_id}{'right'};
					}
					else{

						$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'} = $temp_ref_gene_ref->{$temp_ref_gene_id}{'left'};
					}
					$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'} = $temp_5p;
				}
				#print "$temp_query_transcript_id\t$temp_ref_transcript_id\t$temp_gene_id\t$temp_ref_gene_id\t$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'gene'}\t$temp_gene_5p{$temp_gene_id}{$temp_ref_gene_id}{'transcript'}\n";
			}
		}
	}
	foreach my $temp_query_transcript_id (keys %$temp_query_transcript_ref){
		my $temp_gene_id = $temp_query_transcript_ref->{$temp_query_transcript_id}{'gene_id'};
		my @temp_sorted_gene_id = ();
		if($temp_query_transcript_ref->{$temp_query_transcript_id}{'strand'} eq "-"){
			@temp_sorted_gene_id = sort {($temp_gene_5p{$temp_gene_id}{$b}{'transcript'} <=> $temp_gene_5p{$temp_gene_id}{$a}{'transcript'}) || ($temp_gene_5p{$temp_gene_id}{$b}{'gene'} <=> $temp_gene_5p{$temp_gene_id}{$a}{'gene'}) || ($a cmp $b)} keys %{$temp_gene_5p{$temp_gene_id}};
		}
		else{
			@temp_sorted_gene_id = sort {($temp_gene_5p{$temp_gene_id}{$a}{'transcript'} <=> $temp_gene_5p{$temp_gene_id}{$b}{'transcript'}) || ($temp_gene_5p{$temp_gene_id}{$a}{'gene'} <=> $temp_gene_5p{$temp_gene_id}{$b}{'gene'}) || ($a cmp $b)} keys %{$temp_gene_5p{$temp_gene_id}};
		}
		my $temp_ref_gene_id = join("\|", @temp_sorted_gene_id);
		#print "Sorted: $temp_ref_gene_id\n";
		my @temp_sorted_gene_name = ();
		my @temp_sorted_gene_type = ();
		foreach my $temp_gene_id (@temp_sorted_gene_id){
			if($temp_ref_gene_ref->{$temp_gene_id}{'gene_name'}){
				push(@temp_sorted_gene_name, $temp_ref_gene_ref->{$temp_gene_id}{'gene_name'});
			}
			if($temp_ref_gene_ref->{$temp_gene_id}{'gene_type'}){
				push(@temp_sorted_gene_type, $temp_ref_gene_ref->{$temp_gene_id}{'gene_type'});
			}
		}
		my $temp_ref_gene_name = join("\|", @temp_sorted_gene_name);
		my $temp_ref_gene_type = join("\|", @temp_sorted_gene_type);
		print_exon(
			$temp_fh, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}, 
			$temp_query_transcript_id, 
			$temp_ref_gene_id, 
			$temp_ref_gene_name, 
			$temp_ref_gene_type, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'status'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'asm_gene_id'}, 
			$temp_query_transcript_ref->{$temp_query_transcript_id}{'ref_transcript_id'}, 
			"", 
			"", 
			"", 
			""
		); 
	}
	close $temp_fh;

	return 1;
}

#Print missing exons in the reference
sub print_missing{
	my ($temp_ref_transcript_ref, $temp_missing) = @_;
	open(my $temp_fh, "> $temp_missing") or die "Cannot create $temp_missing!\n";
	foreach my $temp_ref_transcript_id (keys %$temp_ref_transcript_ref){
		if(($temp_ref_transcript_ref->{$temp_ref_transcript_id}{'hit'} == 0) && ($temp_ref_transcript_ref->{$temp_ref_transcript_id}{'novel_gene'} == 0)){
			print_exon(
				$temp_fh, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}, 
				$temp_ref_transcript_id, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_id'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_name'}, 
				$temp_ref_transcript_ref->{$temp_ref_transcript_id}{'gene_type'}, 
				"missing", 
				"", 
				"", 
				"", 
				"", 
				"", 
				""
			);
		}
	}
	close $temp_fh;
}





