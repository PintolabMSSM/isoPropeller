This is a pipeline for isoform identification and annotation using using GENCODE v41 (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz).

## Isoform identification

![NIAP_collapse](https://github.com/PintolabMSSM/isoseq_map_analysis/assets/34832128/7d6f84ca-8968-43a2-a7ad-c74be0cfc77d)

The `isoPropeller` command is used to collapse reads in a .bam file into isoforms in a .gtf file. Reads that are likely intrapriming are filtered out. Monoexonic reads without 5' support are also filtered out. Multiexonic reads with the exact set of exon-exon junction(s) and being mapped in the same direction will be collapsed into one transcript. On the other hand, overlapping monoexonic transcripts will be merged to generate one transcript. In both cases, the longest 5' and 3' ends, among the transcripts to be merged, will be taken as the boundaries of the transcript. The distribution of all 5' and 3' ends will also be recorded.

```
perl isoPropeller -i <prefix>.bam -e -p <prefix> -o <prefix> -g <genome .fa> -f /sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed -t <number of threads>
grep $'\t'"transcript"$'\t' <prefix>.gtf | grep "^chr" | grep -v "^chrM" | grep -v "depth \"1\"" | cut -d'"' -f4 > temp_transcript_id.txt
perl select_gtf_by_attribute_list.pl <prefix>.gtf <prefix>_1more.gtf temp_transcript_id.txt transcript_id
rm temp_transcript_id.txt
```

IsoPropeller, by default, will filter out reads predicted as intrapriming, which can be adjusted using the `-w`, `-r` and `-v` options. Alternatively, considering the accuracy of 3' end in long-read data, it is reasonable to keep an intrapriming read if its 3' end is close to any annotated 3' ends. IsoPropeller provides an option, `-l`, to specify a bed file including regions that are close to any annotated 3' ends. For example, these regions can include any annotated 3' terminal exon plus 500 nucleotides downstream.

```
perl prepare_intrapriming_rescue_bed.pl -i gencode.v41.annotation.gtf -o gencode.v41.annotation_500bp -d 500
perl isoPropeller -i <prefix>.bam -e -p <prefix> -o <prefix> -g <genome .fa> -f /sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed -l gencode.v41.annotation_500bp_terminal.bed -t <number of threads>
```

## Merging multiple samples

![NIAP_merge](https://user-images.githubusercontent.com/34832128/189032861-50eb555b-ff9b-4b69-9bdf-60dd28472765.png)

The script `isoPropeller_merge` is used to merge the isofrom .gtf files of different samples by merging isoforms with the same structure and renaming consistently. Multiexonic transcirpts with the exact set of exon-exon junction(s) and being transcribed from the same direction will be merged into one representative transcript. On the other hand, overlapping monoexonic transcripts will be merged to generate a representative transcript. In both cases, the longest 5' and 3' ends, among the transcripts to be merged, will be taken for the representative transcript. A schematic illustration is shown above.

```
for gtf in `ls *_1more.gtf`; do echo $gtf; done > temp_gtf_list.txt
perl isoPropeller_merge -i temp_gtf_list.txt -o MAP -p MAP -e depth -t <number of threads>
rm temp_gtf_list.txt
```

The input is a plain text file specifying the list of .gtf file(s) to be merged, one for each row. The full path to the file is needed.

A file nameed `MAP.gtf` of the representative transcirpts after merging and a plain text file named `MAP_id.txt`, indicating the source of each representative transcripts, will be generated. The `MAP_id.txt` file is consisted of three columns, where the first two columns indicating the gene ID and transcript ID of the representative transcripts, the third one indicating `","` delineated original IDs before merging.

The IDs of a gene and a transcript after merging as `"MAP_<CHR>_[01]_<LOCI_NUM>_<SUB_LOCI_NUM>"` and `"MAP_<CHR>_[01]_<LOCUS_NUM>.<TRANSCRIPT_NUM>"`, respectively. The `"<CHR>"` value indicates the chromosome name. The `"[01]"` value indicates the strand of the locus, or transcript, with `"0"` meaning the `"+"` strand and `"1"` meaning the `"-"` one. The `"<LOCI_NUM>"` indicates the positional order of a locus at a strand. A locus is defined as a region consisted of overlapping transcripts, with each overlapping region containing at least two transcripts. The `"<SUB_LOCI_NUM>"` indicates sub-loci in on locus containing multiexonic transcripts sharing exon boundaries, or monoexonic transcripts having exonic overlapping with exons from other transcripts. The `"<TRANSCRIPT_NUM>"` is a unique identifier of a transcript in the corresponding locus.

## Clustering of unique splice chain based on 5' and 3' ends

The merged multiexonic isoforms can be further split to different isoforms clusters based on their 5' and 3' ends using the `isoPropeller_merge_split` script. The output files `*_end_dist.txt` from the **Isoform identification** step for all samples are needed, as well as the `MAP_id.txt` file from the **Merging multiple samples**.

```
for ed in `ls *_end_dist.txt`; do echo $ed; done > temp_end_dist_list.txt
perl isoPropeller_merge_split -i temp_end_dist_list.txt -o MAP_clustered -g MAP.gtf -d MAP_id.txt -t <number of threads>
rm temp_end_dist_list.txt
```

Four files will be generated, with `MAP_clustered.gtf` showing all isoforms, `MAP_clustered_exp.txt` showing the read count, `MAP_clustered_tss.bed` showing the TSS regions and `MAP_clustered_tts.bed` showing the TTS regions after the clustering step. All gene_id in the .gtf file will remain the same. A string `_\<cluster ID\>` will be appended to the original transcript_id to indicate which cluster of a unique splice chain the isoform belongs to.

## Defining TSS/TTS regions

Since the TSSs/TTSs are variable according to the data, the regions, instead of single-base positions, can be defined.

```
for ed in `ls *_end_dist.txt`; do echo $ed; done > temp_end_dist_list.txt
perl isoPropeller_end_region -i temp_end_dist_list.txt -o MAP -d MAP_id.txt -t <number of threads>
rm temp_end_dist_list.txt
```

Two files will be generated, with `MAP_tss.bed` and `MAP_tts.bed` indicating the TSS and TTS regions for each isoform, respectively. For considering TSS/TTS support for an isoform, it is recommended to use the TSS/TTS region, instead of a single-base position.

## Merging various metrics of 5' and 3' ends

In the **Isoform identification** step, read depths of all TSS-TTS combination in a unique splice chain are recorded in the `*_end_dist.txt` file. Various metrics can be calculated from this file, including the followings.

```
ends_entropy: The entropy of TSS-TTS combination
num_tss/num_tts:The numer of TSS/TTS
tss_mass_center/tts_mass_center: The depth weighted TSS/TTS position
tss_entropy/tts_entropy: The entropy of TSS/TTS position
tss_sd/tts_sd: The standard deviation of TSS/TTS position
tss_max_diff/tts_max_diff: The range of possible TSS/TTS
```

Each of the above metrics can be combined across all samples into a matrix for downstream analyses.

```
for ed in `ls *_end_dist.txt`
do
  prefix=`basename $ed .txt`
  perl isoPropeller_end_dist -i $ed -o ${ed}_parsed.txt
done

for ed in `ls *_end_dist_parsed.txt`; do echo $ed; done > temp_end_dist_list.txt
perl isoPropeller_merge_end_dist -i temp_end_dist_list.txt -o MAP -m MAP_id.txt

for type in `echo ends_entropy tss_entropy tss_mass_center tss_max_diff tss_sd tts_entropy tts_mass_center tts_max_diff tts_sd`
do
  head -1 MAP_${type}.txt | sed 's/\t/\n/g' | cut -d'/' -f8 | tr '\n' '\t' | sed 's/\t$/\n/' > MAP_${type}_renamed.txt
  tail -n+2 MAP_${type}.txt >> MAP_${type}_renamed.txt
done
```

The output files `MAP_${type}_renamed.txt` are the matrix of all the metrics with proper headers.

## Reference-based annotation

![NIAP_annotate_v1 11](https://github.com/PintolabMSSM/isoseq_map_analysis/assets/34832128/967e8c4d-5c16-42ac-be30-cbbb047eabfb)

In this step, all isoforms in the map will be compared with the reference annotation. `/sc/arion/projects/EPIASD/IsoSeq_ourdata/NIAP_merged/MAP_cleaned_unfiltered_tss.bed` is an optional file specifying the TSS region of each unique isoform for improving the annotation.

The output file `MAP.gtf` will be used in the step. The output file `MAP_gencode.gtf` will contain the annotated information.

```
perl isoPropeller_annotate -q MAP.gtf -r gencode.v41.annotation.gtf -o MAP_gencode.gtf -e /sc/arion/projects/EPIASD/IsoSeq_ourdata/NIAP_merged/MAP_cleaned_unfiltered_tss.bed -t <number of threads>
```

In order to extract the annotation, the block below can be used.

```
for attribute in `echo ref_transcript_id asm_gene_id ref_gene_id gene_type gene_name status`; do echo $attribute; done > temp_attribute_list.txt
perl gtf2summary.pl -i MAP_gencode.gtf -o MAP_gencode -a temp_attribute_list.txt -j merged-intropolis-PEC-GTEX-owndata.SJ.out.tab -r -t <number of threads>
```

