#!/usr/bin/perl

# 28.08.2023 09:01:41 EDT
# Author: Dalila Pinto

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $nFivePrime   = 1000;
my $nThreePrime  = 300;
my $nMinReads    = 2;
my $sOutPrefix   = "output-splicechain-endclust";
my $flVerbose    = 0;
GetOptions("help!"   => \$sHelp,
           "a:i"     => \$nFivePrime,
           "b:i"     => \$nThreePrime,
           "m:i"     => \$nMinReads,
           "o:s"     => \$sOutPrefix,
           "v!"      => \$flVerbose);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -a <5prime-threshold> -b <3prime-threshold> input1 ... inputN
   
   Take a set of NIAP output files consisting of splice chains with 3' and 5'
   terminal exon lengths and their counts, and then split chains by length.
   Each splicechain file is expected to correspond to a single sample.
   
   Arguments:
    -a <integer>
      Maximum length difference at the 5' end to create a new splicechain 
      cluster. Default: $nFivePrime
    -b <integer>
      Maximum length difference at the 3' end to create a new splicechain
      cluster. Default: $nThreePrime
    -m <integer>
      Minimum number of splice chain reads per sample
    -o <string>
      Output file prefix
    -v
      Verbose output with detailed information on individual end-clusters
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read all input files into a single multi-exon splicechain hash, where the splicechain
# serves as the hash (associative array) key. The 5':3' end length combination
# is set as a secondary key for each splice chain (to record all possible terminal
# exon length combinations), and the sample count is recorded in a tertiary key.
# The overall hash structure is as follows:
#
#     $hSCmulti{splicechain}{5len:3len}{sample} = readcount;
#
my %hSCmulti;
foreach my $sFile (@ARGV){
   
   warn "Reading data from '$sFile'\n";
   
   open IN, $sFile or die "Error: can't open '$sFile': $!\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my ($sChr, $sStrand, $sGeneID, $sTranscriptID, $nDepth, $sEndDist, $sSpliceChain, $sEndLengths) = split /\t/, $_, -1;
      
      # Read multi-exons
      if ($sSpliceChain and $sEndLengths and ($sEndLengths ne "end_length") and ($nDepth >= $nMinReads) ){
         my @asEndLengths = split /,/, $sEndLengths;
         foreach my $sEL (@asEndLengths){
            my ($n5p, $n3p, $nRc) = split /:/, $sEL;
            $hSCmulti{$sSpliceChain}{"$n5p:$n3p"}{$sFile} += $nRc;
         }
      }
   }
   close IN;
}
my $nSpliceChainCount = scalar keys %hSCmulti;
warn "Read $nSpliceChainCount unique splicechains with at least $nMinReads reads per sample\n";


# Remove any trailing underscores from the output prefix and open output files for writing
$sOutPrefix =~ s/_$//;
open OUTCLUST, ">${sOutPrefix}_clust.txt"  or die "Error: can't open '${sOutPrefix}_clust.txt': $!\n";
open OUTGTF,   ">${sOutPrefix}.gtf"        or die "Error: can't open '${sOutPrefix}.gtf': $!\n";
open OUTCOUNT, ">${sOutPrefix}_counts.txt" or die "Error: can't open '${sOutPrefix}_counts.txt': $!\n";

# Write headers for cluster and count files
print OUTCLUST join("\t", '#gene_id', 'transcript_id', '5prime_exon_length', '3prime_exon_length', 'splice_chain'), "\n";
print OUTCOUNT join("\t", '#transcript_id', @ARGV), "\n";

###################################
# PROCESS MULTI-EXON SPLICECHAINS #
###################################

# Process each splicechain
warn "Clustering splicechain ends\n";
my %hScClust;
my $nSCcount = 1;
foreach my $sSC (keys %hSCmulti){
   my @asLPsort;
   my @asLPfinal;
   my @asLPraw = keys %{$hSCmulti{$sSC}};
   
   #-------------------------------------------------#
   # Splicechain terminal exon clustering by lengths #
   #-------------------------------------------------#
   
   # First we create an array of arrays of 5' and 3' lengths for sorting.
   # The array has 4 columns, the first two columns store the 5' and 3' length, respectively.
   # The last two columns store the 5' and 3' cluster IDs, respectively. The 5' and 3' cluster
   # IDs are initialized as 1 for all 5' and 3' length pairs (i.e. all pairs in one cluster)
   foreach my $sLP (@asLPraw){
      my ($n5p, $n3p) = split ":", $sLP;
      push @asLPsort, [($n5p, $n3p,1,1)];
   }
   
   # Sort ascending by 5' end lengths and then call 5' end clusters. Each time the
   # the difference in length compared to the previous entry exceeds the 5' length
   # threshold, we increment the cluster number until we processed all records.
   # The new cluster number is updated in @asLPraw and the maximum 5' exon length
   # reached in each cluster is recorded in the %hClust5pMaxLen hash.
   @asLPsort = sort {$a->[0] <=> $b->[0]} (@asLPsort);
   my $nClustID5p = 1;
   my %hClust5pMaxLen;
   $hClust5pMaxLen{$nClustID5p} = $asLPsort[0][0];
   for (my $i=1 ; $i<@asLPsort ; $i++){
      $nClustID5p++ if (  ($asLPsort[$i][0] - $asLPsort[$i-1][0]) >= $nFivePrime );
      $asLPsort[$i][2] = $nClustID5p;
      $hClust5pMaxLen{$nClustID5p} = $asLPsort[$i][0];
   }
   
   # Sort ascending by 3' end lengths and then call 3' end clusters. Each time the
   # the difference in length compared to the previous entry exceeds the 3' length
   # threshold, we increment the cluster number until we processed all records.
   # The new cluster number is updated in @asLPraw and the maximum 3' exon length
   # reached in each cluster is recorded in the %hClust3pMaxLen hash.
   @asLPsort = sort {$a->[1] <=> $b->[1]} (@asLPsort);
   my $nClustID3p = 1;
   my %hClust3pMaxLen;
   $hClust3pMaxLen{$nClustID3p} = $asLPsort[0][1];
   for (my $i=1 ; $i<@asLPsort ; $i++){
      $nClustID3p++ if ( ($asLPsort[$i][1] - $asLPsort[$i-1][1]) >= $nThreePrime );
      $asLPsort[$i][3] = $nClustID3p;
      $hClust3pMaxLen{$nClustID3p} = $asLPsort[$i][1];
   }
   
   # Finally, assign a final cluster number for each distinct combination of 5'
   # and 3' cluster. The @asLPraw hash is sorted by 5' and 3' cluster ID and
   # each time the 5' or 3' cluster ID changes, the final cluster ID is incremented.
   # The final clusters and sample counts per cluster are kept in %hSCclust
   my %hSCclust;
   @asLPsort = sort {$a->[2] <=> $b->[2] or $a->[3] <=> $b->[3]} (@asLPsort);
   my $nLast5pc      = $asLPsort[0][2];
   my $nLast3pc      = $asLPsort[0][3];
   my $nClustIDfinal = 1;
   for (my $i=0 ; $i<@asLPsort ; $i++){
      my ($n5pl, $n3pl, $n5pc, $n3pc) = @{$asLPsort[$i]};
      
      # Increment cluster ID any time the 5' or 3' cluster ID changes
      $nClustIDfinal++ if ($n5pc != $nLast5pc or $n3pc != $nLast3pc);
      ($nLast5pc, $nLast3pc) = ($n5pc, $n3pc);
      
      # Set the 5' and 3' lengths to the max length recorded for each end cluster
      $hSCclust{$nClustIDfinal}{length5p} = $hClust5pMaxLen{$n5pc};
      $hSCclust{$nClustIDfinal}{length3p} = $hClust3pMaxLen{$n3pc};
      
      # Aggregate counts across final clusters for all samples
      foreach my $sFile ( keys %{$hSCmulti{$sSC}{"$n5pl:$n3pl"}} ){
         $hSCclust{$nClustIDfinal}{samples}{$sFile} += $hSCmulti{$sSC}{"$n5pl:$n3pl"}{$sFile};
      }
      
      # Debug all clusters
      print join ("\t", "SC.$nSCcount", $n5pl, $n3pl, $n5pc, $n3pc, $nClustIDfinal), "\n" if ($flVerbose);
   }
   
   #---------------#
   # Write outputs #
   #---------------#
   
   # Write data
   foreach my $nClustID ( sort {$a<=>$b} keys %hSCclust ){
      
      # Create gene and transcript IDs and parse the splicechain into exon chunks
      my $sGeneID         = "SC.$nSCcount";
      my $sTranscriptID   = join('', $sGeneID, '.ec', $nClustID);
      my ($sChr, $sStrand, @anSpliceChain) = split /:/, $sSC;
      my $nFirstExonEnd   = shift @anSpliceChain;
      my $nFirstExonStart = $sStrand eq '+' ? $nFirstExonEnd - $hSCclust{$nClustID}{length5p} + 1 : $nFirstExonEnd - $hSCclust{$nClustID}{length3p} + 1;
      my $nLastExonStart  = pop @anSpliceChain;
      my $nLastExonEnd    = $sStrand eq '+' ? $nLastExonStart + $hSCclust{$nClustID}{length3p} - 1 : $nLastExonStart + $hSCclust{$nClustID}{length5p} - 1;
   
      #----------------------------#
      # Write cluster file entries #
      #----------------------------#
      print OUTCLUST join ("\t", $sGeneID, $sTranscriptID, $hSCclust{$nClustID}{length5p}, $hSCclust{$nClustID}{length3p}, $sSC), "\n";

      #----------------------------#
      # Write GTF file entries     #
      #----------------------------# 
           
      # Write GTF transcript entry
      print OUTGTF join ("\t", $sChr, 'NIAP', 'transcript', $nFirstExonStart, $nLastExonEnd, '.', $sStrand, '.', "gene_id \"$sGeneID\"; transcript_id \"$sTranscriptID\";"), "\n";
      die "Error: transcript end is smaller than start for $sTranscriptID\n" unless ($nFirstExonStart <= $nLastExonEnd);
      
      # Write first GTF exon entry
      print OUTGTF join ("\t", $sChr, 'NIAP', 'exon', $nFirstExonStart, $nFirstExonEnd, '.', $sStrand, '.', "gene_id \"$sGeneID\"; transcript_id \"$sTranscriptID\";"), "\n";
      die "Error: first exon end is smaller than start for $sTranscriptID\n" unless ($nFirstExonStart <= $nFirstExonEnd);
      
      # Write internal GTF exon entries    
      for ( my $i=0 ; $i<@anSpliceChain ; $i+=2 ){
         print OUTGTF join ("\t", $sChr, 'NIAP', 'exon', $anSpliceChain[$i], $anSpliceChain[$i+1], '.', $sStrand, '.', "gene_id \"$sGeneID\"; transcript_id \"$sTranscriptID\";"), "\n";
         die "Error: internal exon end is smaller than start for $sTranscriptID\n" unless ($anSpliceChain[$i] <= $anSpliceChain[$i+1]);
      }
      
      # Write last GTF exon entry
      print OUTGTF join ("\t", $sChr, 'NIAP', 'exon', $nLastExonStart, $nLastExonEnd, '.', $sStrand, '.', "gene_id \"$sGeneID\"; transcript_id \"$sTranscriptID\";"), "\n";
      die "Error: last exon end is smaller than start for $sTranscriptID\n" unless ($nLastExonStart <= $nLastExonEnd);

      #----------------------------#
      # Write count matrix entries #
      #----------------------------# 
      print OUTCOUNT $sTranscriptID;
      foreach my $sFile (@ARGV){
         my $nCount = exists $hSCclust{$nClustID}{samples}{$sFile} ? $hSCclust{$nClustID}{samples}{$sFile} : 0;
         print OUTCOUNT "\t$nCount";
      }
      print OUTCOUNT "\n";
   }
   
   # Increment splicechain counter
   warn "Processed $nSCcount splicechains\n" if $nSCcount % 100000 == 0;
   $nSCcount++;
}
warn "Processed $nSCcount splicechains\n";

# Close the output files
close OUTCLUST;
close OUTGTF;
close OUTCOUNT;
