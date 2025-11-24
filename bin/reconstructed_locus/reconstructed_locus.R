#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

# Read input fusion gene ratios
ratio=read.table(paste(args[1], "_fusion_gene_ratio.txt", sep=""), header=TRUE, sep="\t", check.names=F)

# Read input trackgroups file
trackgroups = read.table(args[2], header=FALSE, sep="\t")
tglist      = c( "all", unique( trackgroups[,2]) )

# Add trackgroup columns to ratio file for further thresholding
# In addition to the separate trackgroups we also add an 'all' group
for ( tg in tglist ){
   if ( tg == "all" ){
      selected = trackgroups[,1]
   } else {
      selected = trackgroups[ trackgroups[,2] == tg, 1]
   }
   ratio[[ paste0( "num_passed_sample_", tg) ]]   = rowSums( ratio[, selected, drop=FALSE ] >   0.5 )
   ratio[[ paste0( "num_failed_sample_",  tg) ]]  = rowSums((ratio[, selected, drop=FALSE] >   0   ) & (ratio[, selected]<=0.5) )
   ratio[[ paste0( "num_missing_sample_", tg) ]]  = rowSums( ratio[, selected, drop=FALSE] ==  0   )
   ratio[[ paste0( "num_empty_sample_",   tg) ]]  = rowSums( ratio[, selected, drop=FALSE] == -1   )
}

# Write ratio outputs for each trackgroup with a threshold applied
for ( tg in tglist ){
   if ( tg == "all" ){
      selected = trackgroups[,1]
   } else {
      selected = trackgroups[ trackgroups[,2] == tg, 1]
   }
   threshold     = ifelse( length(selected)<2, 1, floor(length(selected)/2) )
   row_selection = (ratio$cluster_size==1)  &  (ratio$monoexonic==0)  &  (ratio$num_pc<2)  &  (ratio[[ paste0( "num_passed_sample_", tg) ]]>threshold)
   #if ( any(row_selection) ){
      write.table(ratio$gene_id[ row_selection ], 
                  file      = paste0(args[1], "_reclocus_", tg, "_gene_id.txt"), 
                  sep       = "\n", 
                  col.names = FALSE, 
                  row.names = FALSE, 
                  quote     = FALSE)
   #}
}
