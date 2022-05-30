##############################################################################################
##
## Calculate statistical significance of overlap 
## between apob-M(RNAi) or apob-S(RNAi) differentially expressed transcripts
## and transcripts enriched in lineage-specific subclusters/states in Fincher et al., 2018
##
## using R package "GeneOverlap"
##
## apob expression data from Wong et al., Nat. Comm., 2022
##
## subcluster/state data from Fincher et al., Science, 2018
## download tables at PubMed Central (see below)
##
##############################################################################################


##############################################################################################
##
## load packages
##
##############################################################################################

# install 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GeneOverlap")

# load
library(GeneOverlap)


##############################################################################################
##
## download and pre-process single cell subcluster transcript assignments
##
##############################################################################################

## download Table S2 from Fincher et al., Science, 2018

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6563842/bin/NIHMS1023600-supplement-Table_S2.xlsx

## for each subcluster worksheet,
## delete the top row so that the first row is the header with column names
## save individual subcluster worksheets as tab-delimited .txt files

## then read in the text table

S2_table <- read.delim("TableS2_epidermis.txt", header=TRUE)

## other options:
S2_table <- read.delim("TableS2_cathepsin.txt", header=TRUE)
S2_table <- read.delim("TableS2_intestine.txt", header=TRUE)
S2_table <- read.delim("TableS2_muscle.txt", header=TRUE)
S2_table <- read.delim("TableS2_neural.txt", header=TRUE)
S2_table <- read.delim("TableS2_parenchymal.txt", header=TRUE)
S2_table <- read.delim("TableS2_pharynx.txt", header=TRUE)
S2_table <- read.delim("TableS2_protonephridia.txt", header=TRUE)

## for example, for epidermis

head(S2_table)

# Contig         Best.blast.hit  Accession e.value Organism log.fold.enrichment cluster.enriched.in p.value Adjusted.p.value   AUC Power
# 1  dd_Smed_v4_478_0_1                   <NA>       <NA>      NA     <NA>           0.7623489                   0       0                0 0.746 0.492
# 2  dd_Smed_v4_332_0_1 clone NB.21.11E PROG-1 JX122762.1   0e+00     Smed           0.7609910                   0       0                0 0.771 0.542
# 3  dd_Smed_v4_213_0_1                   <NA>       <NA>      NA     <NA>           0.6634588                   0       0                0 0.751 0.502
# 4 dd_Smed_v4_6912_0_1 nucleobindin 1 (NUCB1) uc002pla.3   2e-30    Human           0.6305754                   0       0                0 0.716 0.432
# 5  dd_Smed_v4_363_0_1                   <NA>       <NA>      NA     <NA>           0.6173929                   0       0                0 0.749 0.498
# 6   dd_Smed_v4_61_0_1                   <NA>       <NA>      NA     <NA>           0.5708968                   0       0                0 0.746 0.492
# etc....

## then edit to enable cross-referencing

# change first column
colnames(S2_table)[1] <- "geneID"
# change geneID to v6
S2_table[1] <- lapply(S2_table[1], gsub, pattern = "_v4_", replacement = "_v6_", fixed = TRUE) #works
# check changes 
head(S2_table)


##############################################################################################
##
## import single cell subcluster state assignments
##
## (order for each subcluster was determined empirically by assessing piwi-1 mRNA expression
## in each subcluster in Fincher et al., Science, 2018, and by cross-referencing with bulk 
## RNA-Seq data from whole planarians 24 hours after lethal irradiation in Zeng et al., 
## Cell, 2018)
##
## note that each subcluster (epidermal, intestine, etc.) is likely to encompass multiple
## cell types/sublineages as well as individual states
##
##############################################################################################

# read in state table (see below for additional options)
state <- read.delim("TableS2_EPIDERMAL_STATE.txt", header = TRUE) 

# other options:
state <- read.delim("TableS2_CATHEPSIN_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_INTESTINE_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_MUSCLE_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_NEURAL_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_PARENCHYMAL_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_PHARYNX_STATE.txt", header = TRUE) 
state <- read.delim("TableS2_PROTONEPHRIDIA_STATE.txt", header = TRUE) 


##############################################################################################
##
## generate lists of transcripts in each state
##
##############################################################################################

# generate list
subcluster_LIST <- list() 
# then run loop
for(i in 1:NROW(state)) {
  temp <- S2_table[grep(paste("^", state$S2_Subcluster[i], "$", sep=""), S2_table$cluster.enriched.in),]
  # then populate a list:
  subcluster_LIST[i] <- list(as.character(temp$geneID))
  # then name the list
  names(subcluster_LIST)[i] <- state$S2_Subcluster[i]
}
  

##############################################################################################
##
## import lists of transcripts dysregulated in control vs. apob-M or apob-S 
##
## logFC > 0 and FDR-adjusted p<.05
##
##############################################################################################

# make blank list
apob_dysreg_LIST <- list()
# now fill first element with pHI:
apob_dysreg_LIST[1] <- read.delim("20210521_apobM_UP_and_DOWN_list.txt", header=FALSE, as.is = TRUE)
# then pLO
apob_dysreg_LIST[2] <- read.delim("20210521_apobS_UP_and_DOWN_list.txt", header=FALSE, as.is = TRUE)
# and rename
names(apob_dysreg_LIST) <- c("apobM_dysreg", "apobS_dysreg")
# print (optional, show list)
print(apob_dysreg_LIST)


##############################################################################################
##
## set transcriptome size 
##
## not shown:
## transcripts with non-zero expression in each subcluster/cell type were determined by 
## normalizing the digital expression matrix in GEP GSE111764 in Seurat as in Fincher et al., 2018
## then extracting individual cell types, and determining the number of transcripts expressed
## in at least 0.5% of all cells within that subcluster/cell type
##
## then overlap was determined between these cell-specific transcripts 
## and those detected in control vs. apob bulk RNA-Seq (Wong et al., 2022)
##
##############################################################################################

# set transcriptome size for subcluster of interest

# transcripts expressed in 0.5% of cells in a cluster and also in control vs. apob bulk RNA-Seq
ts.fincher = 16737 # epidermis 
ts.fincher = 16636 # intestine 
ts.fincher = 15800 # protonephridia 
ts.fincher = 16454 # muscle 
ts.fincher = 16882 # pharynx 
ts.fincher = 16517 # cathepsin 
ts.fincher = 16294 # parenchymal/parapharyngeal 
ts.fincher = 16327 # neural 


##############################################################################################
##
## perform statistical test with GeneOverlap 
##
##############################################################################################

## first, look at the length of each list:

sapply(list_subcluster, length)

# example output for epidermis:
# 2    6    0    1    3    9    7    4    8   12   11   10    5 
# 496   73   53  139  102  261  290  363  409  455  765  140 1484 

sapply(apob_dysreg_LIST, length)

# output:
# apobM_dysreg apobS_dysreg 
# 1981         4507 

## run GeneOverlap, generate plots, output data

# run GeneOverlap test
gom.obj <- newGOM(apob_dysreg_LIST, list_subcluster, 
                  ts.fincher)

# look at matrix
getMatrix(gom.obj, name="pval")

# write pval maatrix to file
df <- as.data.frame(getMatrix(gom.obj, name="pval"))
df <- cbind(Fincher_subcluster = rownames(df), df)
rownames(df) <- NULL
write.table(df, "20210530_GeneOverlap_Matrix_apob_vs_Fincher_subcluster_states.xls", sep="\t", row.names=FALSE)

# draw heat map to file (p values in matrix, colorized by Odds Ratio)
par(mar=c(0.1,0.1,0.1,0.1))
jpeg(filename =  "20210530_GeneOverlap_Matrix_apob_vs_Fincher_subcluster_states_Odds.jpg",
     width = 1600, height = 400, units = 'px') # adjustable
drawHeatmap(gom.obj, ncolused=5, grid.col="Oranges", note.col="black")
graphics.off()

# draw heat map to file (p values in matrix, colorized by Jaccard Index)
par(mar=c(0.1,0.1,0.1,0.1))
jpeg(filename =  "20210530_GeneOverlap_Matrix_apob_vs_Fincher_subcluster_states_Jaccard.jpg",
     width = 1600, height = 400, units = 'px') # adjustable
drawHeatmap(gom.obj, ncolused=5, what="Jaccard", grid.col="Blues", note.col="black")
graphics.off()



##############################################################################################
##
## sessionInfo()
##
##############################################################################################

sessionInfo()

R version 4.1.1 (2021-08-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] GeneOverlap_1.30.0 ggplot2_3.3.6     

loaded via a namespace (and not attached):
  [1] magrittr_2.0.3      munsell_0.5.0       colorspace_2.0-3    R6_2.5.1            rlang_1.0.2         fansi_1.0.3         caTools_1.18.2      tools_4.1.1         grid_4.1.1          gtable_0.3.0       
[11] KernSmooth_2.23-20  utf8_1.2.2          cli_3.3.0           withr_2.5.0         ellipsis_0.3.2      gtools_3.9.2.1      digest_0.6.29       tibble_3.1.7        lifecycle_1.0.1     crayon_1.5.1       
[21] BiocManager_1.30.18 RColorBrewer_1.1-3  farver_2.1.0        bitops_1.0-7        vctrs_0.4.1         glue_1.6.2          labeling_0.4.2      compiler_4.1.1      pillar_1.7.0        gplots_3.1.3       
[31] scales_1.2.0        pkgconfig_2.0.3    

