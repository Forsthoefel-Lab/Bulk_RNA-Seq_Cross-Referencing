##############################################################################################
##
## Calculate statistical significance of overlap 
## between apob-M(RNAi) or apob-S(RNAi) differentially expressed transcripts
## and transcripts enriched in PIWI-HIGH, PIWI-LO, and PIWI-NEG subpopulations
##
## using R package "GeneOverlap"
##
## apob expression data from Wong et al., Nat. Comm., 2022
##
## X1/X2/Xins data from Zeng et al., Cell, 2018
## aligned to dd_Smed_v6 (Brandl et al., Nuc. Acids Res., 2016)
## then analyzed in edgeR (Wong et al., Nat. Comm., 2022)
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

# install
install.packages("VennDiagram")

# load
library(VennDiagram)


##############################################################################################
##
## import apob RNAi data and Zeng X1/X2/Xins and calculate transcriptome size
##
##############################################################################################

## for stringency, set transcriptome size as the total number of transcripts detected in 
## both the apoB data set and the X1/X2/Xins data set

# read in apob data and calculate total number of transcripts detected
edgeR_apob <- edgeR_apob <- read.delim("apob_RNAi_edgeR_DE_expression.txt", header=TRUE)
keep <- !is.na(edgeR_apob$logFC_apobM_vs_egfp)
apob_detected <- edgeR_apob[keep,] 
NROW(apob_detected) # 18436

# read in Zeng X1/X2/Xins data and calculate total number of trans
edgeR_zeng <- read.delim("Zeng_2018_edgeR_DE_expression_PIWI_HI_LO_NEG.txt", header=TRUE)
keep <- !is.na(edgeR_zeng$logFC_piwiHI_vs_piwiNEG)
zeng_PIWI_detected <- edgeR_zeng[keep,] 
NROW(zeng_PIWI_detected) # 14499

# now use Venn Diagram to determine overlap (transcriptome size)
my_venn <- list(apob_detected=apob_detected$geneID, zeng_PIWI_gates_detected=zeng_PIWI_detected$geneID)

# simple Venn
grid.newpage()
venn.plot <- venn.diagram(my_venn, filename = NULL)
grid.draw(venn.plot)

# therefore for apob vs Zeng PIWI gates, transcriptome size = 14012:
ts.zeng <- 14012


##############################################################################################
##
## import PIWI-HI, PIWI-LO, PIWI-NEG and apob-M, apob-S lists
##
##############################################################################################

# Zeng PIWI gates

# make blank list
Zeng_PIWI_gates_LIST <- list()
# now fill first element with pHI:
Zeng_PIWI_gates_LIST[1] <- read.delim("20210521_PIWI_HI_Sig_list.txt", header=FALSE, as.is = TRUE)
# then pLO
Zeng_PIWI_gates_LIST[2] <- read.delim("20210521_PIWI_LO_Sig_list.txt", header=FALSE, as.is = TRUE)
# then pNEG
Zeng_PIWI_gates_LIST[3] <- read.delim("20210521_PIWI_NEG_Sig_list.txt", header=FALSE, as.is = TRUE)
# and rename
names(Zeng_PIWI_gates_LIST) <- c("PIWI_HI_Sig", "PIWI_LO_Sig", "PIWI_NEG_Sig")
# print (optional, show list)
print(Zeng_PIWI_gates_LIST)
# print (optional, show list)
print(Zeng_PIWI_gates_LIST$PIWI_HI_Sig)
# print (optional, show list)
print(Zeng_PIWI_gates_LIST[1])

## Wong apob dysregulated lists

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
# print (optional, show list)
print(apob_dysreg_LIST$apobM_dysreg)
# print (optional, show list)
print(apob_dysreg_LIST[1])


##############################################################################################
##
## perform statistical test with GeneOverlap 
##
##############################################################################################

## first, look at the length of each list:

sapply(Zeng_PIWI_gates_LIST, length)

# output:
# PIWI_HI_Sig  PIWI_LO_Sig PIWI_NEG_Sig 
# 1376          433         3646 

sapply(apob_dysreg_LIST, length)

# output:
# apobM_dysreg apobS_dysreg 
# 1981         4507 

## run GeneOverlap, generate plots, output data

# reminder: transcriptome size 
ts.zeng <- 14012 # from above

# run GeneOverlap test
gom.obj <- newGOM(apob_dysreg_LIST, Zeng_PIWI_gates_LIST, 
                  ts.zeng)

# draw heat map (p values in matrix, colorized by Odds Ratio)
drawHeatmap(gom.obj, ncolused=5, grid.col="Oranges", note.col="black")

# draw heat map (p values in matrix, colorized by Jaccard Index)
drawHeatmap(gom.obj, ncolused=5, what="Jaccard", grid.col="Blues", note.col="black")

# look at matrix
getMatrix(gom.obj, name="pval")

# write pval maatrix to file
df <- as.data.frame(getMatrix(gom.obj, name="pval"))
df <- cbind(Zeng_PIWI_gate = rownames(df), df)
rownames(df) <- NULL
write.table(df, "20210527_GeneOverlap_Matrix_apob_vs_ZENG_PIWI_GATES.xls", sep="\t", row.names=FALSE)


##############################################################################################
##
## sessionInfo()
##
##############################################################################################

sessionInfo()

R version 4.1.1 (2021-08-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] VennDiagram_1.7.3   futile.logger_1.4.3 GeneOverlap_1.30.0  ggplot2_3.3.6      

loaded via a namespace (and not attached):
  [1] magrittr_2.0.3       munsell_0.5.0        colorspace_2.0-3     R6_2.5.1             rlang_1.0.2          fansi_1.0.3          caTools_1.18.2       tools_4.1.1         
[9] gtable_0.3.0         KernSmooth_2.23-20   utf8_1.2.2           cli_3.3.0            lambda.r_1.2.4       withr_2.5.0          ellipsis_0.3.2       gtools_3.9.2.1      
[17] digest_0.6.29        tibble_3.1.7         lifecycle_1.0.1      crayon_1.5.1         formatR_1.12         BiocManager_1.30.18  RColorBrewer_1.1-3   farver_2.1.0        
[25] futile.options_1.0.1 bitops_1.0-7         vctrs_0.4.1          glue_1.6.2           labeling_0.4.2       compiler_4.1.1       pillar_1.7.0         gplots_3.1.3        
[33] scales_1.2.0         pkgconfig_2.0.3    


