##############################################################################################
##
## Cross reference apob-M(RNAi) or apob-S(RNAi) differentially expressed transcripts
## with  PIWI-1-HIGH ("HI")-, PIWI-1-LOW ("LOW")-, and 
## PIWI-1-NEGATIVE ("NEG")-enriched transcripts
##
## apob expression data from Wong et al., Nat. Comm., 2022
##
## PIWI-HI/LO/NEG data from Zeng et al., Cell, 2018
## aligned to dd_Smed_v6 (Brandl et al., Nuc. Acids Res., 2016)
## then analyzed in edgeR (Wong et al., Nat. Comm., 2022)
##
##############################################################################################


##############################################################################################
##
## load packages
##
##############################################################################################

# install.packages("ggplot2")

library(ggplot2)

##############################################################################################
##
## import apob RNAi data and generate DE lists
##
##############################################################################################

edgeR_apob <- read.delim("apob_RNAi_edgeR_DE_expression.txt", header=TRUE)

##
## subset DE apob-M UP and DOWN
##

# limit by fdr<.05
apobM_vs_egfp_fdr05 <- edgeR_apob[edgeR_apob$FDR_apobM_vs_egfp<= 0.05,]
keep <- !is.na(apobM_vs_egfp_fdr05$FDR_apobM_vs_egfp) # be sure no rows with NA
apobM_vs_egfp_fdr05 <- apobM_vs_egfp_fdr05[keep,] 
# limit by logFC
apobM_vs_egfp_fdr05_UP <- apobM_vs_egfp_fdr05[apobM_vs_egfp_fdr05$logFC_apobM_vs_egfp >= 0, ] # 842
apobM_vs_egfp_fdr05_DOWN <- apobM_vs_egfp_fdr05[apobM_vs_egfp_fdr05$logFC_apobM_vs_egfp <= 0, ] # 1139

##
## subset DE apob-S UP and DOWN
##

# limit by fdr<.05
apobS_vs_egfp_fdr05 <- edgeR_apob[edgeR_apob$FDR_apobS_vs_egfp<= 0.05,]
keep <- !is.na(apobS_vs_egfp_fdr05$FDR_apobS_vs_egfp) # be sure no rows with NA
apobS_vs_egfp_fdr05 <- apobS_vs_egfp_fdr05[keep,] 
# limit by logFC
apobS_vs_egfp_fdr05_UP <- apobS_vs_egfp_fdr05[apobS_vs_egfp_fdr05$logFC_apobS_vs_egfp >= 0, ] # 1960
apobS_vs_egfp_fdr05_DOWN <- apobS_vs_egfp_fdr05[apobS_vs_egfp_fdr05$logFC_apobS_vs_egfp <= 0, ] # 2547


##############################################################################################
##
## import PIWI-HI/LO/NEG data and generate DE lists
##
##############################################################################################

edgeR_zeng <- read.delim("Zeng_2018_edgeR_DE_expression_PIWI_HI_LO_NEG.txt", header=TRUE)

##
## subset DE PIWI-HI Signature transcripts
##

## PIWI-HI vs PIWI-NEG, fdr05, up or down
pHI <- edgeR_zeng[edgeR_zeng$FDR_piwiHI_vs_piwiNEG < .05,] 
keep <- !is.na(pHI$FDR_piwiHI_vs_piwiNEG) # be sure no rows with NA
pHI <- pHI[keep,] 
# limit by FC (0)
pHI <- pHI[pHI$logFC_piwiHI_vs_piwiNEG >= 0,] # 4623 up
# limit pHI further vs. pLO:
pHI_Sig <- pHI[pHI$FDR_piwiHI_vs_piwiLO < .05,] # 1514
pHI_Sig <- pHI_Sig[pHI_Sig$logFC_piwiHI_vs_piwiLO >= 0,] # 1376 "signature" PIWI-HI transcripts

##
## subset DE PIWI-LO Signature transcripts
##

## PIWI-LO vs PIWI-NEG, fdr05, up or down
pLO <- edgeR_zeng[edgeR_zeng$FDR_piwiLO_vs_piwiNEG < .05,] 
keep <- !is.na(pLO$FDR_piwiLO_vs_piwiNEG) # be sure no rows with NA
pLO <- pLO[keep,]
# limit by FC (0)
pLO <- pLO[pLO$logFC_piwiLO_vs_piwiNEG >= 0,] # 4039 up
# # limit pLO further vs. pHI:
pLO_Sig <- pLO[pLO$FDR_piwiLO_vs_piwiHI < .05,] # 1553
pLO_Sig <- pLO_Sig[pLO_Sig$logFC_piwiLO_vs_piwiHI >= 0,] # 433 "signature" pLO trancripts

##
## subset DE PIWI-NEG Signature transcripts
##

## PIWI-NEG, transcripts that are logFC<0 (and fdr .05) relative to PIWI-HI **AND** PIWI-LO
# first pHI
pNEG1 <- edgeR_zeng[edgeR_zeng$FDR_piwiHI_vs_piwiNEG < .05,]
keep <- !is.na(pNEG1$FDR_piwiHI_vs_piwiNEG) # get rid of NA'sx
pNEG1 <- pNEG1[keep,] # 9874  
pNEG1 <- pNEG1[pNEG1$logFC_piwiHI_vs_piwiNEG <= 0,] # 5251  
# then pLO
pNEG2 <- pNEG1[pNEG1$FDR_piwiLO_vs_piwiNEG < .05,] #3713
keep <- !is.na(pNEG2$FDR_piwiLO_vs_piwiNEG) # get rid of NA's
pNEG2 <- pNEG2[keep,] # 3713  
pNEG2 <- pNEG2[pNEG2$logFC_piwiLO_vs_piwiNEG <= 0,] # 3646  
# now make sure unique
pNEG_Sig <- unique(pNEG2) # 3646 "signature" pNEG transcripts


##############################################################################################
##
## merge apob-M (mild) and PIWI-HI/LO/NEG subsets and calculate overlap, then plot
##
##############################################################################################

##
## create table with calculations for each gate
##

# create empty matrix (will convert to data frame at end)
zeng_vs_apob <- matrix(, 6, ncol=6) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_vs_apob) <- c("gate", "experiment", "transcript_num_dysreg", "transcripts_in_gate", "percent_dysreg", "mean_logFC")

# calculate percent of dysregulated transcripts and mean logFCs for each gate

# pHI_Sig vs. apobM-UP
pHI_calc <- merge(apobM_vs_egfp_fdr05_UP, pHI_Sig)
zeng_vs_apob[1,1:6] <- c("pHI_Sig", "apobM-UP", as.numeric(NROW(pHI_calc)), as.numeric(NROW(pHI_Sig)), as.numeric(NROW(pHI_calc)/NROW(pHI_Sig)), mean(pHI_calc$logFC_apobM_vs_egfp))
# pHI_Sig vs. apobM-DOWN
pHI_calc <- merge(apobM_vs_egfp_fdr05_DOWN, pHI_Sig)
zeng_vs_apob[2,1:6] <- c("pHI_Sig", "apobM-DOWN", as.numeric(NROW(pHI_calc)), as.numeric(NROW(pHI_Sig)), as.numeric(NROW(pHI_calc)/NROW(pHI_Sig)), mean(pHI_calc$logFC_apobM_vs_egfp))
# pLO_Sig vs. apobM-UP
pLO_calc <- merge(apobM_vs_egfp_fdr05_UP, pLO_Sig)
zeng_vs_apob[3,1:6] <- c("pLO_Sig", "apobM-UP", as.numeric(NROW(pLO_calc)), as.numeric(NROW(pLO_Sig)), as.numeric(NROW(pLO_calc)/NROW(pLO_Sig)), mean(pLO_calc$logFC_apobM_vs_egfp))
# pLO_Sig vs. apobM-DOWN
pLO_calc <- merge(apobM_vs_egfp_fdr05_DOWN, pLO_Sig)
zeng_vs_apob[4,1:6] <- c("pLO_Sig", "apobM-DOWN", as.numeric(NROW(pLO_calc)), as.numeric(NROW(pLO_Sig)), as.numeric(NROW(pLO_calc)/NROW(pLO_Sig)), mean(pLO_calc$logFC_apobM_vs_egfp))
# pLO_Sig vs. apobM-UP
pNEG_calc <- merge(apobM_vs_egfp_fdr05_UP, pNEG_Sig)
zeng_vs_apob[5,1:6] <- c("pNEG_Sig", "apobM-UP", as.numeric(NROW(pNEG_calc)), as.numeric(NROW(pNEG_Sig)), as.numeric(NROW(pNEG_calc)/NROW(pNEG_Sig)), mean(pNEG_calc$logFC_apobM_vs_egfp))
# pNEG_Sig vs. apobM-DOWN
pNEG_calc <- merge(apobM_vs_egfp_fdr05_DOWN, pNEG_Sig)
zeng_vs_apob[6,1:6] <- c("pNEG_Sig", "apobM-DOWN", as.numeric(NROW(pNEG_calc)), as.numeric(NROW(pNEG_Sig)), as.numeric(NROW(pNEG_calc)/NROW(pNEG_Sig)), mean(pNEG_calc$logFC_apobM_vs_egfp))

# convert to data frame
zeng_vs_apob <- data.frame(zeng_vs_apob)
# convert columns to plot to numeric (for plotting)
zeng_vs_apob$percent_dysreg <- as.numeric(as.character(zeng_vs_apob$percent_dysreg)) 
zeng_vs_apob$mean_logFC <- as.numeric(as.character(zeng_vs_apob$mean_logFC)) 
# then reorder levels so apobM-UP plots on top
zeng_vs_apob$experiment = factor(zeng_vs_apob$experiment, levels = c("apobM-UP", "apobM-DOWN"))
# look
zeng_vs_apob

# write output if desired
write.table(zeng_vs_apob, "20220521_PIWI-HI-LO-NEG_vs_apobM.xls", row.names=FALSE, sep="\t")

##
## plot
##

# plot pct_dysregulated transcripts - OPTIONAL: ylim
p <- ggplot(data=zeng_vs_apob, aes(x=gate, y=percent_dysreg, fill=experiment)) +
  geom_col(width = 0.75) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  #                                 panel.grid.minor = element_line(color="gray70"), 
  #                                 panel.grid.major = element_line(color="gray70")) +
  ylim(0,0.7) +
  labs(title="Dysregulated transcripts in apoB(RNAi)-mild", x="gate", y = "percent_dysregulated_transcripts") 
p + scale_fill_manual(values=c('salmon','lightblue')) 

# plot mean logFC
p <- ggplot(data=zeng_vs_apob, aes(x=gate, y=mean_logFC, fill=experiment)) +
  geom_col(width = 0.75, position=position_dodge2()) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  ylim(-0.7, 1.0) +
  labs(title="Mean logFC in apoB(RNAi)-mild", x="gate", y = "mean logFC") 
p + scale_fill_manual(values=c('salmon','lightblue')) 


##############################################################################################
##
## merge apob-S (severe) and PIWI-HI-LO-NEG subsets and calculate overlap, then plot
##
##############################################################################################

##
## create table with calculations for each gate
##

# create empty matrix (will convert to data frame at end)
zeng_vs_apob <- matrix(, 6, ncol=6) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_vs_apob) <- c("gate", "experiment", "transcript_num_dysreg", "transcripts_in_gate", "percent_dysreg", "mean_logFC")

# calculate percent of dysregulated transcripts and mean logFCs for each gate

# pHI_Sig vs. apobS-UP
pHI_calc <- merge(apobS_vs_egfp_fdr05_UP, pHI_Sig)
zeng_vs_apob[1,1:6] <- c("pHI_Sig", "apobS-UP", as.numeric(NROW(pHI_calc)), as.numeric(NROW(pHI_Sig)), as.numeric(NROW(pHI_calc)/NROW(pHI_Sig)), mean(pHI_calc$logFC_apobS_vs_egfp))
# pHI_Sig vs. apobS-DOWN
pHI_calc <- merge(apobS_vs_egfp_fdr05_DOWN, pHI_Sig)
zeng_vs_apob[2,1:6] <- c("pHI_Sig", "apobS-DOWN", as.numeric(NROW(pHI_calc)), as.numeric(NROW(pHI_Sig)), as.numeric(NROW(pHI_calc)/NROW(pHI_Sig)), mean(pHI_calc$logFC_apobS_vs_egfp))
# pLO_Sig vs. apobS-UP
pLO_calc <- merge(apobS_vs_egfp_fdr05_UP, pLO_Sig)
zeng_vs_apob[3,1:6] <- c("pLO_Sig", "apobS-UP", as.numeric(NROW(pLO_calc)), as.numeric(NROW(pLO_Sig)), as.numeric(NROW(pLO_calc)/NROW(pLO_Sig)), mean(pLO_calc$logFC_apobS_vs_egfp))
# pLO_Sig vs. apobS-DOWN
pLO_calc <- merge(apobS_vs_egfp_fdr05_DOWN, pLO_Sig)
zeng_vs_apob[4,1:6] <- c("pLO_Sig", "apobS-DOWN", as.numeric(NROW(pLO_calc)), as.numeric(NROW(pLO_Sig)), as.numeric(NROW(pLO_calc)/NROW(pLO_Sig)), mean(pLO_calc$logFC_apobS_vs_egfp))
# pLO_Sig vs. apobS-UP
pNEG_calc <- merge(apobS_vs_egfp_fdr05_UP, pNEG_Sig)
zeng_vs_apob[5,1:6] <- c("pNEG_Sig", "apobS-UP", as.numeric(NROW(pNEG_calc)), as.numeric(NROW(pNEG_Sig)), as.numeric(NROW(pNEG_calc)/NROW(pNEG_Sig)), mean(pNEG_calc$logFC_apobS_vs_egfp))
# pNEG_Sig vs. apobS-DOWN
pNEG_calc <- merge(apobS_vs_egfp_fdr05_DOWN, pNEG_Sig)
zeng_vs_apob[6,1:6] <- c("pNEG_Sig", "apobS-DOWN", as.numeric(NROW(pNEG_calc)), as.numeric(NROW(pNEG_Sig)), as.numeric(NROW(pNEG_calc)/NROW(pNEG_Sig)), mean(pNEG_calc$logFC_apobS_vs_egfp))

# convert to data frame
zeng_vs_apob <- data.frame(zeng_vs_apob)
# convert columns to plot to numeric (for plotting)
zeng_vs_apob$percent_dysreg <- as.numeric(as.character(zeng_vs_apob$percent_dysreg)) 
zeng_vs_apob$mean_logFC <- as.numeric(as.character(zeng_vs_apob$mean_logFC)) 
# then reorder levels so apobS-UP plots on top
zeng_vs_apob$experiment = factor(zeng_vs_apob$experiment, levels = c("apobS-UP", "apobS-DOWN"))
# look
zeng_vs_apob

# write output if desired
write.table(zeng_vs_apob, "20210522_PIWI-HI-LO-NEG_vs_apobS.xls", row.names=FALSE, sep="\t")

##
## plot
##

# plot pct_dysregulated transcripts - OPTIONAL: ylim
p <- ggplot(data=zeng_vs_apob, aes(x=gate, y=percent_dysreg, fill=experiment)) +
  geom_col(width = 0.75) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
  #                                 panel.grid.minor = element_line(color="gray70"), 
  #                                 panel.grid.major = element_line(color="gray70")) +
  ylim(0,0.7) +
  labs(title="Dysregulated transcripts in apoB(RNAi)-severe", x="gate", y = "percent_dysregulated_transcripts") 
p + scale_fill_manual(values=c('salmon','lightblue')) 

# plot mean logFC
p <- ggplot(data=zeng_vs_apob, aes(x=gate, y=mean_logFC, fill=experiment)) +
  geom_col(width = 0.75, position=position_dodge2()) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  ylim(-0.7, 1.0) +
  labs(title="Mean logFC in apoB(RNAi)-severe", x="gate", y = "mean logFC") 
p + scale_fill_manual(values=c('salmon','lightblue')) 


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
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] ggplot2_3.3.6

loaded via a namespace (and not attached):
  [1] fansi_1.0.3      digest_0.6.29    withr_2.5.0      utf8_1.2.2       crayon_1.5.1     grid_4.1.1       R6_2.5.1         lifecycle_1.0.1  gtable_0.3.0    
[10] magrittr_2.0.3   scales_1.2.0     pillar_1.7.0     rlang_1.0.2      cli_3.3.0        farver_2.1.0     vctrs_0.4.1      ellipsis_0.3.2   labeling_0.4.2  
[19] tools_4.1.1      glue_1.6.2       munsell_0.5.0    compiler_4.1.1   pkgconfig_2.0.3  colorspace_2.0-3 tibble_3.1.7    

