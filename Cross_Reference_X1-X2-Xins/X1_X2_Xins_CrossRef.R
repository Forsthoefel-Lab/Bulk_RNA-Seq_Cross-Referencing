##############################################################################################
##
## Cross reference apob-M(RNAi) or apob-S(RNAi) differentially expressed transcripts
## with X1-, X2-, and Xins-enriched transcripts
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
## import X1/X2/Xins data and generate DE lists
##
##############################################################################################

edgeR_zeng <- read.delim("Zeng_2018_edgeR_DE_expression_X1_X2_Xins.txt", header=TRUE)

##
## subset DE X1 Signature transcripts
##

## X1 vs Xins, fdr05, up or down
X1 <- edgeR_zeng[edgeR_zeng$FDR_X1_vs_Xins < .05,] 
keep <- !is.na(X1$FDR_X1_vs_Xins) # be sure no rows with NA
X1 <- X1[keep,]
# limit by FC (0)
X1 <- X1[X1$logFC_X1_vs_Xins >= 0,] # 4884 up
# then limit also by vs. X2
X1_Sig <- X1[X1$FDR_X1_vs_X2 < .05,] # 3626
X1_Sig <- X1_Sig[X1_Sig$logFC_X1_vs_X2 >= 0,] # 3572 "signature" X1 transcripts

##
## subset DE X2 Signature transcripts
##

## X2_vs_Xins, fdr05, up or down
X2 <- edgeR_zeng[edgeR_zeng$FDR_X2_vs_Xins < .05,] 
keep <- !is.na(X2$FDR_X2_vs_Xins) # be sure no rows with NA
X2 <- X2[keep,] 
# limit by FC (0)
X2 <- X2[X2$logFC_X2_vs_Xins >= 0,] # 3923 up
# further limit X2 to just "progeny" vs. X1:
X2_Sig <- X2[X2$FDR_X2_vs_X1 < .05,] # 2861
X2_Sig <- X2_Sig[X2_Sig$logFC_X2_vs_X1 >= 0,] #819 "signature" X2 transcripts

##
## subset DE Xins Signature transcripts
##

## Xins, genes that are logFC<0 (and fdr .05) relative to X1 **AND** X2
# first X1
Xins1 <- edgeR_zeng[edgeR_zeng$FDR_X1_vs_Xins < .05,]
keep <- !is.na(Xins1$FDR_X1_vs_Xins) # be sure no rows with NA
Xins1 <- Xins1[keep,] # 11285  
Xins1 <- Xins1[Xins1$logFC_X1_vs_Xins <= 0,] # 6401  
# then X2 -- start with Xins1
Xins2 <- Xins1[Xins1$FDR_X2_vs_Xins < .05,] 
keep <- !is.na(Xins2$FDR_X2_vs_Xins) # be sure no rows with NA
Xins2 <- Xins2[keep,] # 4425  
Xins2 <- Xins2[Xins2$logFC_X2_vs_Xins <= 0,] # 3959
# now be sure all unique
Xins_Sig <- unique(Xins2) # 3959 "signature" Xins transcripts


##############################################################################################
##
## merge apob-M (mild) and X1/X2/Xins subsets and calculate overlap, then plot
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

# X1_Sig vs apobM-UP
X1_calc <- merge(apobM_vs_egfp_fdr05_UP, X1_Sig)
zeng_vs_apob[1,1:6] <- c("X1_Sig", "apobM-UP", as.numeric(NROW(X1_calc)), as.numeric(NROW(X1_Sig)), as.numeric(NROW(X1_calc)/NROW(X1_Sig)), mean(X1_calc$logFC_apobM_vs_egfp))
# X1_Sig vs apobM-DOWN
X1_calc <- merge(apobM_vs_egfp_fdr05_DOWN, X1_Sig)
zeng_vs_apob[2,1:6] <- c("X1_Sig", "apobM-DOWN", as.numeric(NROW(X1_calc)), as.numeric(NROW(X1_Sig)), as.numeric(NROW(X1_calc)/NROW(X1_Sig)), mean(X1_calc$logFC_apobM_vs_egfp))
# X2_Sig vs apobM-UP
X2_calc <- merge(apobM_vs_egfp_fdr05_UP, X2_Sig)
zeng_vs_apob[3,1:6] <- c("X2_Sig", "apobM-UP", as.numeric(NROW(X2_calc)), as.numeric(NROW(X2_Sig)), as.numeric(NROW(X2_calc)/NROW(X2_Sig)), mean(X2_calc$logFC_apobM_vs_egfp))
# X2_Sig vs apobM-DOWN
X2_calc <- merge(apobM_vs_egfp_fdr05_DOWN, X2_Sig)
zeng_vs_apob[4,1:6] <- c("X2_Sig", "apobM-DOWN", as.numeric(NROW(X2_calc)), as.numeric(NROW(X2_Sig)), as.numeric(NROW(X2_calc)/NROW(X2_Sig)), mean(X2_calc$logFC_apobM_vs_egfp))
# Xins_Sig vs apobM-UP
Xins_calc <- merge(apobM_vs_egfp_fdr05_UP, Xins_Sig)
zeng_vs_apob[5,1:6] <- c("Xins_Sig", "apobM-UP", as.numeric(NROW(Xins_calc)), as.numeric(NROW(Xins_Sig)), as.numeric(NROW(Xins_calc)/NROW(Xins_Sig)), mean(Xins_calc$logFC_apobM_vs_egfp))
# Xins_Sig vs apobM-DOWN
Xins_calc <- merge(apobM_vs_egfp_fdr05_DOWN, Xins_Sig)
zeng_vs_apob[6,1:6] <- c("Xins_Sig", "apobM-DOWN", as.numeric(NROW(Xins_calc)), as.numeric(NROW(Xins_Sig)), as.numeric(NROW(Xins_calc)/NROW(Xins_Sig)), mean(Xins_calc$logFC_apobM_vs_egfp))

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
write.table(zeng_vs_apob, "20220521_X1-X2-Xins_vs_apobM.xls", row.names=FALSE, sep="\t")

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
  ylim(-0.6, 1.0) +
  labs(title="Mean logFC in apoB(RNAi)-mild", x="gate", y = "mean logFC") 
p + scale_fill_manual(values=c('salmon','lightblue')) 


##############################################################################################
##
## merge apob-S (severe) and X1/X2/Xins subsets and calculate overlap, then plot
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

# X1_Sig vs apobS-UP
X1_calc <- merge(apobS_vs_egfp_fdr05_UP, X1_Sig)
zeng_vs_apob[1,1:6] <- c("X1_Sig", "apobS-UP", as.numeric(NROW(X1_calc)), as.numeric(NROW(X1_Sig)), as.numeric(NROW(X1_calc)/NROW(X1_Sig)), mean(X1_calc$logFC_apobS_vs_egfp))
# X1_Sig vs apobS-DOWN
X1_calc <- merge(apobS_vs_egfp_fdr05_DOWN, X1_Sig)
zeng_vs_apob[2,1:6] <- c("X1_Sig", "apobS-DOWN", as.numeric(NROW(X1_calc)), as.numeric(NROW(X1_Sig)), as.numeric(NROW(X1_calc)/NROW(X1_Sig)), mean(X1_calc$logFC_apobS_vs_egfp))
# X2_Sig vs apobS-UP
X2_calc <- merge(apobS_vs_egfp_fdr05_UP, X2_Sig)
zeng_vs_apob[3,1:6] <- c("X2_Sig", "apobS-UP", as.numeric(NROW(X2_calc)), as.numeric(NROW(X2_Sig)), as.numeric(NROW(X2_calc)/NROW(X2_Sig)), mean(X2_calc$logFC_apobS_vs_egfp))
# X2_Sig vs apobS-DOWN
X2_calc <- merge(apobS_vs_egfp_fdr05_DOWN, X2_Sig)
zeng_vs_apob[4,1:6] <- c("X2_Sig", "apobS-DOWN", as.numeric(NROW(X2_calc)), as.numeric(NROW(X2_Sig)), as.numeric(NROW(X2_calc)/NROW(X2_Sig)), mean(X2_calc$logFC_apobS_vs_egfp))
# Xins_Sig vs apobS-UP
Xins_calc <- merge(apobS_vs_egfp_fdr05_UP, Xins_Sig)
zeng_vs_apob[5,1:6] <- c("Xins_Sig", "apobS-UP", as.numeric(NROW(Xins_calc)), as.numeric(NROW(Xins_Sig)), as.numeric(NROW(Xins_calc)/NROW(Xins_Sig)), mean(Xins_calc$logFC_apobS_vs_egfp))
# Xins_Sig vs apobS-DOWN
Xins_calc <- merge(apobS_vs_egfp_fdr05_DOWN, Xins_Sig)
zeng_vs_apob[6,1:6] <- c("Xins_Sig", "apobS-DOWN", as.numeric(NROW(Xins_calc)), as.numeric(NROW(Xins_Sig)), as.numeric(NROW(Xins_calc)/NROW(Xins_Sig)), mean(Xins_calc$logFC_apobS_vs_egfp))

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
write.table(zeng_vs_apob, "20210522_X1-X2-Xins_vs_apobS.xls", row.names=FALSE, sep="\t")

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
  ylim(-0.6, 1.0) +
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

