##############################################################################################
##
## Cross reference apob-M(RNAi) or apob-S(RNAi) differentially expressed transcripts
## with Fincher single cell subclusters
##
## apob expression data from Wong et al., Nat. Comm., 2022
##
## Fincher subclusters downloaded from open access source (see below)
## Fincher et al., Science, 2018
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
## download and pre-process single cell subcluster transcript assignments
##
##############################################################################################

## download Table S2 from Fincher et al., Science, 2018

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6563842/bin/NIHMS1023600-supplement-Table_S2.xlsx

## for each subcluster worksheet,
## delete the top row so that the first row is the header with column names
## save individual subcluster worksheets as tab-delimited .txt files

## then read in the text table

S2_table<- read.delim("TableS2_epidermis.txt", header=TRUE)

## other options:
S2_table<- read.delim("TableS2_cathepsin.txt", header=TRUE)
S2_table<- read.delim("TableS2_intestine.txt", header=TRUE)
S2_table<- read.delim("TableS2_muscle.txt", header=TRUE)
S2_table<- read.delim("TableS2_neural.txt", header=TRUE)
S2_table<- read.delim("TableS2_parenchymal.txt", header=TRUE)
S2_table<- read.delim("TableS2_pharynx.txt", header=TRUE)
S2_table<- read.delim("TableS2_protonephridia.txt", header=TRUE)

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
## set variables - for subcluster of interest
##
##############################################################################################

# set "subcluster" to table of interest
subcluster = S2_table

# set subcluster of interest: choose one (example, epidermal)
subcluster_name = "S2_epidermal"

# other options:
subcluster_name = "S2_cathepsin"
subcluster_name = "S2_intestine"
subcluster_name = "S2_muscle"
subcluster_name = "S2_neural"
subcluster_name = "S2_parenchymal"
subcluster_name = "S2_pharynx"
subcluster_name = "S2_protonephridia"


##############################################################################################
##
## generate table of dysregulated transcripts per subcluster - apobM (mild)
##
##############################################################################################

##
## Table1 = apobM - UP
##

## apobM - set variables for table generation 

target = apobM_vs_egfp_fdr05_UP
target_name = "apobM_UP"

logFC = "logFC_apobM_vs_egfp"

## MERGE 1 - extract targets from subcluster using merge
S2_sub = merge(target, subcluster, by="geneID") 

# calculate subcluster frequency
table_S2_subclust = as.data.frame(table(S2_sub$cluster.enriched.in, exclude=NULL)) # make table a data frame in the environment
# rename column 1 to reflect subcluster
colnames(table_S2_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
# rename column 2 to reflect target
colnames(table_S2_subclust)[2] = paste("Freq_", target_name, sep="")
# normalize / calculate % of all extracted subclusters... questionable value but leave in for now
table_S2_subclust[,paste("Freq_Pct_subcluster_", target_name, sep="")] = table_S2_subclust$Freq/sum(table_S2_subclust$Freq)

# to normalize to apobM / calculate % of apobM that are dysregulated (use variable "target_name" above to create custom column)
numrow = as.numeric(NROW(as.data.frame(unique(target$geneID))))
table_S2_subclust[,paste("Freq_Pct_", target_name, sep="")] = table_S2_subclust[,2]/numrow # could also name the Frequency with paste but it should always be column2 here

# sort and look
table_S2_subclust <- table_S2_subclust[order(-table_S2_subclust[,4]),]
row.names(table_S2_subclust) <- NULL # resets rownames after sort (have to do last)
# look
table_S2_subclust

## calculate number of transcripts in each subcluster and normalize to that (most appropriate for seeing if a subcluster is affected)

# first create data framae
transcripts_per_subclust = as.data.frame(table(subcluster$cluster.enriched.in, exclude = NULL))
# then can rename column to something more sensible than "Var1" - again use "paste" to create custom colname
colnames(transcripts_per_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
colnames(transcripts_per_subclust)[2] = paste("transcripts_per_subcluster_", target_name, sep="")
# take a look
transcripts_per_subclust

## MERGE 2 - add to table_S2_subclust
transcripts_per_subclust <- merge(table_S2_subclust, transcripts_per_subclust, by.x = colnames(table_S2_subclust)[1], by.y = colnames(transcripts_per_subclust)[1])

# now calculate frequency per subcluster
transcripts_per_subclust[,paste("Freq_Pct_by_subclust_", target_name, sep="")] <- transcripts_per_subclust[,2]/transcripts_per_subclust[,5] # again, could name columns using paste but they should always be in the same position
# take a look: 
transcripts_per_subclust
# can sort if you want: 
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,4]),] # by % of orig
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,6]),] # by % of transcripts in subcluster
row.names(transcripts_per_subclust) <- NULL # resets rownames after sort
# take a look (nothing should be over 1.0): 
transcripts_per_subclust

## now use for loop to calculate logFCs and add to table

# make empty matrix (will convert to data frame at end)
zeng_logFC_by_finch_subclust <- matrix(, nrow=NROW(transcripts_per_subclust), ncol=3) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_logFC_by_finch_subclust) <- c(paste(subcluster_name, "_subcluster", sep=""), paste("transcript_num_", target_name, sep=""), paste("avg_logFC_", target_name, sep=""))

# now run a loop
# this populates the table by grepping out rows, one at a time, based on the entries in the "transcripts_per_subclust" table above
# must use "as.character" to properly print subcluster number since this is a factor in the "transcripts_per_subclust" data frame
# also uses "logFC" variable for column 3, defined above with other variables

for(i in 1:NROW(transcripts_per_subclust)) { # start loop
  S2_sub2 <- S2_sub[grep(paste("^", as.character(transcripts_per_subclust[i,1]), "$", sep=""), S2_sub$cluster.enriched.in),]
  zeng_logFC_by_finch_subclust[i,1:3] <- c(as.character(transcripts_per_subclust[i,1]), nrow(S2_sub2), mean(S2_sub2[,logFC], na.rm = TRUE)) # nice, specify row with i and columns 1:3
}

# convert to data frame
zeng_logFC_by_finch_subclust <- as.data.frame(zeng_logFC_by_finch_subclust)

# look
zeng_logFC_by_finch_subclust 

# merge with counts
Table1 <- merge(transcripts_per_subclust, zeng_logFC_by_finch_subclust) # should not need to specify "by" columns if all above worked - it will use col1, which should be identical

# look, verify that transcript numbers match, and that no Freq_Pct_by_subcluster is >1 (just a QC check)
Table1


##
## Table2 = apobM - DOWN
##

## apobM - set variables for table generation 

target = apobM_vs_egfp_fdr05_DOWN
target_name = "apobM_DOWN"

logFC = "logFC_apobM_vs_egfp"

## MERGE 1 - extract targets from subcluster using merge
S2_sub = merge(target, subcluster, by="geneID") 

# calculate subcluster frequency
table_S2_subclust = as.data.frame(table(S2_sub$cluster.enriched.in, exclude=NULL)) # make table a data frame in the environment
# rename column 1 to reflect subcluster
colnames(table_S2_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
# rename column 2 to reflect target
colnames(table_S2_subclust)[2] = paste("Freq_", target_name, sep="")
# normalize / calculate % of all extracted subclusters... questionable value but leave in for now
table_S2_subclust[,paste("Freq_Pct_subcluster_", target_name, sep="")] = table_S2_subclust$Freq/sum(table_S2_subclust$Freq)

# to normalize to apobM / calculate % of apobM that are dysregulated (use variable "target_name" above to create custom column)
numrow = as.numeric(NROW(as.data.frame(unique(target$geneID))))
table_S2_subclust[,paste("Freq_Pct_", target_name, sep="")] = table_S2_subclust[,2]/numrow # could also name the Frequency with paste but it should always be column2 here

# sort and look
table_S2_subclust <- table_S2_subclust[order(-table_S2_subclust[,4]),]
row.names(table_S2_subclust) <- NULL # resets rownames after sort (have to do last)
# look
table_S2_subclust

## calculate number of transcripts in each subcluster and normalize to that (most appropriate for seeing if a subcluster is affected)

# first create data framae
transcripts_per_subclust = as.data.frame(table(subcluster$cluster.enriched.in, exclude = NULL))
# then can rename column to something more sensible than "Var1" - again use "paste" to create custom colname
colnames(transcripts_per_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
colnames(transcripts_per_subclust)[2] = paste("transcripts_per_subcluster_", target_name, sep="")
# take a look
transcripts_per_subclust

## MERGE 2 - add to table_S2_subclust
transcripts_per_subclust <- merge(table_S2_subclust, transcripts_per_subclust, by.x = colnames(table_S2_subclust)[1], by.y = colnames(transcripts_per_subclust)[1])

# now calculate frequency per subcluster
transcripts_per_subclust[,paste("Freq_Pct_by_subclust_", target_name, sep="")] <- transcripts_per_subclust[,2]/transcripts_per_subclust[,5] # again, could name columns using paste but they should always be in the same position
# take a look: 
transcripts_per_subclust
# can sort if you want: 
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,4]),] # by % of orig
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,6]),] # by % of transcripts in subcluster
row.names(transcripts_per_subclust) <- NULL # resets rownames after sort
# take a look (nothing should be over 1.0): 
transcripts_per_subclust

## now use for loop to calculate logFCs and add to table

# make empty matrix (will convert to data frame at end)
zeng_logFC_by_finch_subclust <- matrix(, nrow=NROW(transcripts_per_subclust), ncol=3) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_logFC_by_finch_subclust) <- c(paste(subcluster_name, "_subcluster", sep=""), paste("transcript_num_", target_name, sep=""), paste("avg_logFC_", target_name, sep=""))

# now run a loop
# this populates the table by grepping out rows, one at a time, based on the entries in the "transcripts_per_subclust" table above
# must use "as.character" to properly print subcluster number since this is a factor in the "transcripts_per_subclust" data frame
# also uses "logFC" variable for column 3, defined above with other variables

for(i in 1:NROW(transcripts_per_subclust)) { # start loop
  S2_sub2 <- S2_sub[grep(paste("^", as.character(transcripts_per_subclust[i,1]), "$", sep=""), S2_sub$cluster.enriched.in),]
  zeng_logFC_by_finch_subclust[i,1:3] <- c(as.character(transcripts_per_subclust[i,1]), nrow(S2_sub2), mean(S2_sub2[,logFC], na.rm = TRUE)) # nice, specify row with i and columns 1:3
}

# convert to data frame
zeng_logFC_by_finch_subclust <- as.data.frame(zeng_logFC_by_finch_subclust)

# look
zeng_logFC_by_finch_subclust 

# merge with counts
Table2 <- merge(transcripts_per_subclust, zeng_logFC_by_finch_subclust) # should not need to specify "by" columns if all above worked - it will use col1, which should be identical

# look, verify that transcript numbers match, and that no Freq_Pct_by_subcluster is >1 (just a QC check)
Table2

##
## merge two tables
##

final_table = merge(Table1, Table2, all.x=TRUE, all.y=TRUE) # need to be TRUE in case some subclusters have NAs


##############################################################################################
##
## merge final_table with subcluster state tables for plotting - apobM (mild)
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

# merge subcluster and state columns
state$Subcluster <- paste(state$S2_Subcluster, " - ", state$State, sep="")
# assign order based on the order of the "state" file - for plotting below
state$order <- rownames(state)
# rename first column to enable merge with final_table
colnames(state)[1] <- paste(subcluster_name, "_subcluster", sep="")

# merge with final_table
final_table = merge(state, final_table, all.x = TRUE, all.y = TRUE) # must be TRUE to include subclusters with NAs

# write table, if desired: apobM-UP/DOWN
write.table(final_table, paste("20210328_apobM-UP-DOWN_extract_from_", subcluster_name, "_BY_SUBCLUSTER.xls", sep=""), row.names=FALSE, sep="\t")


##############################################################################################
##
## plot - apob-M
##
##############################################################################################

# set variables
apob1 = "apobM_UP"
apob2 = "apobM_DOWN"

# OPTION 1: create a "melted" table manually - apobM - ONLY subcluster number
transcript_pct_melt <- data.frame(Subcluster=final_table[,1], Order=as.numeric(final_table[,4]), apob=apob1, Freq_Percent_by_Subcluster=final_table[,9])
numrow = as.numeric(NROW(final_table)) # to calc multiplication factor for Order
transcript_pct_melt2 <- data.frame(Subcluster=final_table[,1], Order=(numrow+as.numeric(final_table[,4])), apob=apob2, Freq_Percent_by_Subcluster=final_table[,16])
transcript_pct_melt <- rbind(transcript_pct_melt,transcript_pct_melt2)
# look
transcript_pct_melt

# OPTION 2: create a "melted" table manually - apobM -subcluster number and states
transcript_pct_melt <- data.frame(Subcluster=paste(final_table[,1], final_table[,2], sep=" - "), Order=as.numeric(final_table[,4]), apob=apob1, Freq_Percent_by_Subcluster=final_table[,9])
numrow = as.numeric(NROW(final_table)) # to calc multiplication factor for Order
transcript_pct_melt2 <- data.frame(Subcluster=paste(final_table[,1], final_table[,2], sep=" - "), Order=(numrow+as.numeric(final_table[,4])), apob=apob2, Freq_Percent_by_Subcluster=final_table[,16])
transcript_pct_melt <- rbind(transcript_pct_melt,transcript_pct_melt2)
# look
transcript_pct_melt

# set plot name
plot_name = "apobM-UP-DOWN"

# create a plot - STACKED, not dodged - OPTIONAL - use ylim to control scale / keep same for all subclusters
p <- ggplot(data=transcript_pct_melt, 
            aes(x=reorder(Subcluster, Order), y=Freq_Percent_by_Subcluster, fill = factor(apob, levels=c("apobM_UP", "apobM_DOWN")))) +
  geom_col(width = 0.75) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=10)) +
  theme(axis.text.y = element_text(size=24)) +
  theme(axis.ticks = element_line(size = 1, color="black")) + # adds tick marks
  theme(axis.ticks.x=element_blank()) + # removes x-axis ticks only
  theme(legend.title=element_text(size=12)) + 
  theme(legend.text=element_text(size=12)) +
  theme(legend.position = 'bottom') +
  ylim(0,0.8) + #80 for apob
  labs(title=paste(subcluster_name, " transcript frequency", sep=""), x="", y = "")
p + scale_fill_manual(values=c('salmon','lightblue')) + labs(fill = plot_name) # change colors and legend title -- uses plot_name variable so you don't need to change it!

# write plot to file:

data_name = "transcript_percent"

# 1. Open a pdf file --  choose one -- set width for consistency based on number of states/subclusters

# epidermis - 13 subclusters, so 13/2 = 6.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4) # creates pdf file with width 6.5 and height 4 
# pharynx - 10 subclusters, so 10/2 = 5.0"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# muscle - 14 subclusters, so 14/2 = 7.0"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# cathepsin - 17 subclusters, so 17/2 = 8.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=8.5, height=4)
# intestine - 9 subclusters, so 9/2 = 4.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# neural - 61 subclusters, so 61/2 = 30.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=30.5, height=4)
# parenchymal - 20 subclusters, so 20/2 = 5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=10, height=4)
# protonephridia - 6 subclusters, so 6/2 = 3"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)

# 2. Create plot from above
p + scale_fill_manual(values=c('salmon','lightblue')) + labs(fill = plot_name) # change colors and legend title -- uses pdf_name variable so you don't need to change it!

# 3. Close the pdf file
dev.off() 


##############################################################################################
##
## generate table of dysregulated transcripts per subcluster - apobS (severe)
##
##############################################################################################

##
## Table1 = apobS - UP
##

## apobS - set variables for table generation 

target = apobS_vs_egfp_fdr05_UP
target_name = "apobS_UP"

logFC = "logFC_apobS_vs_egfp"

## MERGE 1 - extract targets from subcluster using merge
S2_sub = merge(target, subcluster, by="geneID") 

# calculate subcluster frequency
table_S2_subclust = as.data.frame(table(S2_sub$cluster.enriched.in, exclude=NULL)) # make table a data frame in the environment
# rename column 1 to reflect subcluster
colnames(table_S2_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
# rename column 2 to reflect target
colnames(table_S2_subclust)[2] = paste("Freq_", target_name, sep="")
# normalize / calculate % of all extracted subclusters... questionable value but leave in for now
table_S2_subclust[,paste("Freq_Pct_subcluster_", target_name, sep="")] = table_S2_subclust$Freq/sum(table_S2_subclust$Freq)

# to normalize to apobM / calculate % of apobM that are dysregulated (use variable "target_name" above to create custom column)
numrow = as.numeric(NROW(as.data.frame(unique(target$geneID))))
table_S2_subclust[,paste("Freq_Pct_", target_name, sep="")] = table_S2_subclust[,2]/numrow # could also name the Frequency with paste but it should always be column2 here

# sort and look
table_S2_subclust <- table_S2_subclust[order(-table_S2_subclust[,4]),]
row.names(table_S2_subclust) <- NULL # resets rownames after sort (have to do last)
# look
table_S2_subclust

## calculate number of transcripts in each subcluster and normalize to that (most appropriate for seeing if a subcluster is affected)

# first create data framae
transcripts_per_subclust = as.data.frame(table(subcluster$cluster.enriched.in, exclude = NULL))
# then can rename column to something more sensible than "Var1" - again use "paste" to create custom colname
colnames(transcripts_per_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
colnames(transcripts_per_subclust)[2] = paste("transcripts_per_subcluster_", target_name, sep="")
# take a look
transcripts_per_subclust

## MERGE 2 - add to table_S2_subclust
transcripts_per_subclust <- merge(table_S2_subclust, transcripts_per_subclust, by.x = colnames(table_S2_subclust)[1], by.y = colnames(transcripts_per_subclust)[1])

# now calculate frequency per subcluster
transcripts_per_subclust[,paste("Freq_Pct_by_subclust_", target_name, sep="")] <- transcripts_per_subclust[,2]/transcripts_per_subclust[,5] # again, could name columns using paste but they should always be in the same position
# take a look: 
transcripts_per_subclust
# can sort if you want: 
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,4]),] # by % of orig
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,6]),] # by % of transcripts in subcluster
row.names(transcripts_per_subclust) <- NULL # resets rownames after sort
# take a look (nothing should be over 1.0): 
transcripts_per_subclust

## now use for loop to calculate logFCs and add to table

# make empty matrix (will convert to data frame at end)
zeng_logFC_by_finch_subclust <- matrix(, nrow=NROW(transcripts_per_subclust), ncol=3) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_logFC_by_finch_subclust) <- c(paste(subcluster_name, "_subcluster", sep=""), paste("transcript_num_", target_name, sep=""), paste("avg_logFC_", target_name, sep=""))

# now run a loop
# this populates the table by grepping out rows, one at a time, based on the entries in the "transcripts_per_subclust" table above
# must use "as.character" to properly print subcluster number since this is a factor in the "transcripts_per_subclust" data frame
# also uses "logFC" variable for column 3, defined above with other variables

for(i in 1:NROW(transcripts_per_subclust)) { # start loop
  S2_sub2 <- S2_sub[grep(paste("^", as.character(transcripts_per_subclust[i,1]), "$", sep=""), S2_sub$cluster.enriched.in),]
  zeng_logFC_by_finch_subclust[i,1:3] <- c(as.character(transcripts_per_subclust[i,1]), nrow(S2_sub2), mean(S2_sub2[,logFC], na.rm = TRUE)) # nice, specify row with i and columns 1:3
}

# convert to data frame
zeng_logFC_by_finch_subclust <- as.data.frame(zeng_logFC_by_finch_subclust)

# look
zeng_logFC_by_finch_subclust 

# merge with counts
Table1 <- merge(transcripts_per_subclust, zeng_logFC_by_finch_subclust) # should not need to specify "by" columns if all above worked - it will use col1, which should be identical

# look, verify that transcript numbers match, and that no Freq_Pct_by_subcluster is >1 (just a QC check)
Table1


##
## Table2 = apobS - DOWN
##

## apobS - set variables for table generation 

target = apobS_vs_egfp_fdr05_DOWN
target_name = "apobS_DOWN"

logFC = "logFC_apobS_vs_egfp"

## MERGE 1 - extract targets from subcluster using merge
S2_sub = merge(target, subcluster, by="geneID") 

# calculate subcluster frequency
table_S2_subclust = as.data.frame(table(S2_sub$cluster.enriched.in, exclude=NULL)) # make table a data frame in the environment
# rename column 1 to reflect subcluster
colnames(table_S2_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
# rename column 2 to reflect target
colnames(table_S2_subclust)[2] = paste("Freq_", target_name, sep="")
# normalize / calculate % of all extracted subclusters... questionable value but leave in for now
table_S2_subclust[,paste("Freq_Pct_subcluster_", target_name, sep="")] = table_S2_subclust$Freq/sum(table_S2_subclust$Freq)

# to normalize to apobM / calculate % of apobM that are dysregulated (use variable "target_name" above to create custom column)
numrow = as.numeric(NROW(as.data.frame(unique(target$geneID))))
table_S2_subclust[,paste("Freq_Pct_", target_name, sep="")] = table_S2_subclust[,2]/numrow # could also name the Frequency with paste but it should always be column2 here

# sort and look
table_S2_subclust <- table_S2_subclust[order(-table_S2_subclust[,4]),]
row.names(table_S2_subclust) <- NULL # resets rownames after sort (have to do last)
# look
table_S2_subclust

## calculate number of transcripts in each subcluster and normalize to that (most appropriate for seeing if a subcluster is affected)

# first create data framae
transcripts_per_subclust = as.data.frame(table(subcluster$cluster.enriched.in, exclude = NULL))
# then can rename column to something more sensible than "Var1" - again use "paste" to create custom colname
colnames(transcripts_per_subclust)[1] = paste(subcluster_name, "_subcluster", sep="")
colnames(transcripts_per_subclust)[2] = paste("transcripts_per_subcluster_", target_name, sep="")
# take a look
transcripts_per_subclust

## MERGE 2 - add to table_S2_subclust
transcripts_per_subclust <- merge(table_S2_subclust, transcripts_per_subclust, by.x = colnames(table_S2_subclust)[1], by.y = colnames(transcripts_per_subclust)[1])

# now calculate frequency per subcluster
transcripts_per_subclust[,paste("Freq_Pct_by_subclust_", target_name, sep="")] <- transcripts_per_subclust[,2]/transcripts_per_subclust[,5] # again, could name columns using paste but they should always be in the same position
# take a look: 
transcripts_per_subclust
# can sort if you want: 
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,4]),] # by % of orig
transcripts_per_subclust <- transcripts_per_subclust[order(-transcripts_per_subclust[,6]),] # by % of transcripts in subcluster
row.names(transcripts_per_subclust) <- NULL # resets rownames after sort
# take a look (nothing should be over 1.0): 
transcripts_per_subclust

## now use for loop to calculate logFCs and add to table

# make empty matrix (will convert to data frame at end)
zeng_logFC_by_finch_subclust <- matrix(, nrow=NROW(transcripts_per_subclust), ncol=3) # uses NROW to calculate the number of rows needed based on transcripts_per_subclust data frame from above
# name the columns using paste and variables from above
colnames(zeng_logFC_by_finch_subclust) <- c(paste(subcluster_name, "_subcluster", sep=""), paste("transcript_num_", target_name, sep=""), paste("avg_logFC_", target_name, sep=""))

# now run a loop
# this populates the table by grepping out rows, one at a time, based on the entries in the "transcripts_per_subclust" table above
# must use "as.character" to properly print subcluster number since this is a factor in the "transcripts_per_subclust" data frame
# also uses "logFC" variable for column 3, defined above with other variables

for(i in 1:NROW(transcripts_per_subclust)) { # start loop
  S2_sub2 <- S2_sub[grep(paste("^", as.character(transcripts_per_subclust[i,1]), "$", sep=""), S2_sub$cluster.enriched.in),]
  zeng_logFC_by_finch_subclust[i,1:3] <- c(as.character(transcripts_per_subclust[i,1]), nrow(S2_sub2), mean(S2_sub2[,logFC], na.rm = TRUE)) # nice, specify row with i and columns 1:3
}

# convert to data frame
zeng_logFC_by_finch_subclust <- as.data.frame(zeng_logFC_by_finch_subclust)

# look
zeng_logFC_by_finch_subclust 

# merge with counts
Table2 <- merge(transcripts_per_subclust, zeng_logFC_by_finch_subclust) # should not need to specify "by" columns if all above worked - it will use col1, which should be identical

# look, verify that transcript numbers match, and that no Freq_Pct_by_subcluster is >1 (just a QC check)
Table2

##
## merge two tables
##

final_table = merge(Table1, Table2, all.x=TRUE, all.y=TRUE) # need to be TRUE in case some subclusters have NAs


##############################################################################################
##
## merge final_table with subcluster state tables for plotting - apobM (mild)
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

# merge subcluster and state columns
state$Subcluster <- paste(state$S2_Subcluster, " - ", state$State, sep="")
# assign order based on the order of the "state" file - for plotting below
state$order <- rownames(state)
# rename first column to enable merge with final_table
colnames(state)[1] <- paste(subcluster_name, "_subcluster", sep="")

# merge with final_table
final_table = merge(state, final_table, all.x = TRUE, all.y = TRUE) # must be TRUE to include subclusters with NAs

# write table, if desired: apobM-UP/DOWN
write.table(final_table, paste("20210328_apobS-UP-DOWN_extract_from_", subcluster_name, "_BY_SUBCLUSTER.xls", sep=""), row.names=FALSE, sep="\t")


##############################################################################################
##
## plot - apob-S
##
##############################################################################################

# set variables
apob1 = "apobS_UP"
apob2 = "apobS_DOWN"

# OPTION 1: create a "melted" table manually - apobM - ONLY subcluster number
transcript_pct_melt <- data.frame(Subcluster=final_table[,1], Order=as.numeric(final_table[,4]), apob=apob1, Freq_Percent_by_Subcluster=final_table[,9])
numrow = as.numeric(NROW(final_table)) # to calc multiplication factor for Order
transcript_pct_melt2 <- data.frame(Subcluster=final_table[,1], Order=(numrow+as.numeric(final_table[,4])), apob=apob2, Freq_Percent_by_Subcluster=final_table[,16])
transcript_pct_melt <- rbind(transcript_pct_melt,transcript_pct_melt2)
# look
transcript_pct_melt

# OPTION 2: create a "melted" table manually - apobM -subcluster number and states
transcript_pct_melt <- data.frame(Subcluster=paste(final_table[,1], final_table[,2], sep=" - "), Order=as.numeric(final_table[,4]), apob=apob1, Freq_Percent_by_Subcluster=final_table[,9])
numrow = as.numeric(NROW(final_table)) # to calc multiplication factor for Order
transcript_pct_melt2 <- data.frame(Subcluster=paste(final_table[,1], final_table[,2], sep=" - "), Order=(numrow+as.numeric(final_table[,4])), apob=apob2, Freq_Percent_by_Subcluster=final_table[,16])
transcript_pct_melt <- rbind(transcript_pct_melt,transcript_pct_melt2)
# look
transcript_pct_melt

# set plot name
plot_name = "apobS-UP-DOWN"

# create a plot - STACKED, not dodged - OPTIONAL - use ylim to control scale / keep same for all subclusters
p <- ggplot(data=transcript_pct_melt, 
            aes(x=reorder(Subcluster, Order), y=Freq_Percent_by_Subcluster, fill = factor(apob, levels=c("apobS_UP", "apobS_DOWN")))) +
  geom_col(width = 0.75) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size=10)) +
  theme(axis.text.y = element_text(size=24)) +
  theme(axis.ticks = element_line(size = 1, color="black")) + # adds tick marks
  theme(axis.ticks.x=element_blank()) + # removes x-axis ticks only
  theme(legend.title=element_text(size=12)) + 
  theme(legend.text=element_text(size=12)) +
  theme(legend.position = 'bottom') +
  ylim(0,0.8) + #80 for apob
  labs(title=paste(subcluster_name, " transcript frequency", sep=""), x="", y = "")
p + scale_fill_manual(values=c('salmon','lightblue')) + labs(fill = plot_name) # change colors and legend title -- uses plot_name variable so you don't need to change it!

# write plot to file:

data_name = "transcript_percent"

# 1. Open a pdf file --  choose one -- set width for consistency based on number of states/subclusters
# epidermis - 13 subclusters, so 13/2 = 6.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4) # creates pdf file with width 6.5 and height 4 
# pharynx - 10 subclusters, so 10/2 = 5.0"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# muscle - 14 subclusters, so 14/2 = 7.0"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# cathepsin - 17 subclusters, so 17/2 = 8.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=8.5, height=4)
# intestine - 9 subclusters, so 9/2 = 4.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)
# neural - 61 subclusters, so 61/2 = 30.5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=30.5, height=4)
# parenchymal - 20 subclusters, so 20/2 = 5"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=10, height=4)
# protonephridial - 6 subclusters, so 6/2 = 3"  (all by 4" high)
pdf(paste("20211028", subcluster_name, plot_name, data_name, "by_SUBCLUSTER.pdf", sep="_"), width=6.5, height=4)

# 2. Create plot from above
p + scale_fill_manual(values=c('salmon','lightblue')) + labs(fill = plot_name) # change colors and legend title -- uses pdf_name variable so you don't need to change it!

# 3. Close the pdf file
dev.off() 


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
  [1] ggplot2_3.3.6

loaded via a namespace (and not attached):
  [1] fansi_1.0.3         digest_0.6.29       withr_2.5.0         utf8_1.2.2          crayon_1.5.1        grid_4.1.1          R6_2.5.1            lifecycle_1.0.1    
[9] gtable_0.3.0        magrittr_2.0.3      scales_1.2.0        pillar_1.7.0        rlang_1.0.2         cli_3.3.0           farver_2.1.0        vctrs_0.4.1        
[17] ellipsis_0.3.2      labeling_0.4.2      tools_4.1.1         glue_1.6.2          munsell_0.5.0       compiler_4.1.1      pkgconfig_2.0.3     colorspace_2.0-3   
[25] BiocManager_1.30.18 tibble_3.1.7  
