# This script provides an example pipeline to study gene DNA methylation in P.californicus with an annotated genome reference

# load packages


library(data.table)
library(ggplot2)

setwd("/Users/taniachavarriapizarro/Desktop/Methylationfreq")


setwd("/Users/taniachavarriapizarro/Desktop/")

#Plotting methylation average introns vs exons 
data <- fread("P.c_Worker.freq.tsv", h=T)
head(data)


# remove unnecessary columns
data[, c("V9","V10", "V11", "V13","V15"):= NULL]


# namecolumns

names(data) <- c("chrom", "pos", "end", "num_motifs_in_group", "total_bases", "methylated_bases",
                 "methylated_frequency","group_sequence")
head(data)

data <- data[total_bases >= 3]

# remove unnecessary columns

data[, c("num_motifs_in_group","group_sequence","ID"):= NULL]

# add 1 to positions so that they are one-based instead of zero-based
data[, pos:=pos+1]
data

fwrite(data, "Methylfreqqueenexonintrons.tsv", sep="\t")

data <- fread("Methylfreqqueenexonintrons.tsv", h=T)
head(data)

data[, c("ID"):= NULL]
head(data)

####exon vs introns figure

# Load necessary libraries
library(ggplot2)
library(dplyr)

chromosomes_to_plot <- c("LG1", "LG2", "LG3","LG4","LG5","LG6",
                         "LG7","LG8","LG9","LG10","LG11","LG12",
                         "LG13","LG14","LG15","LG16") # Replace with your chromosome names

# Filter the cytoband data and scaffold lengths for only the specified chromosomes
combined_cytoband <- cytoband %>% filter(chromosome %in% chromosomes_to_plot)


# Filter data for "exons" and "introns"
filtered_data <- data %>%
  filter(type %in% c("exon", "intron") & total_bases <= 3000 & chrom %in% chromosomes_to_plot)

filtered_exon_data <- data %>%
  filter(type == "exon" & total_bases <= 3000 & chrom %in% chromosomes_to_plot)

filtered_intron_data <- data %>%
  filter(type == "intron" & total_bases <= 3000 & chrom %in% chromosomes_to_plot
         & methylated_frequency >= 0 & methylated_frequency <= 0.25)



# Load necessary libraries
library(ggplot2)
library(dplyr)
### ("exon" = "#B8DE29FF", "intron" = "darkgreen")

# Create the plot
plot <- ggplot(filtered_data, aes(x = total_bases, y = methylated_frequency, color = type)) +
  geom_point(size = 2, alpha = 0.1) +
  scale_color_manual(values = c("exon" ="#B8DE29FF" , "intron" = "darkgreen")) + # Custom colors
  labs(
    x = "Total bases",
    y = "Methylation Frequency",
    color = "Type"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14), # Improve text readability
    legend.position = "top" # Move legend to the top
  )

# Display the plot
print(plot)


###TEs and intergenetic methylation 
data <- fread("CA4queen.TE.fa.out.tsv", h=F)
head(data)

# remove unnecessary columns
data[, c("V9","V10", "V11"):= NULL]
head(data)

# namecolumns

names(data) <- c("chrom", "pos", "end", "num_motifs_in_group", "total_bases", "methylated_bases",
                 "methylated_frequency","group_sequence","strand","ID")

###namesIntergenetic
names(data) <- c("chrom", "pos", "end", "num_motifs_in_group", "total_bases", "methylated_bases",
                 "methylated_frequency","group_sequence")

head(data)

data <- data[total_bases >= 3]

# remove unnecessary columns

data[, c("num_motifs_in_group","group_sequence"):= NULL]

# add 1 to positions so that they are one-based instead of zero-based
data[, pos:=pos+1]
data

fwrite(data, "MethylationfreqIntergeneticqueen.tsv", sep="\t")

# Install and load the stringr package
install.packages("stringr")
library(stringr)

# Assuming your data frame is named your_data and the column is named your_column
###just TE
data$TE_ID <- str_extract(data$ID, '""(.*?)""')
head(data)
data$type <- str_extract(data$ID, '(.*?)""')
head(data)

data[, c("ID"):= NULL]
head(data)

fwrite(data, "TEqueenIDmethylfreq.tsv", sep="\t")

data <- fread("CA4queen.TE.fa.out.tsv", h=F)
head(data)

data[, c("V8","V9","V10", "V11","V12"):= NULL]
head(data)

# namecolumns

names(data) <- c("chrom", "start", "end", "num_motifs_in_group", "total_bases", "methylated_bases",
                 "methylated_frequency","ID")


head(data)

data <- data[total_bases >= 3]


library(tidyverse)

# Change chromosome names using the mapping
methylation_data <- data %>%
  mutate(chrom = factor(chrom, levels = names(chromosome_name_mapping), labels = chromosome_name_mapping))

head(methylation_data)

###calculation methylation freq average and sd

grouped_data <- methylation_data %>%
  group_by(ID) %>%
  summarise(mean_methylation = mean(methylated_frequency),
            sd_methylation = sd(methylated_frequency))

head(grouped_data)
library(dplyr)

grouped_data <- grouped_data %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

fwrite(grouped_data, "AverageTEMeanMethylperType.tsv", sep="\t")


grouped_data <- grouped_data %>%
  mutate(extracted_name = sub('.*\\"{2}(.*?)\\"{2}.*', '\\1', ID))

head(grouped_data)
grouped_data <- grouped_data %>%
  mutate(name_without_suffix = sub('_(\\d+)$', '', extracted_name))
head(grouped_data)
fwrite(grouped_data, "AverageTEMeanMethylperType.tsv", sep="\t")

TEdata <- fread("summary_TEdataMethyl.tsv", h=T)
library(ggplot2)
library(viridis)

summary_TEdata <- TEdata %>%
  group_by(name_without_suffix) %>%
  summarise(
    mean_value = mean(mean_value),
    sd_value = sd(sd_value),
    .groups = 'drop'
  )

# View the summary dataframe
print(summary_TEdata)
summary_TEdata[is.na(summary_TEdata)] <- 0
fwrite(summary_TEdata, "summary_TEdataMethyl.tsv", sep="\t")
summary_TEdata[is.na(summary_TEdata)] <- 0
summary_TEdata <- summary_TEdata %>%
  mutate(sd_value = ifelse(sd_value < 0, 0, sd_value))


ggplot(summary_TEdata, aes(x = name_without_suffix, y = mean_value, fill = name_without_suffix)) +
  geom_bar(stat = "identity", color = "black") +  # Bar plot with outline
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                width = 0.2, color = "black") +      # Error bars for SD
  scale_fill_viridis(discrete = TRUE) +           # Viridis color scale for discrete categories
  labs(
    title = "Mean and SD by Name Without Suffix",
    x = "Name Without Suffix",
    y = "Mean Value"
  ) +
  theme_minimal()


###intergenetic

grouped_data <- methylation_data %>%
  summarise(mean_methylation = mean(methylated_frequency),
            sd_methylation = sd(methylated_frequency))

head(grouped_data)
fwrite(grouped_data, "MeanInterMethyl_queen.tsv", sep="\t")

grouped_data <- methylation_data %>%
  group_by(chrom) %>%
  summarise(mean_methylation = mean(methylated_frequency),
            sd_methylation = sd(methylated_frequency))

head(grouped_data)
fwrite(grouped_data, "AveInterMethylperChromH1.tsv", sep="\t")

###TE
ggplot(grouped_data, aes(x = TE_ID, y = mean_methylation, fill = TE_ID)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean_methylation - sd_methylation, ymax = mean_methylation + sd_methylation),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "Average Methylation by TE_ID Pleometrotic P.californicus",
       x = "TE_ID",
       y = "Methylation Frequency") +
  scale_fill_viridis_d(direction = -1)+
  theme_minimal()

###Intergenetic

ggplot(grouped_data, aes(x = chrom, y = mean_methylation, fill = chrom)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean_methylation - sd_methylation, ymax = mean_methylation + sd_methylation),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "Average Methylation Intergenetic by chromosom Haplometrotic P.californicus",
       x = "Chromoso",
       y = "Methylation Frequency") +
  scale_fill_viridis_d(direction = -1)+
  theme(legend.position = "none")


####Finally Methylation mean by type

data <- fread("MeanMethylperTypeH1.tsv", h=T)
head(data)

ggplot(data, aes(x = type, y = mean_methylation, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean_methylation - sd_methylation, ymax = mean_methylation + sd_methylation),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "Average Methylation by type Haplometrotic P.californicus",
       x = "type",
       y = "Methylation Frequency") +
  scale_fill_viridis_d(direction = -1)+
  theme_minimal()


# This script provides an example pipeline to study gene DNA methylation in P.californicus with an annotated genome reference

# load packages
library(data.table)
library(ggplot2)

setwd("/Users/taniachavarriapizarro/Desktop/")

# data format must be one line per cytosine:
# chromosome	position	strand(+/-)	mc_class=context(CG/CHG/CHH)	methylated_bases(number_methylated_C_reads)	total_bases(total_number_reads_mapped)
# read in the data
data <- fread("P.c_Worker.freq.tsv", h=T)
data

#5952667

###totalnumber calledsites
gc_calledsites <- data$called_sites
total_called_sites <- sum(gc_calledsites) #497298407.  ##queen 491381680

####totalnumbermethylatedcites
gc_methylatedsites <- data$called_sites_methylated
total_called_sites <- sum(gc_methylatedsites) #13316233.  3%  ###13124213 $3%


##6076070. 5952667  56%
###totalcpgmethylated

#3180956. 3357346 



# Only keep positions with 3 or more reads mapped
data <- data[called_sites >= 3] #3352854 


# read in the ONT data (one line per cytosine)
# before reading in the ONT file, make sure you replace spaces by tabs in the file
# in bash: tr " " "\t" < input_file > input_file_tables

# please note that ONT output files often change format, make sure the columns are correctly named below

names(data) <- c("chrom", "pos", "end", "num_motifs_in_group", "total_bases", "methylated_bases", 
                 " methylated_frequency", "group_sequence")


data[, c("end","num_motifs_in_group", "group_sequence"):= NULL]

# add 1 to positions so that they are one-based instead of zero-based
data[, pos:=pos+1]
data


#### this part of the script runs a binomial test on every gene to infer its methylation state

gff3 <- fread("/Users/taniachavarriapizarro/Desktop/Methylation Pogo/Pcal3.1.clean.submitted.polished2.gff3", h=F)
gff3
gff3[, c("V2", "V6","V8") := NULL] # delete unnecessary columns
gff3

colnames(gff3) <- c("chr", "type", "start", "stop", "strand", "ID")
gff3
gff3 <- gff3[type == "gene"] # retain only gene information
gff3[, GeneName := strsplit(ID, 'ID=', perl=T)[[1]][2], by=1:nrow(gff3)]# extract gene name from ID column
fwrite(gff3, "genes.gff3", sep="\t")
gff3[, c("ID") := NULL] # delete unnecessary columns
gff3
fwrite(gff3, "genes.gff3", sep="\t")
gff3 <- fread("/Users/taniachavarriapizarro/Desktop/Methylation Pogo/genes.gff3", h=T)


# Create output file for Gene cytosine methylation counts
header <- data.table("Gene" = character(), "Chromosome" = character(), "start" = numeric(), "stop" = numeric(), "strand" = character(), 
                     "Number_CG_total_reads" = numeric(), "Number_methylated_CG_reads" = numeric())
fwrite(header, file = "gene_methylationP.c_Worker.txt", sep="\t", append=F)

# loop on each gene to count number of methylated and unmethylated cytosines in each context 
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentstart <- gff3$start[i]
  Currentstop <- gff3$stop[i]
  Currentstrand <- gff3$strand[i]
  
  # extract gene methylation for current gene
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentstart)&(pos<=Currentstop)]
  
  # check there is something
  if (nrow(GeneMethylation)>0) {
    # Sum the number of methylated and total reads
    sumMeth <- sum(GeneMethylation$methylated_bases)
    sumtotal <- sum(GeneMethylation$total_bases)
    GeneData <- data.table(as.character(CurrentGene), as.character(CurrentChr), as.numeric(Currentstart), as.numeric(Currentstop), as.character(Currentstrand), sumtotal, sumMeth)
    
    # add gene info to output file
    fwrite(GeneData, file="gene_methylationP.c_Worker.txt", sep="\t", append=TRUE)
  }
}


# load Gene data
GeneData <- fread(file="gene_methylationP.c_Worker.txt", h=T)
GeneData

# computation of average CG gene methylation across all genes
proba_methylated_CG <- sum(GeneData$Number_methylated_CG_reads, na.rm=T)/sum(GeneData$Number_CG_total_reads, na.rm=T)

# binom test on each gene to test for significant CG methylation compared to other genes
GeneData[,  CG_Binom_pvalue_greater := if (Number_CG_total_reads > 0) {binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "greater")$p.value} else {as.numeric(NA)}, by = 1:nrow(GeneData)]
GeneData[,  CG_Binom_pvalue_less := if (Number_CG_total_reads > 0) {binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "less")$p.value} else {as.numeric(NA)}, by = 1:nrow(GeneData)]
GeneData

# adjust p-values for multiple tests with Benjamini & Hochberg (1995) ("BH")
GeneData[,  CG_Binom_pvalue_greater_adjusted := p.adjust(CG_Binom_pvalue_greater, "BH")]
GeneData[,  CG_Binom_pvalue_less_adjusted := p.adjust(CG_Binom_pvalue_less, "BH")]
GeneData

# Classify gene methylation
Methylation <- function(u,x) {
  # u the adjusted p-value for CG methylation to be higher than expected by chance
  # x the adjusted p-value for CG methylation to be lower than expected by chance
  # outputs mCG if gene is significantly CG body methylated (u < 0.05)
  if ((!is.na(u))&(u < 0.05)) {
    # CG body methylated gene
    return("gbM")
  } else if ((!is.na(x))&(x < 0.05)) {
    # CG unmethylated gene
    return("UM")
  } else if ((!is.na(x))&(x > 0.05)&(!is.na(u))&(u > 0.05)) {
    # CG intermediately methylated gene
    return("IM")
  } else {
    return("NA")
  }
}
GeneData[, GeneMethylationCall := mapply(Methylation, CG_Binom_pvalue_greater_adjusted, CG_Binom_pvalue_less_adjusted)]
GeneData
table(GeneData$GeneMethylationCall)
GeneData <-  fwrite(GeneData, file="gene_methylationCallP.c_Worker.txt", sep="\t", append=TRUE)
# Note that in some papers intermediately methylated genes 'IM' are merged together with unmethylated genes 'UM' 



###################################################################################################
#          Plotting Gene methylation metaplot and studying promoter methylation
###################################################################################################
# loop on all genes to compute methylation proportions on windows upstream of the gene, downstream and inside the gene.
header <- list("SumMethylatedReads", "TotalReads", "window", "GeneName")
fwrite(header, file="Gene_windows_P.c_Worker.txt", sep="\t", append=FALSE)
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases))]
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows_P.c_Worker.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows <- fread("Gene_windows_P.c_Worker.txt", h=T)
Windows

# Summary over all genes
SummaryWindows <- Windows[TotalReads > 0, .(sum(SumMethylatedReads)/sum(TotalReads), binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][1], binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][2]), by=c("window")]
names(SummaryWindows) <- c("window", "Average_Methylation", "CI_low", "CI_high")
fwrite(SummaryWindows, file="Summarywindows_PcqWorker.txt", sep="\t", append=TRUE)


# plot metaplot
jpeg(f="windows_weighteg_methylation_Genes.jpeg", width=8, height=4, units = "in", res=500, pointsize = 18, quality=100)
ggplot(data=SummaryWindows, aes(x = window, y = Average_Methylation)) + 
  theme_bw(base_size=20) +
  geom_point(aes(x = window, y = Average_Methylation), size = 1) +  
  geom_line(aes(x = window, y = Average_Methylation), linewidth=0.5) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), linetype=1, size=0.2, alpha=0.1, show.legend = FALSE) +
  labs(x=paste("position", sep=''), y="Weighted methylation") + 
  guides(col=guide_legend(title="Methylation type:")) +
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000'))

dev.off()





###comparing methods Nanopore and WGBS

SummaryWindows <- fread("/Users/taniachavarriapizarro/Desktop/New fig methylome/Summarywindows_Q_WGBS6vNano10.txt", h=T)

ggplot(data = SummaryWindows, aes(x = window, y = Average_Methylation, group = type, color = type, fill = type)) +
  theme_bw(base_size = 20) +
  geom_point(size = 1, shape = 20) +  
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), linetype = 1, size = 0.6, alpha = 0.3) +
  labs(x=paste("position", sep=''), y="Weighted methylation") + 
  scale_y_continuous(limits = c(0, 0.15)) +
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000')) +
  scale_color_manual(values = c("#31a354","#756bb1")) +
  scale_fill_manual(values = c(alpha("#31a354", 0.2),alpha("#756bb1", 0.2))) +
  theme(legend.position = "right")


###comparing life stages and castes haplometrotic 

SummaryWindows <- fread("Summarywindows_5L_4P_W_Q.txt", h=T)

ggplot(data = SummaryWindows, aes(x = window, y = Average_Methylation, group = type, color = type, fill = type)) +
  theme_bw(base_size = 20) +
  geom_point(size = 1, shape = 20) +  
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), linetype = 1, size = 0.2, alpha = 0.1) +
  labs(x=paste("position", sep=''), y="Weighted methylation") + 
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000')) +
  scale_color_manual(values = c("#ffcc33","#993300","#ff0000","#ff6600")) +
  scale_fill_manual(values = c(alpha("#ffcc33", 0.1),alpha("#993300", 0.1),alpha("#ff0000", 0.1),alpha("#ff6600", 0.1))) +
  theme(legend.position = "right")+
  ylim(0, 0.11)


# gene promoter methylation (use window 20 which corresponds to 200bp upstream of Transcription start site)
# Nota Bene: In animals, promoters are usually longer than the typical plant 200pb promoter, 
# you can take more windows than just window 20 (for example, window 15 to 20)
Windows <- Windows[window==20]
Windows[, window := NULL]

# compute promoter average proportion of methylated CG, CHG and CHH sites
promoter_proba_methylated_CG <- sum(Windows$SumMethylatedReads, na.rm=T)/sum(Windows$TotalReads, na.rm=T)

# compare promoter methylation to genome average using a binomial test
BinomTest_greater <- function(x, n, p){ if(n<1) { return(as.double(NA))} else { return(binom.test(x, n, p, alternative = "greater")$p.value) }}
Binom_pvalues_greater_CG <- Windows[, BinomTest_greater(SumMethylatedReads, TotalReads, promoter_proba_methylated_CG), by = 1:nrow(Windows)]

# adjust p-values
Windows[, Padjusted_greater_CG :=  p.adjust(Binom_pvalues_greater_CG$V1, "BH")]

# inference of CG promoter methylated genes
Windows[, promoter_CG_methylated := if((!is.na(Padjusted_greater_CG))&(Padjusted_greater_CG <= 0.05)){"yes"}else{"no"}, by = 1:nrow(Windows)]
fwrite(Windows, "Promoter_methylation_PcqWorker.txt", sep="\t")


####Venn plot comparing methylated genes

# Install and load the necessary package
install.packages("VennDiagram")
library(VennDiagram)

# Read the gene lists from text files
genes_larvae <- read.table("GBM_5L.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_pupae <- read.table("GBM_4P.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_queen <- read.table("GBM_Q.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_workers <- read.table("GBM_workers2.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# Load the VennDiagram library
library(VennDiagram)


##Find unique genes pupae and queen

library(data.table)


# Find unique genes for each category
unique_larvae <- setdiff(genes_larvae, union(union(genes_pupae, genes_queen), genes_workers))
unique_pupae <- setdiff(genes_pupae, union(union(genes_larvae, genes_queen), genes_workers))
unique_queen <- setdiff(genes_queen, union(union(genes_larvae, genes_pupae), genes_workers))
unique_workers <- setdiff(genes_workers, union(union(genes_larvae, genes_pupae), genes_queen))

# Save unique genes to files
write.table(unique_larvae, "unique_genes_larvae.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique_pupae, "unique_genes_pupae.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique_queen, "unique_genes_queen.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique_workers, "unique_genes_workers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Output a summary
cat("Unique genes for Larvae:", length(unique_larvae), "\n")
cat("Unique genes for Pupae:", length(unique_pupae), "\n")
cat("Unique genes for Queen:", length(unique_queen), "\n")
cat("Unique genes for Workers:", length(unique_workers), "\n")

# Load necessary library
library(VennDiagram)

# Read the gene lists from text files
genes_larvae <- read.table("GBM_5L.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_pupae <- read.table("GBM_4P.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_queen <- read.table("GBM_Q.txt", header = FALSE, stringsAsFactors = FALSE)$V1
genes_workers <- read.table("GBM_workers2.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# Load the VennDiagram library
library(VennDiagram)

# Calculate pairwise, triple, and quadruple intersections
n12 <- length(intersect(genes_larvae, genes_pupae))
n13 <- length(intersect(genes_larvae, genes_queen))
n14 <- length(intersect(genes_larvae, genes_workers))
n23 <- length(intersect(genes_pupae, genes_queen))
n24 <- length(intersect(genes_pupae, genes_workers))
n34 <- length(intersect(genes_queen, genes_workers))

n123 <- length(Reduce(intersect, list(genes_larvae, genes_pupae, genes_queen)))
n124 <- length(Reduce(intersect, list(genes_larvae, genes_pupae, genes_workers)))
n134 <- length(Reduce(intersect, list(genes_larvae, genes_queen, genes_workers)))
n234 <- length(Reduce(intersect, list(genes_pupae, genes_queen, genes_workers)))

n1234 <- length(Reduce(intersect, list(genes_larvae, genes_pupae, genes_queen, genes_workers)))

# Create the Venn diagram
venn.plot <- draw.quad.venn(
  area1 = length(genes_larvae),
  area2 = length(genes_pupae),
  area3 = length(genes_queen),
  area4 = length(genes_workers),
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = c("Larvae", "Pupae", "Queen", "Workers"),
  fill = c("#ffcc33","#993300","#ff0000","#ff6600"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#ffcc33","#993300","#ff0000","#ff6600")
)

# Save the Venn diagram as an image file (optional)
# png("venn_diagram_larvae_pupae_queen_workers.png")
# grid.draw(venn.plot)
# dev.off()


#!/usr/bin/env Rscript

# load packages
library(ggplot2)
library(data.table)
library(glmmTMB)
library(DHARMa)
library(MuMIn)


library(reader)
library(readxlsb)
library(tidyverse)
library(rvest)
library(dplyr)
library(data.table)



#=================================================================================
# Plotting link between DNA methylation and gene expression in P.californicus
#=================================================================================
# load data


Col0 <- fread("gene_methylationCallH.txt", h=T)
Col0

table(Col0$GeneMethylationCall)
# Remove the specified substring from the "GeneName" column
Col0[, Gene := gsub(";Name=[^;]+", "", Gene)]

#mydata
datamethy <- read.delim("/Users/taniachavarriapizarro/Desktop/datamethy.tsv", header = T, sep = "\t")
head(datamethy)
###datamethy <- datamethy %>%
###mutate(log = as.numeric(log))

# Assuming you have a dataframe 'df' with columns 'gene', 'raw_counts', 'gene_length_bp' for each sample
datamethy <- datamethy %>%
  mutate(gene_length_kb = gene_length_bp / 1000)  # Convert gene length to kb

# Calculate FPKM
datamethy <- datamethy %>%
  group_by(Geneid) %>%
  mutate(total_mapped_reads = sum(counts),    # Total counts for each sample
         FPKM = (counts * 1e9) / (total_mapped_reads * gene_length_kb)) %>%
  ungroup()

datamethy <- datamethy %>%
  mutate(total_mapped_reads = sum(counts),    # Total counts for each sample
         FPKM = (counts * 1e9) / (total_mapped_reads * gene_length_kb)) 


setwd("/Users/taniachavarriapizarro/Desktop/New fig methylome")

head(datamethy)

fwrite(datamethy, "/Users/taniachavarriapizarro/Desktop/datamethyFPKM.tsv", sep="\t")

datamethy <- read.delim("/Users/taniachavarriapizarro/Desktop/New fig methylome/datamethyFPKM.tsv", header = T, sep = "\t")
head(datamethy)

# Assuming your dataframe is called 'df'
# Columns: gene, FPKM, and sample

# Normalize FPKM
datamethy <- datamethy %>%
  mutate(total_FPKM = sum(FPKM),              # Calculate total FPKM per sample
         normalized_FPKM = FPKM / total_FPKM,
         log_normalized_FPKM = log10(normalized_FPKM + 1)) 

datamethy <- datamethy %>%
  mutate(log_normalized_FPKM = log10(FPKM + 1))



# View the normalized dataframe
head(datamethy_normalized)

fwrite(datamethy_normalized, "/Users/taniachavarriapizarro/Desktop/datamethy_normalizedFPKM.tsv", sep="\t")


# Example data structure
# Assuming 'df' contains 'FPKM' and 'MethylationCall' columns
range(datamethy$log_normalized_FPKM)


# Boxplot with customized colors and y-axis limit
ggplot(datamethy, aes(x = GeneMethylationCall, y = log_normalized_FPKM, fill = GeneMethylationCall)) +
  theme_bw(base_size=16) +
  geom_boxplot(outlier.colour=NA) +
  coord_cartesian(ylim=c(0, 10)) +
  labs(fill = "Methylation\nstate", y="Normalized Expression Level") +
  scale_fill_manual(values=c("magenta", "cyan"), labels=c("gbM", "UM")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Genes", y="log_normalized_FPKM")

ggplot(datamethy, aes(x = GeneMethylationCall, y = FPKM, fill = GeneMethylationCall)) +
  theme_bw(base_size=16) +
  geom_boxplot(outlier.colour=NA) +
  coord_cartesian(ylim=c(0, 130000)) +
  labs(fill = "Methylation\nstate", y="Normalized Expression Level") +
  scale_fill_manual(values=c("magenta", "cyan"), labels=c("gbM", "UM")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Genes", y="FPKM")



# Statistical comparison of gbM and UM genes

# Assuming 'df_normalized' contains 'log_normalized_FPKM' and 'MethylationCall'

# Subset the data into two groups
group1 <- datamethy %>% filter(GeneMethylationCall == "gbM") %>% pull(log_normalized_FPKM)
group2 <- datamethy %>% filter(GeneMethylationCall == "UM") %>% pull(log_normalized_FPKM)

# Perform Wilcoxon test
wilcox_test <- wilcox.test(group1, group2, exact = FALSE)

# Display the results
print(wilcox_test)


# Subset the data into two groups
group1 <- datamethy %>% filter(GeneMethylationCall == "gbM") %>% pull(FPKM)
group2 <- datamethy %>% filter(GeneMethylationCall == "UM") %>% pull(FPKM)

# Perform Wilcoxon test
wilcox_test <- wilcox.test(group1, group2, exact = FALSE)

# Display the results
print(wilcox_test)


# Example dataset

# Calculate methylation frequency
datamethy$methylation_frequency <- datamethy$Number_methylated_CG_reads / datamethy$Number_CG_total_reads

# Perform correlation analysis
correlation_result <- cor.test(datamethy$methylation_frequency, datamethy$log_normalized_FPKM)

# Extract correlation coefficient and p-value
correlation_value <- round(correlation_result$estimate, 3)
p_value <- signif(correlation_result$p.value, 3)

head(datamethy)
range(datamethy$methylation_frequency)
range(datamethy$FPKM)



# Create the scatter plot with annotation and custom colors
ggplot(datamethy, aes(x = methylation_frequency, y = log_normalized_FPKM , color = GeneMethylationCall)) +
  geom_point(size = 1) + # Smaller points
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(
    title = "Methylation Frequency vs Gene Expression",
    x = "Methylation Frequency",
    y = "Expression Counts",
    color = "Category"
  ) +
  scale_color_manual(
    values = c("magenta", "cyan"),
    labels = c("gbM", "UM")
  ) +
  theme_classic()



###heatmap methylation frequency per genes and stages

# Step 1: Read the TSV file

data <- fread("/Users/taniachavarriapizarro/Desktop/gene_methylationP.c_Worker.tsv", h=T)
head(data)


# Step 2: Calculate methylation frequency and add as a new column
# Assuming 'total_bases' is the column with total bases and 'methylated_bases' is the column with methylated bases
data$methylation_frequency <- data$Number_methylated_CG_reads / data$Number_CG_total_reads

# Step 3: Overwrite the original file with the updated data (including the new column)
fwrite(data, "Gene_methylationPc_worker.tsv", sep=",")

data <- fread("Gene_methylationPc_worker.tsv", h=T, sep=",")


# Read each file
queen <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/queengenes.tsv", h=T, sep="\t")
head(queen)
worker <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/workersgenes2.tsv", header=T, sep="\t")
head(worker)
larvae <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/larvaegenes.tsv", h=T, sep="\t")
head(larvae)
pupae <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/pupaegenes.tsv", h=T, sep="\t")
head(pupae)


setwd("/Users/taniachavarriapizarro/Desktop/Methylationfreq")
data <- fread("common_genes_four_stages.tsv", h=T, sep="\t")
head(data)

data_matrix <- as.matrix(data[, -1])
# Convert methylation frequencies to percentages
data_matrix_percentage <- data_matrix * 100


rownames(data_matrix_percentage) <- data$GeneName
head(data_matrix_percentage)


# Install necessary packages if you haven't already
install.packages("pheatmap")
install.packages("viridis")

# Load necessary libraries
library(pheatmap)
library(viridis)

# Set a seed for reproducibility (optional)
set.seed(789)  # You can change the seed number for different random samples

# Randomly select 100 genes (rows)
selected_genes <- data_matrix_percentage[sample(1:nrow(data_matrix_percentage), 100), ]
selected_genes_clean <- selected_genes[complete.cases(selected_genes), ]


selected_genes[is.na(selected_genes)] <- 0
selected_genes[is.nan(selected_genes)] <- 0
selected_genes[is.infinite(selected_genes)] <- 0
# Check for any NA, NaN, or Inf values again
any(is.na(selected_genes_clean))
any(is.nan(selected_genes_clean))
any(is.infinite(selected_genes_clean))


###just with candidate genes
data <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/common_genes_four_stages.tsv", h=T, sep="\t")
head(data)

# Define the list of candidate genes
candidate_genes <- c("PCAL_000191", "PCAL_000651", "PCAL_000835","PCAL_000844","PCAL_000955","PCAL_001345", 
                     "PCAL_001375", "PCAL_001744", "PCAL_001766","PCAL_001788","PCAL_002551","PCAL_002632",
                     "PCAL_002649", "PCAL_002745","PCAL_002747", "PCAL_002885", "PCAL_003303","PCAL_003443",
                     "PCAL_003745","PCAL_007693","PCAL_006796")


# Read in the main data file
data <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/common_genes_four_stages.tsv", h=T, sep="\t")

# Read in the additional file with gene names and corresponding protein names
gene_function <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/Candidategeneslistmethyl2.tsv", h=T, sep="\t")

# Assuming 'GeneName' is the column with gene names and 'ProteinName' is the column with protein names
# Merge the two datasets based on the 'GeneName' column
merged_data <- merge(data, gene_function, by = "GeneName")
head(merged_data)
fwrite(merged_data, "/Users/taniachavarriapizarro/Desktop/Methylationfreq/candidategenesfunctionfourstages.tsv", sep=",")
# Combine 'GeneName' and 'GeneFunction' into a single string for each row
merged_data$CombinedName <- paste(merged_data$GeneName, merged_data$GeneFunction, sep = " - ")
head(merged_data)
fwrite(merged_data, "/Users/taniachavarriapizarro/Desktop/Methylationfreq/merge_datafourstages.tsv", sep=",")

data <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/merge_datafourstages.tsv", h=T, sep="\t")
head(data)

# Remove Gene and Protein columns and convert to matrix
# Exclude GeneName and ProteinName columns
data_matrix <- as.matrix(data[, -1])

rownames(data_matrix) <- data$CombinedName  # Use the protein names as row names
head(data_matrix)
fwrite(data_matrix, "/Users/taniachavarriapizarro/Desktop/Methylationfreq/data_matrixfourstages.tsv", sep=",")

data <- fread("/Users/taniachavarriapizarro/Desktop/Methylationfreq/merge_datafourstagesDNMT.tsv", h=T, sep="\t")
head(data)
data_matrix <- as.matrix(data[, -1])
# Convert methylation frequencies to percentages
data_matrix <- data_matrix * 100

# Create heatmap
pheatmap(
  data_matrix,
  color = viridis(10),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  main = "Differential Methylation Frequency for Candidate Genes"
)




###Plot GENES methylated at promoter vs Gene expression

datamethy <- read.delim("/Users/taniachavarriapizarro/Desktop/PomoterGenemethylexpr_H.csv", header = T, sep = ";")
head(datamethy)
datamethy <- datamethy %>%
  mutate(mean2 = as.numeric(counts))

ggplot(data=datamethy, aes(x=promoter_CG_methylated, y=log(mean2+1),fill=promoter_CG_methylated)) +
  theme_bw(base_size=16) +
  geom_boxplot(outlier.colour=NA) +
  coord_cartesian(ylim=c(0,20)) +
  labs(fill = "Methylation\nstate", y="Normalized Expression Level") +
  scale_fill_manual(values=c("cyan","magenta"), labels=c("no","yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Promoter Methylation", y="Gene expression level log(mean+1)")



###methylation comparation bisulfite and nanopore


setwd("/Users/taniachavarriapizarro/Desktop/")


setwd("/Users/taniachavarriapizarro/Desktop/")

#Plotting methylation average introns vs exons 
data <- fread("Pcp_PV_L7_1.methylation_compa.tsv", h=F)
head(data)


# remove unnecessary columns
data[, c("V8", "V9","V10", "V11"):= NULL]


# namecolumns

names(data) <- c("chrom", "pos", "end", "num_motifs_in_group", "total_bases", "methylated_bases",
                 "methylated_frequency1","methylated_frequencybisulper","methylbasesbisul","totalmethylbisul")
head(data)

data <- data[total_bases >= 3]

data <- data[totalmethylbisul >= 3]


# add 1 to positions so that they are one-based instead of zero-based
data[, pos:=pos+1]
data




fwrite(data, "MethylationfreqcomparaPcpPV7.tsv", sep="\t")


library(tidyverse)

# Change chromosome names using the mapping
methylation_data <- data %>%
  mutate(chrom = factor(chrom, levels = names(chromosome_name_mapping), labels = chromosome_name_mapping))

head(methylation_data)
###calculation methylation freq average and sd

# Calculate methylated frequencies for each row
methylated_frequencies <- numeric(nrow(methylation_data))
for (i in 1:nrow(methylation_data)) {
  methylated_frequencies[i] <- methylation_data$methylated_frequencybisulper[i] /100
}

# Add the calculated frequencies as a new column
methylation_data$methylated_frequency2 <- methylated_frequencies

head(methylation_data)
fwrite(methylation_data, "MethylationfreqcompaPcpPV7.tsv", sep="\t")

filtered_data <- methylation_data[methylation_data$methylated_bases != 0 & methylation_data$methylbasesbisul != 0, ]

filtered_data <- na.omit(methylation_data[, c("methylated_frequency1", "methylated_frequency2")])


fwrite(filtered_data, "MethylationfreqcompaPcpPV7_filtered_data.tsv", sep="\t")

# Calculate the correlation
c <- cor(filtered_data$methylated_frequency1, filtered_data$methylated_frequency2)



library(ggplot2)
library(RColorBrewer)


# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)


title <- sprintf("N = %d r = %.3f", nrow(filtered_data), 0.82)
ggplot(filtered_data, aes(methylated_frequency1, methylated_frequency2)) +
  geom_bin2d(bins=25) + scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Bisulfite Methylation Frequency") +
  ylab("Nanopolish Methylation Frequency") +
  theme_bw(base_size=20) +
  ggtitle(title)

###regular correlation plot

title <- sprintf("N = %d r = %.3f", nrow(filtered_data), 0.82)
ggplot(filtered_data, aes(methylated_frequency1, methylated_frequency2)) +
  geom_point(size = 1) +                                # Add points to the plot
  geom_smooth(method = "lm", col = "red") +     # Add the regression line
  labs(x = "Bisulfite Methylation Frequency",
       y = "Nanopolish Methylation Frequency") +
  theme_minimal()   

dev.off()



####Now with WGBS


setwd("/Users/taniachavarriapizarro/Desktop/MethylationfreqWGBS")

# data format must be one line per cytosine:
# chromosome	position	strand(+/-)	mc_class=context(CG/CHG/CHH)	methylated_bases(number_methylated_C_reads)	total_bases(total_number_reads_mapped)
# read in the data
data <- fread("Pcq_6b.CpG_report.txt", h=F)
data

names(data) <- c("chrom", "pos","strand", "methylated_bases", "total_bases", "motifs_in_group", "group_sequence")
data

data <- data[order(data$chrom), ]
data


#5952667

###totalnumber calledsites
gc_calledsites <- data$total_bases       #pupa wgbs 194804760  #queen 
total_called_sites <- sum(gc_calledsites) #497298407.  ##queen 491381680

####totalnumbermethylatedcites
gc_methylatedsites <- data$methylated_bases ##5623387 pupa 2.8%
total_called_sites <- sum(gc_methylatedsites) #13316233.  3%  ###13124213 $3%



# Only keep positions with 3 or more reads mapped
data <- data[total_bases >= 3] #3352854 


# read in the ONT data (one line per cytosine)
# before reading in the ONT file, make sure you replace spaces by tabs in the file
# in bash: tr " " "\t" < input_file > input_file_tables

# please note that ONT output files often change format, make sure the columns are correctly named below

data[, c("motifs_in_group", "group_sequence"):= NULL]

# add 1 to positions so that they are one-based instead of zero-based
data[, pos:=pos+1]
data


#### this part of the script runs a binomial test on every gene to infer its methylation state

gff3 <- fread("/Users/taniachavarriapizarro/Desktop/Methylation Pogo/Pcal3.1.clean.submitted.polished2.gff3", h=F)
gff3
gff3[, c("V2", "V6","V8") := NULL] # delete unnecessary columns
gff3

colnames(gff3) <- c("chr", "type", "start", "stop", "strand", "ID")
gff3
gff3 <- gff3[type == "gene"] # retain only gene information
gff3[, GeneName := strsplit(ID, 'ID=', perl=T)[[1]][2], by=1:nrow(gff3)]# extract gene name from ID column
fwrite(gff3, "genes.gff3", sep="\t")
gff3[, c("ID") := NULL] # delete unnecessary columns
gff3
fwrite(gff3, "genes.gff3", sep="\t")
gff3 <- fread("/Users/taniachavarriapizarro/Desktop/Methylation Pogo/genes.gff3", h=T)


# Create output file for Gene cytosine methylation counts
header <- data.table("Gene" = character(), "Chromosome" = character(), "start" = numeric(), "stop" = numeric(), "strand" = character(), 
                     "Number_CG_total_reads" = numeric(), "Number_methylated_CG_reads" = numeric())
fwrite(header, file = "gene_methylationPcq_6b.CpG_report.txt", sep="\t", append=F)

# loop on each gene to count number of methylated and unmethylated cytosines in each context 
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentstart <- gff3$start[i]
  Currentstop <- gff3$stop[i]
  Currentstrand <- gff3$strand[i]
  
  # extract gene methylation for current gene
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentstart)&(pos<=Currentstop)]
  
  # check there is something
  if (nrow(GeneMethylation)>0) {
    # Sum the number of methylated and total reads
    sumMeth <- sum(GeneMethylation$methylated_bases)
    sumtotal <- sum(GeneMethylation$total_bases)
    GeneData <- data.table(as.character(CurrentGene), as.character(CurrentChr), as.numeric(Currentstart), as.numeric(Currentstop), as.character(Currentstrand), sumtotal, sumMeth)
    
    # add gene info to output file
    fwrite(GeneData, file="gene_methylationPcq_6b.CpG_report.txt", sep="\t", append=TRUE)
  }
}


# load Gene data
GeneData <- fread("gene_methylationPcq_6b.CpG_report.txt", h=T)
GeneData
head(GeneData)

# computation of average CG gene methylation across all genes
proba_methylated_CG <- sum(GeneData$Number_methylated_CG_reads, na.rm=T)/sum(GeneData$Number_CG_total_reads, na.rm=T)

# binom test on each gene to test for significant CG methylation compared to other genes
###GeneData[,  CG_Binom_pvalue_greater := if (Number_CG_total_reads > 0) {binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "greater")$p.value} else {as.numeric(NA)}, by = 1:nrow(GeneData)]

GeneData[, CG_Binom_pvalue_greater := 
           if (Number_CG_total_reads > 0 & Number_CG_total_reads >= Number_methylated_CG_reads) {
             binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "greater")$p.value
           } else {
             as.numeric(NA)
           }, 
         by = 1:nrow(GeneData)
]



###GeneData[,  CG_Binom_pvalue_less := if (Number_CG_total_reads > 0) {binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "less")$p.value} else {as.numeric(NA)}, by = 1:nrow(GeneData)]

GeneData[, CG_Binom_pvalue_less := 
           if (Number_CG_total_reads > 0 & Number_CG_total_reads >= Number_methylated_CG_reads) {
             binom.test(Number_methylated_CG_reads, Number_CG_total_reads, proba_methylated_CG, alternative = "less")$p.value
           } else {
             as.numeric(NA)
           }, 
         by = 1:nrow(GeneData)
]


GeneData

# adjust p-values for multiple tests with Benjamini & Hochberg (1995) ("BH")
GeneData[,  CG_Binom_pvalue_greater_adjusted := p.adjust(CG_Binom_pvalue_greater, "BH")]
GeneData[,  CG_Binom_pvalue_less_adjusted := p.adjust(CG_Binom_pvalue_less, "BH")]
GeneData

# Classify gene methylation
Methylation <- function(u,x) {
  # u the adjusted p-value for CG methylation to be higher than expected by chance
  # x the adjusted p-value for CG methylation to be lower than expected by chance
  # outputs mCG if gene is significantly CG body methylated (u < 0.05)
  if ((!is.na(u))&(u < 0.05)) {
    # CG body methylated gene
    return("gbM")
  } else if ((!is.na(x))&(x < 0.05)) {
    # CG unmethylated gene
    return("UM")
  } else if ((!is.na(x))&(x > 0.05)&(!is.na(u))&(u > 0.05)) {
    # CG intermediately methylated gene
    return("IM")
  } else {
    return("NA")
  }
}
GeneData[, GeneMethylationCall := mapply(Methylation, CG_Binom_pvalue_greater_adjusted, CG_Binom_pvalue_less_adjusted)]
GeneData
table(GeneData$GeneMethylationCall)
GeneData <-  fwrite(GeneData, file="gene_CallmethylationPcq_6b.CpG_report.txt", sep="\t", append=TRUE)
# Note that in some papers intermediately methylated genes 'IM' are merged together with unmethylated genes 'UM' 



###################################################################################################
#          Plotting Gene methylation metaplot and studying promoter methylation
###################################################################################################
# loop on all genes to compute methylation proportions on windows upstream of the gene, downstream and inside the gene.
header <- list("SumMethylatedReads", "TotalReads", "window", "GeneName")
fwrite(header, file="Gene_windows_Pcq_6b.CpG_report.txt", sep="\t", append=FALSE)
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases))]
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows_Pcq_6b.CpG_report.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows <- fread("Gene_windows_Pcq_12b.CpG_report.txt", h=T)
Windows

# Summary over all genes
SummaryWindows <- Windows[TotalReads > 0, .(sum(SumMethylatedReads)/sum(TotalReads), binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][1], binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][2]), by=c("window")]

names(SummaryWindows) <- c("window", "Average_Methylation", "CI_low", "CI_high")
fwrite(SummaryWindows, file="SummaryWindows_Pcq_6b.CpG_report.txt", sep="\t", append=TRUE)


# plot metaplot
jpeg(f="windows_weighteg_methylation_Genes.jpeg", width=8, height=4, units = "in", res=500, pointsize = 18, quality=100)
ggplot(data=SummaryWindows, aes(x = window, y = Average_Methylation)) + 
  theme_bw(base_size=20) +
  geom_point(aes(x = window, y = Average_Methylation), size = 1) +  
  geom_line(aes(x = window, y = Average_Methylation), linewidth=0.5) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), linetype=1, size=0.2, alpha=0.1, show.legend = FALSE) +
  labs(x=paste("position", sep=''), y="Weighted methylation") + 
  guides(col=guide_legend(title="Methylation type:")) +
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000'))

dev.off()



#####comparinglarva_pupa_queen_worker

SummaryWindows <- fread("Summarywindows_5L_4P_Q.txt", h=T)

ggplot(data = SummaryWindows, aes(x = window, y = Average_Methylation, group = type, color = type, fill = type)) +
  theme_bw(base_size = 20) +
  geom_point(size = 1, shape = 20) +  
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), linetype = 1, size = 0.2, alpha = 0.1) +
  labs(x=paste("position", sep=''), y="Weighted methylation") + 
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000')) +
  scale_color_manual(values = c("#ffcc33","#993300","#ff0000","#ff6600")) +
  scale_fill_manual(values = c(alpha("#ffcc33", 0.1),alpha("#993300", 0.1),alpha("#ff0000", 0.1),alpha("#ff6600", 0.1))) +
  theme(legend.position = "right")



# gene promoter methylation (use window 20 which corresponds to 200bp upstream of Transcription start site)
# Nota Bene: In animals, promoters are usually longer than the typical plant 200pb promoter, 
# you can take more windows than just window 20 (for example, window 15 to 20)
Windows <- Windows[window==20]
Windows[, window := NULL]

# compute promoter average proportion of methylated CG, CHG and CHH sites
promoter_proba_methylated_CG <- sum(Windows$SumMethylatedReads, na.rm=T)/sum(Windows$TotalReads, na.rm=T)

# compare promoter methylation to genome average using a binomial test
###BinomTest_greater <- function(x, n, p){ if(n<1) { return(as.double(NA))} else { return(binom.test(x, n, p, alternative = "greater")$p.value) }}

# Define a function that handles the checks and ensures numeric output
BinomTest_greater <- function(x, n, p) {
  if (is.numeric(x) && is.numeric(n) && n > 0 && n >= x) {
    return(as.numeric(binom.test(x, n, p, alternative = "greater")$p.value))
  } else {
    return(as.numeric(NA))
  }
}

# Apply the function to the data
Binom_pvalues_greater_CG <- Windows[, .(
  BinomPValueGreater = BinomTest_greater(SumMethylatedReads, TotalReads, promoter_proba_methylated_CG)
), by = 1:nrow(Windows)]




Binom_pvalues_greater_CG <- Windows[, BinomTest_greater(SumMethylatedReads, TotalReads, promoter_proba_methylated_CG), by = 1:nrow(Windows)]

# adjust p-values
Windows[, Padjusted_greater_CG :=  p.adjust(Binom_pvalues_greater_CG$V1, "BH")]

# inference of CG promoter methylated genes
Windows[, promoter_CG_methylated := if((!is.na(Padjusted_greater_CG))&(Padjusted_greater_CG <= 0.05)){"yes"}else{"no"}, by = 1:nrow(Windows)]
fwrite(Windows, "Promoter_methylation_PcC1_1_.CpG.txt", sep="\t")






setwd("/Users/taniachavarriapizarro/Desktop/New fig methylome")

# data format must be one line per cytosine:
# chromosome	position	strand(+/-)	mc_class=context(CG/CHG/CHH)	methylated_bases(number_methylated_C_reads)	total_bases(total_number_reads_mapped)
# read in the data
methylation_data2 <- fread("Pcq_8w_bedGraph.tsv", h=F)
methylation_data2

names(methylation_data2) <- c("chromosome", "start","end", "methylated_frequency")
methylation_data2

methylation_data2 <- methylation_data2[order(methylation_data2$chromosome), ]
methylation_data2

# Change chromosome names using the mapping
methylation_data2 <- methylation_data2 %>%
  mutate(chrom = factor(chromosome, levels = names(chromosome_name_mapping), labels = chromosome_name_mapping))

# List of specific chromosomes to include in the single plot
chromosomes_to_plot <- c("1","2", "3","4","5", "6","7","8","9","10","11","12","13","14","15","16") # Replace with your chromosome names

# Filter the cytoband data and scaffold lengths for only the specified chromosomes

combined_methylation2 <- methylation_data2 %>% filter(chrom %in% chromosomes_to_plot)

fwrite(combined_methylation2, "/Users/taniachavarriapizarro/Desktop/WGBSdataFig1.txt", sep="\t")

combined_methylation2 <- fread("/Users/taniachavarriapizarro/Desktop/WGBSdataFig1.txt", sep="\t")

####loading direct data for plot

combined_methylation <- fread("/Users/taniachavarriapizarro/Desktop/NanoporedataFig1.txt", sep="\t")
combined_methylation2 <- fread("/Users/taniachavarriapizarro/Desktop/WGBSdataFig1.txt", sep="\t")
te_islands <- fread("/Users/taniachavarriapizarro/Desktop/TEdataFig2.txt", sep="\t") 

# Load WGBS methylation data
wbgs_data <- combined_methylation2
head(wbgs_data)

wbgs_data$methylated_frequency <- as.numeric(wbgs_data$methylated_frequency)
wbgs_data <- wbgs_data %>%
  mutate(methylated_frequency = methylated_frequency/100) 


# Verify the range to ensure filtering
print(range(wbgs_data$methylated_frequency, na.rm = TRUE))

wgbs_data <- wbgs_data %>%
  mutate(start = pmin(start, end),
         end = pmax(start, end))

nanopore_data <- combined_methylation


library(dplyr)
library(ggplot2)

# Define the window size
window_size <- 200000

# Combine datasets with an indicator column for the method
combined_data <- bind_rows(
  nanopore_data %>%
    mutate(method = "Nanopore"),
  wgbs_data %>%
    mutate(method = "WGBS")
)



#####No mirror

combined_binned <- combined_data %>%
  mutate(
    window_start = (start %/% window_size) * window_size,
    window_end = window_start + window_size - 1
  ) %>%
  group_by(method, chrom, window_start, window_end) %>%
  summarize(
    average_methylated_frequency = mean(methylated_frequency, na.rm = TRUE),
    CI_low = mean(methylated_frequency, na.rm = TRUE) - 1.96 * sd(methylated_frequency, na.rm = TRUE) / sqrt(n()),
    CI_high = mean(methylated_frequency, na.rm = TRUE) + 1.96 * sd(methylated_frequency, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()



plot <- ggplot(combined_binned, aes(x = window_start, y = average_methylated_frequency, group = interaction(chrom, method))) +
  geom_line(aes(color = method)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = method), alpha = 0.2) +
  labs(
    x = "Genomic Window Start (bp)",
    y = "Average Methylation Frequency",
    title = "Methylation Frequency Across Methods",
    color = "Method",
    fill = "Method"
  ) +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1")) +
  scale_fill_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1"))+
  facet_wrap(~ chrom, scales = "free_x")


####chromosome linked
head(combined_binned)

# Step 1: Group and summarize to get max_window
chrom_lengths_step1 <- combined_binned %>%
  group_by(chrom) %>%
  summarize(max_window = max(window_end, na.rm = TRUE))

# Step 2: Add cumulative start positions
chrom_lengths_step2 <- chrom_lengths_step1 %>%
  mutate(cum_start = lag(cumsum(max_window), default = 0))

# Step 3: Select desired columns
chrom_lengths <- chrom_lengths_step2 %>%
  dplyr::select(chrom, cum_start)



# Merge cumulative starts with the original data
combined_binned <- combined_binned %>%
  left_join(chrom_lengths, by = "chrom") %>%
  mutate(
    cum_window_start = window_start + cum_start,
    cum_window_end = window_end + cum_start
  )

# Step 2: Calculate chromosome label positions
chrom_labels <- chrom_lengths_step2 %>%
  mutate(label_position = cum_start + max_window / 2)

# Step 3: Plot with continuous x-axis
plot <- ggplot(combined_binned, aes(x = cum_window_start, y = average_methylated_frequency, group = interaction(chrom, method))) +
  geom_line(aes(color = method), size = 1.2) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = method), alpha = 0.2) +
  labs(
    x = "Genomic Position (bp)",
    y = "Average Methylation Frequency",
    title = "Methylation Frequency Across Methods",
    color = "Method",
    fill = "Method"
  ) +
  theme_minimal() +
  scale_x_continuous(
    labels = scales::comma,
    breaks = chrom_labels$cum_start, # Add breaks at chromosome boundaries
    sec.axis = dup_axis(
      labels = chrom_labels$chrom, # Add chromosome names as secondary labels
      breaks = chrom_labels$label_position
    )
  ) +
  scale_color_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1")) +
  scale_fill_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1"))

# Print the plot
print(plot)

####just methylation frequency with Nanopore

# Filter the data for Nanopore only
nanopore_data <- combined_binned %>%
  filter(method == "Nanopore")

# Create the scatterplot
nanopore_plot <- ggplot(nanopore_data, aes(x = cum_window_start, y = average_methylated_frequency)) +
  geom_point(color = "#31A354", size = 2, alpha = 0.6) + # Scatter points for Nanopore
  labs(
    x = "Genomic Position (bp)",
    y = "Average Methylation Frequency",
    title = "Methylation Frequency for Nanopore"
  ) +
  theme_minimal() +
  scale_x_continuous(
    labels = scales::comma,
    breaks = chrom_labels$cum_start, # Add breaks at chromosome boundaries
    sec.axis = dup_axis(
      labels = chrom_labels$chrom, # Add chromosome names as secondary labels
      breaks = chrom_labels$label_position
    )
  )

# Print the plot
print(nanopore_plot)

####now with TE islands


head(chrom_lengths)
head(te_islands)

te_islands <- te_islands %>%
  mutate(joiner = as.integer(joiner))
unique(te_islands$joiner)
unique(chrom_lengths$chrom)

te_islands <- te_islands %>%
  mutate(joiner = as.character(joiner))
chrom_lengths <- chrom_lengths %>%
  mutate(chrom = as.character(chrom))

te_islands <- te_islands %>%
  left_join(chrom_lengths, by = c("joiner" = "chrom"))

head(te_islands)


# Step 1: Map scaffold and joiner to chromosome
te_islands <- te_islands %>%
  left_join(chrom_lengths, by = c("joiner" = "chrom")) %>%
  mutate(
    cum_window_start = start + cum_start,
    cum_window_end = stop + cum_start
  )

# Step 2: Merge TE island information with Nanopore data (optional if overlay only)
nanopore_data <- combined_binned %>%
  filter(method == "Nanopore")



###now just TE islands

# Filter TE islands data to include only relevant categories
te_islands_filtered <- te_islands %>%
  filter(TEisland == "TEisland") # Retain only "TEisland" rows

# Create the combined plot
combined_plot <- ggplot() +
  # TE islands as bar-like areas
  geom_bar(
    data = te_islands_filtered,
    aes(x = cum_window_start, y = value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Methylation frequency as scatter points
  geom_point(
    data = nanopore_data,
    aes(x = cum_window_start, y = average_methylated_frequency),
    color = "#31A354", size = 2, alpha = 0.6
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  labs(
    x = "Genomic Position (Mb)",
    y = "Methylation Frequency / TE Coverage",
    fill = "TE Island",
    title = "TE Island Data and Methylation Frequency"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines")
  ) +
  scale_fill_manual(values = c("TEisland" = "orange"), na.value = "gray")

# Print the plot
print(combined_plot)

####
# Filter TE islands to only include relevant data
te_islands_filtered <- te_islands %>%
  filter(TEisland == "TEisland") %>%
  mutate(value = -value) # Keep TE islands below 0

# Create the plot with shared origin at zero
aligned_plot <- ggplot() +
  # TE islands as bar-like areas below 0
  geom_bar(
    data = te_islands_filtered,
    aes(x = cum_window_start, y = value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Methylation frequency as scatter points above 0
  geom_point(
    data = nanopore_data,
    aes(x = cum_window_start, y = average_methylated_frequency),
    color = "#31A354", size = 2, alpha = 0.6
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  scale_y_continuous(
    name = "Methylation Frequency (+) / TE Coverage (-)",
    labels = abs, # Show absolute values on the y-axis
    expand = expansion(mult = c(0.1, 0.1)) # Add space above and below the plot
  ) +
  labs(
    x = "Genomic Position (Mb)",
    fill = "TE Island",
    title = "TE Islands and Methylation Frequency with Shared Zero Origin"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines")
  ) +
  scale_fill_manual(values = c("TEisland" = "orange"), na.value = "gray")

# Print the plot
print(aligned_plot)


####

# Filter TE islands to only include relevant data
te_islands_filtered <- te_islands %>%
  filter(TEisland == "TEisland") # Retain only "TEisland" rows

# Create the plot with both datasets as bar-like areas
aligned_plot <- ggplot() +
  # TE islands as bar-like areas
  geom_bar(
    data = te_islands_filtered,
    aes(x = cum_window_start, y = value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Methylation frequency as bar-like areas
  geom_bar(
    data = nanopore_data,
    aes(x = cum_window_start, y = average_methylated_frequency, fill = "Methylation"),
    stat = "identity",
    alpha = 0.5,
    color = NA
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  scale_y_continuous(
    name = "TE Coverage and Methylation Frequency",
    expand = expansion(mult = c(0.1, 0.1)) # Add space above the plot
  ) +
  labs(
    x = "Genomic Position (Mb)",
    fill = "Legend",
    title = "TE Islands and Methylation Frequency as Bar-Like Areas"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines")
  ) +
  scale_fill_manual(
    values = c("TEisland" = "orange", "Methylation" = "#31A354"),
    na.value = "gray"
  )

# Print the plot
print(aligned_plot)

####mirrorplot

# Filter TE islands to only include relevant data
te_islands_filtered <- te_islands %>%
  filter(TEisland == "TEisland") # Retain only "TEisland" rows

# Invert methylation frequency for the negative side
nanopore_data_mirrored <- nanopore_data %>%
  mutate(average_methylated_frequency = -average_methylated_frequency) # Make methylation frequency negative

# Create the mirrored plot
mirrored_plot <- ggplot() +
  # TE islands as bar-like areas on the positive side
  geom_bar(
    data = te_islands_filtered,
    aes(x = cum_window_start, y = value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Methylation frequency as bar-like areas on the negative side
  geom_bar(
    data = nanopore_data_mirrored,
    aes(x = cum_window_start, y = average_methylated_frequency, fill = "Methylation"),
    stat = "identity",
    alpha = 0.5,
    color = NA
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  scale_y_continuous(
    name = "TE Coverage (+) / Methylation Frequency (-)",
    labels = abs, # Show absolute values on the y-axis
    expand = expansion(mult = c(0.1, 0.1)) # Add space above and below the plot
  ) +
  labs(
    x = "Genomic Position (Mb)",
    fill = "Legend",
    title = "Mirrored TE Islands and Methylation Frequency"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines")
  ) +
  scale_fill_manual(
    values = c("TEisland" = "orange", "Methylation" = "#31A354"),
    na.value = "gray"
  )

# Print the plot
print(mirrored_plot)

####scaling better
# Scaling factors
methylation_scaling_factor <- 5  # Adjust for methylation frequency
te_island_scaling_factor <- 2     # Adjust for TE islands

# Scale TE islands and methylation frequency
te_islands_scaled <- te_islands %>%
  filter(TEisland == "TEisland") %>%
  mutate(scaled_value = value * te_island_scaling_factor) # Scale TE islands

nanopore_data_scaled <- nanopore_data %>%
  mutate(
    scaled_methylation = -average_methylated_frequency * methylation_scaling_factor # Scale and invert methylation
  )

# Create the mirrored plot with both datasets scaled
mirrored_plot <- ggplot() +
  # Scaled TE islands as bar-like areas on the positive side
  geom_bar(
    data = te_islands_scaled,
    aes(x = cum_window_start, y = scaled_value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Scaled methylation frequency as bar-like areas on the negative side
  geom_bar(
    data = nanopore_data_scaled,
    aes(x = cum_window_start, y = scaled_methylation, fill = "Methylation"),
    stat = "identity",
    alpha = 0.5,
    color = NA
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  scale_y_continuous(
    name = "TE Coverage (+) / Methylation Frequency (-)",
    sec.axis = sec_axis(~ . / methylation_scaling_factor, name = "Methylation Frequency"), # Secondary axis for methylation
    labels = abs, # Show absolute values on the y-axis
    expand = expansion(mult = c(0.1, 0.1)) # Add space above and below the plot
  ) +
  labs(
    x = "Genomic Position (Mb)",
    fill = "Legend",
    title = "Mirrored TE Islands and Scaled Methylation Frequency"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines"),
    plot.margin = margin(10, 10, 10, 10), # Add space around the plot
    aspect.ratio = 0.3 # Make the plot wider
  ) +
  scale_fill_manual(
    values = c("TEisland" = "orange", "Methylation" = "#31A354"),
    na.value = "gray"
  )

# Print the plot
print(mirrored_plot)

####now with TE as bars and methylation and points

# Scaling factors
methylation_scaling_factor <- 5  # Adjust for methylation frequency
te_island_scaling_factor <- 2    # Adjust for TE islands

# Scale TE islands and methylation frequency
te_islands_scaled <- te_islands %>%
  filter(TEisland == "TEisland") %>%
  mutate(scaled_value = value * te_island_scaling_factor) # Scale TE islands

nanopore_data_scaled <- nanopore_data %>%
  mutate(
    scaled_methylation = -average_methylated_frequency * methylation_scaling_factor # Scale and invert methylation
  )

# Create the mirrored plot
mirrored_plot <- ggplot() +
  # Scaled TE islands as bar-like areas on the positive side
  geom_bar(
    data = te_islands_scaled,
    aes(x = cum_window_start, y = scaled_value, fill = TEisland),
    stat = "identity",
    alpha = 0.7,
    color = NA
  ) +
  # Scaled methylation frequency as scatter points on the negative side
  geom_point(
    data = nanopore_data_scaled,
    aes(x = cum_window_start, y = scaled_methylation, color = "Methylation"),
    size = 2,
    alpha = 0.6
  ) +
  scale_x_continuous(
    labels = function(x) x / 1000000, # Genomic positions in Mb
    breaks = chrom_labels$cum_start,
    sec.axis = dup_axis(
      labels = chrom_labels$chrom,
      breaks = chrom_labels$label_position
    )
  ) +
  scale_y_continuous(
    name = "TE Coverage (+) / Methylation Frequency (-)",
    sec.axis = sec_axis(~ . / methylation_scaling_factor, name = "Methylation Frequency"), # Secondary axis for methylation
    labels = abs, # Show absolute values on the y-axis
    expand = expansion(mult = c(0.1, 0.1)) # Add space above and below the plot
  ) +
  labs(
    x = "Genomic Position (Mb)",
    fill = "Legend",
    color = "Legend",
    title = "Mirrored TE Islands and Methylation Frequency (Scatter)"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(-1, "lines"),
    plot.margin = margin(10, 10, 10, 10), # Add space around the plot
    aspect.ratio = 0.3 # Make the plot wider
  ) +
  scale_fill_manual(
    values = c("TEisland" = "orange"),
    na.value = "gray"
  ) +
  scale_color_manual(
    values = c("Methylation" = "#31A354")
  )

# Print the plot
print(mirrored_plot)


#####now with TE families

methylation_data3 <- fread("/Users/taniachavarriapizarro/Desktop/TEqueenIDmethylfreq.tsv", h=F)
methylation_data3

names(methylation_data3) <- c("chromosome", "start","end","called_sites",
                              "called_sites_methylated", "methylated_frequency", "strand")
methylation_data3

methylation_data3 <- methylation_data3[order(methylation_data3$chromosome), ]
methylation_data3

# Change chromosome names using the mapping
methylation_data3 <- methylation_data3 %>%
  mutate(chrom = factor(chromosome, levels = names(chromosome_name_mapping), labels = chromosome_name_mapping))

# List of specific chromosomes to include in the single plot
chromosomes_to_plot <- c("1","2", "3","4","5", "6","7","8","9","10","11","12","13","14","15","16") # Replace with your chromosome names

# Filter the cytoband data and scaffold lengths for only the specified chromosomes

combined_methylation3 <- methylation_data3 %>% filter(chrom %in% chromosomes_to_plot)



# Load TE methylation data
TE_data <- combined_methylation3
head(TE_data)

TE_data$methylated_frequency <- as.numeric(TE_data$methylated_frequency)

TE_data$chrom <- as.numeric(TE_data$chrom)

# Verify the range to ensure filtering
print(range(TE_data$methylated_frequency, na.rm = TRUE))

TE_data <- TE_data %>%
  mutate(start = pmin(start, end),
         end = pmax(start, end))

nanopore_data <- combined_methylation

library(dplyr)
library(ggplot2)

# Define the window size
window_size <- 200000

# Combine datasets with an indicator column for the method
combined_data2 <- bind_rows(
  nanopore_data %>%
    mutate(method = "Nanopore"),
  TE_data %>%
    mutate(method = "TE")
)


###methylation depth

# Load ggplot2 library
library(ggplot2)

# Calculate methylation depth for each row

# Calculate methylation coverage
total_called_sites <- sum(methylation_data$called_sites > 0)
methylated_sites <- sum(methylation_data$called_sites_methylated > 0)


# Calculate mean methylation depth across all rows

methylation_data$methylation_depth <- methylation_data$methylated_frequency/methylation_data$called_sites_methylated
methylation_data <- methylation_data %>%
  mutate(methylation_depth = ifelse(is.nan(methylation_depth), 0, methylation_depth))
range(methylation_data$methylation_depth)
range(methylation_data$called_sites_methylated)

###methylation frequency = methylation coverage × methylation depth.

# Step 2: Plot with a wider confidence interval
ggplot(methylation_data, aes(x =methylated_frequency, y = methylation_depth)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#31A354", 
              fill = "#31A354", alpha = 0.2, level = 0.99, linetype = "dashed") +
  labs(
    title = "Methylation depth vs methylation coverage",
    x = "methylation_coverage",
    y = "methylation_depth"
  ) +
  theme_bw() +
  theme(text = element_text(size = 16))


ggplot(methylation_data, aes(x =methylated_frequency, y = methylation_depth)) +
  geom_point(alpha = 0.3, color = "#756BB1")+
  labs(
    title = "Methylation depth vs methylation coverage",
    x = "methylation_coverage",
    y = "methylation_depth"
  ) +
  theme_bw() +
  theme(text = element_text(size = 16))




###load WGBS

methylation_data2 <- fread("Pcq_8w_bedGraph.tsv", h=F)
methylation_data2


names(methylation_data2) <- c("chromosome", "start","end", "methylated_frequency")
methylation_data2

methylation_data2 <- methylation_data2[order(methylation_data2$chromosome), ]
methylation_data2


# Load WGBS methylation data
wbgs_data <- methylation_data2
head(wbgs_data)


# If WGBS 
wbgs_data <- wbgs_data %>%
  mutate(methylated_frequency_scaled = methylated_frequency/100) %>%
  filter(methylated_frequency_scaled > 0.25)

# Verify the range to ensure filtering
print(range(wbgs_data$methylated_frequency_scaled, na.rm = TRUE))

methylation_data3 <- fread("Pcq_8w.CpG_report.tsv", h=F)
methylation_data3


names(methylation_data3) <- c("chromosome", "start","strand", "called_sites_methylated","called_sites",
                              "group","group_sequence")
methylation_data3

methylation_data3 <- methylation_data3[order(methylation_data3$chromosome), ]
methylation_data3

methylation_data3 <- methylation_data3 %>%
  mutate(end = start + 1)

methylation_data3 <- methylation_data3 %>%
  select(chromosome, start, end, everything())

head(methylation_data3)
methylation_data3 <- methylation_data3 %>%
  mutate(
    end = start + 1,
    methylated_frequency = ifelse(is.nan(called_sites_methylated / called_sites), 0, called_sites_methylated / called_sites)
  ) %>%
  select(chromosome, start, end, -strand, called_sites_methylated, called_sites, methylated_frequency, -group, -group_sequence)


head(methylation_data3)



fwrite(methylation_data3, "/Users/taniachavarriapizarro/Desktop/Pcq_8w.CpG_totalfreq.txt", sep="\t")


methylation_data <- fread("/Users/taniachavarriapizarro/Desktop/CNA6freq.tsv", h=T)
methylation_data <- subset(methylation_data, called_sites >= 3)

head(methylation_data)
str(methylation_data)


# List of specific chromosomes to include in the single plot
chromosomes_to_plot <- c("LG1", "LG2", "LG3","LG4","LG5","LG6",
                         "LG7","LG8","LG9","LG10","LG11","LG12",
                         "LG13","LG14","LG15","LG16")


methylation_data <- methylation_data %>%
  filter(chromosome %in% chromosomes_to_plot)


###WGBS
methylation_data2 <- fread("/Users/taniachavarriapizarro/Desktop/Pcq_8w.CpG_totalfreq.tsv", h=T)
methylation_data2 <- subset(methylation_data2, called_sites >= 3)

head(methylation_data2)
str(methylation_data2)


# List of specific chromosomes to include in the single plot
chromosomes_to_plot <- c("LG1", "LG2", "LG3","LG4","LG5","LG6",
                         "LG7","LG8","LG9","LG10","LG11","LG12",
                         "LG13","LG14","LG15","LG16")


methylation_data2 <- methylation_data2 %>%
  filter(chromosome %in% chromosomes_to_plot)

head(methylation_data2)

library(dplyr)
library(ggplot2)
library(mgcv)  # For GAM smoothing

# Example data for Nanopore and WGBS datasets
nanopore_data <- methylation_data
nanopore_data <- nanopore_data %>%
  mutate(method = "Nanopore")
head(nanopore_data)

wgbs_data <- methylation_data2
wgbs_data <- wgbs_data %>%
  mutate(method = "WGBS")
head(wgbs_data)


# Combine the datasets
combined_data <- bind_rows(nanopore_data, wgbs_data)
fwrite(combined_data, "/Users/taniachavarriapizarro/Desktop/combined_data.txt", sep="\t")


# Step 1: Calculate mean and standard deviation by coverage for each dataset
methylation_stats <- combined_data %>%
  group_by(method, called_sites) %>%
  summarize(
    mean_methylated_frequency = mean(methylated_frequency, na.rm = TRUE),
    sd_methylated_frequency = sd(methylated_frequency, na.rm = TRUE)
  )

# Step 2: Plot with a wider confidence interval for both datasets
ggplot(combined_data, aes(x = called_sites, y = methylated_frequency, color = method, fill = method)) +
  geom_smooth(
    method = "lm",
    alpha = 0.2,
    level = 0.99,
    linetype = "dashed"
  ) +
  labs(
    title = "Mean Methylated Frequency by Coverage for Nanopore vs. WGBS",
    x = "Coverage (X)",
    y = "Methylated Frequency"
  ) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_color_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1")) +
  scale_fill_manual(values = c("Nanopore" = "#31A354", "WGBS" = "#756BB1")) +
  scale_x_continuous(limits = c(0, NA)) +   # Ensures x-axis starts at 0
  scale_y_continuous(limits = c(0, NA))     # Ensures y-axis starts at 0

pdf(file = "/Users/taniachavarriapizarro/Desktop/Combined_Circos_plot_Nanopore_vs_WGBS.pdf", 
    width = 6, height = 6,compress = TRUE)

dev.off()


#####TEannotation

library(data.table)
library(ggplot2)


setwd("/Users/taniachavarriapizarro/Desktop/")

#Plotting methylation average introns vs exons 
data <- fread("Pcal3.1.TE.clean.final.tsv", header=FALSE, fill=TRUE)

head(data)


# remove unnecessary columns
data[, c("V1","V2", "V3", "V4","V12","V13","V14","V15","V16"):= NULL]
data[, c("V9"):= NULL]

# namecolumns

names(data) <- c("chrom","pos", "end", "left", "matchingrepeat", "repeatclass/family")
head(data)





fwrite(data, "Pcal.3.1.TEcleanfinal.tsv", sep="\t")


data <- fread("PcCA4queen_TE.tsv", header=F)
head(data)

# remove unnecessary columns
data[, c("V8","V9", "V10", "V11"):= NULL]


# namecolumns

names(data) <- c("chrom","pos", "end", "motif","called_sites", "methylated_sites","methylated_frequency","left","TE_ID" ,"TEfamily")
head(data)



library(tidyverse)

# Change chromosome names using the mapping
methylation_data <- data %>%
  mutate(chrom = factor(chrom, levels = names(chromosome_name_mapping), labels = chromosome_name_mapping))

head(methylation_data)
###calculation methylation freq average and sd

grouped_data <- methylation_data %>%
  group_by(TEfamily) %>%
  summarise(mean_methylation = mean(methylated_frequency),
            se_methylation = sd(methylated_frequency) / sqrt(n()))

head(grouped_data)
fwrite(grouped_data, "MeanMethylbytypeTE.tsv", sep="\t")


grouped_data <- fread("MeanMethylbytypeTE.tsv", h=T)
head(grouped_data)

ggplot(grouped_data, aes(x = type, y = mean_methylation, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean_methylation - se_methylation, ymax = mean_methylation + se_methylation),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "Average Methylation by Region queen P.californicus",
       x = "Region",
       y = "Methylation Frequency") +
  scale_fill_viridis_d(direction = -1)+
  theme_minimal()


####phylogenetic tree DNMT1

install.packages(c("ape", "phangorn", "Biostrings"))
library(ape)
library(phangorn)
library(Biostrings)

# Install required packages if not already installed
install.packages(c("ape", "phangorn", "Biostrings", "msa"))

# Load necessary libraries
library(ape)
library(phangorn)
library(Biostrings)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")
library(msa)

# Define DNA sequences for the six species
sequences <- DNAStringSet(c(
  "Apis mellifera" = paste(
    "ATGAAGTCTAATACATCAAAATTTAAACAAACATCTATAATTGATATGTTTTCTCAACAACTTGCTAAAAAAAAGCAAATGAACAATGATGCAACTAATAAAGAAGACATACAAAAATTACAAAATATAAAAACAGAAAA",
    "TAATATATTGCAAAAATGTGAACAAAATAAAGAAAATCATGATCCAAAAACAGATGAGGAATATTTAGATAAGCACATTGAAAAAAAAGCTAAACTTACAGATGTTAATATATTACCAACACCTCAAATTAAAAATATAA",
    sep = ""),
  "Bombus impatiens" = paste(
    "ATGAAGTCTGTTGCATCCAAATCTAATTATCAACCATCTATAGTTGACATGTTTTATAAACAACTTGCTATGGGAAAACAAATAAGTTATGATGCAAATAACGATAGAGAATTACAGAATTCACAAAACATTAAAACAGAA",
    "AATAACATATCTTTGAAAAGTGAGCAAGATAAAGAAAACCATTATCCAGGAACAAATGAGGACTATTCTGATAAACATGTTGAAAAAAAAGTCAAGTTTAAAGATACTACAATTTTACCAACATCTAAAATTAAAAGTACA",
    sep = ""),
  "Camponotus floridianus" = paste(
    "ATGATCGATAAATACGCGCGCCACAACCCACATTCTCTCGTAATCTGTAAAAAACCGGACGACGTTTCATCGATCACGACGTCATTCACGAGAATTATGCCACCCACTGTGCCAGAGGACGAACCAGGTTCCAGCAATAGA",
    "AGAATGTCTCGTAGCAGAAAATTTGTTGATATAGAAGAATCAATTGAAGCAAATGGTTCGACGGCTAAAGCAACTCCATCCGTATTGGATTTATTTGCTATACAATCAGCAAAAAGAAAAAAGTTAAACAATGAAACAAAC",
    sep = ""),
  "Oocera biroi" = paste(
    "ATGTCATATTGCGGAATTATGCCACCTACTGTACCAAGCGAACCAGACTCCAGTAGCAGAAGGATGTCCCGTAGTCGGAAGCATGCCAGCATGGAGGAATCCGTTGAAGCGGGTGATCTAACAGCGGACACAACTCCGTCC",
    "GCGCTGGACTCATTTTCCAAGCAACTGGTGAAAAGGAAGAAGATGAATACCGAGGCGAACGGCGAGGCAAAATCCGTGTCTGATCCTGAGTACGAAGTAAAAATGGAGTCGGACGACAATACCACCATAGAGTTCGACAAG",
    sep = ""),
  "Pogonomyrmex barbatus" = paste(
    "ATGAATCGATACGCGGATCCCTTTTCACATTTCCCAACATCTTTGCTAGCCGCGATGGCTCGTGGTACAGATCGTGCGATCGAACGTAACATTGCTAGCCGGCGAGTAAGAGATGAGTCACGAGAAAATGAACCAGGATGC",
    "AGTAAGAAGAAAATACTTCGGAGTATGAAACGTAACAAAGAGTCAGATGAAGTAAAAGAACAAGTGGAATTAAATTCCAAAGTAAACTTAGATACATCCGTATTAGATAAAACTCCAAAGCAATCAGCAAAGAAAAACATGA",
    sep = ""),
  "Pogonomyrmex californicus" = paste(
    "ATGACCGAGAACGAACCGGGTTGCAGTAAAAAAAAGATCCTTCGATCCATGAAGCGGAATAAAGAAAGTGATGAGGTAAAGGAGCAGGTCGAACTGAACAGCAAAGTTAACCTCGATATCTCCGTCCTAGATAAGACGCCC",
    "AAACAATCGGCGAAAAAGAACATGACGACGGAGACAAACAAAGAGTCTGTAGTTTCCGATTGCCATGAGATTGAGATAAAAACGGATGGGAACGCTATACCCGAATTCGATAAAGAGAACCATGACCCCGAGATCAATCAGG",
    sep = "")
))

#####DNMT3

sequences <- DNAStringSet(c(
  "Apis mellifera" = gsub("\n", "",
                          "ATGTTCGGTCTCGCTCGGTTTCATAGCCTCCCGCCGATCACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGGCAACGGCAGCGACGGTAACGACCACGACGACACCGACGACGACGACGACGGCGACGAGGATGACGAGAATGACATCCCAGGAGGAAATGAAAGCAGCGAGCACGAAGGATACGTTGGACCAGTGGGCGGGAAAAAGTGTTTGGTCGTTGGTGTGCGCGATGGAATTGGGCCTGATATCGGATACGACCGCCTCTCGCGGCAACGTCGAGGTCTCGTCTCGCTCCGAGGTCGCGGCCAAATCGCCGCGAAACGCATCGTCATCGCCCCCCTTACGGATCTCGCCACGGTTGAGAACCGTAACGCCGCGGTTTGGAAAGAAACGTGCCAACGTAGGACGAAGAGGAAGAAGAAGAGGGAAACGAGGGTCGAAATCTGGGAAGAAAGTGGGGAGGAAGGTGTACGCGAACGATTCAACGATCACTGTGGATCAATTTTGGAGATCTCGCAAGAGAAAACTCGGCCGGCCGAAGAAGCAAAACCGTTCCGCCGATCTGCTGTCCTCCTCGTTCAGTCGTCCGACTTTGTCCCACGGCTCGTGTCGAACGACTCGGTCGCGCTCGGTCAGCAATAACAATAACAATAACAACAACAACAACAACGACCTTCTCGCCAACAGCACGGACGCTCCACGCTTCGAGCGACGAGTGTCTCGCGACAATAACGGCGCGTGGAAAAATTCGTCGAAAACTGTGCTCCACGAAAAAACGAACGACGAAAATACGAACGAGAACGATGTTAACGAAGGGAAGGGATTCGAGTCGAAGGATAAAGGGAAACGGGACTCGACGGAATCGACGAAATTCGAGGAGATTCACACGCTTCCTTTGATGAAGAGAATTTTTCGCGAGAGTGGTGTTGTTGTTTACGGTGTTGAGAACGAGAGAAAAGATGGCCGACGAGATACGGACGAGAAGTCGACGACGATTTGCGAGGCTGGGAAGAAAGATTGCTCCGTGAAAAAAGCTTGCTCCTCGGAAGGGAGGAAGGGGGTGGTTCGAGCGAAAGCGTTGAAAGGTAACGAGGAAAATCAGGAGGTGAGAAGAGTGACGCGAGGATCGATGAAAATTGGGAAGGATCTGTCCGTCGGCAAATTGGTGTGGGGCTATTGCGCGGGATGGTGGCCAGCGTTAATTATCGATGCCGACCATGTAGGAATGTTGTCCGAGGAAGGGAAATTATGGGTTTACTGGATCGGAGAAGCTCGTATATCATTATTGAACGAAAAGACGCAGATCGAACCGTTCTCGTGTAACCTCAAGGCTAGACTGACGCAGAATTTAAACGTACCGCGAATTCGTGCCATTGACGCAACTATGCAGATGTTGCGCAAAAAATTGGGCGGTACTTTGACCAAGCCGTATTTCACGTGGATCGAAAGCAATTTCCCCAAAAATATGATCGAAATGTTGGACGAAATCAAGTTTTATCCGTATCCGGTAAAGATGCAACAAAGATTGGACCATCTTAGGGAAAAAAACGCAAAAGTAACGGAAAGATATTTATTGGATCAAAAACGGGAAAATCAAGAGAAAAAATTGGCTGAGAAGTCGAAAGATTCGCCGCAAAAAGTCAACGTAGATTTGACGCTTTTGCCGTTGAAAGAGCAGAAACCGGGAATTATAGCCTGGGCTAAAATTGCCGGGCACAATTGGTGGCCAGCGATGATTATCGATTATCGCGACTGTTGCATGCGAGAACCGACGTTCGGTTGCCAATGGATCATGTGGTACGGCGACTACAAACTGTCGGAGGTGCATCATCAATTGTTCTTGAGGTTCGACA",
                          "AGGGGATGGAGAAAATGCGCGACTACACGAGCAACACGAAGAAGCATATCTACCTCGTAGGAGTTCTCCAAGCCTCCAAGGATTATTGCTCCCGTCTAGGATTCGATACCTCCAACTGGACTTTGGACGACGCTTTCGAATATTTTTCCAAGCCTAATCATTACGATTACGCGTCGTCTGCGAACACTTGGAGGAGAGAAGACTCGGTCAAGATCTACGACAAGTACTCGGCCCGCATCGCGGAAAAATTGAACGAGCTGAAGGACAATCCGAACGTGGACGACCAACGGGCCAACGATATAAACAACAGCGATGACCTGCGATCGGCGATAAAAGGGGAAATCTCGTTCGACTCGTTGTGCCTCAAGTGCCTTCGAGTTTCCAACGACGAAATGGACATTCATCCGTTCTTCGAGGGATCTTTGTGCAAAGATTGTTCGGTAAAGAACAACTTTTTCTCTTTTCTCGAATCCTTGCTCCTTCGAGTATTCTCTAACCAAAAAAATCTTTGTTCCAGGAACGATATAAACCCTTTTTATTGCACCGTGTGCGCCGCCTCTGGGATGGTGATCATCTGCGACAAAGAAGATTGTCCAAGGGTTTACTGCACCGCCTGCATGAAACACCTCCTGTGCCCGACGACTTACGAGCAAGTGCTCCAGGAGGATCCTTGGGAATGCTTCCTCTGCAAGTCAAGATCTTTCACCACAGACACCATCGTCAGGCCTCGAGCTAATTGGAAGGACAAGATAATCAATATGTTCCGTACCAGTTGCGACTCGAATGTGGAACACCTGGTTGCCAAGCATAATTCAGAGAAGAGAAAAATACGCGTTCTCTCCCTGTTCGACGGCCTTGGCACAGGTTTGCTCGTACTTCTGAAACTGGGCTTCATCGTGGACGCGTATTACGCGAGCGAGATAGACCAAGACGCGTTGATGGTTACCGCCTCGCACTTCGGCGACCGTATCCTCCAATTGGGCAACGTGAAGGACATCACGTGTAACACGATCAAAGAGATAGCTCCCATCGATCTTCTGATCGGTGGATCGCCGTGCAACGATCTCAGCTTGGCCAATCCTGCCCGCCTCGGCCTCCACGATCCCAGAGGAACGGGTGTACTGTTCTTCGAGTACCGTCGAATTCTGAAACTGGTGAGGAAGCTTAATAACGAGCGGCATTTGTTCTGGTTGTACGAGAACGTCGCCTCGATGCCCAGCGAGTACAGATTGGAGATCAACAAACATCTAGGACAGGAACCGGACGTGATAGACTCGGCAGACTTCTCGCCCCAGCACAGATTGAGACTGTACTGGCACAACTTTCCCATAGAGCCACGTCTGTTATCGTCCCAAAGGGAGCAGGACGTGCAGGATATACTTACGCCGCATTGCCAACGTTACTCGTTGGTCAAAAAGATACGCACAGTCACTACGAAAGTGAACTCGTTGAAGCAAGGCAAACTTGCGTTGAAACCAATTTTGATGAAGGACGAAAGCGACTCGTTGTGGATAACAGAGCTGGAGGAAATATTCGGGTTCCCGCGTCATTACACGGACGTAAAGAATTTATCTGCGACGAAAAGGCAGAGGTTAATAGGTAAATCTTGGAGCGTGCAAACGTTGACAGCCATTTTCGAATCTCTTTGTCCCTTTTTCGAGCGCGATATCGTGGAAATCGAGGGATGA"
  ),
  "Bombus impatiens" = gsub("\n", "",
                            "ATGTTCGCGCGCCATCCAAGCCACGATCGGAGAACCGTAGTCGAAGACTCGGGGCTCAACTTCGAGTCGGAGTCGCCGAATAGACAGGCAGAGGAGGAAAGGTGGGAGGTCTTTGAGCGAACGAATCTACCCAGCAAGGATGCGTTGCAACGATGGGCGGGGAACAGCGTTTGGTCCATGGTATGCGCCATGGAACTGGGCTTGATATCAGACGTTTCCCGAAACAACGCGCCGTTCCACACGGAACGACCGTCGGATGAACGACCGTCCGTCCAGACTAACGCCACCGTTTCCAAGCCCTCGTTGCCCTTGGAAAATGACTGCATCGGTGCCGCGAGTGAGACGACGAAGAGATCGAGATCGTCCACCAGGAGGAAACGAAAGTGCAAATTTCCATCTAGATCCAAAGCCAAATACTCGTCCATCGGGAAAGCGGTAAAGAAGAGCTCGCGCAAAGATAGAGGAGATGATTCGACCACCGCGAACAACGACTCCAACACCGTAGATAAGTTTTGGAGGTCCCGAAAGGCGGGTCTCGACAGATCGGAGAAACGAACCAACCAAGAGTCCGAAAGTTCGCAACATTTTCTTCCTGACGAAGCTGGAAACCAGCACCACCATCGAGACTCCATCGTAGCCGAACGCGACGACGATGATAATTCGGCGAATCGATCGAATCACGCGGCAAACGGCCTCGCCATACAAATTGACGTTCCGCAGACTCGTGCGTCCGATCAGCAGGTGTCTTACGAAAACGACGTGGCTAACTTCGTTTCACCGGACTATAATCGCACGGAGGAAAGCGATCTCTGTGTGATCACGGCGCAGATCTTGGAAACCAACGAACTCGACGCTTTTCAGACCTCGACGCTTGTCCTCGACGAGCAGAAGGAGCAAATGTTGGTCGACAACGACGAGGAACGTCGTTCCTCGAAAGAAGAATGTAACAACGAGGGTACCTCGACGACGAACAACCACGTCGAAGATAACGCGTCGATGGCCGATCGAATGGAAACTGAGGATTCGTCGGAGCGACCACCGGAGAACGACATTGCCCTCGAACGGCTCGACCAAGTCGAATCGACGAGTTTCGACGATACGGATACGCTTCCTCTGGCGAAAAGAATCTCGTCCGAAAGTCGAACGTCCAGCGAAGGCGAACGAAACAACGAAACGACCAAGCCGGTGAAAAGAGCGGCGAGAAACAAAGGATCCAGGTGGAGCGAAGAATATTACGAGGTCAGAAGGGTTACGCGAGGATCGATGAAAATCAGCAAGGATCTTTTCGTCGGCAAGCTGGTGTGGGGCTGTTGTTCGGGATGGTGGCCAGCGTTAATTATCGATGCCGACCATGTAGGAATGTTGTCCGAAGCAGGCAAATCATGGGTTTACTGGATTGGAGAAGCTCGTATATCATTATTGAATGAAAAGACACAGATCGAACCGTTCTCGTGCAACCTCAAGAGCCGACTAGCGCATAATTCAAACGAAGCGCGAATGCGTGTCATCGACGCAACTATGCAGATGCTGCGGAAATTATTGGGCGGCACTTTAACCAAGCCTTACTTCACCTGGATCGAGAACAATTTACCAAACACGATCGAAACGTTAGACGAAATCAGGTTCTACCCATATCCAGACAAGATACAACAGAGATTAAACAGGCTTAGGGAGAAGAACGCAAAAGTTACGGAGAAATTCTTGCTAGACCAAAAACGTATGTTGTTACTTTCTTCGTTCTTTTTTTTTTCCTCTTCTGTTTCTCCTCCGATGCGAATATCAATTGCTGGATTTCTTTCGATAGGCGAAAGTCAAGGGAAAAGGTTGGCTGAAAAAGCGAAAGATTCGCTGCAACGAGGCAACGTAGATT",
                            "TAACGCATTTGCCGTTGAAAGAACAGAATCCGGGGATTATAACGTGGGCTAAAATCGCCGGACACAATTGGTGGCCAGCGATGATTATCGATTATCGTGACTGTTGCATGCGTGAACCAAGCTTTGGTTGCCAGTGGATTATGTGGTACGGTGATTACCAGCTCTCAGAGGTGCACCATCAGTTGTTCTTGAGATTCGACAAGGGCATGGAAAAAATGCGCAACTACATAAATAATACGAAAAAGCACATATACCTCGTTGGAGTTCTCCAAGCCTCCAAGGATTATTGCTCTCGTCTGGGATTCGAAACCAAGGACTGGACTGTAAACGACGCGTTGAGATATTTCGCAAGGAACAAAGATTCCAAAGAATCGTGCGCTCAAAGGAAAGAAGATTCCGTCAAAATATACGATAAATATTCAGCTTGCATAGCCAAGAAATTGAACGAGTTGAAGAACAATACAAATGTAGACGATGACCGGACGAATGACATAAAGAACAGTGATGACTTACGTTCTGCGATGAACGGAAATATCGCGTTCGAATCGTTGTGCCTAAAGTGTCTACGAGTTGCTGTGGGTAAGACGGAAATTCATCCATTCTTCGAGGGATCTTTGTGCAAAGATTGCTCGGATCGCTATAAGCCCTGCATGTTCCTTTTCGGCAATGACGACAAGTGCTTTTACTGCACGGTGTGCGCAGCCTCTGGCATGATGATCATATGTGACAAAGAGGACTGCCCAAGAGTGTATTGTACCGCCTGCATGAAACACCTGTTGTGCCCAACGACCTATGACCAAGTTCTTCAGGAAGATCCTTGGGAGTGTTTCCTCTGCGCAACTAACAATAAGCCTCGACCGTCCTTCGATTCGATCATCAAGCCACGAGCAAATTGGAAAGACAAGATGATCAATATGTTCCGTACGAAATCTAACCCGCAACAGTTAGTGAAACGTAATCGGGAGAAGAAAAAGATACGCGTGCTGTCCTTATTCGATGGCCTCGGTACAGGTTTACTGGTACTTCTCAAGTTAGGTTTAGCCGTGGACACGTATTACGCGAGCGAGATCGACCCGGATGCGTTAATGGTGACCACTGCGCATTTTGGTGATCGAATCGTTCATTTGGGTAACGTGCAAGATATCACGAAGGAAAAGATTCAAGAAATAGCGCCTATCGATTTATTAATCGGTGGATCGCCGTGCAATGATCTCAGCTTGGCTAATCCAGCGCGACTCGGCCTCTATGATCCAAAAGGGACTGGTATCCTATTCTTCGAGTACTGTCGTATCAAGGATCTGCTAAAGGAAGTCAACAATGAATGCCATTTGTTTTGGTTATATGAAAATGTCGCTTCAATGCCGACCGAATATAGATTGGAAATTAACAGACATTTAGGACAAGAACCTGACACGATAGACTCGGCAGATTTTTCTCCTCAGCACAGATTGAGGTTATACTGGCATAATTTACCGATAGAAATACACTCATTGTCTACACAACGAGAGCAGGACGTTCAGGATATACTTACGCCGCATTGTCAACGTTATTCGCTTGTTAAAAAGATTCGGACTGTTACTACGAGAGTAAACTCGCTAAAGCAAGGCAAATTGGCTTTGAAACCAATTCTTATGAAGGACGAAAGTGATTCCTTGTGGATAACCGAACTCGAAGAGATATTCGGATTCCCGCGTCACTACACGGACGTGAAGAATTTATCGGCGACAAAAAGACAAGGGTTGATAGGAAAGTCTTGGAGCGTTCAAACACTAACTGCCATTTTTAGATCTCTTTGTCCTTTTTTCGAGTGTAATACGGTTGAAACGAACTAA"
  ),
  "Camponotus floridianus" = gsub("\n", "",
                                  "ATGTTACGGTTTAACGACGAAATGGAGGAAGCAGTGATCGACAATCCACCTGGGCGCGATTATTCGGACGCGATGCTCCGGAGGCCGACGTCGGCGGAGGGAGACCGCCCGCTGCTCCAAAACAACGACGATCGTCGAGACGAAGCGGATGTTGCGGGACTTGAGCGTGCGGTTGTGCAGCAACAGAAATACCAAGAGTCGGATCGCGTTAAGAGGAAGTCATCGAGATCGGACAAGGTAGTGAGCGCTTGGTCCATAGTCTGCGCCATAGATCTGGGTCTGGTCCAGAAATCATCCCTGAAGGATGTCACTATCGCGAGTTACTTGGACTTGGATTTTCTGCGTACCACGAAGAAAGATGCGCAGAGATTCTCCGATTTCGAGGCGAATAATCGCGAGGACGTCGTATCCGATCGCACGTCCGCGAATATCTACGATGGCTCCTGGCCGTTCATTGACCCGGCGCGAATTACGCGAATCACTTCGAGAAAGTCCGTCCGAAAACGCGGCAGGAGGAGGAGATGGTCGACGAAACGGCAAAGGAAACTAAACGTAGATAAGACGCGTTGGGATATCAACGACGCGACTACGAATCTTTCGATCAAGCGAGCCGCATCGTCGTCGCCCGAAGACATCCGACGTGAAAATCTCGGGGTATCGAATGATGGCCCGACGGGATTATCGGAGAACAAGTTTTCGTTGACGCAACTTCACGCTCGGTGGTCCGATCGCGATGAAGACGAGAATCGCGAGCGTCTATCCGGAAGTTACGGATCGACCGACGAAGACGAAGGAGTTTGTCGAGCAACGTCGGAGACGAGCTCGAGCGAGATCGTTCAGATTGATTATTCCGGTGCTGGATATTCCAAATTCTACACGTCGAATTGCCAAAAGAAAGAATCGCAAAATCGAACGGAATTGTCGGAATATCCCGATTACTTTTCAAACGTCTTGCGTATAGACGACAACAATACCGAGGACGATGATTTTAATGCAACGCCAGATTGGGATAATTACAACACGCGAAGCAAATTTACAGAAATTAATACGCGTACATCCGTTAAGAACGGAAACACTTTAGATAACATTAGAAGGAGGTCGTCGAGAAAAAAGAGTCATCGTTTAAGAAGTAGCGAGAACCATTTGAAGAAACCTCTAAAGGAAAGAAGCAACAAATTGATGCGACGGAAAAAAGAATCTGGAAATTCTGTAACAAATAATGACAATTACAATTTGAGAATTTCCCGCGATTTGGAGGAAAACTCTCCGGAATGCAAATCAGACGCCAGAGACAGCTTCGCTGCTGACACGCGTGATACCGTAAAAATTGTGCGGGGACGTTGCGATTGTAATTTTACAAATGCCGCCTGTACATGTAAGAGGAATTTAAAGAGCACAGCTACCGCGATTATGCGAGATATTCAGGCGCATATTCGTAGCGCGGAAATTGATGTCGAAAATATTCTGGATTTTACGATTTCCGCTGATAGTCTGGGCGAAAAGCCAGGTGATAAGGATCTTAATCGGGATATTTCCATGGAATACGATTGGATCGGATTGGAGGACTACGAGGACCTTAACAATGATGAAATATCGCTCCACAAGTCGAAACAATCATTCGATCTGGTGAAGGCCGATCCAACGGAATCATTTATTGACATGAGTATGAAACGCAGGCGTCGTAATAGACGATTAAACAGATTTACGATCGACTGTAAATTTGTCGACAGAACGGAAACGGCTTATCGAAGTGCCTTCCATGAATGCGTCGATTTTGACTCTCCGGTAACGACATCCTCTCTGAAAATAAACATGGATGAGAGCTGGAGATCGGAAGACGACGACAGCGGCGTAGACATGATCCATTCATCTTGGCGCCGATCGGTATCCTTCATTAAAGATTGCGATAACGAATTTGATAGCGCGCTGTTACTGGCACGAAAAATTATTAATGACGATTCGGCGAAGGATGAAAAATATTACGATGACGAGCACGAGACATTGGAAAACACGTGTTACAAATATCACGAAACAGCGGATGTTTATGCGATTAGTAGAGATACCTCTCTTCATGATTCATGCTCGGATGAGGAGGACGCAAACTTTGTAAACCAACACGATAAAACAAAGAATATAATCACATATAAAAAAAAGGGTAGGGATATTTGGAGATTGTCTTCAATGGATCGCAATTTCGTCAAAAATCGAAGCGATAAGGTTCAAGTTACGAGAAACAAGAAAGTCGATAGCAAGGTACGTTTAATGCGAAGAGATCGAGATGAAAAGGGTACCGATTTGGAATATGATAGTGATAAAAAGGAATACGTTCCGCGAAGAATAACGAGAGCTTTGGTAAAGACTACTAAGGAGAATCATGCGGGCAAGTTAGTATGGGGATACTGCTCAGGATGGTGGCCAGCTCTGATTGTCGATGCCGAACATGCAGGAATGGCTCCTATTTCCGGCAAATTGTGGGTCTACTGGATCGGAGAATCGCAAATATCATTGCTCAGGGAAAGCACGCAGATACAACCATTCTCGAGGAACCTG",
                                  "GAGCACCGATTGACGCAACCATCTCGTAAAAGTATTCGTTCGCGCGCCATTGACGCGACTATACAGATGCTTCGCCAACGTTTCGGCTGTACTCTGACAAAGCCTTATTATTATTGGATACAACGGAATATAAATGGATTCGAGACGTTGGATGATCTAATATTTTATCCATACCCGAAAAATATACAAGAAAGACTGGATACCTTGAGAGAAAAAAATATCAAAGCTACTGAGAGATTCATATCAAACCAAAGAAGTTCACCAGAGACATCACCGACAAAGAAACAGGCTGTTGAAAAATGCAAAGATACAAATACAAGTTGGAAGCAAATCGAAGATGAGCGTTTACCTTTACAAAACCAAAATCCTGGTGTAATAGCTTGGGCGAAAATCGCTGGACATTGTTGGTGGCCCGCAATGATTATCGATTATCGCGACTGCTGCTTGAAGGAACCGAGTTTTGGTTGTCAATGGATTATGTGGTACGGCGATTACAAAGTTTCAGAGGTACGTCATTTGGAATTTTTGAAATTTTACAAGGGACTGGAGAAGATGCGCGACTACATTCAAAATACAGTTAAGCAGTGTTATCTCGACGGCGTGCTTCAAGCTTCGAAGGATTATTGTTCGCGTTTAGGATGCTCCACGGATAACTGGACTTTAGATAATGTGTTTGAATATTTTTCGAATATGAATAATATTCACGTACCATATAACCAATTGCAAGTTTCAGACTCTAACAAGATATATGATAAATATTCCGACGAGATAGTAAAAAAGATCAACGAATTTAAATCCAAGCCGAATGTAGATGCTGAGCGAAAGAATGACATAAAAACGAGTGATGCTCTGCATCGCGTAATATCAGGAGAATGCACAGTGGAAAAATTATGTTTGAAATGTTTAAGGTTTTCTAAAGATAAAATGGAGGAACATCCGTTCTTTATAGGATCTTTATGTAAAGAATGTTCGTATGCATTCAAACCGTGTATGTTTGTGCATGGAAACGATGGAAAATGCTTTTATTGCACAGTTTGCGCCGCCACCGGAACTGTTCTTATTTGTGACACAGACGATTGTCCACGTGTCTATTGTACTGCCTGCTTGAAATTTCTGATATGTCCAAAAGCGTACGACGATATGCTATCGGAAGATCCGTGGAAATGTTTTCTCTGCAGAGATGAATCTAAACAACCTGCGAATATGTTATTGCATCCGCGATTGGATTGGAAAGAAAAAGTTTCTACAATGTTCCGCACAGCCTCTAATCCTGCCTCCAAAGACATAAATTTTGAAAATTATAAGAATCAAAAGAAACCCATACGTGTATTGTCACTTTTCGATGGTTTGAGTACTGGATTTTTGGTGCTTTTAAATTTGGGAATTGTCGTGGATGTTTATTACGCGAGTGAAATTGATAAGAATGCCTTAACGATTAGCTCTGCTCATTTTGGCGATCGAATCACTTACTTGGGAGATGTGAGAGATATAACGAAGGAAAAGATTCAAGAAATTGCACCAATCGATTTGTTAATCGGCGGATCACCGTGCAATGATTTGAGTCTGGTCAATCCGGCGCGATTAGGTCTCCACGATCCCAAAGGAACGGGAATTCTCTTTTTCGAGTACTGCCGAATTAAAAAACTAATGAAGAAAGCCAATAAAGATCGCCACTTGTTCTGGTTGTACGAAAATGTTGCTTCAATGCCGAGCGAGTACAGATTAGAAATTAACAAGCATTTGGGACGGGAACCCGACGTCATAGACGCAGCGGATTTCTCGCCTCAGCATAGATTGCGACTTTACTGGCACAATTTCCCCTTTAATCCGTATATGCCGTTGTTTCAAAATCAACAGGACGTTCAAGATAAGCTTACACCAAATTTAAATCGTAAAGCTCTTTGTAAGAAACTGCGCACGGTTACATGCAGGACAGGTTCTCTGTTACAAGGAAAAGCGGAAGTAAAACCAATCATGATGAAGGGTGAAAGTGATAGATTGTGGATCACCGAACTCGAAGAGATCTTTGGTTTTCCGCGACATTTTACGGACGTGAAGAATTTGTCGGCCACAAATCGGCAGAAATTGCTCGGCAAGTCTTGGAGCGTACAAACACTTACTGCCATATTGCGACCTTTATGTTTCTACTTTAAATGCAAAGAAGACGAGACGAGCAACAATATTTCTGCTTTACATAAAGACATCTCACTTCTGTATCGCGGTAAACAATTTTAG"
  ),
  "Oocera biroi" = gsub("\n", "",
                        "ATGCCTTTAGGAACGAGTGCGGTCGAATCAAGTTGGCAGTACTGGCGCTTCCCGCGTCGAGTGTACGCATGGCCGCGATGCGTGAGTGATCGGCCGTTAAAAAAATCGTTTCGGAGTAGTGTCGCGATGGAAAGAATGGCGATTAACAAATCCTACGAAGACGATTACGGGAATCCGATGGCCTGGACGTCGGCACCTACCGAGATGAGCGACGCGTGGCTCCAAAACAGCGATCCTCCAAACGAAGCGAACTTTCCGGAACAGCGGCCGGCTACAGCGCAGCAACAGAAACATCAAAATTTGGACCGACACAAAAGCGGGCCCAGGAAAAAATCCACAAAAGATCCAAAAAAGATAGCGAGCGTGTGGTCCGTAGAATGTGCCATGCAATTGGGTCTGCTCCCAGAATCATCCCTGAAAAACATCACCATCGTGAATACCTTGGACTTCGATTCTATGCGTAACGACAGGAGCGCAGAGAGCTTTGTCGATCTCGAGACGAGCGATCACGAGACCGTGCTCGAGACGAGACTGTCCAAGTTCACGAACAGGTTGAATGCCACGTGGCCCTTGATCATCCCGGTGTCCGATCCCTCGTCGAGGCCGATTCGCAAGCGAAGCAGGCGGAAAAAGAGATGGACGATAACACGACGTAGGAAACCGAGGAACAGCGTAGATGAGGAACTCAACGACGCGACGGCGGACTCGACGATCAAGAGATCCGCGTCGTCGCCGCTTCGACGTAAGCCGAAACGGGGGAATCTCGGCGTGTCGACTGATGCTTGGGAGAGATTATCAGACGACAACTTTTCGTCGACGCAACTTCACGCTCGGTCACGCGATCATGATGAGGATGACGACGAGGATCGCGAGCTGTCCAAGGATTACGGATTGGCCGGCAGAGAGAAAGTTCGTCGCGCGACGTCGACGAGAAACGCGACGAGAGAGATCACGCAGCCCGATTACTCCGAGTCCGAATATTCCAGATCCTACACGTCGAGTCGGCAGAAAGAGGAGTCGAGGGAGTTGTCGGAATACCCCGAATACTTTATGAACGTCCTGCGTATAGATAACAATACCGGGGAGGATGATTCCGACACTTCGAAGGACACGAAGCCAGATTGTGATACTTCCGTGGACGATGCTCGTGGCAGAAGGTCCGCGAAACTCTCGGAACGCGCTCCGCGCGCGATTAATGTGCGTAGGAACGAAAGAACTTTGGCCGTGGACAATAGGGTGCTATCGAAGAGGAGGAATCGTCATTCTTCAATGGACAAGAAACCTGCCAGGAACTTGCAGAATAGGAGCATTAGAAGAAGAAGAAGAAGAAGAAGACGGGCCAACGAATGGAAGAAGGAAGACTCCGCGACGGATAGTAATTACAACTCAAATTCGCAGGATGACTCCTCGGCATGCAGATCCGACGGCAGGAAAGGGGATGCGCGCGGAACCGCGAGAAACGGACGTAAACCTTGCGACTTGTGTAACACGCGACGACACTCCGGCAATCTTACGAATAACGCCTGCACGTGCAAGAGGAATCTGCGGAACACAGCCAGCGCGATCATGCGGGACATTCAGGAGCAGATACGTGACACGGAAGTCACGCGGGATGCCGAAAATGCTGGCAACGCGATTCCTGACATCTCCGTCGATGATCCAAGCCAGGAACATTCAAAAGAGAAGCTCGATCGAGACGATGTTTCCGTGATGCGCGGTGATGAGGAAGACAACGATAATCGACGGGATAACGTGGAGGACCTGTTATTCCAAGAGTCGAGACCATCCTCTAAAAACGTCGAGCCAAAGATTCAATCGCCAAATTTCACAGTGAACAGCGGTAAACGTAGGAAAAAGAGATTAAATAGATGGATGTCATCTATCAATTGTGATTTAGTCGACAGAACGCAGACAACTTACCACAGCACTTTCCACGAGTGCGATCTCGAGTCTTCAACAACGTCGTCGCTGAAGATAAATATGACGGAGAACTGGAAACCGGAAGACGGTGACAACGCCGGAATGATCCATTCGTCATTCCGGTCTTCCTTATCCGAAGGATGCAGCGACGAGCATGATAGCGCGTTGTTATTGGCTCGAAAGATTATTAGCGACGATTTGTCGAAGGACGAAAAGACCGACGATGAAAACAAGACGATGGGAAGCTCATGTTACAAATATGAGGATAAATGTGACGAAACATTGAATAGTCAGCAGACGGAACAGGGACGTAGCAGGGATAATTTTATTTACAATTCGTCTTCGGAAGATGAAGACACGAAGGAGGATACGGAATTTGAAAATCAATATGACAATCATTTAATGAAAGACATAGTTACGTATAGAAATAAGAATAAAGTTACCCGAAGATCATCGGTCGATCAGAATCTCGTCAAGGATGAAAGGGATGATACCAAGGTCCGAACGGAAAAAGGTATGAGAAAAACAGAAATTCGTAGCAGAACACATCTAGTCCACAGAAACGCAAGTGGCAAAGATAACAACTTGGAAGATGACAGTGATAAAAAAGAATATGTTCCACAAAGAATAACAAGGGCTTTGATAAAGCAGAAAAAGGATAAGGGGAATCATGTTGGCAAGTTAGTGTGGGGATACTGCTCTGGATGGTGGCCAGCTCTCATTGTCGATGCAGCGCATGTAGGAATGGTTCCGACTTCCGGCAAGTTATGGGTATACTGGATCGGAGAATCGCAAATATCATTGCTCAACGAAAAGACACAAATACAACCATTCTCGAGGAACCTGGAACATCGGTTGGCCCAACCTCGTAAAAATATTCAATCACGCGCCATTGACGCGACTATACAGATACTTCGCCAATATTTCGACTGTACTCTGACAAAGCCTTATTATTTTTGGATACAACGGA",
                        "ATATTTCTGATATGGAGACACTGGACGACCTGACGTTTTACCCATATCCGAAAAATATACAAGAAAGACTAGATGTCTTAAAGGAAAAGAACGCCAAAGTTACGGAAAAATTTATATCAAGCCAGAGAAGTTCACCGGAAACAACACCGCCAAGGCAAAGCGCAGCGAAGAAGCAGATTGCTGATAAATTAAAAAATATAAATACGAGTTGGAAACAGATCGAGGATGAGCGTTTACCGTTGCAACATCAAAAACCTGGCGTAATTGCGTGGGCGAAAATTGCTGGACATTCTTGGTGGCCTGCAATGATTATTGATTACCGCGATTGCTGCTTGAAGGAACCAAGTTTTGGTTGTCAATGGATCATGTGGTATGGGGATTATAAAGTTTCAGAGGTACGTCACCTAGAATTTTTAAAATTTCACAAAGGAATAGACAAAATGAAGGACTATATACAAAATACATCGAAACAATCGTTCCTTGAGGGTGTGCTCCAAGCATCAAAGGATTATTGTTCGCGTTTAGGATGCAAGACGGATAACTGGACTTTAACTAATGTCATTGAATATTTTTCGAATATGAATAATATTCACGTACCATATAACCAATTGCAAGTTTCAGAGTCTAACAAGATATATGATAAGTATTCCGACGAGATAATAAAAAAGATCAACGAATTTAGATCGAAGCCAACTGTAGATGATGAACGAAAGAGTGACATAAAAGAGAGTGATGCACTGCGGAGTGTAATATCAGGAGAATCCACCATGGAGAAATTATGTTTGAAATGTTTGAAGTTTTCTAAAAGTAAAATGGAGCAACATCCGTTCTTTGTTGGATCTTTGTGTAAAGAATGCTCGAACGAGTTTAAGCCATGCATGTTCGTGTTTGGGAACGATGGCAAGTGCTTTTATTGTACCATTTGTGCTAGCATGGGAACGGTCATTATTTGCGACAGGGAAGATTGCCCTCGTGTTTACTGCACAGCTTGTTTGAAATATTTGATATGTCCCAAAACGTATGACGACATGTTATTGGAAGACCCGTGGGAATGCTTTCTCTGTAGAGATGAGTCTAGGCAACCTGCGAACATGTTGTTGCAGCCACGACCGGATTGGAAGAGCAAAGTCACCATCATGTTTCGCACCTGTAACTCCACATCGGATAATGTGGATCTCAAAAGTTACCAAAAGAAAAAGCGACCCATACGCGTGCTGTCGCTTTTCGATGGTTTGAGTACCGGTTTGCTGGTGCTTCTGAAATTAGGACTCGACGTGGAAGTTTATTACGCGAGCGAGATCGACGATGATGCTCTAATGGTCAGCTCCGCTCATTTTGGCGATCGTATCACGTATCTGGGAGACGTAAGAGGCATAACGAAGGAAAAGATCCAAGAAATTGCACCGATCGACTTGCTAATCGGTGGGTCACCCTGCAATGATCTGAGTTTGGTCAATCCGGCACGAATGGGTCTCCACAATCCTAAAGGAACGGGTATTTTGTTCTTCGATTATTGTCGAATCAAGAAACTGTTGAACAAGGCGAACAAAGGTCGTCATTTGTTCTGGCTGTTCGAGAACGTCGCTTCCATGCCCACCGAGTACCGATTAGAAATAAACAAGCATCTGGGGCAAGAACCTGATGTTATAGACTCGGCAGACTTCTCACCGCAGCACAGACTGAGGCTCTACTGGCACAATCTTCCTCTCAATCCATACATGCCATTATTCCAAAAGCAACAAGACGTTCAGGATATACTGACGCCGCACTGTAATCGTTACGCTTTAGTGAAGAAAATACGCACGGTCACAACCAGAACGAATTCCCTAAGACAAGGAAGCAGAGCAGAACTGAAGCCAATTATGATGAGGGGTGAATGCGACACGTTATGGATCACAGAGCTTGAGGAGGTCTTTGGTTTTCCGCGACATTACACGGATGTGAGAAATTTGTCGGCCACAAACCGACAAAAGTTGATCGGCAGGTCTTGGAGTGTACAGACACTTACTGCCATATTGCAACCTTTATGCTCCTATTTTAAGTGCAATGAGAACGAGACGAGCTGATCCGATTGA"
  ),
  "Pogonomyrmex barbatus" = gsub("\n", "",
                                 "ATGTACAATGCGCACAACCACGAACGCGAGGATGTGTTGCGACGGAGCAAGCGCGAGAAAGAGACAACGGAGAAATCGCGCGAATGCCAAAAGCCGAATCAACTTAAAAGGAGGTCATCAAAATTATCTAAGGTAGCAAACGCTTGGTCCTTATGTTGCGCCTTGGAATTGGGTTTGCTTCGCGAATCATCTCTCAAGGAAGTTGCCATCGCGAGTTACTCGGATCTGGATTCTCTGTCTACCGTGGAAGATGCGACATATCGCGACAATCTCGAAACAAATGATATCGAAGATGATACATCCAAGTCTGCGAGTATCTTCGGTAACTCGTGGCCGTTCATCGATCCGGCGCAAATTACGCGAGTCTCCCGCACGAAAATCACCCGAAGGCGGAGGAGAAGACCGAGATGGTCGATGAGACGGCAAAGAAAACCGAGAATCGACAACACGACCGAGAGAACCGTTTCGCTACAAGGAACTGAGCAAGAAAATCTTGGCGTATCTAACGATCGGTCAGTAAGATTATCATTAAGGAAGAAGTCTTTGCGATCTCCTCCCCGATCGCTCGCTCGCGACAAGCGTTTATCCGAAAGCTACAGATCAACTGACGATGACGACGAAGAAGAGGAAGAAGTAATTAGTCAAACGGCATCGGAAATAAGCTCAGAAGAGATCGCTCCAACTGATTATTCCGGTTCTGGGTACTCCAGATTCTACTCACCGAGTCAAGCAAAGGAAGACCTTCAAAATCGGACGGAATTGTCGGAATATCCCGATTACTTTATGAACGTCTTACGCATTGGTAATAGTAACAATGATCTCGATACCTCGGAGGATGTAACGTCGGTTTCTGATACTTGGAGGGATCGTTATTCCTCGCAAACCGACAAGTCGATGAAACCTGTGAGACGCCATATTTCAAACAAGATCCTCGAGAATGGAAACTCTCCGCGTAGAATGAAAAAATCGTCATCGAGAAAGAGAATTATGCAGGATATCCAAAAGTGCATAAATACGGAAGTCACGAAGATCGGGGAAGCCGATGAAGTCTGTGAGGTCAGTGAGAATACCATTTTCGATGTCGAGAAGAATATGCTCGAAATTCTCGGTTATAAAGATCTCAGTCAGGATATCTCTATGGAATACGATTGGGACAATGACGAGCTTAACGAACGGGATTTATCATCCCAGGACGTGAAATGCTCACTCGTCAAACCGACAAATGTCGAGGAATCTTCGATAAAAGTTAGTACAATTGGCAGACGTCGTAAGAGGAGGAACAGATGGCATACGATCGATTGCAAATTCGACGAGAAAACGGAGACGGCTTACCGAAGCGTTTTCCACGATTGCGTCGATTTCGATACTCCGAGATCGACGCCGTCCCTGAAGATAAACATGGAGGAGAACTGGAAATCAGAGGACAGCAGTGACAACGACGATGGGATCCACTCATCCTGGCGCCGATCGACCTTAGGATCCTTAACTGAAGAATGCGGTAACGAAGTTGATGGCGCGCTGCTATTGGCTCGGAAGATTGTTAGCGACGATTCGACGAAGGACGAAAAATGTGATGACGAGCACGAGACGTTGAGTAACATGTGTTACAAGTATCAGGAAACGTCAGATACTTTCCCAGGAATGACGACAAAAGATTCTTTTTACGATTCGGATGAGGACATGAAAGACATGTTAATAAGTCAACGTGACGAAACGAAGGATATAATTACGTATAAAAATAAAAATAGATCAGTTAGGAGATCAGCATCAATCAACTGCAGCGTAATCAAGGATGAAGATGACGACACTAAGGTTCAAGATAAAATTACGAAGACCACAAAAGCTAGTAGTACTATAAACGTAATGCGACAATATCCGCGTGACAAGAGTAATAATTTGGAGGATGAAAGTGATAAAATAGAAAATGGTCCACGAAGAGTGACCAGAGCCTTAATAAAAACTAAAAATCATAAAGAGTATCATGCCGGCAAACTCGTATGGGGATTCTTCGCAGGATGGTGGCCAGCATTGATCGTCAAAGCTGAACATGTGGGAATGGCTTCTGAACCCGGGAAATTGTGGGTCTACTGGCTCGGAGAATCGCAAATATCGTTGCTCAAAGAAGAAACACAGATACAACCATTCTCGAGGAACTTGGAGGATCGATTGTCGCAACAGTCTCATAAAAGTATTCGTTCGCGCGCCATTAACGCGACTATACAGATGCTTCGCCAACGTTTCAACTGTACTCTGACAAAGCCTTATTATTATTGGATACGACGGAACATAGCTGATGTGGAAAAATTGGATGATTTAACATTTTATCCGTACCCGGAAAATATACAAGAAAGATTGGATGCTTTGAAGGAAAAGAATGCCA",
                                 "AAGCTACAAAGAAATATATATTAAGTCAAAGAAGTTCACCAGAGTCACCACCGAAAAAACAGAATACGGAAAAAAACAAAGATGTAAATACAAGTTGGGAACAGATTGACAATGAGCGTTTACCTTTACAAAATCAACAGCCTGGTGTAATAGCCTGGGCAAAAATCCCAGGACACTGTTGGTGGCCCGCGATGATTATTGATTATCGCGATTGCTGTTTGAAAGAACCGAGTTTTGGCTGTCAATGGATCATGTGGTATGGCGATTATCAAGTTTCGGAGGTGCGTCATCTCGAATTTCTGAAATTTCACAAGGGAATAGAAAAAATGCGCGAATATATTCAAAATTCGAGTAGGCAGACATATCTTGAAGGCGTACTTCAAGCTTCTAAGGATTATTGTTCGCGTTTAGGATGCAACGTGGATAGCTGGACTTTAGATAATGTCTTTGAATATTTTTCGAATGATATTCACATACCATGTAACCAATTGCAAGTTTCAGACTCTAACAAGGTATATGATAAGTATTCCGACGAGATAGTAAAAAAAATCAACGAATTTAAATTCAAATCAGATGTAGATGATGAACGAAAACGTGACATAAAAGCGAGCGATGATTTACGTCGCGTAATGTCAGGAGAATGCGAAGTGGAAAAATTGTGTTTGAAATGTTTGAAGGTTCCTAAAGGTAGAAAAGATGAACATCCATTCTTTATAGGATCTTTATGTAAAGAATGCTCGAACGAGTTTAAACCGTGTATGTTCGTGCATGGAAACGATGGCAAGTGCTTTTATTGCACAATTTGCGCTGCTACTGGAACTGTTCTTATTTGTGATTCAGAAGATTGTCCACGTGTCTACTGTACAGCTTGTTTGAAATTTTTGATTTGTCCAAAATCATACGACGATATGCTGTTGGAAGATCCATGGGAATGTTTCCTCTGTAGAGATGAAGCCAAACAGCCTATGAATATGATTTTGCGACCGCGGTTAAATTGGAAAGAGAAATTTACAAGAATGTTCCACACAGCCTCTAATATTACATCTAATGTAGACATTGTAAATTATAAGAAAGATAATAGAGCTGTGCGTGTATTGTCGCTTTTCGACGGCTTGAGCACCGGTTTTTTGGTTCTTCAAAAATTGGGACTTGTTGTGGAAGTTTATTATGCTAGTGAAATCGATGTAAATGCTTTAACAATTAGCTCCGCCCATTTCGGTGATCGTATTACTTACTTAGGAGATGTAAGAGGCATAACAAAGGAGAAAATCCAAGAAATTGCACCCATTGATTTGTTAATCGGTGGATCACCATGCAATGATTTAAGTTTGGTCAATCCGGCGCGATTAGGTCTCCACGATCCAAAGGGGACAGGTGTTCTCTTTTTCGAATATTGCCGAATTAAAAAACTAGTAGAAAAAGCTAACAAAAGTCGGCATATGTTTTGGTTATTCGAGAACGTAGCTTCTATGCCAAGCGAATATAGATTAGAAATTAACAAACATTTGGGACGAGAACCTGATGTTATAGATTCAGCAGACTTTTCACCGCAACACAGACTGAGACTTTATTGGCATAATCTTCCCTTTAATCCATATATGCCGTTGTTTCAAAGTCGACAGGACGTTCAAGATAAGCTTACACCAAACTTAAACCGTAAGGCTCTCTGTAAGAAGCTTCACACAGTTACCGGCAGAACAGGTTCTTTATTGCAAGGAAAAGCGGAATTGAAACCAATCATAATGAAAGGCAAAACTGATACGATATGGATCACTGAGCTCGAGGAAATCTTCGGTTTTCCGCGACATTATACGGATGTAAAAAATTTATCAGCTACAAATCGGCAAAAATTACTTGGCAAGTCATGGAGCGTACAAACGCTTACTGCCATATTACGACCTTTATGCTTATATTTTAAGTGCAATGAAGATGAAACGAGAAAGTCAACTTCTTGA"
  ),
  "Pogonomyrmex californicus" = gsub("\n", "",
                                     "ATGTACCGAGTCGACCACCTAAAGCGTGATGAGTTTTTGTTTCAGAAACAATCTACTACGGAGCAGAAGCGGCTCGACAAAAGCTACCTTAGCGCTAGACTGTTAATGAGTGCTAGATTTCGGGTGATTAGGA
TAGCGGCCCCCCTCCAGTTAATTGTAGCAAGAGTAAGTAATACGCATGCAGTCCACGAGTCAGAGTGCCAGAAACCTAACCAGCTCAAGCGACGGTCCTCAAAATTATCCAAGGTAGCTAATGCGTGGTCCTT
ATGTTGTGCACTAGAGCTAGGATTGCTCCGCGAATCTAGCCTGAAAGAAGTCGCTATAGCTAGTTATAGTGATTTAGACTCTCTCTCGACGGTGGAAGACGCAACTTATCGGGACAACCTAGAGACAAATGAC
ATCGAAGACGATACAAGTAAGTCTGCATCCATATTCGGCAACTCGTGGCCCTTTATTGATCCGGCCCAAATTACCCGGGTAAGCCGCACCAAAATCACACGGAGGCGGAGACGCAGGCCTCGATGGAGTATGC
GTCGCCAACGTAAGTTGAGAATAGACAACACAACTGAGCGGACAGTTTCGCTACAAGGTACCGAACAAGAGAATCTGGGCGTTTCTAACGATCGATCTGTTCGGCTTTCCCTGCGCAAGAAAAGCAGTCGCAG
CCCCCCACGCTCACTAGCACGTGATAAGCGTTTAAGCGAATCGTACCGGTCCACTGATGACGACGACGAGGAGGAAGAAGAAGTCATATCACAAACAGCGTCCGAGATTTCCTCCGAAGAGATAGCTCCCACC
GATTACTCTGGCTCAGGGTATTCTAGGTTCTACAGCCCCTCGCAGGCAAAGGAAGACTTACAGAACCGGACGGAGCTTAGCGAGTATCCCGACTACTTTATGAATGTGCTTCGCATAGGAAATTCGAATAACG
ATCTAGACACATCCGAAGACGTTACTTCCGTTAGCGATACCTGGAGAGATAGATACTCCTCCCAAACCGACAAAAGCATGAAGCCGGTCCGCCGACACATCTCTAACAAGATACTAGAAAATGGTAACAGCCC
AAGGAGGATGAAGAAAAGTAGCTCCCGTAAACGCAATACAGAGGTGGGTAACATCCTGTCTGGTGAACACAGTATAATGCCGAAAAAGCGACGCAACCGGCTTGCGCGGCAAAAGGAGGAAATCAGAAACATC
CCTGTAAACAATGACTACAATCGCCGACGCATATCTTACGACTCCCGAGAAAATTTTCCCAAGAACAAAAGTGACCTCAAGGATGACGTGGCGGCTGACCTCCGTGACTTGCGTGATGCCGCGGAAGTCCCCG
TTCAGAAATGCTGTAACTGTAAAACGCAACAGTCTTTTCATAATTTTACTAACACCGCATGCCCTTGTAAGAGGAATTCTAAAAACACCGCTGACGCTATTATGCAAGATATTCAGAAATGCATTAATACAGA
AGTCACCAAAGTAGAGGAAGCTGACGAAGTTTGTGAGGTCAGTGAGAACACAATCTTTGATGTCGAAAAAAATATGCTGGAAATACTTGGCTACAAAGATCTGAGTCAGGACATATCCATGGAGTATGATTGG
GATAACGATGAGCTTAACGAGCGAGACCTAAGTTCGCAGGATGTTAAATGTTCATTGGTAAAACCGACGAACGTGGAGGAGAGCTCTATTAAGGTCTCCACGATAGGGCGACGTCGAAAGAGACGAAATCGCT
GGCATACGATTGATTGCAAATTTGACGAGAAGACTGAAACAGCATACCGGTCGGTATTCCATGACTGCGTGGACTTCGACACCCCGAGATCTACGCCCTCTTTGAAGATCAATATGGAGGAGAACTGGAAGAG
TGAGGATTCTTCGGACAACGATGACGGTATTCATAGTTCTTGGCGTCGCTCCACTCTGGGGTCCTTAACCGAGGAATGCGGCAACGAAGTAGACGGGGCTTTGCTACTCGCTCGTAAAATCGTAAGCGACGAC
TCAACCAAGGATGAGAAATGCGACGATGAGCACGAAACTCTCAGCAACATGTGTTACAAGTACCAGGAAACATCCGACACGTTTCCAGGTATGACTACCAAGGACTCGTTTTACGATTCGGATGAAGACATGAA
AGACATGCTTATTTCCCAGCGAGACGAAACGAAGGATATTATTACATACAAGAACAAAAACAGGAGTGTGCGGCGCTCAGCGAGTATCAACTGCAGCGTAATCAAAGATGAAGATGACGACACAAAGGTGCAAG
ATAAGATCACTAAAACGACAAAGGCCTCATCAACTATCAACGTTATGAGACAGTACCCTAGGGACAAGTCTAACAACCTAGAGGACGAATCGGACAAAATTGAGAACGGGCCGCGACGAGTTACAAGGGCTCTT
ATAAAGACAAAGAACCATAAGGAGTATCATGCAGGAAAACTGGTATGGGGTTTCTTCGCAGGATGGTGGCCGGCCCTTATTGTCAAGGCGGAACATGTAGGAATGGCCTCTGAGCCAGGCAAATTATGGGTGTA
CTGGTTGGGCGAATCGCAGATATCGCTGTTGAAAGAGGAAACTCAAATCCAACCCTTCAGCCGGAACTTAGAGGACCGTCTCAGTCAACAATCTCACAAAAGTATAAGGAGTCGCGCGATCAACGCAACTATAC",
                                     "AGATGCTTCGCCAGAGGTTTAATTGTACTTTGACAAAACCGTACTACTATTGGATACGGCGTAACATAGCGGACGTCGAGAAACTTGACGATCTCACTTTTTACCCATATCCTGAAAATATCCAGGAACGGCTG
GATGCTCTCAAGGAAAAAAATGCAAAGGCAACCAAGAAGTACATACTGTCTCAGAGATCGTCGCCCGAATCTCCTCCCAAGAAACAGAACACGGAAAAGAACAAAGACGTAAATACATCGTGGGAACAGATCGA
TAATGAGCGTCTCCCCCTACAAAATCAGCAGCCAGGAGTAATAGCCTGGGCGAAGATACCAGGACATTGTTGGTGGCCCGCCATGATCATTGACTATAGGGATTGTTGTCTAAAGGAGCCTTCGTTCGGATGTC
AGTGGATCATGTGGTACGGGGACTATCAGGTTTCCGAGGTGCGCCATTTAGAGTTCCTCAAGTTTCATAAGGGAATAGAAAAAATGCGGGAATATATCCAGAACTCGTCGCGGCAGACCTATTTGGAAGGAGTC
CTTCAAGCTTCAAAGGATTACTGCTCCAGGTTAGGGTGCAACGTTGACTCGTGGACACTAGATAACGTGTTCGAGTACTTTTCGAATGACATACATATACCATGTAATCAGCTACAGGTATCCGATTCAAACAA
AGTGTACGATAAGTATAGTGATGAAATAGTGAAAAAGATCAATGAATTTAAATTCAAGAGCGACGTCGACGATGAACGGAAGCGTGATATAAAAGCATCGGACGACCTCAGGCGCGTCATGTCTGGGGAGTGTG
AGGTCGAAAAGCTGTGTCTGAAATGCCTGAAAATTCCGAAAGGTCGGAAGGATGAGCATCCATTCTTTATTGGTAGTCTGTGCAAAGAATGTTCCAATGAGTTCAAGCCTTGCATGTTTGTGCACGGTAACGAT
GGAAAATGCTTTTACTGTACAATTTGCGCTGCCACGGGAACAGTTTTAATCTGCGATAGCGAAGATTGTCCGAGAGTTTATTGTACGGCCTGCCTTAAATTTTTAATTTGCCCCAAAAGTTATGATGACATGTT
GTTAGAAGACCCGTGGGAGTGTTTCCTATGTCGCGATGAAGCAAAGCAACCAATGAACATGATCCTACGGCCCCGCTTAAACTGGAAGGAGAAGTTCACACGAATGTTTCATACAGCCTCGAATATAACTTCGA
ATGTAGACATCGTAAACTACAAAAAAGACAATCGGGCCGTAAGAGTTCTGTCGCTCTTTGACGGATTGTCGACGGGATTCCTAGTCTTACAGAAACTCGGATTGGTGGTTGAGGTGTATTATGCAAGCGAAATA
GACGTAAATGCCTTAACCATCAGCTCTGCACACTTCGGGGACCGGATTACGTACCTCGGAGACGTTAGGGGCATTACAAAGGAAAAGATCCAAGAAATCGCTCCGATAGACCTCTTAATTGGCGGTTCCCCCTG
TAACGATCTAAGTCTGGTGAATCCGGCTCGATTGGGCCTCCATGATCCGAAGGGAACAGGTGTGTTGTTTTTCGAGTACTGTCGCATAAAGAAGCTAGTCGAAAAGGCTAACAAGTCGCGTCACATGTTCTGGC
TATTTGAGAATGTAGCGTCGATGCCCTCCGAATACAGACTCGAAATAAACAAACACTTGGGGCGCGAGCCCGACGTCATTGACTCTGCCGACTTTAGTCCTCAACATCGGCTTAGATTATACTGGCATAATCTG
CCCTTCAATCCCTATATGCCACTGTTTCAATCTCGCCAGGATGTACAAGACAAGCTAACGCCCAATTTGAATAGAAAGGCCTTATGCAAAAAGCTGCATACTGTTACTGGCCGGACGGGGAGCCTTCTCCAGGG
AAAGGCGGAACTTAAACCTATAATCATGAAGGGTAAGACGGATACCATATGGATCACTGAACTAGAGGAGATTTTTGGCTTTCCCCGCCACTATACGGAC
GTCAAAAATCTATCTGCAACGAATCGACAAAAACTATTAGGAAAATCATGGTCGGTGCAAACTTTAACCGCCATCCTACGACCTTTATGCCTTTATTTCAAATGTAATGAAGACGAAACGCGCAAGTCGACTTCC")
))



# Perform multiple sequence alignment (ClustalW)
alignment <- msa(sequences, method = "ClustalW")

# Convert alignment to DNA format for phylogenetic analysis
aligned_seqs <- as.DNAbin(alignment)

# Create a distance matrix using Kimura 2-parameter (K80) model
dist_matrix <- dist.dna(aligned_seqs, model = "K80")

# Build the Neighbor-Joining (NJ) Tree
nj_tree <- nj(dist_matrix)

# Plot the NJ tree
plot(nj_tree, main = "Neighbor-Joining Tree for Gene Sequences")

# Convert sequences to phyDat format for Maximum Likelihood analysis
phyDat_seqs <- phyDat(as.character(alignment), type = "DNA")

class(nj_tree)

nj_tree <- as.phylo(nj_tree)

phyDat_seqs <- phyDat(aligned_seqs, type = "DNA")

tree_labels <- nj_tree$tip.label
data_labels <- names(sequences)

setdiff(tree_labels, data_labels)  # Show mismatched names (if any)
setdiff(data_labels, tree_labels)  # Show missing labels in tree

nj_tree$tip.label <- data_labels  # Rename tree labels to match sequences

fit <- pml(nj_tree, data = phyDat_seqs)
fit <- optim.pml(fit, model = "GTR", optGamma = TRUE)
plot(fit$tree, main = "Maximum Likelihood Tree")


