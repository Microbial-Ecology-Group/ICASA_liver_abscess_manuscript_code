
# Start Here, General Notes -----------------------------------------------
#Project Name: ICASA 
#Script Generation by JDY Started on 5Aug24 Finished on 9Dec24 
#Samples were generated from a subset of cattle (bos taurus origin) enrolled on a RCBD Trial
#Sample types include rumen, ileum, spiral colon, and fecal swabs
#Goal: The major goal of this project was to evaluate the differences in microbial community compositions of animals with different liver outcomes at harvest
#Major note: All data imported into this script was sourced directly from a QIIME2 pipeline completed on the TAMU HPRC (Grace portal) on 3Aug24 and were not altered in any way before importing into this R script
#Note: Functions in the homemade functions tab are a collection of code sourced from collaborators in the VERO group, and from JDY
#To generate matching images, figures were built with the legend, then the get_legend function was used to extract the legend into it's own object, then the image was re-ran with the leged omitted and the two objects were latter re-merged into one image

# Library -----------------------------------------------------------------

setwd("/scratch/user/jdyoung/ICASA/ICASA_R")

library(phyloseq)
library(ggplot2)
library(GUniFrac)
library(stringr)
library(dplyr)
library(metagMisc)
library(metagenomeSeq)
library(vegan)
library(ggdendro)
library(pairwiseAdonis)
library(randomcoloR)
library(btools)
library(multcompView)
library(gghalves)
library(ANCOMBC)
library(tidyr)
library(cowplot)



# Homemade functions ------------------------------------------------------
#change Silva names
changeSILVAtaxa <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Kingdom = str_replace(x[,1], "k__",""),
                          Phylum = str_replace(x[,2], "p__",""),
                          Class = str_replace(x[,3], "c__",""),
                          Order = str_replace(x[,4], "o__",""),
                          Family = str_replace(x[,5], "f__",""),
                          Genus = str_replace(x[,6], "g__",""),
                          Species = str_replace(x[,7],"s_",""),
                          stringsAsFactors = FALSE)}

# Example of a function that we can use (and re-use) to remove unwanted taxa
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

### function to calculate generalized UniFrac using GUniFrac package
library(GUniFrac, quietly = T)

gunifrac <- function(x) {
  y <- GUniFrac(as(t(otu_table(x)), "matrix"), phy_tree(x), alpha = c(0.5))$unifracs
  z <- y[, ,"d_0.5"]
  z.d <- as.dist(z)
  
}

merge_low_abundance <- function(data, threshold=1){
  transformed <- transform_sample_counts(data, function(x) {x/sum(x)}*100)
  otu.table <- as.data.frame(otu_table(transformed))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "zzzOther"
      tax_table(merged)[i,1:7] <- "zzzOther "}
  }
  return(merged)
}

###merge anything less than top code
merge_less_than_top <- function(data, top=20){
  transformed <- transform_sample_counts(data, function(x) {x/sum(x)}*100)
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "zzzOther"
      tax_table(merged)[i,1:7] <- "zzzOther"} # change to represent the number of levels
  }
  return(merged)
}


# Function to extract legend
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


calculate_rarefaction_curves <- function(otu_matrix, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(otu_matrix, measures, depth) {
    if(max(sample_sums(otu_matrix)) < depth) return()
    otu_matrix <- prune_samples(sample_sums(otu_matrix) >= depth, otu_matrix)
    
    rarified_otu_matrix <- rarefy_even_depth(otu_matrix, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_otu_matrix, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, otu_matrix = otu_matrix, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}


# Color Palettes ----------------------------------------------------------

trt_palette <- c("yellow1", "red2", "green4", "dodgerblue2")
liver_palette <- c("grey0", "maroon1","goldenrod", "darkviolet", "darkorange1", "#fff3af", "#8FAD48")

combo_palette <- c("grey0", "maroon1","goldenrod", "darkviolet", "darkorange1", "#fff3af", "#8FAD48","green","blue", "red", "yellow")

binom_palette <- c("grey0", "red")

#used for pub figures
sample_site_palette <- c( "red2","gold2",  "dodgerblue2","green3")

simple_palette <- c("darkseagreen","lightskyblue", "mediumpurple3", "darkorange3")
git_simple_palette <- c("darkseagreen", "darkorange3")

git_simple_dendro_palette <- c("darkseagreen", "darkorange3", "darkolivegreen2","deepskyblue3","yellow1", "firebrick3")

fecal_dendro_palette <- c("darkseagreen","lightskyblue", "mediumpurple3", "darkorange3", "darkolivegreen2","deepskyblue3", "yellow1", "firebrick3")

treatment_palette <- c("darkolivegreen2", "deepskyblue3", "yellow1","firebrick3")

full_dendro_palette <- c("darkseagreen","lightskyblue", "mediumpurple3", "darkorange3", "darkolivegreen2","deepskyblue3", "yellow1", "firebrick3","red2","gold2",  "dodgerblue2","green3")

li_ra_palette <- c("#1d6d00",	"#db13a6",	"#efc74f",	"#2381c4",	"#c28ee2",	"#c6ffff",	"#f9687c",	"#e2e06c",	"#18ed09",	"#1cbc16",	"#bee88b",	"#a9fccd",	"#d80816",	"#c916ae",	"#0a7204",	"#fcec99",	"#6c2999",	"#b89ee5",	"#ebed89",	"#ccf9ff",	"#9e270f",	"#18f4f1",	"#76ed86",	"#860fc1",	"#7cbbff",	"#1fd394",	"#fff1b7",	"#f45d94",	"#08e08a",	"#55dd2c",	"#f2b252",	"#13256d",	"#ea8272",	"#df34f9",	"#bceff4",	"#a8f95c",	"#22ad40",	"#ffc6c1",	"#e5127b",	"#bdeafc",	"#b73a45",	"#FFF0F5",	"#a6d60a",	"#8B8B83")

all_samples_ra_palette <- c("#1d6d00",	"#db13a6",	"#efc74f",	"#2381c4",	"#c28ee2",	"#F0E68C",	"#c6ffff",	"#f9687c",	"#e2e06c",	"#18ed09",	"#1cbc16",	"#bee88b",	"#a9fccd",	"#d80816",	"#c916ae",	"#fcec99",	"#6c2999",	"#b89ee5",	"#ebed89",	"#ccf9ff",	"#9e270f",	"#860fc1",	"#7cbbff",	"#1fd394",	"#fff1b7",	"#f45d94",	"#08e08a",	"#55dd2c",	"#f2b252",	"#13256d",	"#ea8272",	"#df34f9",	"#bceff4",	"#a8f95c",	"#22ad40",	"#ffc6c1",	"#bdeafc",	"#b73a45",	"#FFF0F5",	"#CD661D",	"#8B8B83")

ru_ra_palette <- c("#1d6d00",	"#db13a6",	"#2381c4",	"#c28ee2",	"#f9687c",	"#e2e06c",	"#1cbc16",	"#bee88b",	"#c916ae",	"#4B1AEC",	"#0a7204",	"#6c2999",	"#ebed89",	"#ccf9ff",	"#D74F75",	"#9e270f",	"#75855A",	"#4F06A7",	"#860fc1",	"#7cbbff",	"#1fd394",	"#fff1b7",	"#f45d94",	"#FFB640",	"#55dd2c",	"#f2b252",	"#13256d",	"#ea8272",	"#df34f9",	"#bceff4",	"#a8f95c",	"#22ad40",	"#ffc6c1",	"#3DC26B",	"#b73a45",	"#D5F12D",	"#435A39",	"#682169",	"#FFF0F5",	"#CD661D",	"#a6d60a",	"#8B8B83")

si_ra_palette <- c("#1d6d00",	"#db13a6",	"#efc74f",	"#2381c4",	"#c28ee2",	"#F0E68C",	"#c6ffff",	"#e2e06c",	"#1cbc16",	"#bee88b",	"#d80816",	"#c916ae",	"#4B1AEC",	"#fcec99",	"#6c2999",	"#b89ee5",	"#33F959",	"#ebed89",	"#ccf9ff",	"#D74F75",	"#860fc1",	"#7cbbff",	"#1fd394",	"#fff1b7",	"#f45d94",	"#55dd2c",	"#13256d",	"#df34f9",	"#bceff4",	"#a8f95c",	"#ffc6c1",	"#b73a45",	"#9F5FE4",	"#FFF0F5",	"#CD661D",	"#8B8B83")

fec_ra_palette <- c("#1d6d00",	"#db13a6",	"#efc74f",	"#2381c4",	"#c28ee2",	"#F0E68C",	"#c6ffff",	"#f9687c",	"#e2e06c",	"#18ed09",	"#1cbc16",	"#bee88b",	"#a9fccd",	"#d80816",	"#c916ae",	"#fcec99",	"#6c2999",	"#ebed89",	"#ccf9ff",	"#9e270f",	"#860fc1",	"#7cbbff",	"#1fd394",	"#fff1b7",	"#f45d94",	"#08e08a",	"#55dd2c",	"#f2b252",	"#13256d",	"#ea8272",	"#df34f9",	"#bceff4",	"#22ad40",	"#bdeafc",	"#b73a45",	"#FFF0F5",	"#8B8B83")

# Import Data -------------------------------------------------------------

data <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
data #705 samples 68786 taxa

#write sample names
write.csv(sample_names(data), "sample_names.csv")
#read meta data
metadata <- import_qiime_sample_data("sample_names.txt")

#merge metadata
data1 <- merge_phyloseq(data,metadata)
data1 #705 samples 68786 taxa


# Explore Data ------------------------------------------------------------
###Fix names
rank_names(data1)
colnames(tax_table(data1)) <- c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
rank_names(data1)

#change Silva names
tax.data <- data.frame(tax_table(data1))
tax.data.names <- changeSILVAtaxa(tax.data)

###change NA to better name 
#convert to charcter
for(i in 1:7){tax.data.names[,i] <- as.character(tax.data.names[,i])}
#replace with empty string
tax.data.names[is.na(tax.data.names)] <- ""

head(tax.data.names)

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}

#check headers re-insert matrix and check tail
head(tax.data.names)
tax_table(data1) <- as.matrix(tax.data.names)
tail(tax_table(data1))

#check for Eukaryota
data2 <- subset_taxa(data1, Domain!= "Eukaryota")#none present

###Split samples and controls
samples <- subset_samples(data1, Sample_type == "Sample")
samples <- prune_taxa(taxa_sums(samples)>0, samples)
samples #68642 taxa and 686 samples

controls <- subset_samples(data1, Sample_type == "Control")
controls <- prune_taxa(taxa_sums(controls)>0, controls)
controls #3712 taxa and 19 samples
sort(sample_sums(controls)) #the weird UDP0290 samples are contributing heavily these are not mine, also need to look at EB7 and EB5

#investigate EBs with reads
controls_genus <- tax_glom(controls, taxrank = "Genus", NArm = F) %>% 
  psmelt()
length(unique(controls_genus$Genus))#369 genera
ggplot(controls_genus, aes(x=Sample_name, y = Abundance, fill = Genus))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")
##Note : UDP samples are really skewing this need to find out what these are.

#look at number of reads per sample and sample distributions
sample_sum_df <- data.frame(sum= sample_sums(samples))
ggplot(sample_sum_df, aes(x= sum))+
  geom_histogram(color= "black", fill = "indianred", binwidth = 5000)+
  ggtitle("Distribution of sample sequencing depth")+ 
  xlab("Read Counts")+
  theme(axis.title.y = element_blank())

mean(sample_sums(samples)) #613341
median(sample_sums(samples)) #621846
sort(sample_sums(samples)) #if we cut at 300k we would loose 17 samples all SI and LI except one fecal 6611
min(sample_sums(samples))#2
max(sample_sums(samples))# 1083504
IQR(sample_sums(samples))# 140419.2

#look at removed samples for min and max
samples_rm <- prune_samples(sample_sums(samples)<250000, samples)
samples_rm # 14 samples
samples_rm <- prune_taxa(taxa_sums(samples)>0, samples_rm)

min(sample_sums(samples_rm))#min 2
max(sample_sums(samples_rm))#225258


##cut at 250,000 
samples <- prune_samples(sample_sums(samples)>250000, samples)
samples #cuts down to 672 samples
samples <- prune_taxa(taxa_sums(samples)>0, samples)
samples #68284 taxa and 672 samples

min(sample_sums(samples))#273,678
max(sample_sums(samples)) #1,083,504
mean(sample_sums(samples)) #625,022.9
median(sample_sums(samples))#625442.5
IQR(sample_sums(samples))#137055.8

# Checking classification efficiency  -------------------------------------
# First we have to extract the tax_table and convert it to a "data.frame"
taxa.df <- as.data.frame(tax_table(samples))
taxa.df
#KINGDOM
# Now let's search for which taxa start with "unclassified"
unclassified_kingdom.df <- taxa.df %>% filter(grepl('Unassigned', Kingdom))
# Pull out the unique 
unclassified_kingdom <- row.names(unclassified_kingdom.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_kingdom.ps = pop_taxa(samples, unclassified_kingdom)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_kingdom.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_kingdom.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_kingdom.ps)) / sum(sample_sums(samples)) * 100

###Kingdom 99.99423%

#Phylum
# Now let's search for which taxa start with "unclassified"
unclassified_phy.df <- taxa.df %>% filter(grepl('unclassified', Phylum))
# Pull out the unique 
unclassified_phy <- row.names(unclassified_phy.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_phy.ps = pop_taxa(samples, unclassified_phy)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_phy.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_phy.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_phy.ps)) / sum(sample_sums(samples)) * 100

###Phylum 99.94697%

#CLASS
# Now let's search for which taxa start with "unclassified"
unclassified_class.df <- taxa.df %>% filter(grepl('unclassified', Class))
# Pull out the unique 
unclassified_class <- row.names(unclassified_class.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_class.ps = pop_taxa(samples, unclassified_class)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_class.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_class.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_class.ps)) / sum(sample_sums(samples)) * 100

###Class 99.90943%

#ORDER
# Now let's search for which taxa start with "unclassified"
unclassified_order.df <- taxa.df %>% filter(grepl('unclassified', Order))
# Pull out the unique 
unclassified_order <- row.names(unclassified_order.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_order.ps = pop_taxa(samples, unclassified_order)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_order.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_order.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_order.ps)) / sum(sample_sums(samples)) * 100

###Order 99.87953%

#FAMILY
# Now let's search for which taxa start with "unclassified"
unclassified_family.df <- taxa.df %>% filter(grepl('unclassified', Family))
# Pull out the unique 
unclassified_family <- row.names(unclassified_family.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_family.ps = pop_taxa(samples, unclassified_family)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_family.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_family.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_family.ps)) / sum(sample_sums(samples)) * 100

###Family 99.69113%

#GENUS
# Now let's search for which taxa start with "unclassified"
unclassified_genus.df <- taxa.df %>% filter(grepl('unclassified', Genus))
# Pull out the unique 
unclassified_genus <- row.names(unclassified_genus.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_genus.ps = pop_taxa(samples, unclassified_genus)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_genus.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_genus.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_genus.ps)) / sum(sample_sums(samples)) * 100

###Genus 58.81937% ??? Need to ask about this


# Rarefaction curves -------------------------------------------------------

#subset ASV table
otu_matrix <- otu_table(samples)
otu_matrix <- as.matrix(otu_matrix)


#Calculate rarefaction curve
set.seed(42)

#calculate observed and Shannon at different depths (caps at 1 mil)
rarefaction_curve_data <- calculate_rarefaction_curves(otu_matrix, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary1 <- rarefaction_curve_data_summary %>% mutate(Sample = gsub("\\.","-",Sample))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary1, data.frame(sample_data(samples)), by.x = 'Sample', by.y = 'row.names')


View(rarefaction_curve_data_summary_verbose)

#subset only observed dont need shannons right now
#subset observed data
sub_rare_data <- subset(rarefaction_curve_data_summary_verbose, Measure == 'Observed')

#order variable for plot
sub_rare_data$Sample_site <- factor(sub_rare_data$Sample_site, levels = c("RU", "SI","LI","Fec"))
sub_rare_data$Liver_Simple <- factor(sub_rare_data$Liver_Simple, levels = c("Edible","A-","A","A+"))

#only plot observed
rare_plot <- ggplot(data = sub_rare_data, 
                    mapping = aes(
                      x = Depth,
                      y = Alpha_diversity_mean,
                      ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                      ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                      colour = Liver_Simple, 
                      group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = simple_palette, name= "Liver Score" )+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

rare_plot

rare_legend <- get_legend(rare_plot)

# Rarefaction plots split by sample type ----------------------------------
#Rumen
ru_rare <- subset(sub_rare_data, Sample_site=="RU")

#Order liver variable
ru_rare$Liver_Simple <- factor(ru_rare$Liver_Simple, levels = c("Edible","A+"))

ru_rare_plot<- ggplot(data = ru_rare, 
          mapping = aes(
            x = Depth,
            y = Alpha_diversity_mean,
            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
            colour = Liver_Simple, 
            group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = git_simple_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)),
                     labels = c("0","","","", "100,000","","","","200,000","","","","300,000","","","","400,000",
                                "","","","500,000","","","","600,000","","","","700,000","","","","800,000",
                                "","","","900,000","","","","1,000,000"),
                     expand = c(0.01,0.01,0,0))+
  ylim(0,4500)+
  coord_cartesian(xlim = c(0,825000))+
  labs(y = "Species Richness", x = "Sample Depth", title = "A.Rumen")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black", angle = 45,hjust = 1),
        axis.title.x = element_text(size = 30, vjust = 2.),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

ru_rare_plot

#Small intestine
si_rare <- subset(sub_rare_data, Sample_site== "SI")

#order liver for plot
si_rare$Liver_Simple <- factor(si_rare$Liver_Simple, levels = c("Edible","A+"))

si_rare_plot<- ggplot(data = si_rare, 
          mapping = aes(
            x = Depth,
            y = Alpha_diversity_mean,
            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
            colour = Liver_Simple, 
            group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = git_simple_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)),
                     labels = c("0","","","", "100,000","","","","200,000","","","","300,000","","","","400,000",
                                "","","","500,000","","","","600,000","","","","700,000","","","","800,000",
                                "","","","900,000","","","","1,000,000"),
                     expand = c(0.01,0.01,0,0))+
  ylim(0,4500)+
  coord_cartesian(xlim = c(0,825000))+
  labs(y = "Species Richness", x = "Sample Depth", title = "B.Small Intestine")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black", angle = 45,hjust = 1),
        axis.title.x = element_text(size = 30, vjust = 2.),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

si_rare_plot

#Large Intestine
li_rare <- subset(sub_rare_data, Sample_site=="LI")

#Order liver
li_rare$Liver_Simple <- factor(li_rare$Liver_Simple, levels = c("Edible","A+"))

li_rare_plot<- ggplot(data = li_rare, 
          mapping = aes(
            x = Depth,
            y = Alpha_diversity_mean,
            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
            colour = Liver_Simple, 
            group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = git_simple_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)),
                     labels = c("0","","","", "100,000","","","","200,000","","","","300,000","","","","400,000",
                                "","","","500,000","","","","600,000","","","","700,000","","","","800,000",
                                "","","","900,000","","","","1,000,000"),
                     expand = c(0.01,0.01,0,0))+
  ylim(0,4500)+
  coord_cartesian(xlim = c(0,825000))+
  labs(y = "Species Richness", x = "Sample Depth", title = "C.Large Intestine")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black", angle = 45,hjust = 1),
        axis.title.x = element_text(size = 30, vjust = 2.),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
li_rare_plot

#Fecal
fec_rare <- subset(sub_rare_data, Sample_site=="Fec")

#order liver variable
fec_rare$Liver_Simple <- factor(fec_rare$Liver_Simple, levels = c("Edible","A-","A","A+"))

fec_rare_plot<- ggplot(data = fec_rare, 
          mapping = aes(
            x = Depth,
            y = Alpha_diversity_mean,
            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
            colour = Liver_Simple, 
            group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = simple_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)),
                     labels = c("0","","","", "100,000","","","","200,000","","","","300,000","","","","400,000",
                                "","","","500,000","","","","600,000","","","","700,000","","","","800,000",
                                "","","","900,000","","","","1,000,000"),
                     expand = c(0.01,0.01,0,0))+
  ylim(0,4500)+
  coord_cartesian(xlim = c(0,825000))+
  labs(y = "Species Richness", x = "Sample Depth", title = "D.Fecal")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black", angle = 45,hjust = 1),
        axis.title.x = element_text(size = 30, vjust = 2.),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
fec_rare_plot

#join plots

rare_plot_final <- plot_grid(plot_grid(ru_rare_plot, si_rare_plot, li_rare_plot, fec_rare_plot, axis = "btrl", align = "hv"),
                             rare_legend, rel_widths = c(1,.5))

rare_plot_final

# Comparing sampling depth between samples --------------------------------

sample_sum_df_new <- data_frame(ASV_count = sample_sums(samples))
metadata.df <- as(sample_data(samples), "data.frame")
seqDepth_metadata <- cbind(metadata.df, sample_sum_df_new)

#sample type
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Sample_site)# p = 0.0001168
pairwise.wilcox.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Sample_site)

depth_plot <- ggplot(seqDepth_metadata, aes(x = Sample_site, y = ASV_count))+
                     geom_boxplot()
depth_plot

#comparing trt
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Treatment) #p = 0.7532

#comparing block
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Block) #p = 0.08481

#liver score 
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Liver_Score) #p = 0.08579

#liver simple
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Liver_Simple) #p = 0.1274

#liver binom
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Liver_Binom) #p = 0.1441

#rumen score
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Rumen_Score) #p = 0.5369


# alpha diversity -------------------------------------------------------------
alphadiv1 <- estimate_richness(samples, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
alphadiv2 <- estimate_pd(samples)

alpha_div <- cbind(alphadiv1, alphadiv2)
alpha_div

alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div_meta <- cbind(metadata.df,alpha_div)
alpha_div_meta

##Look at sample site first
sample_site_anov_observed <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Sample_site, p.adjust.method = "BH")
sample_site_anov_observed #everything different except LI and fecal

sample_site_anov_shannon <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Sample_site, p.adjust.method = "BH")
sample_site_anov_observed #everything different except LI and fecal

##subset fecal samples
fecal_alpha <- subset(alpha_div_meta, Sample_site== "Fec")

#only fecal by trt
trt_observed_anov <- kruskal.test(fecal_alpha$Observed, fecal_alpha$Treatment)
trt_observed_anov# differences between trt p = 0.0000057
pairwise.wilcox.test(fecal_alpha$Observed, fecal_alpha$Treatment, p.adjust.method = "BH")#diet effect

trt_shan_anov <- kruskal.test(fecal_alpha$Shannon, fecal_alpha$Treatment)
trt_shan_anov #no difference p = 0.2428

#only fecal by liver score
ls_observed_anova <- pairwise.wilcox.test(fecal_alpha$Observed, fecal_alpha$Liver_Score, p.adjust.method = "BH")
ls_observed_anova # NS

ls_shannon_anova <- pairwise.wilcox.test(fecal_alpha$Shannon, fecal_alpha$Liver_Score, p.adjust.method = "BH")
ls_shannon_anova #no difference

#only fecal by liver simple
kruskal.test(fecal_alpha$Observed, fecal_alpha$Liver_Simple)# p = 0.08
LS_observed_anova <- pairwise.wilcox.test(fecal_alpha$Observed, fecal_alpha$Liver_Simple, p.adjust.method = "BH")
LS_observed_anova #no difference

kruskal.test(fecal_alpha$Shannon, fecal_alpha$Liver_Simple) #p=0.23
LS_shannon_anova <- pairwise.wilcox.test(fecal_alpha$Shannon, fecal_alpha$Liver_Simple, p.adjust.method = "BH")
LS_shannon_anova #no difference

#only fecal binomial 
lsb_observed_anova <- pairwise.wilcox.test(fecal_alpha$Observed, fecal_alpha$Liver_Binom, p.adjust.method = "BH")
lsb_observed_anova #no difference

lsb_shan_anova <- pairwise.wilcox.test(fecal_alpha$Shannon, fecal_alpha$Liver_Binom, p.adjust.method = "BH")
lsb_shan_anova #no difference

##Look at rumen
##subset rumen samples
ru_alpha <- subset(alpha_div_meta, Sample_site== "RU")

#liver score
LS_ru_observed_anova <- pairwise.wilcox.test(ru_alpha$Observed, ru_alpha$Liver_Simple, p.adjust.method = "BH")
LS_ru_observed_anova # no difference

LS_ru_shann_anova <- pairwise.wilcox.test(ru_alpha$Shannon, ru_alpha$Liver_Simple, p.adjust.methods = "BH")
LS_ru_shann_anova # no difference

#treatment
kruskal.test(ru_alpha$Observed, ru_alpha$Treatment) # p = 0.0001
pairwise.wilcox.test(ru_alpha$Observed, ru_alpha$Treatment, p.adjust.method = "BH")# complex but mainly diet affect?

kruskal.test(ru_alpha$Shannon, ru_alpha$Treatment) #p = 0.000002
pairwise.wilcox.test(ru_alpha$Shannon, ru_alpha$Treatment, p.adjust.method = "BH")#diet affect

##Look at small intestine
##subset si
si_alpha <- subset(alpha_div_meta, Sample_site== "SI")

#Liver score
LS_si_observed_anova <- pairwise.wilcox.test(si_alpha$Observed, si_alpha$Liver_Simple, p.adjust.method = "BH")
LS_si_observed_anova # was a difference p = 0.016

LS_si_shan_anova <- pairwise.wilcox.test(si_alpha$Shannon, si_alpha$Liver_Simple, p.adjust.method = "BH")
LS_si_shan_anova #p = 0.02

#treatment
kruskal.test(si_alpha$Observed, si_alpha$Treatment) # p = 0.6526

kruskal.test(si_alpha$Shannon, si_alpha$Treatment) #p = 0.0123
pairwise.wilcox.test(si_alpha$Shannon, si_alpha$Treatment, p.adjust.method = "BH")# diet affect


##Look at Large Intestine
#subset
li_alpha <- subset(alpha_div_meta, Sample_site=="LI")

#liver score
LS_li_observed_anova <- pairwise.wilcox.test(li_alpha$Observed, li_alpha$Liver_Simple, p.adjust.method = "BH")
LS_li_observed_anova #no difference

LS_li_shan_anova <- pairwise.wilcox.test(li_alpha$Shannon, li_alpha$Liver_Simple, p.adjust.method = "BH")
LS_li_shan_anova #no difference

#treatment
kruskal.test(li_alpha$Observed, li_alpha$Treatment)# p = 0.3548

kruskal.test(li_alpha$Shannon, li_alpha$Treatment) #p = 0.8856


# Alpha plots  Exploratory------------------------------------------------------------
#sample site observed
sample_site_richness_boxplot <- ggplot(alpha_div_meta, aes(x= Sample_site, y= Observed, fill = Sample_site)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  scale_x_discrete(limits=c("RU", "SI", "LI", "Fec"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values =sample_site_palette)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
sample_site_richness_boxplot

#sample site shannon
sample_site_shan_boxplot <- ggplot(alpha_div_meta, aes(x= Sample_site, y= Shannon, fill = Sample_site)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  scale_x_discrete(limits=c("RU","SI", "LI", "Fec"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_site_palette)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
sample_site_shan_boxplot

#fecal trt observed
fecal_trt_observed_boxplot <- ggplot(fecal_alpha, aes(x= Treatment, y= Observed, fill = Treatment)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = trt_palette) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
fecal_trt_observed_boxplot

#fecal trt shannon
fecal_trt_shannon_boxplot <- ggplot(fecal_alpha, aes(x= Treatment, y= Shannon, fill = Treatment)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = trt_palette) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
fecal_trt_shannon_boxplot

#Liver score observed
liver_score_observed_boxplot <- ggplot(fecal_alpha, aes(x= Liver_Score, y= Observed, fill = Liver_Score)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = randomColor(7)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
liver_score_observed_boxplot

#liver score shannon
liver_score_shan_boxplot <- ggplot(fecal_alpha, aes(x= Liver_Score, y= Shannon, fill = Liver_Score)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = liver_palette) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
liver_score_shan_boxplot

#Liver score binom observed
liver_score_binom_observed_boxplot <- ggplot(fecal_alpha, aes(x= Liver_Binom, y= Observed, fill = Liver_Binom)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = randomColor(7)) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
liver_score_binom_observed_boxplot

#liver score shannon
liver_score_binom_shan_boxplot <- ggplot(fecal_alpha, aes(x= Liver_Binom, y= Shannon, fill = Liver_Binom)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = liver_palette) +
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
liver_score_binom_shan_boxplot




# Alpha plots for pubs Liver score ----------------------------------------------------
#set order of Liver simple to match the plot
fecal_alpha$Liver_Simple <- factor(fecal_alpha$Liver_Simple, levels = c("Edible","A-","A","A+"))

#Liver simple observed
liver_simple_observed_boxplot <- ggplot(fecal_alpha, aes(x= Liver_Simple, y= Observed, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "G.Fecal") +
  scale_x_discrete(limits=c("Edible", "A-", "A", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = simple_palette) +
  scale_color_manual(values = simple_palette)+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 400, label=(expression(paste(italic("P"), " = 0.08"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
    plot.margin = unit(c(1,1,1,1),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))

liver_simple_observed_boxplot


# Extract the legend
#Note: run the plot one time with the legend then remove the legend 
# and run again after you extract the legend to a new object to make final plot
legend_alpha <- get_legend(liver_simple_observed_boxplot)

#liver simple shannon
liver_simple_shan_boxplot <-ggplot(fecal_alpha, aes(x= Liver_Simple, y= Shannon, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "H.Fecal") +
  scale_x_discrete(limits=c("Edible", "A-", "A", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = simple_palette) +
  scale_color_manual(values = simple_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 2, label=(expression(paste(italic("P"), " = 0.23"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_shan_boxplot

#rumen
#set order
ru_alpha$Liver_Simple <- factor(ru_alpha$Liver_Simple, levels = c("Edible", "A+"))
levels(ru_alpha$Liver_Simple)

liver_simple_ru_observed_boxplot <- ggplot(ru_alpha, aes(x= Liver_Simple, y= Observed, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "A.Rumen") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 400, label=(expression(paste(italic("P"), " = 0.61"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 3, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_ru_observed_boxplot


#liver simple shannon
liver_simple_ru_shan_boxplot <-ggplot(ru_alpha, aes(x= Liver_Simple, y= Shannon, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "B.Rumen") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 2, label=(expression(paste(italic("P"), " = 0.66"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_ru_shan_boxplot

##small intestine
#set order
si_alpha$Liver_Simple <- factor(si_alpha$Liver_Simple, levels = c("Edible", "A+"))
levels(si_alpha$Liver_Simple)

liver_simple_si_observed_boxplot <- ggplot(si_alpha, aes(x= Liver_Simple, y= Observed, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "C.Small Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 400, label=(expression(paste(italic("P"), " = 0.02"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 3, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_si_observed_boxplot


#liver simple shannon
liver_simple_si_shan_boxplot <-ggplot(si_alpha, aes(x= Liver_Simple, y= Shannon, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "D.Small Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 2, label=(expression(paste(italic("P"), " = 0.02"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_si_shan_boxplot


##large intestine
#set order
li_alpha$Liver_Simple <- factor(li_alpha$Liver_Simple, levels = c("Edible", "A+"))
levels(li_alpha$Liver_Simple)

liver_simple_li_observed_boxplot <- ggplot(li_alpha, aes(x= Liver_Simple, y= Observed, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "E.Large Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 400, label=(expression(paste(italic("P"), " = 0.25"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_li_observed_boxplot


#liver simple shannon
liver_simple_li_shan_boxplot <-ggplot(li_alpha, aes(x= Liver_Simple, y= Shannon, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "F.Large Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,2.2))+
  annotate("text", x = 2.35, y = 2, label=(expression(paste(italic("P"), " = 0.10"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
liver_simple_li_shan_boxplot

#combine plots
alpha_plot_combo <- plot_grid(plot_grid(liver_simple_ru_observed_boxplot, liver_simple_ru_shan_boxplot, 
                              liver_simple_si_observed_boxplot, liver_simple_si_shan_boxplot, 
                              liver_simple_li_observed_boxplot, liver_simple_li_shan_boxplot,
                              liver_simple_observed_boxplot, liver_simple_shan_boxplot,
                              ncol = 2,align = "hv", axis = "btrl"), 
                              legend_alpha, align = "hv", axis = "btrl", rel_widths = c(1,0.5) )
alpha_plot_combo

alpha_plot_wide <- plot_grid(plot_grid(liver_simple_ru_observed_boxplot, liver_simple_si_observed_boxplot, liver_simple_li_observed_boxplot,liver_simple_observed_boxplot,
                                       liver_simple_ru_shan_boxplot,liver_simple_si_shan_boxplot, liver_simple_li_shan_boxplot,
                                         liver_simple_shan_boxplot,
                                        nrow =  2,align = "hv", axis = "btrl"), 
                              legend_alpha, align = "hv", axis = "btrl", rel_widths = c(1,0.5) )

alpha_plot_wide
# Alpha plots for pub Treatment -------------------------------------------
###Fecal
#set order of Treatment to match the plot
fecal_alpha$Treatment <- factor(fecal_alpha$Treatment, levels = c("68","68T","67","67T"))

#Treatment observed
trt_fec_observed_boxplot <- ggplot(fecal_alpha, aes(x= Treatment, y= Observed, fill = Treatment)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "G.Fecal") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR")) +
  scale_color_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR"))+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,4.2))+
  annotate("text", x = 4.4, y = 400, label=(expression(paste(italic("P"), " < 0.001"))), size = 5, color = "black")+
  annotate("rect", xmin = 4.0, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  annotate("text", x = 1, y = 4500, label = "a", size = 5, color = "black")+
  annotate("text", x = 2, y = 4500, label = "a", size = 5, color = "black")+
  annotate("text", x = 3, y = 4500, label = "b", size = 5, color = "black")+
  annotate("text", x = 4, y = 4500, label = "b", size = 5, color = "black")+
    theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_fec_observed_boxplot

alpha_trt_legend <- get_legend(trt_fec_observed_boxplot)

trt_fec_shann_boxplot<- ggplot(fecal_alpha, aes(x= Treatment, y= Shannon, fill = Treatment)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "H.Fecal") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 2, label=(expression(paste(italic("P"), " = 0.24"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_fec_shann_boxplot

###Rumen
#set order
ru_alpha$Treatment <- factor(ru_alpha$Treatment, levels = c("68","68T","67","67T"))

trt_ru_observed_boxplot<- ggplot(ru_alpha, aes(x= Treatment, y= Observed, fill = Treatment)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "A.Rumen") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR")) +
  scale_color_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR"))+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,4.2))+
  annotate("text", x = 4.4, y = 400, label=(expression(paste(italic("P"), " < 0.001"))), size = 5, color = "black")+
  annotate("rect", xmin = 4.0, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  annotate("text", x = 1, y = 3500, label = "a", size = 5, color = "black")+
  annotate("text", x = 2, y = 3500, label = "ab", size = 5, color = "black")+
  annotate("text", x = 3, y = 3500, label = "bc", size = 5, color = "black")+
  annotate("text", x = 4, y = 3500, label = "c", size = 5, color = "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_ru_observed_boxplot

trt_ru_shann_boxplot<- ggplot(ru_alpha, aes(x= Treatment, y= Shannon, fill = Treatment)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "B.Rumen") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 2, label=(expression(paste(italic("P"), " < 0.001"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  annotate("text", x = 1, y = 6.5, label = "a", size = 5, color = "black")+
  annotate("text", x = 2, y = 6.5, label = "a", size = 5, color = "black")+
  annotate("text", x = 3, y = 6.5, label = "b", size = 5, color = "black")+
  annotate("text", x = 4, y = 6.5, label = "b", size = 5, color = "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_ru_shann_boxplot

###Small intestine
#set order
si_alpha$Treatment <- factor(si_alpha$Treatment, levels = c("68","68T","67","67T"))


trt_si_observed_boxplot<- ggplot(si_alpha, aes(x= Treatment, y= Observed, fill = Treatment)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "C.Small Intestine") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR")) +
  scale_color_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR"))+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,4.2))+
  annotate("text", x = 4.4, y = 400, label=(expression(paste(italic("P"), " = 0.65"))), size = 5, color = "black")+
  annotate("rect", xmin = 4.0, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_si_observed_boxplot

trt_si_shann_boxplot<- ggplot(si_alpha, aes(x= Treatment, y= Shannon, fill = Treatment)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "D.Small Intestine") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 2, label=(expression(paste(italic("P"), " = 0.01"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  annotate("text", x = 1, y = 6.5, label = "a", size = 5, color = "black")+
  annotate("text", x = 2, y = 6.5, label = "a", size = 5, color = "black")+
  annotate("text", x = 3, y = 6.5, label = "b", size = 5, color = "black")+
  annotate("text", x = 4, y = 6.5, label = "b", size = 5, color = "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_si_shann_boxplot


###Large intestine
#set order
li_alpha$Treatment <- factor(li_alpha$Treatment, levels = c("68","68T","67","67T"))

trt_li_observed<- ggplot(li_alpha, aes(x= Treatment, y= Observed, fill = Treatment)) +
  theme_bw() + labs(y= "Observed ASVs", x= "", title = "E.Large Intestine") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR")) +
  scale_color_manual(values = treatment_palette,labels = c("CONREG", "CONERR", "HOTREG","HOTERR"))+
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000),
                     labels = c("500","1,000","1,500","2,000","2,500","3,000","3,500","4,000","4,500","5,000"))+
  coord_cartesian(ylim = c(400,4500), xlim = c(1,4.2))+
  annotate("text", x = 4.4, y = 400, label=(expression(paste(italic("P"), " = 0.35"))), size = 5, color = "black")+
  annotate("rect", xmin = 4.0, xmax = 6, ymin = 0, ymax = 500, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_li_observed

trt_li_shann_boxplot<- ggplot(li_alpha, aes(x= Treatment, y= Shannon, fill = Treatment)) +
  theme_bw() + labs(y= "Shannon's Diversity Index", x= "", title = "F.Large Intestine") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_continuous(breaks = c(2,3,4,5,6,7),
                     labels = c("2","3","4","5","6","7"))+
  coord_cartesian(ylim = c(2,7), xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 2, label=(expression(paste(italic("P"), " = 0.89"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 2.2, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_li_shann_boxplot


#merge plots

trt_alpha <- plot_grid(plot_grid(trt_ru_observed_boxplot, trt_ru_shann_boxplot, 
                                 trt_si_observed_boxplot, trt_si_shann_boxplot, 
                                 trt_li_observed, trt_li_shann_boxplot,
                                 trt_fec_observed_boxplot, trt_fec_shann_boxplot,
                                 ncol = 2,align = "hv", axis = "btrl"), 
                       alpha_trt_legend, align = "hv", axis = "btrl", rel_widths = c(1,0.5) )

trt_alpha

# Beta Diversity ----------------------------------------------------------

####CSS Transform
any(taxa_sums(samples)==0)
samples.css <- phyloseq_transform_css(samples, log = F)
samples.css.df <- as(sample_data(samples.css), "data.frame")


###Split by sample site
#rumen
RU.css <- subset_samples(samples.css, Sample_site== "RU")
RU.css <- prune_taxa(taxa_sums(RU.css)>0, RU.css)
RU.css #84 samples 15,158 taxa
RU.css.df <- as(sample_data(RU.css), "data.frame")

#SI
SI.css <- subset_samples(samples.css, Sample_site== "SI")
SI.css <- prune_taxa(taxa_sums(SI.css)>0, SI.css)
SI.css #73 samples 22,144 taxa
SI.css.df <- as(sample_data(SI.css), "data.frame")

#LI
LI.css <- subset_samples(samples.css, Sample_site== "LI")
LI.css <- prune_taxa(taxa_sums(LI.css)>0, LI.css)
LI.css #79 samples 31,291 taxa
LI.css.df <- as(sample_data(LI.css), "data.frame")

#FEcal
Fec.css <- subset_samples(samples.css, Sample_site== "Fec")
Fec.css <- prune_taxa(taxa_sums(Fec.css)>0, Fec.css)
Fec.css #436 samples 43,209 taxa
Fec.css.df <- as(sample_data(Fec.css), "data.frame")


# Unifrac and Ordination --------------------------------------------------

#all samples
samples.dist <- gunifrac(samples.css)
samples.ord <- ordinate(samples.css, method = "NMDS", distance = samples.dist)

samples.css@sam_data$Sample_site <- factor(samples.css@sam_data$Sample_site, levels = c("RU", "SI", "LI", "Fec"))

#sample type
sample_site_ord <- plot_ordination(samples.css, samples.ord, type = "samples", color = "Sample_site")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Sample_site), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = sample_site_palette, labels = c("Rumen", "Small Intestine", "Large Intestine", "Fecal"))+
  scale_fill_manual(values = sample_site_palette, labels = c("Rumen", "Small Intestine", "Large Intestine", "Fecal"))+
  guides(fill=guide_legend(title = "Sample Location"), color= guide_legend(title = "Sample Location"))+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

sample_site_ord

##treatment
#set order
samples.css@sam_data$Treatment <- factor(samples.css@sam_data$Treatment, levels = c("68","68T","67","67T"))
#calculate segments to make spider plot 
trt_gunifrac.ord <- metaMDS(comm = samples.dist, try = 20, trymax = 500, autotransform = F)
trt_gunifrac.plot <- ordiplot(trt_gunifrac.ord$points)
trt.gunifrac.scrs <- scores(trt_gunifrac.plot, display = "sites")
trt_gunifrac.scrs <- cbind(as.data.frame(trt.gunifrac.scrs), treatment= samples.css.df$Treatment)
trt_gunifrac.cent <- aggregate(cbind(MDS1,MDS2) ~ treatment, data = trt_gunifrac.scrs, FUN = mean)
trt_gunifrac.segs <- merge(trt_gunifrac.scrs, setNames(trt_gunifrac.cent, c("treatment","cMDS1","cMDS2")), by = 'treatment', sort = F)
trt_gunifrac.segs

trt_gunifrac.segs <- factor(trt_gunifrac.segs$treatment)

#treatment normal
treatment_ord <- plot_ordination(samples.css, samples.ord, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= treatment), alpha = 0.1, lty= 2, level = 0.95, size= 1)+
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 26, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
treatment_ord
#spider lines
ggplot(trt_gunifrac.segs, aes(fill = treatment, colour = treatment)) + theme_bw() +
  labs(x="NMDS1", y="NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(aes(x=MDS1,y=MDS2, colour= treatment), alpha = 0.5, shape = 18, size = 5) +
  stat_ellipse(geom = "polygon", aes(x=MDS1,y=MDS2), alpha = c(0.05), lty=2, level = 0.95, linewidth = 0.2) +
  geom_segment(aes(x=MDS1,y=MDS2, xend=cMDS1, yend=cMDS2), alpha = 0.1, linewidth = 1) +
  geom_point(aes(x=cMDS1, y=cMDS2), size = 18, shape = 18) +
  #geom_text(aes(x=cMDS1, y=cMDS2, label = treatment), colour = "white", size = 5) +
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

#liver score
liver_score_ord <- plot_ordination(samples.css, samples.ord, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = simple_palette)+
  scale_fill_manual(values = simple_palette)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))
liver_score_ord



# Subset/make Ordination plots for pub -----------------------------------------------
#Rumen
ru.sub <- subset_samples(samples.css, Sample_site %in% c("RU"))
ru.dist <- gunifrac(ru.sub)
ru.ord <- ordinate(ru.sub, method = "NMDS", distance = ru.dist)

#set order
RU.css@sam_data$Liver_Simple <- factor(RU.css@sam_data$Liver_Simple, levels = c("Edible","A+"))

ru_liver_score_ord <- plot_ordination(RU.css, ru.ord, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "A.Rumen", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = git_simple_palette)+
  scale_fill_manual(values = git_simple_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
ru_liver_score_ord

#set correct trt order
RU.css@sam_data$Treatment <- factor(RU.css@sam_data$Treatment, levels = c("68", "68T", "67", "67T"))

ru_trt_ord <- plot_ordination(RU.css, ru.ord, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "A.Rumen", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Treatment), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 26, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
ru_trt_ord


#SI
si.sub <- subset_samples(samples.css, Sample_site %in% c("SI"))
si.dist <- gunifrac(si.sub)
si.ord <- ordinate(si.sub, method = "NMDS", distance = si.dist)

#set order
SI.css@sam_data$Liver_Simple <- factor(SI.css@sam_data$Liver_Simple, levels = c("Edible", "A+"))

si_liver_score_ord <- plot_ordination(SI.css, si.ord, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "B.Small Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = git_simple_palette)+
  scale_fill_manual(values = git_simple_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
si_liver_score_ord

#set order
SI.css@sam_data$Treatment <- factor(SI.css@sam_data$Treatment, levels = c("68","68T","67","67T"))

si_trt_ord <- plot_ordination(SI.css, si.ord, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "B.Small Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Treatment), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 26, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
si_trt_ord

#LI
li.sub <- subset_samples(samples.css, Sample_site %in% c("LI"))
li.dist <- gunifrac(li.sub)
li.ord <- ordinate(li.sub, method = "NMDS", distance = li.dist)


#set order
LI.css@sam_data$Liver_Simple <- factor(LI.css@sam_data$Liver_Simple, levels = c("Edible","A+"))
#plot as is with outliers included
li_liver_score_ord <- plot_ordination(LI.css, li.ord, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "C.Large Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = LI.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = .1, lty= 1, size= 1)+
  scale_colour_manual(values = git_simple_palette)+
  scale_fill_manual(values = git_simple_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
li_liver_score_ord

#set order
LI.css@sam_data$Treatment <- factor(LI.css@sam_data$Treatment, levels = c("68","68T", "67","67T"))

li_trt_w_out <- plot_ordination(LI.css, li.ord, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "C.Large Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = LI.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Treatment), alpha = .1, lty= 1, size= 1)+
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 26, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
li_trt_w_out

#remove outlier
#remove 6768
li.out.sub <- subset_samples(li.sub, Sample_name != "ICASA-6768-LI")
li.dist.outrm <- gunifrac(li.out.sub)
li.ord.outrm <- ordinate(li.out.sub, method = "NMDS", distance = li.dist.outrm)


#plot as is with outliers included
li_liver_score_ord_outrm <- plot_ordination(LI.css, li.ord.outrm, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "C.Large Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = LI.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = .1, lty= 1, size= 1)+
  scale_colour_manual(values = git_simple_palette)+
  scale_fill_manual(values = git_simple_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
li_liver_score_ord_outrm


li_trt_ord_outrm <- plot_ordination(LI.css, li.ord.outrm, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "C.Large Intestine", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = LI.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Treatment), alpha = .1, lty= 1, size= 1)+
  scale_colour_manual(values = treatment_palette)+
  scale_fill_manual(values = treatment_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
li_trt_ord_outrm

#Fecal
fec.sub <- subset_samples(samples.css, Sample_site %in% c("Fec"))
fec.dist <- gunifrac(fec.sub)
fec.ord <- ordinate(fec.sub, method = "NMDS", distance = fec.dist)

Fec.css@sam_data$Liver_Simple <- factor(Fec.css@sam_data$Liver_Simple, levels = c("Edible","A-", "A", "A+"))

fec_liver_score_ord <- plot_ordination(Fec.css, fec.ord, type = "samples", color = "Liver_Simple")+
  theme_bw()+
  labs(title = "D.Fecal", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = Fec.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Simple), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = simple_palette)+
  scale_fill_manual(values = simple_palette)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, colour = "black"),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

fec_liver_score_ord

#extract legend from fec
ord_legend <- get_legend(fec_liver_score_ord)

#set order
Fec.css@sam_data$Treatment <- factor(Fec.css@sam_data$Treatment, levels = c("68","68T","67","67T"))

fec_trt_ord <- plot_ordination(Fec.css, fec.ord, type = "samples", color = "Treatment")+
  theme_bw()+
  labs(title = "D.Fecal", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = Fec.css@sam_data$Sample_name)+
  stat_ellipse(geom = "polygon", aes(fill= Treatment), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  scale_fill_manual(values = treatment_palette, labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  guides(fill=guide_legend(title = "Treatment"), color= guide_legend(title = "Treatment"))+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 26, colour = "black"),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
fec_trt_ord

trt_ord_legend <- get_legend(fec_trt_ord)

#join th plots liver
ord_plot_final_liver <- plot_grid(plot_grid(ru_liver_score_ord, si_liver_score_ord, li_liver_score_ord, fec_liver_score_ord,
                                                 ncol = 2,align = "hv", axis = "btrl"), 
                                       ord_legend, align = "hv", axis = "btrl", rel_widths = c(2,1))

ord_plot_final_trt

#remove outlier liver
ord_plot_final_no_outlier_liver <- plot_grid(plot_grid(ru_liver_score_ord, si_liver_score_ord, li_liver_score_ord_outrm, fec_liver_score_ord,
                                      ncol = 2,align = "hv", axis = "btrl"), 
                            ord_legend, align = "hv", axis = "btrl", rel_widths = c(1,0.5))

ord_plot_final_no_outlier_liver

#join th plots trt
ord_plot_final_trt <- plot_grid(plot_grid(ru_trt_ord, si_trt_ord, li_trt_w_out, fec_trt_ord,
                                      ncol = 2,align = "hv", axis = "btrl"), 
                            trt_ord_legend, align = "hv", axis = "btrl", rel_widths = c(2,1))

ord_plot_final_trt

#remove outlier liver
ord_plot_final_no_outlier_trt <- plot_grid(plot_grid(ru_trt_ord, si_trt_ord, li_trt_ord_outrm, fec_trt_ord,
                                                 ncol = 2,align = "hv", axis = "btrl"), 
                                       trt_ord_legend, align = "hv", axis = "btrl", rel_widths = c(1,0.5))

ord_plot_final_no_outlier_trt



# Relative Abundance ------------------------------------------------------

rel_abund <- transform_sample_counts(samples.css, function(x) {x/sum(x)} * 100)
ra_family <- tax_glom(rel_abund, taxrank = "Family", NArm = F)
ra_family_filt <- merge_low_abundance(ra_family, threshold = 0.001) 
ra_family_filt #142 taxa 672 samples


# Adonis for all the things -----------------------------------------------
##Start with all the samples

#create data frame
microbiome_df = as.data.frame(as(sample_data(samples.css),"matrix"))

#look at multiple factors
adonis_all <- adonis2(samples.dist ~ Sample_site + Liver_Simple + Treatment, by = "margin", microbiome_df)
adonis_all



#bodysite
adonis_microbiome <- adonis2(samples.dist ~Sample_site, microbiome_df)
adonis_microbiome # p = 0.001 r = 0.35981
#pairwise
adonis_microbiome_body_site <- pairwise.adonis2(samples.dist ~ Sample_site, microbiome_df)
# view the results
adonis_microbiome_body_site
#check dispersion
permutest(betadisper(samples.dist, microbiome_df$Sample_site), pairwise = TRUE) #p = 0.001


#liver score
adonis_liver_score <- adonis2(samples.dist~ Liver_Score, microbiome_df)
adonis_liver_score #p= 0.001
pairwise_liver_score <- pairwise.adonis2(samples.dist ~Liver_Score, microbiome_df)
pairwise_liver_score #edible is different than most things
permutest(betadisper(samples.dist, microbiome_df$Liver_Score), pairwise = T) # p = 0.001

#liver simple
adonis_liver_simple <- adonis2(samples.dist~Liver_Simple, microbiome_df, strata = microbiome_df$Treatment)
adonis_liver_simple #p = 0.001 r = 0.03687
pairwise_liver_simple <- pairwise.adonis2(samples.dist ~ Liver_Simple, microbiome_df)
pairwise_liver_simple # all classes are different from each other
permutest(betadisper(samples.dist, microbiome_df$Liver_Simple), pairwise = T) # p = 0.001

#trt
adonis_trt <- adonis2(samples.dist~Treatment, microbiome_df)
adonis_trt #p = 0.001 r = 0.02385
pairwise.adonis2(samples.dist~ Treatment, microbiome_df)# diet affect between samples
permutest(betadisper(samples.dist, microbiome_df$Treatment), pairwise = T)# p = 0.005

#block
adonis_block <- adonis2(samples.dist~Block, microbiome_df)
adonis_block # p = 0.001
pairwise.adonis2(samples.dist~Block, microbiome_df)
permutest(betadisper(samples.dist, microbiome_df$Block), pairwise = T)# p = 0.005

#rumen score
adonis_rumen <- adonis2(samples.dist~Rumen_Score, microbiome_df)
adonis_rumen # p 0.014
pairwise.adonis2(samples.dist~Rumen_Score, microbiome_df)
permutest(betadisper(samples.dist, microbiome_df$Rumen_Score), pairwise = T) #p = 0.001

# Adonis for rumen -------------------------------------------------------
#make data fram
ru_df <- as.data.frame(as(sample_data(RU.css), "matrix"))

#liver score
adonis_ru_liver_simple <- adonis2(ru.dist~Liver_Simple, ru_df, strata = ru_df$Treatment)
adonis_ru_liver_simple #p = 0.699 r = 0.00852
pairwise_ru_liver_simple <- pairwise.adonis2(ru.dist~Liver_Simple, ru_df)
permutest(betadisper(ru.dist, ru_df$Liver_Simple), pairwise = T) #p = 0.399


#treatment
adonis_ru_trt <- adonis2(ru.dist~Treatment, ru_df, strata = ru_df$Block)# p = 0.001 r =0.10449
adonis_ru_trt
trt_ru_pairwise <- pairwise.adonis2(ru.dist~Treatment, ru_df)
trt_ru_pairwise #diet effect
permutest(betadisper(ru.dist, ru_df$Treatment), pairwise = T)#p = 0.026


# Adonis for Small Intestine ----------------------------------------------
#make data frame
si_df <- as.data.frame(as(sample_data(SI.css), "matrix"))

#Liver Simple
adonis_si_liver_simple <- adonis2(si.dist~Liver_Simple, si_df, strata = si_df$Treatment)
adonis_si_liver_simple # p = 0.017 r = 0.03037 
permutest(betadisper(si.dist, si_df$Liver_Simple), pairwise = T)# p = 0.139

#treatment
adonis_si_trt <- adonis2(si.dist~Treatment, si_df, strata = si_df$Block) #p = 0.052 r = 0.05989
adonis_si_trt
pairwise_trt_si <- pairwise.adonis2(si.dist~Treatment, si_df)
pairwise_trt_si #lots of differences; hard to explain
permutest(betadisper(si.dist, si_df$Treatment), pairwise = T)# p = 0.038


# Adonis for Large Intestine ----------------------------------------------
#make data fram
li_df <- as.data.frame(as(sample_data(LI.css), "matrix"))

#liver simple
adonis_li_liver_simple <- adonis2(li.dist~Liver_Simple, li_df, strata = li_df$Treatment) #p = 0.439 r = 0.01206
adonis_li_liver_simple
permutest(betadisper(li.dist,li_df$Liver_Simple), pairwise = T)# p = 0.151

#treatment
adonis_li_trt <- adonis2(li.dist~Treatment, li_df, strata = li_df$Block) #p = 0.001 r = 0.06681
adonis_li_trt
pairwise_trt_li <- pairwise.adonis2(li.dist~ Treatment,data = li_df, strata = 'Block')
pairwise_trt_li #diet effect
permutest(betadisper(li.dist, li_df$Treatment), pairwise = T)# p = 0.403


# Adonis for Fecal --------------------------------------------------------

 
#make data frame
fecal_df <- as.data.frame(as(sample_data(Fec.css), "matrix"))

#liver score
adonis_fec_liver_score <- adonis2(fec.dist~Liver_Score, fecal_df)
adonis_fec_liver_score #P = 0.051
pairwise_fec_liver_score <- pairwise.adonis2(fec.dist~Liver_Score, fecal_df)
pairwise_fec_liver_score #no difference between A+ categories should probably used the simple score to add power
permutest(betadisper(fec.dist, fecal_df$Liver_Score), pairwise = T) #P = 0.222

#Liver simple
adonis_fec_liver_simple <- adonis2(fec.dist~Liver_Simple, fecal_df, strata = fecal_df$Treatment)
adonis_fec_liver_simple # p = 0.031 r = 0.01201
pairwise_fec_liver_simple <- pairwise.adonis2(fec.dist~Liver_Simple, fecal_df, strata = 'Treatment')
pairwise_fec_liver_simple # every comparison is sig except Ed vs A-
permutest(betadisper(fec.dist, fecal_df$Liver_Simple), pairwise = T) #p = 0.298

#liver binomial
adonis_fec_liver_binom <- adonis2(fec.dist~Liver_Binom, fecal_df)
adonis_fec_liver_binom #p = 0.058
permutest(betadisper(fec.dist, fecal_df$Liver_Binom), pairwise = T) # p = 0.041

#treatment
adonis_fec_trt <- adonis2(fec.dist~Treatment, fecal_df, strata = fecal_df$Block)
adonis_fec_trt #p = 0.001 r = 0.05507
pairwise.adonis2(fec.dist~Treatment, fecal_df, strata = 'Block')# diet affect
permutest(betadisper(fec.dist, fecal_df$Treatment), pairwise = T) #P = 0.024

#block
adonis_fec_block <- adonis2(fec.dist~Block, fecal_df)
adonis_fec_block #p =0.001
pairwise.adonis2(fec.dist~Block, fecal_df)
permutest(betadisper(fec.dist, fecal_df$Block), pairwise = T)# P = 0.21

#rumen score
adonis_fec_rumen <- adonis2(fec.dist~Rumen_Score, fecal_df)
adonis_fec_rumen #p =0.021
pairwise.adonis2(fec.dist~Rumen_Score, fecal_df)
permutest(betadisper(fec.dist, fecal_df$Rumen_Score),pairwise = T)# p = 0.002

# Check individual sample RA for weird things -----------------------------
#Note need to build dendrograms first to make these final plots for the alignment to match: use the respecive dendro.order files
#Rumen
ra_RU <- subset_samples(rel_abund, Sample_site =="RU")
ra_RU <- prune_taxa(taxa_sums(ra_RU) > 0, ra_RU)
ra_RU #15158 taza and 84 samples
RU_phy <- tax_glom(ra_RU, taxrank = "Phylum", NArm = F)
RU_class <- tax_glom(ra_RU, taxrank = "Class", NArm = F)
RU_family <- tax_glom(ra_RU, taxrank = "Family", NArm = F)
RU_genus <- tax_glom(ra_RU, taxrank = "Genus", NArm = F)

RU_family_filt <- merge_low_abundance(RU_family, threshold = 0.1)
RU_family_filt  #42 taxa
RU_family_melt <- psmelt(RU_family_filt)

#right tax table to make colors match
write.csv(tax_table(RU_family_filt), "Ru_tax_table.csv")

#calculate top ten
ru_fam_rank <- RU_family_melt %>% group_by(Family) %>% 
  summarise(mean = mean(Abundance, na.rm = T))


RUPlot <- ggplot(RU_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal()+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = ru.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(values = ru_ra_palette, labels = function(x) {x <- gsub("_"," ",x)
  gsub("^zzz","",x)}) +
  coord_cartesian(clip = "off")+
  theme(#legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 45, hjust = 0.5))
RUPlot

#subset legend 
ru_ra_legend <- get_legend(RUPlot)


#SI
ra_SI <- subset_samples(rel_abund, Sample_site =="SI")
ra_SI <- prune_taxa(taxa_sums(ra_SI) > 0, ra_SI)
ra_SI #22144 taza and 73 samples
SI_phy <- tax_glom(ra_SI, taxrank = "Phylum", NArm = F)
SI_class <- tax_glom(ra_SI, taxrank = "Class", NArm = F)
SI_family <- tax_glom(ra_SI, taxrank = "Family", NArm = F)
SI_genus <- tax_glom(ra_SI, taxrank = "Genus", NArm = F)

SI_family_filt <- merge_low_abundance(SI_family, threshold = 0.1)
SI_family_filt  #36 taxa
SI_family_melt <- psmelt(SI_family_filt)

#write tax table to make colors match
write.csv(tax_table(SI_family_filt),"SI_tax_table.csv")

#calculate top 10
si_fam_rank <- SI_family_melt %>% group_by(Family) %>% 
  summarise(mean = mean(Abundance, na.rm = T))


SIPlot <- ggplot(SI_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal()+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = si.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(values = si_ra_palette, labels = function(x) {x <- gsub("_"," ",x)
  gsub("^zzz","",x)}) +
  coord_cartesian(clip = "off")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 45, hjust = 0.5))
SIPlot
#Extract legend
si_ra_legend <- get_legend(SIPlot)

#LI
ra_LI <- subset_samples(rel_abund, Sample_site =="LI")
ra_LI <- prune_taxa(taxa_sums(ra_LI) > 0, ra_LI)
ra_LI #31291 taxa and 79 samples
LI_phy <- tax_glom(ra_LI, taxrank = "Phylum", NArm = F)
LI_class <- tax_glom(ra_LI, taxrank = "Class", NArm = F)
LI_family <- tax_glom(ra_LI, taxrank = "Family", NArm = F)
LI_genus <- tax_glom(ra_LI, taxrank = "Genus", NArm = F)

LI_family_filt <- merge_low_abundance(LI_family, threshold = 0.1)
LI_family_filt  #45 taxa
LI_family_melt <- psmelt(LI_family_filt)

LI_genus_filt <- merge_low_abundance(LI_genus, threshold = 0.1) %>% 
  psmelt()
#starting_palette <- randomColor(45)
write.csv(starting_palette,"LI_colors.csv")
write.csv(tax_table(LI_family_filt), "LI_family_tax_table.csv")

li_fam_rank <- LI_family_melt %>% group_by(Family) %>% 
  summarise(mean = mean(Abundance, na.rm = T))

LIPlot <- ggplot(LI_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal()+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = li.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(values = li_ra_palette, labels = function(x) {x <- gsub("_"," ",x)
  gsub("^zzz","",x)}) +
  coord_cartesian(clip = "off")+
  theme(legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 45, hjust = 0.5))
LIPlot

li_ra_legend <- get_legend(LIPlot)

#Fecal
ra_fec <- subset_samples(rel_abund, Sample_site =="Fec")
ra_fec <- prune_taxa(taxa_sums(ra_fec) > 0, ra_fec)
ra_fec #43209 taxa and 436 samples
Fec_phy <- tax_glom(ra_fec, taxrank = "Phylum", NArm = F)
Fec_class <- tax_glom(ra_fec, taxrank = "Class", NArm = F)
Fec_family <- tax_glom(ra_fec, taxrank = "Family", NArm = F)
Fec_genus <- tax_glom(ra_fec, taxrank = "Genus", NArm = F)

Fec_family_filt <- merge_low_abundance(Fec_family, threshold = 0.1)
Fec_family_filt  #38 taxa
Fec_family_melt <- psmelt(Fec_family_filt)

#write tax table to make colors match
write.csv(tax_table(Fec_family_filt),"Fec_tax_table.csv")

fec_fam_rank <- Fec_family_melt %>% group_by(Family) %>% 
  summarise(mean = mean(Abundance, na.rm = T))

FecPlot <- ggplot(Fec_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal()+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = fec.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(values = fec_ra_palette, labels = function(x) {x <- gsub("_"," ",x)
  gsub("^zzz","",x)}) +
  coord_cartesian(clip = "off")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 45, hjust = 0.5))
FecPlot

#get legend
fec_ra_legend <- get_legend(FecPlot)


# Merge things together ---------------------------------------------------
##Rumen
#merge legends
ru_legends <- plot_grid(plot_grid(liver_legend, treatment_legend, rel_widths = c(1,1,0.5)), ru_ra_legend,
                        ncol = 3)
ru_legends
#merge plots
ru_ra_plots <- plot_grid(ru_dendro_plot, RUPlot, ncol = 1, align = "hv", axis = "btrl")
ru_ra_plots

#merge both
ru_combo <- plot_grid(ru_ra_plots, ru_legends, rel_widths = c(1,0.5), ncol = 1, align = "v")
ru_combo

##Small intestine
si_legend <- plot_grid(plot_grid(liver_legend, treatment_legend, rel_widths = c(1,1)), si_ra_legend,
                       ncol = 1, rel_heights = c(1,0.5))
si_legend

si_ra_plots <- plot_grid(si_dendro_plot, SIPlot, ncol = 1, align = "hv", axis = "btrl")
si_ra_plots

si_combo <- plot_grid(si_ra_plots, si_legend, rel_widths = c(1,0.5))
si_combo

##Large Intestine
li_legends <- plot_grid(plot_grid(liver_legend, treatment_legend,rel_widths = c(1,1)), li_ra_legend,
                        ncol = 1, rel_heights= c(1,0.5))
li_legends

li_ra_plot <- plot_grid(li_dendro_plot,LIPlot, ncol = 1,
                         align = "hv", axis = "btrl")
li_ra_plot

#merge
li_ra_combo <- plot_grid(li_ra_plot, li_legends, rel_widths = c(1,.5))
li_ra_combo

##Fecal
fec_legends <- plot_grid(plot_grid(fec_liver_legend, treatment_legend,rel_widths = c(1,1)), fec_ra_legend,
                         ncol = 1, rel_heights= c(1,0.5))

fec_ra_plot <- plot_grid(fec_liver_score_dendro, FecPlot, ncol = 1, align = "hv", axis = "btrl")
fec_ra_plot

fec_combo <- plot_grid(fec_ra_plot,fec_legends, rel_widths = c(1,0.5))
fec_combo

# RA for everything  ------------------------------------------------------
#calculate abundance
ra_samples <- prune_taxa(taxa_sums(samples.css) > 0, samples.css)
ra_samples#68284 taxa, 672 samples
samples_phy <- tax_glom(ra_samples, taxrank = "Phylum", NArm = F)
samples_class <- tax_glom(ra_samples, taxrank = "Class", NArm = F)
samples_order <- tax_glom(ra_samples, taxrank = "Order", NArm = F)
samples_family <- tax_glom(ra_samples, taxrank = "Family", NArm = F)
samples_genus <- tax_glom(ra_samples, taxrank = "Genus", NArm = F)

#filter at family level
samples_family_filt <- merge_low_abundance(samples_family, threshold = 0.1)
samples_family_filt  #42 taxa
samples_family_melt <- psmelt(samples_family_filt)

#write tax table to match colors
write.csv(tax_table(samples_family_filt),"all_samples_family.csv")


#note to make final RA plots need to make the dendrogram first to have the dendro.order file to set order of plot.
#plot
samples_ra_Plot <- ggplot(samples_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(values = all_samples_ra_palette, labels = function(x) {x <- gsub("_"," ",x)
  gsub("^zzz","",x)}) +
  coord_cartesian(clip = "off")+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 45, hjust = 0.5))
samples_ra_Plot

all_samples_ra_legend <- get_legend(samples_ra_Plot)

# Dendrogram ------------------------------------------------

#cluster
samples.hclust <- hclust(samples.dist, method = "ward.D2")
samples.dendro.plot <- plot(samples.hclust)
samples.dendro <- as.dendrogram(samples.hclust)
samples.dendro.data <- dendro_data(samples.dendro, type = "rectangle")
metadata_for_dendro.samples <- as_tibble(samples.css@sam_data)
samples.dendro.data$labels <- samples.dendro.data$labels %>% 
  left_join(metadata_for_dendro.samples, by = c("label" = "Sample_name"))
samples.dendro.data$labels
samples.dendro.order <- samples.dendro.data$labels$label
samples.dendro.order

#set label order
samples.dendro.data$labels$Liver_Simple <- factor(samples.dendro.data$labels$Liver_Simple, levels = c("Edible", "A-", "A", "A+"))
samples.dendro.data$labels$Treatment <- factor(samples.dendro.data$labels$Treatment, levels = c("68", "68T", "67", "67T"))
samples.dendro.data$labels$Sample_site <- factor(samples.dendro.data$labels$Sample_site, levels = c("RU", "SI", "LI", "Fec"))

#liver and body site
microbiome_dendro_plot <- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, fill = Treatment,), color = "black", size = 7, shape = 21, position = position_nudge(y=-.05))+
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, fill = Sample_site), color = "black", size = 7, shape = 24, position = position_nudge(y= -0.15))+
  scale_x_discrete(limits = li.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Ledgend", values = full_dendro_palette, guide = guide_legend(order = 1))+
  coord_cartesian(clip = "off")+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
microbiome_dendro_plot

#make legend for sample site since this is only place they are all together
#sample site legend
site_legend_plot <- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, fill = Sample_site), color = "black", size = 9, shape = 24, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0))+
  scale_fill_manual(name = "Sample Type", values = sample_site_palette, guide = guide_legend(order = 1), labels = c("Rumen", "Small Intestine", "Large Intestine", "Fecal"))+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2, "cm"))

site_legend_plot

site_legend <- get_legend(site_legend_plot)


# Merge all samples dendro and RA -----------------------------------------
#merge plots
all_samples_plot <- plot_grid(microbiome_dendro_plot,samples_ra_Plot, ncol = 1, align = "hv", axis = "btrl")
all_samples_plot

#merge legends
all_samples_legend <- plot_grid(plot_grid(liver_legend, treatment_legend, site_legend, rel_widths = c(.5,.5,.5), ncol = 3),
                                all_samples_ra_legend,ncol = 1, rel_heights= c(1,0.5))
all_samples_legend

#merge legend with plot
all_samples_paired_ra_final <- plot_grid(all_samples_plot,all_samples_legend, rel_widths = c(1,0.35))
all_samples_paired_ra_final



# Random Dendros just to see ----------------------------------------------


#just liver
liver_dendro <- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, colour = Liver_Score), size = 9, shape = 15, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = liver_palette)+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
liver_dendro

#sample Site
sample_site_dendro<- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, colour = Sample_site), size = 9, shape = 15, position = position_nudge(y=-.2)) +
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = sample_site_palette)+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
sample_site_dendro

#treatment
treatment_dendro<- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, colour = Treatment), size = 9, shape = 15, position = position_nudge(y=-.2)) +
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = trt_palette)+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
treatment_dendro

#liver binom
fec_liver_binom_dendro <- ggplot(fec.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, colour = Liver_Binom), size = 9, shape = 15, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = fec.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = binom_palette)+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
fec_liver_binom_dendro

#set random variables as factors
fec.dendro.data$labels$Block <- as.factor(fec.dendro.data$labels$Block)
fec.dendro.data$labels$Rumen_Score <- as.factor(fec.dendro.data$labels$Rumen_Score)

#look at paired things to see if we can explain anything
fec_random_dendro <- ggplot(fec.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, colour = Liver_Simple), size = 9, shape = 15, position = position_nudge())+
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, colour = Rumen_Score), size = 9, shape = 15, position = position_nudge(y=-0.5))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = fec.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = randomColor(19))+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
fec_random_dendro



# Subset for dendro by region --------------------------------------------------

###Rumen
ru.hclust <- hclust(ru.dist, method = "ward.D2")
ru.dendro.plot <- plot(ru.hclust)
ru.dendro <- as.dendrogram(ru.hclust)
ru.dendro.data <- dendro_data(ru.dendro, type = "rectangle")
ru_metadata_for_dendro.samples <- as_tibble(RU.css@sam_data)
ru.dendro.data$labels <- ru.dendro.data$labels %>% 
  left_join(ru_metadata_for_dendro.samples, by = c("label" = "Sample_name"))
ru.dendro.data$labels
ru.dendro.order <- ru.dendro.data$labels$label
ru.dendro.order

#set order of factors for the legend
ru.dendro.data$labels$Liver_Simple <- factor(ru.dendro.data$labels$Liver_Simple, levels = c("Edible", "A-", "A", "A+"))
ru.dendro.data$labels$Treatment <- factor(ru.dendro.data$labels$Treatment, levels = c("68", "68T", "67", "67T"))

#plot

ru_dendro_plot<- ggplot(ru.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = ru.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  geom_point(data = ru.dendro.data$labels, aes(x=x, y=y, fill = Treatment,), color = "black", size = 7, shape = 21, position = position_nudge(y=-.05))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = ru.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Ledgend", values = git_simple_dendro_palette, guide = guide_legend(order = 1))+
  guides(fill = guide_legend(overide.aes = list(shape = c(22,21), title = "Liver Score")))+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
ru_dendro_plot

##Make the legends Note: cant make the legend very easy and keep the shapes i want since i need to use fill for both
##Solution: make each legend and extract from a plot I dont want to use then merge the legends at the end for the final figure

#liver legend 
liver_legend_plot <- ggplot(ru.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = ru.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = ru.dendro.order, expand = c(0.035,0,0,0))+
  scale_fill_manual(name = "Liver Score", values = git_simple_dendro_palette, guide = guide_legend(order = 1))+
  guides(fill = guide_legend(overide.aes = list(shape = c(22,21), title = "Liver Score")))+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2, "cm"))
liver_legend_plot

liver_legend <- get_legend(liver_legend_plot)

#treatment legend
treat_legend_plot <- ggplot(ru.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = ru.dendro.data$labels, aes(x=x, y=y, fill = Treatment), color = "black", size = 9, shape = 21, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = ru.dendro.order, expand = c(0.035,0,0,0))+
  scale_fill_manual(name = "Treatment", values = treatment_palette, guide = guide_legend(order = 1), labels = c("CONREG", "CONERR", "HOTREG", "HOTERR"))+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2, "cm"))
treat_legend_plot

treatment_legend <- get_legend(treat_legend_plot)



###Small intestine
#cluster
si.hclust <- hclust(si.dist, method = "ward.D2")
si.dendro.plot <- plot(si.hclust)
si.dendro <- as.dendrogram(si.hclust)
si.dendro.data <- dendro_data(si.dendro, type = "rectangle")
si_metadata_for_dendro.samples <- as_tibble(SI.css@sam_data)
si.dendro.data$labels <- si.dendro.data$labels %>% 
  left_join(si_metadata_for_dendro.samples, by = c("label" = "Sample_name"))
si.dendro.data$labels
si.dendro.order <- si.dendro.data$labels$label
si.dendro.order

#set order of factors for the legend
si.dendro.data$labels$Liver_Simple <- factor(si.dendro.data$labels$Liver_Simple, levels = c("Edible", "A-", "A", "A+"))
si.dendro.data$labels$Treatment <- factor(si.dendro.data$labels$Treatment, levels = c("68", "68T", "67", "67T"))

#plot

si_dendro_plot<- ggplot(si.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = si.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  geom_point(data = si.dendro.data$labels, aes(x=x, y=y, fill = Treatment,), color = "black", size = 7, shape = 21, position = position_nudge(y=-.05))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = si.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Ledgend", values = git_simple_dendro_palette, guide = guide_legend(order = 1))+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
si_dendro_plot

###Large Intestine
#cluster
li.hclust <- hclust(li.dist, method = "ward.D2")
li.dendro.plot <- plot(li.hclust)
li.dendro <- as.dendrogram(li.hclust)
li.dendro.data <- dendro_data(li.dendro, type = "rectangle")
li_metadata_for_dendro.samples <- as_tibble(LI.css@sam_data)
li.dendro.data$labels <- li.dendro.data$labels %>% 
  left_join(li_metadata_for_dendro.samples, by = c("label" = "Sample_name"))
li.dendro.data$labels
li.dendro.order <- li.dendro.data$labels$label
li.dendro.order

#set order of factors for the legend
li.dendro.data$labels$Liver_Simple <- factor(li.dendro.data$labels$Liver_Simple, levels = c("Edible", "A-", "A", "A+"))
li.dendro.data$labels$Treatment <- factor(li.dendro.data$labels$Treatment, levels = c("68", "68T", "67", "67T"))

#plot
li_dendro_plot<- ggplot(li.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = li.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  geom_point(data = li.dendro.data$labels, aes(x=x, y=y, fill = Treatment,), color = "black", size = 7, shape = 21, position = position_nudge(y=-.05))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = li.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Ledgend", values = git_simple_dendro_palette, guide = guide_legend(order = 1))+
  coord_cartesian(clip = "off")+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,.2),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
li_dendro_plot

###Fecal
#cluster
fec.hclust <- hclust(fec.dist, method = "ward.D2")
fec.dendro.plot <- plot(fec.hclust)
fec.dendro <- as.dendrogram(fec.hclust)
fec.dendro.data <- dendro_data(fec.dendro, type = "rectangle")
fec_metadata_for_dendro.samples <- as_tibble(Fec.css@sam_data)
fec.dendro.data$labels <- fec.dendro.data$labels %>% 
  left_join(fec_metadata_for_dendro.samples, by = c("label" = "Sample_name"))
fec.dendro.data$labels
fec.dendro.order <- fec.dendro.data$labels$label
fec.dendro.order

#set order of factors for the legend
fec.dendro.data$labels$Liver_Simple <- factor(fec.dendro.data$labels$Liver_Simple, levels = c("Edible", "A-", "A", "A+"))
fec.dendro.data$labels$Treatment <- factor(fec.dendro.data$labels$Treatment, levels = c("68", "68T", "67", "67T"))


#liver score and treatment
fec_liver_score_dendro <- ggplot(fec.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Wards Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, fill = Treatment,), color = "black", size = 7, shape = 21, position = position_nudge(y=-.05))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = fec.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Ledgend", values = fecal_dendro_palette, guide = guide_legend(order = 1))+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
fec_liver_score_dendro

#only liver simple to extract legend
fec_liver_simple_dendro <- ggplot(fec.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = fec.dendro.data$labels, aes(x=x, y=y, fill = Liver_Simple), color = "black", size = 9, shape = 22, position = position_nudge())+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = fec.dendro.order, expand = c(0.01,0,0,0))+
  scale_fill_manual(name = "Liver Score", values = fecal_dendro_palette, guide = guide_legend(order = 1))+
  theme(plot.title = element_text(size = 20),
         #legend.position = "none",
         plot.margin = unit(c(0,0,0,0),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line.y = element_line(size = 0.75, colour = "black"),
         axis.title.y = element_text(size = 24),
         axis.text.y = element_text(size = 12, colour = "black"),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.y = element_line(colour = "black", size = 0.7),
         legend.text = element_text(size = 40),
         legend.title = element_text(size = 45, hjust = 0.5),
         legend.key.height = unit(2, "cm"))
fec_liver_simple_dendro

fec_liver_legend <- get_legend(fec_liver_simple_dendro)



# Firmicutes to Bacteroidetes Ratio ---------------------------------------

#Subset frimicutes
firm_bac_sub <- subset_taxa(rel_abund, Phylum %in% c("Bacteroidota", "Firmicutes"))

#subset fecal
firm_bac_fec <- subset_samples(firm_bac_sub, Sample_site == "Fec")

#calculate average Firm
firm_fec <- subset_taxa(firm_bac_fec, Phylum == "Firmicutes")

firm_avg <- firm_fec %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_firm = sum(Abundance, na.rm = T),
            sd_firm = sd(Abundance, na.rm = T)) 
#calculate average Bac
bac_fec <- subset_taxa(firm_bac_fec, Phylum == "Bacteroidota")

bac_avg <- bac_fec %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_bac = sum(Abundance, na.rm = T),
            sd_bac = sd(Abundance, na.rm = T))
#subset metadata
fec_metadata <- subset(metadata.df,Sample_site=="Fec")


#join the bac average and firm average
f_to_b_fecal <- fec_metadata %>% left_join(firm_avg, by = "Sample_name") %>% left_join(bac_avg, by = "Sample_name")

#calculate F:B
f_to_b_fecal$ratio <- (f_to_b_fecal$sum_firm/f_to_b_fecal$sum_bac)

#test for a difference in liver
kruskal.test(f_to_b_fecal$ratio,f_to_b_fecal$Liver_Simple) #p = 0.4011

#test for difference in treatment
kruskal.test(f_to_b_fecal$ratio,f_to_b_fecal$Treatment) # p < 0.001
pairwise.wilcox.test(f_to_b_fecal$ratio,f_to_b_fecal$Treatment, p.adjust.method = "BH")


#set order of variable
f_to_b_fecal$Liver_Simple <- factor(f_to_b_fecal$Liver_Simple, levels = c("Edible","A-","A","A+"))

#Liver simple F:B
F_to_b_fec_boxplot <- ggplot(f_to_b_fecal, aes(x= Liver_Simple, y= ratio, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "D.Fecal") +
  scale_x_discrete(limits=c("Edible", "A-", "A", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = simple_palette) +
  scale_color_manual(values = simple_palette)+
  scale_y_log10()+
  coord_cartesian(xlim = c(1,4))+
  annotate("text", x = 4.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.40"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 5, ymin = 0, ymax = 0.75, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

F_to_b_fec_boxplot

firm_ratio_legend <- get_legend(F_to_b_fec_boxplot)

#treatment
#set order
f_to_b_fecal$Treatment <- factor(f_to_b_fecal$Treatment, levels = c("68","68T","67","67T"))

trt_f_to_b_fec <- ggplot(f_to_b_fecal, aes(x= Treatment, y= ratio, fill = Treatment)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "D.Fecal") +
  scale_x_discrete(limits=c("68", "68T", "67", "67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 0.5, label=(expression(paste(italic("P"), " < 0.001"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax =.75, color = "black", alpha = 0)+
  annotate("text", x = 1, y = 80, label = "a", size = 5, color = "black")+
  annotate("text", x = 2, y = 80, label = "a", size = 5, color = "black")+
  annotate("text", x = 3, y = 80, label = "b", size = 5, color = "black")+
  annotate("text", x = 4, y = 80, label = "b", size = 5, color = "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_f_to_b_fec

###Rumen

#subset rumen
firm_bac_ru <- subset_samples(firm_bac_sub, Sample_site == "RU")

#calculate average Firm
firm_ru <- subset_taxa(firm_bac_ru, Phylum == "Firmicutes")

firm_avg_ru <- firm_ru %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_firm = sum(Abundance, na.rm = T),
            sd_firm = sd(Abundance, na.rm = T)) 
#calculate average Bac
bac_ru <- subset_taxa(firm_bac_ru, Phylum == "Bacteroidota")

bac_avg_ru <- bac_ru %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_bac = sum(Abundance, na.rm = T),
            sd_bac = sd(Abundance, na.rm = T))
#subset metadata
ru_metadata <- subset(metadata.df,Sample_site=="RU")


#join the bac average and firm average
f_to_b_ru <- ru_metadata %>% left_join(firm_avg_ru, by = "Sample_name") %>% left_join(bac_avg_ru, by = "Sample_name")

#calculate F:B
f_to_b_ru$ratio <- (f_to_b_ru$sum_firm/f_to_b_ru$sum_bac)

#test for a difference liver
kruskal.test(f_to_b_ru$ratio,f_to_b_ru$Liver_Simple) #p = 0.99

#test for difference trt
kruskal.test(f_to_b_ru$ratio, f_to_b_ru$Treatment)#p = 0.4188

#set order of variable
f_to_b_ru$Liver_Simple <- factor(f_to_b_ru$Liver_Simple, levels = c("Edible","A-","A","A+"))

#Liver simple F:B
F_to_b_ru_boxplot <- ggplot(f_to_b_ru, aes(x= Liver_Simple, y= ratio, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "A.Rumen") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,2))+
  annotate("text", x = 2.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.99"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 6, ymin = 0, ymax = 0.53, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
F_to_b_ru_boxplot

#treatment
#set order 
f_to_b_ru$Treatment <- factor(f_to_b_ru$Treatment, levels = c("68","68T","67","67T"))

trt_f_to_b_ru <- ggplot(f_to_b_ru, aes(x= Treatment, y= ratio, fill = Treatment)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "A.Rumen") +
  scale_x_discrete(limits=c("68", "68T", "67", "67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.42"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 0.53, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_f_to_b_ru

###Small intestine
#subset small intestine
firm_bac_si <- subset_samples(firm_bac_sub, Sample_site == "SI")

#calculate average Firm
firm_si <- subset_taxa(firm_bac_si, Phylum == "Firmicutes")

firm_avg_si <- firm_si %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_firm = sum(Abundance, na.rm = T),
            sd_firm = sd(Abundance, na.rm = T)) 
#subset bac
bac_si <- subset_taxa(firm_bac_si, Phylum == "Bacteroidota")
#calculate avg for bac
bac_avg_si <- bac_si %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_bac = sum(Abundance, na.rm = T),
            sd_bac = sd(Abundance, na.rm = T))
#subset metadata
si_metadata <- subset(metadata.df,Sample_site=="SI")


#join the bac average and firm average
f_to_b_si <- si_metadata %>% left_join(firm_avg_si, by = "Sample_name") %>% left_join(bac_avg_si, by = "Sample_name")

#calculate F:B
f_to_b_si$ratio <- (f_to_b_si$sum_firm/f_to_b_si$sum_bac)

#test for a difference liver score
kruskal.test(f_to_b_si$ratio,f_to_b_si$Liver_Simple) #p = 0.004

#test for difference treatment
kruskal.test(f_to_b_si$ratio, f_to_b_si$Treatment)#p = 0.763

#set order of variable
f_to_b_si$Liver_Simple <- factor(f_to_b_si$Liver_Simple, levels = c("Edible","A-","A","A+"))

#Liver simple F:B
F_to_b_si_boxplot <- ggplot(f_to_b_si, aes(x= Liver_Simple, y= ratio, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "B.Small Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,2))+
  annotate("text", x = 2.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.004"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 5, ymin = 0, ymax = 0.75, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

F_to_b_si_boxplot

#treatment
#set order
f_to_b_si$Treatment <- factor(f_to_b_si$Treatment, levels = c("68","68T","67","67T"))

trt_f_to_b_si <- ggplot(f_to_b_si, aes(x= Treatment, y= ratio, fill = Treatment)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "B.Small Intestine") +
  scale_x_discrete(limits=c("68","68T","67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = .5, label=(expression(paste(italic("P"), " = 0.76"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 0.75, color = "black", alpha = 0)+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_f_to_b_si

###Large intestine
#subset rumen
firm_bac_li <- subset_samples(firm_bac_sub, Sample_site == "LI")

#calculate average Firm
firm_li <- subset_taxa(firm_bac_li, Phylum == "Firmicutes")

firm_avg_li <- firm_li %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_firm = sum(Abundance, na.rm = T),
            sd_firm = sd(Abundance, na.rm = T)) 
#calculate average Bac
bac_li <- subset_taxa(firm_bac_li, Phylum == "Bacteroidota")

bac_avg_li <- bac_li %>% psmelt() %>% group_by(Sample_name) %>% 
  summarise(sum_bac = sum(Abundance, na.rm = T),
            sd_bac = sd(Abundance, na.rm = T))
#subset metadata
li_metadata <- subset(metadata.df,Sample_site=="LI")


#join the bac average and firm average
f_to_b_li <- li_metadata %>% left_join(firm_avg_li, by = "Sample_name") %>% left_join(bac_avg_li, by = "Sample_name")

#calculate F:B
f_to_b_li$ratio <- (f_to_b_li$sum_firm/f_to_b_li$sum_bac)

#test for a difference liver score
kruskal.test(f_to_b_li$ratio,f_to_b_li$Liver_Simple) #p = 0.50

#test for difference treatment
kruskal.test(f_to_b_li$ratio, f_to_b_li$Treatment)# p = 0.02801
pairwise.wilcox.test(f_to_b_li$ratio, f_to_b_li$Treatment, p.adjust.method = "BH")


#set order of variable
f_to_b_li$Liver_Simple <- factor(f_to_b_li$Liver_Simple, levels = c("Edible","A-","A","A+"))

#Liver simple F:B
F_to_b_li_boxplot <- ggplot(f_to_b_li, aes(x= Liver_Simple, y= ratio, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "C. Large Intestine") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette) +
  scale_color_manual(values = git_simple_palette)+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,2))+
  annotate("text", x = 2.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.50"))), size = 5, color = "black")+
  annotate("rect", xmin = 2, xmax = 5, ymin = 0, ymax = 0.75, color = "black", alpha = 0)+
  guides(fill=guide_legend(title = "Liver Score"), color= guide_legend(title = "Liver Score"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
F_to_b_li_boxplot

#Treatment F:B
#set order
f_to_b_li$Treatment <- factor(f_to_b_li$Treatment, levels = c("68","68T","67","67T"))

trt_f_to_b_li <- ggplot(f_to_b_li, aes(x= Treatment, y= ratio, fill = Treatment)) +
  theme_bw() + labs(y= "Firmicutes:Bacteroidota Ratio", x= "", title = "C. Large Intestine") +
  scale_x_discrete(limits=c("68", "68T", "67","67T"))+
  geom_boxplot(aes(color= Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Treatment), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = treatment_palette) +
  scale_color_manual(values = treatment_palette)+
  scale_y_continuous(breaks = c(0,150,300,450),
                     labels = c("0","150","300","450"))+
  scale_y_log10()+
  coord_cartesian( xlim = c(1,4.2))+
  annotate("text", x = 4.35, y = 0.5, label=(expression(paste(italic("P"), " = 0.03"))), size = 5, color = "black")+
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 0.75, color = "black", alpha = 0)+
  annotate("text", x =1, y = 225, label= "ab",size = 5, color = "black")+
  annotate("text", x = 2, y = 225, label = "a", size = 5, color = "black")+
  annotate("text", x = 3, y = 225, label = "b", size = 5, color = "black")+
  annotate("text", x = 4, y = 225, label = "ab", size = 5, color = "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 20, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
trt_f_to_b_li

#join plots into one figure

f_to_b_plot <- plot_grid(plot_grid(F_to_b_ru_boxplot, F_to_b_si_boxplot, F_to_b_li_boxplot, F_to_b_fec_boxplot,
                                   ncol = 2,align = "hv", axis = "btrl" ),
                         firm_ratio_legend, ncol = 2,align = "hv", axis = "btrl", rel_widths = c(1,0.5))

f_to_b_plot

trt_f_to_b <- plot_grid(plot_grid(trt_f_to_b_ru, trt_f_to_b_si, trt_f_to_b_li, trt_f_to_b_fec,
                                  ncol = 2,align = "hv", axis = "btrl" ),
                        alpha_trt_legend, ncol = 2,align = "hv", axis = "btrl", rel_widths = c(1,0.5))


trt_f_to_b
# ANCOM for all -----------------------------------------------------------
ancom_samples <- samples
ancom_samples #672 taxa 

#reorder liver simple to make edible the comparator
ancom_samples@sam_data$Liver_Simple <- factor(ancom_samples@sam_data$Liver_Simple, levels = c("Edible", "A-","A", "A+"))
ancom_samples@sam_data$Liver_Simple

ancom_samples@sam_data

##running ancombc on the variable of interest (liver simple) with score and trt and Sample site as fixed no random
##dunnets turned on
ancombc_output_all_samples <-ancombc2(data= ancom_samples, assay_name = "counts", tax_level = "Family",
                                           fix_formula = "Liver_Simple + Treatment + Sample_site", rand_formula = NULL, 
                                           p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                           group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                                           global = TRUE, pairwise = TRUE,
                                           dunnet = TRUE, trend = FALSE,
                                           iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                           em_control = list(tol = 1e-5, max_iter = 100),
                                           lme_control = lme4::lmerControl(),
                                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                           trend_control = TRUE)

#extract results
ancom_res <- ancombc_output_all_samples$res


# ANCOM for Fecal-------------------------------------------------------------------
#subset only fecal samples
ancom_fec <- subset_samples(samples, Sample_site== "Fec")
#Subset family counts must have untransformed counts
ancom_fec_family <- tax_glom(ancom_fec, taxrank = "Family", NArm = F)
#reorder liver simple to make edible the comparator
ancom_fec_family@sam_data$Liver_Simple <- factor(ancom_fec_family@sam_data$Liver_Simple, levels = c("Edible", "A-","A", "A+"))
ancom_fec_family@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_fec_family <-ancombc2(data= ancom_fec_family, assay_name = "counts", tax_level = "Family",
                                 fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                 p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                 group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                 global = TRUE, pairwise = TRUE,
                                 dunnet = TRUE, trend = FALSE,
                                 iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                 trend_control = TRUE)

# extract results from pairwise comparisons 
fec_res_pair <- ancombc_output_fec_family$res_pair ##lfc is log fold change

#### melt into long format for plotting 
ancom_liver_simple_fec <- fec_res_pair %>% 
  filter (`diff_Liver_SimpleA-`== 1 | diff_Liver_SimpleA== 1 | `diff_Liver_SimpleA+`== 1)%>% #filtering the significant differences(1 is significant)
  mutate(across(starts_with("lfc_Liver_Simple"), 
                ~ifelse(get(gsub("lfc_", "diff_", cur_column()))== 1, .x, 0)),
         across(starts_with("lfc_Liver_Simple"), 
                ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), 
                ~ifelse(.x==1 & get(gsub("diff_", "passed_ss_", cur_column()))==1,
                        "darkred","darkgrey"),.names = "{.col}_color"))%>%
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_") %>%
  pivot_longer(cols=contains("_color"), names_to = "color_group", values_to ="color", names_prefix = "diff_")%>%
  filter(gsub("_rounded", "", group)== gsub ("_color", "",color_group))%>%
  select(-color_group) %>%
  arrange(taxon)
ancom_liver_simple_fec

#Ensure the "group" column matches between LFC values and colors for direct comparison 
## we're getting rid of _rounded suffix using sub command
ancom_liver_simple_fec$group ##want to get rid of "_rounded"
ancom_liver_simple_fec$group<- sub("_rounded", "", ancom_liver_simple_fec$group) 
ancom_liver_simple_fec$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_liver_simple_fec <- ancom_liver_simple_fec %>%
  mutate(group= case_when(
    group == "Liver_SimpleA-" ~ "Edible vs A-",
    group == "Liver_SimpleA" ~ "Edible vs A",
    group == "Liver_SimpleA+" ~ "Edible vs A+",
    TRUE ~ group ##keeps original name for groups not specified (individual)
  ))
ancom_liver_simple_fec$group ##Now the group names are shorter and more manageable
ancom_liver_simple_fec<- ancom_liver_simple_fec %>%
  rename(Family = taxon) ##This ancombc was done at the family level, changing the column name to that

##Dunns
#extract dunns
fec_res_dunn <- ancombc_output_fec_family$res_dunn

#### melt into long format for plotting Dunns
ancom_liver_simple_fec_dunn <- fec_res_dunn %>% 
  filter (`diff_Liver_SimpleA-`== 1 | diff_Liver_SimpleA== 1 | `diff_Liver_SimpleA+`== 1)%>% #filtering the significant differences(1 is significant)
  mutate(across(starts_with("lfc_Liver_Simple"), 
                ~ifelse(get(gsub("lfc_", "diff_", cur_column()))== 1, .x, 0)),
         across(starts_with("lfc_Liver_Simple"), 
                ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), 
                ~ifelse(.x==1 & get(gsub("diff_", "passed_ss_", cur_column()))==1,
                        "darkred","darkgrey"),.names = "{.col}_color"))%>%
  pivot_longer(cols=contains("_rounded"), names_to = "group", values_to = "value", names_prefix= "lfc_") %>%
  pivot_longer(cols=contains("_color"), names_to = "color_group", values_to ="color", names_prefix = "diff_")%>%
  filter(gsub("_rounded", "", group)== gsub ("_color", "",color_group))%>%
  select(-color_group) %>%
  arrange(taxon)
ancom_liver_simple_fec_dunn

#Ensure the "group" column matches between LFC values and colors for direct comparison 
## we're getting rid of _rounded suffix using sub command
ancom_liver_simple_fec_dunn$group ##want to get rid of "_rounded"
ancom_liver_simple_fec_dunn$group<- sub("_rounded", "", ancom_liver_simple_fec_dunn$group) 
ancom_liver_simple_fec_dunn$group #names don't have "grounded" anymore

##rework our group names so they're shorter and more manageable 
ancom_liver_simple_fec_dunn <- ancom_liver_simple_fec_dunn %>%
  mutate(group= case_when(
    group == "Liver_SimpleA-" ~ "Edible vs A-",
    group == "Liver_SimpleA" ~ "Edible vs A",
    group == "Liver_SimpleA+" ~ "Edible vs A+",
    TRUE ~ group ##keeps original name for groups not specified (individual)
  ))
ancom_liver_simple_fec_dunn$group ##Now the group names are shorter and more manageable
ancom_liver_simple_fec_dunn<- ancom_liver_simple_fec_dunn %>%
  rename(Family = taxon) ##This ancombc was done at the family level, changing the column name to that


# Try Using Maggies Code but using the ANCOM analysis from above ----------
#subset data from dunns test
new_fec_ancom <- ancombc_output_fec_family$res_dunn

##Data prep for heatmap
df_fec_liver_simple_dunn <- new_fec_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA-`==1|
                  diff_Liver_SimpleA==1|
                  `diff_Liver_SimpleA+`==1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA-`= ifelse(`diff_Liver_SimpleA-`==1,round(`lfc_Liver_SimpleA-`, 2), `lfc_Liver_SimpleA-`),
                lfc_Liver_SimpleA= ifelse(diff_Liver_SimpleA==1, round(lfc_Liver_SimpleA, 2), lfc_Liver_SimpleA),
                `lfc_Liver_SimpleA+`= ifelse(`diff_Liver_SimpleA+`== 1, round(`lfc_Liver_SimpleA+`, 2), `lfc_Liver_SimpleA+`)) %>% 
  dplyr::transmute(taxon,
            "Liver_SimpleA-"= round(`lfc_Liver_SimpleA-`),
            "Liver_SimpleA"= round(lfc_Liver_SimpleA),
            "Liver_SimpleA+"= round(`lfc_Liver_SimpleA+`)) %>% 
  tidyr::pivot_longer(cols = "Liver_SimpleA-":"Liver_SimpleA":"Liver_SimpleA+", names_to = "group", values_to = "value") %>% 
  dplyr::arrange(taxon)

#sensititvity testing
df_ses_fec_Liver_simple <- new_fec_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA-`== 1 |
                  diff_Liver_SimpleA == 1 |
                  `diff_Liver_SimpleA+` == 1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA-`= ifelse(`passed_ss_Liver_SimpleA-`==1 & `diff_Liver_SimpleA-`==1, "white", "black"),
                lfc_Liver_SimpleA= ifelse(passed_ss_Liver_SimpleA==1 & diff_Liver_SimpleA ==1, "white", "black"),
                `lfc_Liver_SimpleA+`= ifelse(`passed_ss_Liver_SimpleA+`==1 & `diff_Liver_SimpleA+`==1, "white", "black"),) %>% 
  dplyr::transmute(taxon,
                  "Liver_SimpleA-"= `lfc_Liver_SimpleA-`,
                  "Liver_SimpleA"= lfc_Liver_SimpleA,
                  "Liver_SimpleA+"= `lfc_Liver_SimpleA+`) %>% 
  tidyr::pivot_longer(cols = "Liver_SimpleA-":"Liver_SimpleA":"Liver_SimpleA+", names_to = "group", values_to = "color") %>% 
  dplyr::arrange(taxon)

#join sensitivity with dunn
df_fec_liver_simple_dunn1 <- df_fec_liver_simple_dunn %>% 
  dplyr::left_join(df_ses_fec_Liver_simple, by = c("taxon", "group"))

#change order of X-value 
df_fec_liver_simple_dunn1$group <- factor(df_fec_liver_simple_dunn1, levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"))


#Assign your map colors to values
lo = floor(min(df_fec_liver_simple_dunn1$value))
up = ceiling(max(df_fec_liver_simple_dunn1$value))
mid = 0

fec_liver_simple_heatmap = df_fec_liver_simple_dunn1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "black", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold change compared to Edible") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fec_liver_simple_heatmap


##Tornado plot 

#do a tornado
volcano_fec_ancom <- ancombc_output_fec_family$res_dunn
View(volcano_fec_ancom)

#clean data
df_fec_long_volcano<- volcano_fec_ancom %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>%
pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
select(-color_group) %>%
arrange(taxon)%>%
mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fec_long_volcano$group <- sub("_rounded", "", df_fec_long_volcano$group)
#Just double check that your groups turned out correct
df_fec_long_volcano$group
unique(df_fec_long_volcano$group)

#Filter the results to only include significant features
significant_features = df_fec_long_volcano$taxon

# Extract the names of the significant features
significant_feature_names = as.character(significant_features)

# Extract the bias_correct_log_table from the ancombc output
log_table = ancombc_output_fec_family$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant = log_table[rownames(log_table) %in% significant_feature_names, ]

# Convert the row names into a column
log_table_significant$Feature = rownames(log_table_significant)

# Reshape the data into a format suitable for ggplot2
df_fec = reshape2::melt(log_table_significant, id.vars = "Feature")

# Calculating the average value for each Feature
median_values <- df_fec %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
median_values

df_fec$value

# Merge the df data frame with the sample data
fec_merged_df = merge(median_values, df_fec_long_volcano, by.x = "Feature", by.y = "taxon")

colnames(fec_merged_df)

fec_merged_df$Feature <- factor(fec_merged_df$Feature, levels = unique(fec_merged_df$Feature))

# Update the tornado_plot with new specifications
fec_tornado_plot <- ggplot(fec_merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "Fecal") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-1, 1)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code
  

fec_tornado_plot


# Fecal ANCOM Filtered ----------------------------------------------------

#filter fecal to average of 5 counts per sample
##436 samples so if it occurs 5 times per sample would be a sum of 2,180 counts

#Look at whats in unfiltered object
ancom_fec_family #492 taxa 436 samples

ancom_fec_family_5_min <- prune_taxa(taxa_sums(ancom_fec_family)>2180, ancom_fec_family)
ancom_fec_family_5_min #129 taxa 436 samples

#reorder liver simple to make edible the comparator
ancom_fec_family_5_min@sam_data$Liver_Simple <- factor(ancom_fec_family_5_min@sam_data$Liver_Simple, levels = c("Edible", "A-","A", "A+"))
ancom_fec_family_5_min@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_fec_family_5_min <-ancombc2(data= ancom_fec_family_5_min, assay_name = "counts", tax_level = "Family",
                                     fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                     group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE,
                                     iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                     em_control = list(tol = 1e-5, max_iter = 100),
                                     lme_control = lme4::lmerControl(),
                                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                     trend_control = TRUE)


#do a tornado
volcano_fec_ancom_5_min <- ancombc_output_fec_family_5_min$res_dunn

#clean data
df_fec_long_volcano_5_min<- volcano_fec_ancom_5_min %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fec_long_volcano_5_min$group <- sub("_rounded", "", df_fec_long_volcano_5_min$group)
#Just double check that your groups turned out correct
df_fec_long_volcano_5_min$group
unique(df_fec_long_volcano_5_min$group)

#Filter the results to only include significant features
significant_features_5_min = df_fec_long_volcano_5_min$taxon

# Extract the names of the significant features
significant_feature_names_5_min = as.character(significant_features_5_min)

# Extract the bias_correct_log_table from the ancombc output
log_table_5_min = ancombc_output_fec_family_5_min$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant_5_min = log_table[rownames(log_table_5_min) %in% significant_feature_names_5_min, ]

# Convert the row names into a column
log_table_significant_5_min$Feature = rownames(log_table_significant_5_min)

# Reshape the data into a format suitable for ggplot2
df_fec_5_min = reshape2::melt(log_table_significant_5_min, id.vars = "Feature")

# Calculating the average value for each Feature
median_values_5_min <- df_fec_5_min %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
median_values_5_min

df_fec_5_min$value

# Merge the df data frame with the sample data
fec_merged_df_5_min = merge(median_values_5_min, df_fec_long_volcano_5_min, by.x = "Feature", by.y = "taxon")

colnames(fec_merged_df_5_min)

fec_merged_df_5_min$Feature<- factor(fec_merged_df_5_min$Feature, levels = unique(fec_merged_df_5_min$Feature))

# Update the tornado_plot with new specifications
fec_tornado_plot_5_min <- ggplot(fec_merged_df_5_min, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "Fecal") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-1, 1)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

fec_tornado_plot_5_min

#filter fecal to average of 10 counts per sample
##436 samples so if it occurs 10 times per sample would be a sum of 4,360 counts

#Look at whats in unfiltered object
ancom_fec_family #492 taxa 436 samples

ancom_fec_family_10_min <- prune_taxa(taxa_sums(ancom_fec_family)>4360, ancom_fec_family)
ancom_fec_family_10_min #117 taxa 436 samples

#reorder liver simple to make edible the comparator
ancom_fec_family_10_min@sam_data$Liver_Simple <- factor(ancom_fec_family_10_min@sam_data$Liver_Simple, levels = c("Edible", "A-","A", "A+"))
ancom_fec_family_10_min@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_fec_family_10_min <-ancombc2(data= ancom_fec_family_10_min, assay_name = "counts", tax_level = "Family",
                                           fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                           p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                           group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                                           global = TRUE, pairwise = TRUE,
                                           dunnet = TRUE, trend = FALSE,
                                           iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                           em_control = list(tol = 1e-5, max_iter = 100),
                                           lme_control = lme4::lmerControl(),
                                           mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                           trend_control = TRUE)


#do a tornado
volcano_fec_ancom_10_min <- ancombc_output_fec_family_10_min$res_dunn

#clean data
df_fec_long_volcano_10_min<- volcano_fec_ancom_10_min %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>%
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fec_long_volcano_10_min$group <- sub("_rounded", "", df_fec_long_volcano_10_min$group)
#Just double check that your groups turned out correct
df_fec_long_volcano_10_min$group
unique(df_fec_long_volcano_10_min$group)

#Filter the results to only include significant features
significant_features_10_min = df_fec_long_volcano_10_min$taxon

# Extract the names of the significant features
significant_feature_names_10_min = as.character(significant_features_10_min)

# Extract the bias_correct_log_table from the ancombc output
log_table_10_min = ancombc_output_fec_family_10_min$bias_correct_log_table

# Subset the table to include only the significant features
log_table_significant_10_min = log_table[rownames(log_table_10_min) %in% significant_feature_names_10_min, ]

# Convert the row names into a column
log_table_significant_10_min$Feature = rownames(log_table_significant_10_min)

# Reshape the data into a format suitable for ggplot2
df_fec_10_min = reshape2::melt(log_table_significant_10_min, id.vars = "Feature")

# Calculating the average value for each Feature
median_values_10_min <- df_fec_10_min %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
median_values_10_min

df_fec_10_min$value

# Merge the df data frame with the sample data
fec_merged_df_10_min = merge(median_values_10_min, df_fec_long_volcano_10_min, by.x = "Feature", by.y = "taxon")

colnames(fec_merged_df_10_min)

fec_merged_df_10_min$Feature<- factor(fec_merged_df_10_min$Feature, levels = unique(fec_merged_df_10_min$Feature))

# Update the tornado_plot with new specifications
fec_tornado_plot_10_min <- ggplot(fec_merged_df_10_min, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "A. Edible vs A-", "Liver_SimpleA" = "B. Edible vs A", "Liver_SimpleA+"= "C. Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "Fecal") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-1, 1)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

fec_tornado_plot_10_min


# Ancom Fecal Binomial ----------------------------------------------------
#use the filtered at 10 min data set
#check order of variable
ancom_fec_family_10_min@sam_data$Liver_Binom <- factor(ancom_fec_family_10_min@sam_data$Liver_Binom, levels = c("No", "Yes"))
ancom_fec_family_10_min@sam_data$Liver_Binom


##running ancombc on the variable of interest (liver binom) with score and trt as fixed no random
##dunnets turned on
ancombc_output_fec_family_binom <-ancombc2(data= ancom_fec_family_10_min, assay_name = "counts", tax_level = "Family",
                                            fix_formula = "Liver_Binom + Treatment", rand_formula = NULL, 
                                            p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                            group = "Liver_Binom", struc_zero = TRUE, neg_lb = TRUE,
                                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                                            global = TRUE, pairwise = TRUE,
                                            dunnet = TRUE, trend = FALSE,
                                            iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                            em_control = list(tol = 1e-5, max_iter = 100),
                                            lme_control = lme4::lmerControl(),
                                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                            trend_control = TRUE)

#do a tornado
volcano_fec_ancom_binom <- ancombc_output_fec_family_binom$res
View(volcano_fec_ancom_binom)

#need to come back to this lets wait and see if it matters
df_fec_long_volcano_binom<- volcano_fec_ancom_binom %>%
  mutate(across(starts_with("lfc_Liver_Binom"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Binom"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Binom"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% # maggies original code was not filtering non different taxa i.e. if the pass sensitivity but werent diffeernt they were still included so i added another filter level
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_fec_long_volcano_binom$group <- sub("_rounded", "", df_fec_long_volcano_binom$group)
#Just double check that your groups turned out correct
df_fec_long_volcano_binom$group


#Filter the results to only include significant features
fec_significant_features_binom = df_fec_long_volcano_binom$taxon

# Extract the names of the significant features
fec_significant_feature_names_binom = as.character(fec_significant_features_binom)

# Extract the bias_correct_log_table from the ancombc output
fec_log_table_binom = ancombc_output_fec_family_binom$bias_correct_log_table

# Subset the table to include only the significant features
fec_log_table_significant_binom = log_table[rownames(fec_log_table_binom) %in% fec_significant_feature_names_binom, ]

# Convert the row names into a column
fec_log_table_significant_binom$Feature = rownames(fec_log_table_significant_binom)

# Reshape the data into a format suitable for ggplot2
df_fec_binom = reshape2::melt(fec_log_table_significant_binom, id.vars = "Feature")

# Calculating the average value for each Feature
fec_median_values_binom <- df_fec_binom %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
fec_median_values_binom

df_fec_binom$value

# Merge the df data frame with the sample data
fec_merged_df_binom = merge(fec_median_values_binom, df_fec_long_volcano_binom, by.x = "Feature", by.y = "taxon")

colnames(fec_merged_df_binom)

fec_merged_df_binom$Feature <- factor(fec_merged_df_binom$Feature, levels = unique(fec_merged_df_binom$Feature))

# Update the tornado_plot with new specifications
fec_tornado_plot_binom <- ggplot(fec_merged_df_binom, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_BinomYes", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_BinomYes" = "Edible vs Abscess", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-3, 3)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code
fec_tornado_plot_binom





# ANCOM for Rumen ---------------------------------------------------------

#subset only rumen samples
ancom_ru <- subset_samples(samples, Sample_site== "RU")
#Subset family counts must have untransformed counts
ancom_ru_family <- tax_glom(ancom_ru, taxrank = "Family", NArm = F)
ancom_ru_family@sam_data$Liver_Simple
#reorder liver simple to make edible the comparator, but only have A+ and Edible so dont need the other scores
ancom_ru_family@sam_data$Liver_Simple <- factor(ancom_ru_family@sam_data$Liver_Simple, levels = c("Edible", "A+"))
ancom_ru_family@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_ru_family <-ancombc2(data= ancom_ru_family, assay_name = "counts", tax_level = "Family",
                                     fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                     group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE,
                                     iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                     em_control = list(tol = 1e-5, max_iter = 100),
                                     lme_control = lme4::lmerControl(),
                                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                     trend_control = TRUE)

##CLean up the data

#subset data from dunns test
ru_ancom <- ancombc_output_ru_family$res

View(ru_ancom)

##Data prep for heatmap
df_ru_liver_simple_dunn <- ru_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+`==1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`diff_Liver_SimpleA+`== 1, round(`lfc_Liver_SimpleA+`, 2), `lfc_Liver_SimpleA+`)) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= round(`lfc_Liver_SimpleA+`)) %>% 
  tidyr::pivot_longer(cols = "Liver_SimpleA+", names_to = "group", values_to = "value") %>% 
  dplyr::arrange(taxon)

#sensititvity testing
df_ses_ru_Liver_simple <- ru_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+` == 1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`passed_ss_Liver_SimpleA+`==1 & `diff_Liver_SimpleA+`==1, "white", "black"),) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= `lfc_Liver_SimpleA+`) %>% 
  tidyr::pivot_longer(cols ="Liver_SimpleA+", names_to = "group", values_to = "color") %>% 
  dplyr::arrange(taxon)

#join sensitivity with dunn
df_ru_liver_simple_dunn1 <- df_ru_liver_simple_dunn %>% 
  dplyr::left_join(df_ses_ru_Liver_simple, by = c("taxon", "group"))

#change order of X-value 
df_ru_liver_simple_dunn1$group <- factor(df_ru_liver_simple_dunn1, levels = c("Liver_SimpleA+"))


#Assign your map colors to values
lo = floor(min(df_ru_liver_simple_dunn1$value))
up = ceiling(max(df_ru_liver_simple_dunn1$value))
mid = 0

ru_liver_simple_heatmap = df_ru_liver_simple_dunn1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "black", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold change compared to Edible") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ru_liver_simple_heatmap


#do a tornado
volcano_ru_ancom <- ancombc_output_ru_family$res
View(volcano_ru_ancom)
#need to come back to this lets wait and see if it matters
df_ru_long_volcano<- volcano_ru_ancom %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% # maggies original code was not filtering non different taxa i.e. if the pass sensitivity but werent diffeernt they were still included so i added another filter level
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_ru_long_volcano$group <- sub("_rounded", "", df_ru_long_volcano$group)
#Just double check that your groups turned out correct
df_ru_long_volcano$group


#Filter the results to only include significant features
ru_significant_features = df_ru_long_volcano$taxon

# Extract the names of the significant features
ru_significant_feature_names = as.character(ru_significant_features)

# Extract the bias_correct_log_table from the ancombc output
ru_log_table = ancombc_output_ru_family$bias_correct_log_table

# Subset the table to include only the significant features
ru_log_table_significant = log_table[rownames(log_table) %in% ru_significant_feature_names, ]

# Convert the row names into a column
ru_log_table_significant$Feature = rownames(ru_log_table_significant)

# Reshape the data into a format suitable for ggplot2
df_ru = reshape2::melt(ru_log_table_significant, id.vars = "Feature")

# Calculating the average value for each Feature
ru_median_values <- df_ru %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
ru_median_values

df_ru$value

# Merge the df data frame with the sample data
ru_merged_df = merge(ru_median_values, df_ru_long_volcano, by.x = "Feature", by.y = "taxon")

colnames(ru_merged_df)

ru_merged_df$Feature <- factor(ru_merged_df$Feature, levels = unique(ru_merged_df$Feature))

# Update the tornado_plot with new specifications
ru_tornado_plot <- ggplot(ru_merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "A.Rumen") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-1, 1)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

  
ru_tornado_plot 

# ANCOM RU Filtered -------------------------------------------------------

#look at total for unfiltered
ancom_ru_family #492 taxa 84 samples

#filter
#84 samples so if 10 occurrences per sample then sum = 840
ancom_ru_family_10_min <- prune_taxa(taxa_sums(ancom_ru_family)>840, ancom_ru_family)
ancom_ru_family_10_min # 123 taxa

#check level
ancom_ru_family_10_min@sam_data$Liver_Simple
#level is already ordered


##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_ru_family_10_min <-ancombc2(data= ancom_ru_family_10_min, assay_name = "counts", tax_level = "Family",
                                    fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                    group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                    alpha = 0.05, n_cl = 2, verbose = TRUE,
                                    global = TRUE, pairwise = TRUE,
                                    dunnet = TRUE, trend = FALSE,
                                    iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                    em_control = list(tol = 1e-5, max_iter = 100),
                                    lme_control = lme4::lmerControl(),
                                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                    trend_control = TRUE)


#do a tornado
volcano_ru_ancom_10_min <- ancombc_output_ru_family_10_min$res

#need to come back to this lets wait and see if it matters
df_ru_long_volcano_10_min<- volcano_ru_ancom_10_min %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% # maggies original code was not filtering non different taxa i.e. if the pass sensitivity but werent diffeernt they were still included so i added another filter level
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_ru_long_volcano_10_min$group <- sub("_rounded", "", df_ru_long_volcano_10_min$group)
#Just double check that your groups turned out correct
df_ru_long_volcano_10_min$group


#Filter the results to only include significant features
ru_significant_features_10_min = df_ru_long_volcano_10_min$taxon

# Extract the names of the significant features
ru_significant_feature_names_10_min = as.character(ru_significant_features_10_min)

# Extract the bias_correct_log_table from the ancombc output
ru_log_table_10_min = ancombc_output_ru_family_10_min$bias_correct_log_table

# Subset the table to include only the significant features
ru_log_table_significant_10_min = log_table[rownames(ru_log_table_10_min) %in% ru_significant_feature_names_10_min, ]

# Convert the row names into a column
ru_log_table_significant_10_min$Feature = rownames(ru_log_table_significant_10_min)

# Reshape the data into a format suitable for ggplot2
df_ru_10_min = reshape2::melt(ru_log_table_significant_10_min, id.vars = "Feature")

# Calculating the average value for each Feature
ru_median_values_10_min <- df_ru_10_min %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
ru_median_values_10_min

df_ru_10_min$value

# Merge the df data frame with the sample data
ru_merged_df_10_min = merge(ru_median_values_10_min, df_ru_long_volcano_10_min, by.x = "Feature", by.y = "taxon")

colnames(ru_merged_df_10_min)

ru_merged_df_10_min$Feature <- factor(ru_merged_df_10_min$Feature, levels = unique(ru_merged_df_10_min$Feature))

# Update the tornado_plot with new specifications
ru_tornado_plot_10_min <- ggplot(ru_merged_df_10_min, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "A.Rumen") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-1, 1)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code
ru_tornado_plot_10_min



# ANCOM for SI ------------------------------------------------------------

#subset only rumen samples
ancom_si <- subset_samples(samples, Sample_site== "SI")
#Subset family counts must have untransformed counts
ancom_si_family <- tax_glom(ancom_si, taxrank = "Family", NArm = F)
ancom_si_family@sam_data$Liver_Simple
#reorder liver simple to make edible the comparator, but only have A+ and Edible so dont need the other scores
ancom_si_family@sam_data$Liver_Simple <- factor(ancom_si_family@sam_data$Liver_Simple, levels = c("Edible", "A+"))
ancom_si_family@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_si_family <-ancombc2(data= ancom_si_family, assay_name = "counts", tax_level = "Family",
                                    fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                    group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                    alpha = 0.05, n_cl = 2, verbose = TRUE,
                                    global = TRUE, pairwise = TRUE,
                                    dunnet = TRUE, trend = FALSE,
                                    iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                    em_control = list(tol = 1e-5, max_iter = 100),
                                    lme_control = lme4::lmerControl(),
                                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                    trend_control = TRUE)
##CLean up the data

#subset data from dunns test
si_ancom <- ancombc_output_si_family$res

##Data prep for heatmap
df_si_liver_simple_dunn <- si_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+`==1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`diff_Liver_SimpleA+`== 1, round(`lfc_Liver_SimpleA+`, 2), `lfc_Liver_SimpleA+`)) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= round(`lfc_Liver_SimpleA+`)) %>% 
  tidyr::pivot_longer(cols = "Liver_SimpleA+", names_to = "group", values_to = "value") %>% 
  dplyr::arrange(taxon)

#sensititvity testing
df_ses_si_Liver_simple <- si_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+` == 1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`passed_ss_Liver_SimpleA+`==1 & `diff_Liver_SimpleA+`==1, "white", "black"),) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= `lfc_Liver_SimpleA+`) %>% 
  tidyr::pivot_longer(cols ="Liver_SimpleA+", names_to = "group", values_to = "color") %>% 
  dplyr::arrange(taxon)

#join sensitivity with dunn
df_si_liver_simple_dunn1 <- df_si_liver_simple_dunn %>% 
  dplyr::left_join(df_ses_si_Liver_simple, by = c("taxon", "group"))

#change order of X-value 
df_si_liver_simple_dunn1$group <- factor(df_si_liver_simple_dunn1, levels = c("Liver_SimpleA+"))


#Assign your map colors to values
lo = floor(min(df_si_liver_simple_dunn1$value))
up = ceiling(max(df_si_liver_simple_dunn1$value))
mid = 0

si_liver_simple_heatmap = df_si_liver_simple_dunn1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "black", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold change compared to Edible") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
si_liver_simple_heatmap

#do a tornado
volcano_si_ancom <- ancombc_output_si_family$res

#need to come back to this lets wait and see if it matters
df_si_long_volcano<- volcano_si_ancom %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% 
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_si_long_volcano$group <- sub("_rounded", "", df_si_long_volcano$group)
#Just double check that your groups turned out correct
df_si_long_volcano$group


#Filter the results to only include significant features
si_significant_features = df_si_long_volcano$taxon

# Extract the names of the significant features
si_significant_feature_names = as.character(si_significant_features)

# Extract the bias_correct_log_table from the ancombc output
si_log_table = ancombc_output_si_family$bias_correct_log_table

# Subset the table to include only the significant features
si_log_table_significant = log_table[rownames(log_table) %in% si_significant_feature_names, ]

# Convert the row names into a column
si_log_table_significant$Feature = rownames(si_log_table_significant)

# Reshape the data into a format suitable for ggplot2
df_si = reshape2::melt(si_log_table_significant, id.vars = "Feature")

# Calculating the average value for each Feature
si_median_values <- df_si %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
si_median_values

df_si$value

# Merge the df data frame with the sample data
si_merged_df = merge(si_median_values, df_si_long_volcano, by.x = "Feature", by.y = "taxon")

colnames(si_merged_df)

si_merged_df$Feature <- factor(si_merged_df$Feature, levels = unique(si_merged_df$Feature))

# Update the tornado_plot with new specifications
si_tornado_plot <- ggplot(si_merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "", title = "B.Small Intestine") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-2, 2)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

si_tornado_plot 


# ANCOM for SI filtered ---------------------------------------------------
#look at the unfiltered data
ancom_si_family #492 taxa 73 samples
  
#if 73 samples and want 10 per sum = 730
ancom_si_family_10_min <- prune_taxa(taxa_sums(ancom_si_family)>730, ancom_si_family)
ancom_si_family_10_min #127 taxa

#check order of the variable
ancom_si_family@sam_data$Liver_Simple
  
##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_si_family_10_min <-ancombc2(data= ancom_si_family_10_min, assay_name = "counts", tax_level = "Family",
                                      fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                      p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                      group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                      alpha = 0.05, n_cl = 2, verbose = TRUE,
                                      global = TRUE, pairwise = TRUE,
                                      dunnet = TRUE, trend = FALSE,
                                      iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                      em_control = list(tol = 1e-5, max_iter = 100),
                                      lme_control = lme4::lmerControl(),
                                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                      trend_control = TRUE)
##CLean up the data
#do a tornado
volcano_si_ancom_10_min <- ancombc_output_si_family_10_min$res
  
#need to come back to this lets wait and see if it matters
df_si_long_volcano_10_min<- volcano_si_ancom_10_min %>%
    mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
           across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
           across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                  ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                          ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                  .names = "{.col}_color")
    ) %>% 
    pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
    pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
    filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
    select(-color_group) %>%
    arrange(taxon)%>%
    mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")
  
# Ensure 'group' column matches between LFC values and colors for direct comparison
df_si_long_volcano_10_min$group <- sub("_rounded", "", df_si_long_volcano_10_min$group)
#Just double check that your groups turned out correct
df_si_long_volcano_10_min$group
  
  
#Filter the results to only include significant features
si_significant_features_10_min = df_si_long_volcano_10_min$taxon
  
# Extract the names of the significant features
si_significant_feature_names_10_min = as.character(si_significant_features_10_min)
  
# Extract the bias_correct_log_table from the ancombc output
si_log_table_10_min = ancombc_output_si_family_10_min$bias_correct_log_table
  
# Subset the table to include only the significant features
si_log_table_significant_10_min = log_table[rownames(si_log_table_10_min) %in% si_significant_feature_names_10_min, ]
  
# Convert the row names into a column
si_log_table_significant_10_min$Feature = rownames(si_log_table_significant_10_min)
  
# Reshape the data into a format suitable for ggplot2
df_si_10_min = reshape2::melt(si_log_table_significant_10_min, id.vars = "Feature")
  
# Calculating the average value for each Feature
si_median_values_10_min <- df_si_10_min %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
si_median_values_10_min
  
df_si_10_min$value
  
# Merge the df data frame with the sample data
si_merged_df_10_min = merge(si_median_values_10_min, df_si_long_volcano_10_min, by.x = "Feature", by.y = "taxon")
  
colnames(si_merged_df_10_min)
  
si_merged_df_10_min$Feature <- factor(si_merged_df_10_min$Feature, levels = unique(si_merged_df_10_min$Feature))
  
# Update the tornado_plot with new specifications
si_tornado_plot_10_min <- ggplot(si_merged_df_10_min, aes(x = value, y = Feature)) +
    geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
    scale_color_identity() +
    facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                        labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
               scales = "free_y") +
    scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
    labs(x = "LogFC Value", y = "", title = "B.Small Intestine") +
    scale_y_discrete(labels = function(x) gsub("_"," ", x))+
    theme_minimal() +
    theme(
      panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
      panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
      axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
      axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
      axis.title.y = element_text(size = 26),  # Removes y-axis label
      axis.title.x = element_text(size = 26),
      plot.title = element_text(size = 30, hjust = 0.5),
      panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
      strip.background = element_blank(),  # Optional: removes background from panel labels
      strip.text.x = element_text(face = "bold"))+ # Bold text for group names
    scale_x_continuous(limits = c(-2, 2)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code
  
si_tornado_plot_10_min 
  

# ANCOM LI ----------------------------------------------------------------

#subset only rumen samples
ancom_li <- subset_samples(samples, Sample_site== "LI")
#Subset family counts must have untransformed counts
ancom_li_family <- tax_glom(ancom_li, taxrank = "Family", NArm = F)
ancom_li_family@sam_data$Liver_Simple
#reorder liver simple to make edible the comparator, but only have A+ and Edible so dont need the other scores
ancom_li_family@sam_data$Liver_Simple <- factor(ancom_li_family@sam_data$Liver_Simple, levels = c("Edible", "A+"))
ancom_li_family@sam_data$Liver_Simple

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_li_family <-ancombc2(data= ancom_li_family, assay_name = "counts", tax_level = "Family",
                                    fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                    group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                    alpha = 0.05, n_cl = 2, verbose = TRUE,
                                    global = TRUE, pairwise = TRUE,
                                    dunnet = TRUE, trend = FALSE,
                                    iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                    em_control = list(tol = 1e-5, max_iter = 100),
                                    lme_control = lme4::lmerControl(),
                                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                    trend_control = TRUE)
##CLean up the data

#subset data from dunns test
li_ancom <- ancombc_output_li_family$res

##Data prep for heatmap
df_li_liver_simple_dunn <- li_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+`==1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`diff_Liver_SimpleA+`== 1, round(`lfc_Liver_SimpleA+`, 2), `lfc_Liver_SimpleA+`)) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= round(`lfc_Liver_SimpleA+`)) %>% 
  tidyr::pivot_longer(cols = "Liver_SimpleA+", names_to = "group", values_to = "value") %>% 
  dplyr::arrange(taxon)

#sensititvity testing
df_ses_li_Liver_simple <- li_ancom %>% 
  dplyr::filter(`diff_Liver_SimpleA+` == 1) %>% 
  dplyr::mutate(`lfc_Liver_SimpleA+`= ifelse(`passed_ss_Liver_SimpleA+`==1 & `diff_Liver_SimpleA+`==1, "white", "black"),) %>% 
  dplyr::transmute(taxon, "Liver_SimpleA+"= `lfc_Liver_SimpleA+`) %>% 
  tidyr::pivot_longer(cols ="Liver_SimpleA+", names_to = "group", values_to = "color") %>% 
  dplyr::arrange(taxon)

#join sensitivity with dunn
df_li_liver_simple_dunn1 <- df_li_liver_simple_dunn %>% 
  dplyr::left_join(df_ses_li_Liver_simple, by = c("taxon", "group"))


#Assign your map colors to values
lo = floor(min(df_li_liver_simple_dunn1$value))
up = ceiling(max(df_li_liver_simple_dunn1$value))
mid = 0

li_liver_simple_heatmap = df_li_liver_simple_dunn1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "black", high = "red", mid = "blue", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold change compared to Edible") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
li_liver_simple_heatmap

#do a tornado
volcano_li_ancom <- ancombc_output_li_family$res

#need to come back to this lets wait and see if it matters
df_li_long_volcano<- volcano_li_ancom %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% 
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_li_long_volcano$group <- sub("_rounded", "", df_li_long_volcano$group)
#Just double check that your groups turned out correct
df_li_long_volcano$group


#Filter the results to only include significant features
li_significant_features = df_li_long_volcano$taxon

# Extract the names of the significant features
li_significant_feature_names = as.character(li_significant_features)

# Extract the bias_correct_log_table from the ancombc output
li_log_table = ancombc_output_li_family$bias_correct_log_table

# Subset the table to include only the significant features
li_log_table_significant = log_table[rownames(log_table) %in% li_significant_feature_names, ]

# Convert the row names into a column
li_log_table_significant$Feature = rownames(li_log_table_significant)

# Reshape the data into a format suitable for ggplot2
df_li = reshape2::melt(li_log_table_significant, id.vars = "Feature")

# Calculating the average value for each Feature
li_median_values <- df_li %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
li_median_values

df_li$value

# Merge the df data frame with the sample data
li_merged_df = merge(li_median_values, df_li_long_volcano, by.x = "Feature", by.y = "taxon")

colnames(li_merged_df)

li_merged_df$Feature <- factor(li_merged_df$Feature, levels = unique(li_merged_df$Feature))

# Update the tornado_plot with new specifications
li_tornado_plot <- ggplot(li_merged_df, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "Feature", title = "C.Large Intestine") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-2, 2)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

li_tornado_plot 


# ANCOM for LI filtered ---------------------------------------------------
#look at unfiltered
ancom_li_family #492 taxa 79 samples

#if want 10 per sample then sum = 790
ancom_li_family_10_min <- prune_taxa(taxa_sums(ancom_li_family)>790, ancom_li_family)
ancom_li_family_10_min #119 taxa

#check order of variable
ancom_li_family@sam_data$Liver_Simple 

##running ancombc on the variable of interest (liver simple) with score and trt as fixed no random
##dunnets turned on
ancombc_output_li_family_10_min <-ancombc2(data= ancom_li_family_10_min, assay_name = "counts", tax_level = "Family",
                                    fix_formula = "Liver_Simple + Treatment", rand_formula = NULL, 
                                    p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                    group = "Liver_Simple", struc_zero = TRUE, neg_lb = TRUE,
                                    alpha = 0.05, n_cl = 2, verbose = TRUE,
                                    global = TRUE, pairwise = TRUE,
                                    dunnet = TRUE, trend = FALSE,
                                    iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE), #You can increase the iterations here to 100 if needed
                                    em_control = list(tol = 1e-5, max_iter = 100),
                                    lme_control = lme4::lmerControl(),
                                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                    trend_control = TRUE)

#do a tornado
volcano_li_ancom_10_min <- ancombc_output_li_family_10_min$res

#need to come back to this lets wait and see if it matters
df_li_long_volcano_10_min<- volcano_li_ancom_10_min %>%
  mutate(across(starts_with("lfc_Liver_Simple"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x, .x)), #FOR THE LOVE OF GOD PUT THIS AS .X OR IT JUST TAKE ALL THINGS TO 0, it won't show log fold changes 
         across(starts_with("lfc_Liver_Simple"), ~round(.x, 2), .names = "{.col}_rounded"),
         across(starts_with("diff_Liver_Simple"), #This whole string makes sure that it applies NA to the values that don't pass sensitivity testing, thank you Valeria
                ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "red3",
                        ifelse(.x == 1 | get(gsub("diff_", "passed_ss_", cur_column())) == 0, NA, "black")), 
                .names = "{.col}_color")
  ) %>% 
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_",values_drop_na = FALSE) %>%
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_",values_drop_na = FALSE) %>%
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  select(-color_group) %>%
  arrange(taxon)%>%
  mutate(value = ifelse(is.na(color), NA, value)) #here is where it actually makes sure the value AND color say "NA")

# Ensure 'group' column matches between LFC values and colors for direct comparison
df_li_long_volcano_10_min$group <- sub("_rounded", "", df_li_long_volcano_10_min$group)
#Just double check that your groups turned out correct
df_li_long_volcano_10_min$group


#Filter the results to only include significant features
li_significant_features_10_min = df_li_long_volcano_10_min$taxon

# Extract the names of the significant features
li_significant_feature_names_10_min = as.character(li_significant_features_10_min)

# Extract the bias_correct_log_table from the ancombc output
li_log_table_10_min = ancombc_output_li_family_10_min$bias_correct_log_table

# Subset the table to include only the significant features
li_log_table_significant_10_min = log_table[rownames(li_log_table_10_min) %in% li_significant_feature_names_10_min, ]

# Convert the row names into a column
li_log_table_significant_10_min$Feature = rownames(li_log_table_significant_10_min)

# Reshape the data into a format suitable for ggplot2
df_li_10_min = reshape2::melt(li_log_table_significant_10_min, id.vars = "Feature")

# Calculating the average value for each Feature
li_median_values_10_min <- df_li_10_min %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) # na.rm = TRUE to handle any NA values, REMOVE NA OR ELSE IT WILL ALL BE NA!
li_median_values_10_min

df_li_10_min$value

# Merge the df data frame with the sample data
li_merged_df_10_min = merge(li_median_values_10_min, df_li_long_volcano_10_min, by.x = "Feature", by.y = "taxon")

colnames(li_merged_df_10_min)

li_merged_df_10_min$Feature <- factor(li_merged_df_10_min$Feature, levels = unique(li_merged_df_10_min$Feature))

# Update the tornado_plot with new specifications
li_tornado_plot_10_min <- ggplot(li_merged_df_10_min, aes(x = value, y = Feature)) +
  geom_point(aes(size =log2(abs(median_value)+15), color = color)) +# Shape 21 is a filled circle
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Add a dashed line at x=0
  scale_color_identity() +
  facet_grid(~ factor(group,levels = c("Liver_SimpleA-", "Liver_SimpleA", "Liver_SimpleA+"), 
                      labels = c("Liver_SimpleA-" = "Edible vs A-", "Liver_SimpleA" = "Edible vs A", "Liver_SimpleA+"= "Edible vs A+")), 
             scales = "free_y") +
  scale_size_continuous(range = c(2, 4), guide = "none") +  # Adjust the size range as needed
  labs(x = "LogFC Value", y = "", title = "C.Large Intestine") +
  scale_y_discrete(labels = function(x) gsub("_"," ", x))+
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),  # Remove minor vertical grid lines
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),  # Keep horizontal grid lines
    axis.text.y = element_text(color = "black",lineheight = 1.5, size = 12),  # Uniform color for y-axis labels
    axis.ticks.y = element_line(color = "black"),  # Uniform color for y-axis ticks
    axis.title.y = element_text(size = 26),  # Removes y-axis label
    axis.title.x = element_text(size = 26),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(1, "lines"),  # Adjust spacing between panels
    strip.background = element_blank(),  # Optional: removes background from panel labels
    strip.text.x = element_text(face = "bold"))+ # Bold text for group names
  scale_x_continuous(limits = c(-2, 2)) #I manually adjusted my x-axis values so the dots didn't get cut off, if you don't care, use following line of code

li_tornado_plot_10_min 



# merge plots for ancom ---------------------------------------------------

ancom_git <- plot_grid(ru_tornado_plot, si_tornado_plot, li_tornado_plot,align = "hv", axis = "btrl", nrow = 1 )
ancom_git

ancom_git_10 <- plot_grid(ru_tornado_plot_10_min, si_tornado_plot_10_min, li_tornado_plot_10_min, align =  "hv", axis = "btrl", nrow = 1)
ancom_git_10

# Subset taxa just to see -------------------------------------------------
#subset Prevotellaceae
prev_sub <- subset_taxa(ra_family, Family == "Prevotellaceae") 
prev_melt <- psmelt(prev_sub)
prev_sub_si <- subset(prev_melt, Sample_site=="SI") 

#order variable
prev_sub_si$Liver_Simple <- factor(prev_sub_si$Liver_Simple, levels = c("Edible","A+"))

prev_boxplot <- ggplot(prev_sub_si, aes(x= Liver_Simple, y= Abundance, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "Prevotellaceae") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_color_manual(values = git_simple_palette, name = "Liver Score")+
  scale_fill_manual(values = git_simple_palette, name = "Liver Score")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30, face = "italic"),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
prev_boxplot

diff_legend <- get_legend(prev_boxplot)

#Synergistaceae
syn_sub <- subset_taxa(ra_family, Family=="Synergistaceae")
syn_melt <- psmelt(syn_sub)
syn_sub_si <- subset(syn_melt, Sample_site=="SI")

#set order
syn_sub_si$Liver_Simple <- factor(syn_sub_si$Liver_Simple, levels = c("Edible", "A+"))

syn_boxplot <- ggplot(syn_sub_si, aes(x= Liver_Simple, y= Abundance, fill = Liver_Simple)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "Synergistaceae") +
  scale_x_discrete(limits=c("Edible", "A+"))+
  geom_boxplot(aes(color= Liver_Simple), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color= Liver_Simple), size = 2, position = position_jitter(seed = 1, w=.2))+
  scale_fill_manual(values = git_simple_palette, name= "Liver Score") +
  scale_color_manual(values = git_simple_palette, name = "Liver Score")+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30, face = "italic"),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
syn_boxplot

syn_legend <- get_legend(syn_boxplot)

plot_grid(plot_grid(prev_boxplot,syn_boxplot,align = "hv", axis = "tbrl", ncol = 1), syn_legend, rel_widths = c(1,.5))


