library(phyloseq); packageVersion("phyloseq")
library(ggplot2) #Package to generate plots
library(vegan) #community ecology package
library(scales) #To modify scale of plot's axis
library(metagMisc) #Package to manipulate phyloseq objects
library(dplyr) #package used for data manipulation 
library(reshape2) #package used for data manipulation 
library(data.table) #package used for data manipulation 

#Upload metadata and format it
#We are changing the factors order just to make plots more coherent
#with biological gradients Bulk>Rhizo>Endo, and Short>Medium>Tall
samples <- read.delim2("../metadata_per_sample.tsv", header = TRUE)
samples$Compartment <- as.factor(samples$Compartment)
samples$Spartina <- as.factor(samples$Spartina)
samples$Spartina = factor(samples$Spartina,levels(samples$Spartina)[c(2,1,3)])
samples$Compartment = factor(samples$Compartment,levels(samples$Compartment)[c(1,3,2)])
rownames(samples) <- samples$Seq_ID
samples$Transect <- as.factor(samples$Transect)
samples$Point <- as.factor(samples$Point)

#Create a phyloseq object phy using the ASVs table, metadata, and taxonomy table:
phy <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samples), 
                tax_table(taxa))

phy
#Explore number of sequence reads per sample 
summary(sample_sums(phy))

###Remove Chloroplast, Mitochondria, Eukaryots, and Prokaryotes with unknown Phylum
#First, explore identified taxa at the Kingdom level
table(tax_table(phy)[, "Kingdom"], exclude = NULL)

#remove chloroplast 16S rRNA gene
phy.raw = subset_taxa(phy, Order != "Chloroplast" | is.na(Order))
#remove mitochondrial 16S rRNA gene
phy.raw = subset_taxa(phy.raw, Family != "Mitochondria" | is.na(Family))
phy.raw
table(tax_table(phy.raw)[, "Kingdom"], exclude = NULL)

#Remove ASVs from unknown Phylum, and Eukaryotic ASVs (if any)
phy.filt = subset_taxa(phy.raw, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized") & !Kingdom %in% c("", "Eukaryota"))
phy.filt
table(tax_table(phy.filt)[, "Kingdom"], exclude = NULL)

#Explore number of sequence reads per sample after taxonomy filtration
summary(sample_sums(phy.filt))

############ALPHA DIV

#Calculate the number of counts of each sample
counts_all <- as.data.frame(sample_sums(phy.filt))
names(counts_all) <- "Counts"

#Retrieve metadata
samples_all <- data.frame(sample_data(phy.filt))
samples_all$Seq_ID <- as.character(samples_all$Seq_ID)

counts_all$Seq_ID <- samples_all$Seq_ID

##Estimate alpha diversity indexes per sample
alpha_div <- estimate_richness(phy.filt, split = TRUE,
                               measures =c("Observed", "Shannon", "Simpson")) 

#In order to merge the two databases, you need to have a variable with the same name in both dfs.
#In this case, it will be the Seq_ID
head(samples_all$Seq_ID)
alpha_div$Seq_ID <- samples_all$Seq_ID

#Merge metadata with alpha diversity indexes
alpha_div <- merge(samples_all, alpha_div)
alpha_div <- merge(alpha_div, counts_all)

#Plot Shannon

#Parameters for ggplot2. This is just to make the plot look prettier
theme_alpha <- theme_bw() + theme(axis.title.x = element_blank(),
                                  axis.title.y = element_text(size = 18),
                                  axis.text = element_text(size = 16),
                                  axis.text.x = element_text(angle = 15),
                                  panel.grid = element_blank(),
                                  legend.position = "top", legend.title = element_blank(),
                                  legend.text = element_text(size = 16))


fig_richness <- ggplot(alpha_div, aes (x = Compartment, y = Observed, fill = Spartina)) +
  geom_boxplot() + theme_alpha + labs(y = "Richness (ASVs)", x = "") + expand_limits(y = 0)

fig_Shannon <- ggplot(alpha_div, aes (x = Compartment, y = Shannon, fill = Spartina)) +
  geom_boxplot() + theme_alpha + labs(y = "Shannon index", x = "") 

fig_Simpson <- ggplot(alpha_div, aes (x = Compartment, y = Simpson, fill = Spartina)) +
  geom_boxplot() + theme_alpha + labs(y = "Simpson index", x = "") 

fig_richness
fig_Shannon
fig_Simpson

# Q9: What trend do you see in alpha diversity by plant compartment and marsh microhabitat?
# Did you expect these results? How do you explain them?


##Now, we will construct a rank-abundance plot.
#We will group samples by plant compartment, to make the figure readable 

#First, we will transform the counts of each ASV per sample into relative abundance values (0%-100%)
phy_ra <- transform_sample_counts(phy.filt, function(x) x / sum(x) * 100)

#Will generate a list with relative abundance, tax, and metadata as data.frames
otu_table(phy_ra) = t(otu_table(phy_ra))
data <- list(abund = as.data.frame(otu_table(phy_ra)@.Data),
             tax = data.frame(tax_table(phy_ra)@.Data, ASV = rownames(tax_table(phy_ra))),
             sample = suppressWarnings(as.data.frame(as.matrix(sample_data(phy_ra)))))

##The next steps are just data manipulation to calculate the rank abundance plot:
abund <- data[["abund"]]
tax <- data[["tax"]]
tax$ASV <- paste(rep("ASV", length(tax$ASV)), 1:length(tax$ASV), sep = "_")
sample <- data[["sample"]]

#Generate a long dataframe with ASVs, Sample name, and relative abundance as vectors 
abund3 <- cbind.data.frame(Display = tax[, "ASV"], abund) %>% 
  reshape2::melt(id.var = "Display", value.name = "Abundance", variable.name = "Sample")
head(abund3)

#Sort abund3 dataframe
abund3 <- data.table(abund3)[, `:=`(Abundance, sum(Abundance)),
                             by = list(Display, Sample)] %>% setkey(Display, Sample) %>%
  unique() %>% as.data.frame()
head(abund3)

#Generate a dataframe to group samples (and find mean values for the rank-abundance plot)
grp <- data.frame(Sample = rownames(sample), Group = sample[, "Compartment"])

abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
abund5 <- abund3
head(abund5)

#Sort abundance dataframe by decreasing abundance, per sample
abund5 <- abund5[order(abund5$Abundance, decreasing = T),]
abund5_ord <- abund5[order(abund5$Sample),]
head(abund5_ord)

#Create a new vector in the dataframe, with the Rank number
abund5_ord$Rank <- rep(1:length(unique(abund5_ord$Display)), length(unique(abund5_ord$Sample)))
head(abund5_ord)
tail(abund5_ord)

#Calculate the mean relative abundance per group and per rank value
temp3 <- dplyr::group_by(abund5_ord, Group, Rank) %>% dplyr::summarise(Mean = mean(Abundance))
head(temp3)
tail(temp3)

#Remove rank values that have Mean = 0 
temp3 <- temp3[temp3$Mean > 0, ]
head(temp3)
tail(temp3)

#Calculate the cumulative relative abundance:
TotalCounts <- temp3[with(temp3, order(-Mean)), ] %>%
  group_by(Group) %>% mutate(dummy = 1) %>% mutate(Cumsum = cumsum(Mean),
                                                   Rank = cumsum(dummy)) %>% as.data.frame()

head(TotalCounts) #Note that Mean and cumsum have different values for Rank != 1
tail(TotalCounts)

TotalCounts$Group <- as.factor(TotalCounts$Group)
TotalCounts$Group <- factor(TotalCounts$Group, levels(TotalCounts$Group)[c(1,3,2)])

#Plot the rank-abundance plot!

rank_abund_plot <- ggplot(data = TotalCounts, aes(x = Rank, y = Cumsum, color = Group)) + 
  geom_point(size = 2) + geom_line(size = 1) + ylim(0, 100) + scale_x_log10() +
  theme(legend.position = c(0.88,0.12), axis.title = element_text(size = 18),
        axis.text = element_text(size = 16), legend.title = element_blank(),
        legend.text = element_text(size = 16), panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "ASV Rank", y = "Cumulative rel. abundance (%)")

rank_abund_plot
#Note the log scale in the X axis
# Q10 How do you interpret the Rank-Abundance plot? Is it coherent with your 
# other alpha diversity metrics? Please explain

#Make a final plot summarizing all information
#ggarrange from the ggpubr package, allows you to make composite figures with ggplot2 objects
#It is a nice package for making figures for your own manuscripts/work/presentations
alpha_fig <- ggpubr::ggarrange(ggpubr::ggarrange(fig_richness, fig_Shannon, fig_Simpson, ncol = 3,
                                                 common.legend = TRUE, legend = "top"),
                               rank_abund_plot, ncol = 1)

#Save a jpeg file in your Workshop folder. You can include the figure in your report
jpeg(filename = "../fig_alpha.jpeg", res = 300, width = 3250, height = 2500)
alpha_fig
dev.off()

####We finished assessing alpha diversity. Now, we will proceed to analyze our communities'
#beta diversity.

#First,we will perform a multivariate ordination.
#Specifically, we will use Non-metric multidimensional scaling
#with Bray-Curtis dissimilarity index

#As before, this is only a piece of code to make our ggplots pretty:
theme_ord <- theme(axis.title.x = element_text(size = 18),
                   axis.title.y = element_text(size = 18),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 16),
                   plot.title = element_blank(),
                   legend.position = "bottom",
                   legend.title = element_blank(),
                   legend.text = element_text(size = 16)) 

#We will use 3 axis in our nmds ordination (parameter k)
ord.nmds.bray <- ordinate(phy.filt, method="NMDS", distance="bray", k = 3)
ord.nmds.bray
# Q11 What is the stress value? Is it good enough? How does it change if you modify the k value?
#To have reproducible figures, leave k = 3

#Now we will plot the NMDS ordination. We will make different plots.
#The only think that will change is the coloring pattern based on different factors

fig_nmds_bray_Spartina <- plot_ordination(phy.filt, ord.nmds.bray, color="Spartina",
                                          shape = "Spartina") + 
  geom_point(size =4) + coord_fixed(expand= TRUE) + theme_bw() + theme_ord
fig_nmds_bray_Spartina

fig_nmds_bray_Compartment <- plot_ordination(phy.filt, ord.nmds.bray, color="Compartment",
                                             shape = "Compartment") + 
  geom_point(size =4) + coord_fixed(expand= TRUE) + theme_bw() + theme_ord
fig_nmds_bray_Compartment

fig_nmds_bray_Location <- plot_ordination(phy.filt, ord.nmds.bray, color="Location",
                                          shape = "Location", title="Bray NMDS") + 
  geom_point(size =4) + coord_fixed(expand= TRUE) + theme_bw() + theme_ord
fig_nmds_bray_Location

#NMDS 1 and 2 does not discriminate by Location.
#What would happen if we try to plot now NMDS axis 1 and 3?

fig_nmds_bray_Location <- plot_ordination(phy.filt, ord.nmds.bray, color="Location",
                                             shape = "Location", title="Bray NMDS", axes = c(1,3)) + 
  geom_point(size =4) + coord_fixed(expand= TRUE) +
  theme_bw() + theme_ord
fig_nmds_bray_Location

fig_nmds <- ggpubr::ggarrange(fig_nmds_bray_Compartment, fig_nmds_bray_Spartina, 
                              fig_nmds_bray_Location, ncol = 2, nrow = 2, common.legend = FALSE,
                              legend = "top")

#Make jpeg figure in the Workshop folder. You can use it for your report
jpeg("../nmds_fig.jpeg", width = 3000, height = 3000, res = 300)
fig_nmds
dev.off()

# Q12: How do you interpret the NMDS plots? Which factors are explained by which NMDS axes?



######PERMANOVA
#You will perform a PERMANOVA analysis in order to partition distance matrices (Bray-Curtis) 
#of our microbial communities 
#among sources of variation (factors: i.e., plant compartment, Spartina microhabitat, etc)

#Check the documentation of the adonis function in R
?adonis

distance.bray <- phyloseq::distance(phy.filt, method="bray")
sample_metadata <- data.frame(sample_data(phy.filt))

adonis(distance.bray ~ Compartment + Spartina + Location, data = sample_metadata)

# Q13: How do you interpret the results from the PERMANOVA analysis?
# What is the R2 of the different factors? How do you interpret them?

##################### RELATIVE ABUNDANCE ANALYSIS ######################

##Relative abundance analysis of putative sulfur oxidizers, sulfate reducers, and nitrifiers 
#We will do two types of plots
#1. The relative abundance of the most abundant genus from three putative functional groups:
#Nitrifiers, sulfur oxidizers, and sulfate reducers
#2. The relative abundance of all Genus from putative nitrifiers, S oxidizers, and sulfate reducers

#Open file containing list of Genus potentially associated to S oxidation and sulfate reduction 
function_taxa <- read.csv("../taxonomy_function.csv", stringsAsFactors = FALSE)
S_Ox <- function_taxa$Sulfur_oxidizer
S_Re <- function_taxa$Sulfur_reducers

#Inspect Genus
S_Ox
S_Re

#Theme for ggplot2 plots
theme_ra <- theme_bw() + theme(axis.title.x = element_blank(),
                               axis.title.y = element_text(size = 18),
                               axis.text = element_text(size = 16),
                               legend.position = "top",
                               legend.text = element_text(size = 16),
                               legend.title = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               strip.text.x = element_text(size = 14),
                               strip.background = element_rect(fill = "grey90"))

#Function to scale the axis of plots with square root transformations
mysqrt_trans <- function() {
  trans_new("mysqrt", 
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}

#Transform reads to relative abundance
phy.filt.ra  = transform_sample_counts(phy.filt, function(x) x / sum(x))

#Group ASVs into Genus level, summing their abundances
Genus_ra <- tax_glom(phy.filt.ra, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#Change names of Compartment factors for plotting purposes
sample_data(Genus_ra)$Compartment <- gsub("Bulk sediment", "Bulk",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- gsub("Rhizosphere", "Rhizo",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- gsub("Endosphere", "Endo",  sample_data(Genus_ra)$Compartment)
sample_data(Genus_ra)$Compartment <- as.factor(sample_data(Genus_ra)$Compartment)

#Melt phyloseq object into a large data.frame object
Genus_ra_melt  <- psmelt(Genus_ra)

####Subset taxa of putative nitrifiers (Genus names starting with "Nitro")
nitri_genus_phy <- subset_taxa(Genus_ra, grepl("Nitro", Genus))

#Select most abundant genus, you can play with the abund.thr value to select less or more genera
nitri_genus_phy_filt <- phyloseq_filter_prevalence(nitri_genus_phy, prev.trh = 0, 
                                                   abund.trh = 0.01, abund.type = "total",
                                                   threshold_condition = "AND")

#Melt phyloseq object into a large data.frame object
nitri_genus_melt <- psmelt(nitri_genus_phy_filt)

#Include the Order name in the Genus (for plotting purposes)
nitri_genus_melt$Genus <- paste(nitri_genus_melt$Order, nitri_genus_melt$Genus, sep = ", ")

#Subset only some columns
nitri_genus_agg <- nitri_genus_melt[,c("Seq_ID", "Spartina", "Compartment", 
                            "Location", "Genus", "Abundance")]

#Order compartment levels for plotting purposes
nitri_genus_agg$Compartment <- as.factor(nitri_genus_agg$Compartment)
nitri_genus_agg$Compartment <- factor(nitri_genus_agg$Compartment,
                                      levels(nitri_genus_agg$Compartment)[c(1,3,2)])

nitri_genus_agg$Spartina <- as.factor(nitri_genus_agg$Spartina)

#Plot most important nitrifier Genus
fig_nitri <- ggplot(nitri_genus_agg, aes (x = Compartment, y = Abundance*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 4.5), breaks = c(0,0.4,1.7, 4)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18), legend.position = "top") + 
  labs(y = "Relative abundance (%)", x = "")

fig_nitri

##Relative abundance of all nitrifiers

nitri_genus_melt_all <- psmelt(nitri_genus_phy)

#Sum the relative abundance of all nitrifier Genus by 
nitri_ra_sum <- aggregate(nitri_genus_melt_all$Abundance,
                          by=list(nitri_genus_melt=nitri_genus_melt_all$Seq_ID, 
                                  Spartina=nitri_genus_melt_all$Spartina,
                                  Compartment=nitri_genus_melt_all$Compartment, 
                                  Location=nitri_genus_melt_all$Location),
                          FUN=sum)

names(nitri_ra_sum)[5] <- "nitri_ra_sum"

#To get a percentage value, multiply by 100
nitri_ra_sum$nitri_ra_sum <- nitri_ra_sum$nitri_ra_sum*100

##Plotting
nitri_ra_sum$Compartment <- as.factor(nitri_ra_sum$Compartment)
nitri_ra_sum$Compartment <- factor(nitri_ra_sum$Compartment,
                                   levels(nitri_ra_sum$Compartment)[c(1,3,2)])

nitri_ra_sum$Spartina <- as.factor(nitri_ra_sum$Spartina)


g <- ggplot(data = nitri_ra_sum, aes(x = Compartment, y = nitri_ra_sum,
                                     fill = factor(Spartina))) 
fig_nitri_all <- g + geom_boxplot() + 
  theme_ra + labs(y = "Nitrifiers (%)") + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 5.5), breaks = c(0,0.3, 1.3, 3 ,5))

fig_nitri_all

# Q14: What trend do you see in the relative abundance of putative nitrifiers 
# between plant compartments and microhabitat of the marsh? How do you explain this trend?

#Analysis of Sulfur oxidizers

#Subset putative S oxidizers
S_Ox_genus_phy <- subset_taxa(Genus_ra, Genus %in% S_Ox)

#Select most abundant Genus, again you can play with the abund.thr value
S_Ox_genus_phy_filt <- phyloseq_filter_prevalence(S_Ox_genus_phy, prev.trh = 0, 
                                                  abund.trh = 0.8, abund.type = "total",
                                                  threshold_condition = "AND")

#Melt phyloseq object into a long dataframe object
S_Ox_genus_melt <- psmelt(S_Ox_genus_phy_filt)
#Add Order name in Genus
S_Ox_genus_melt$Genus <- paste(S_Ox_genus_melt$Order, S_Ox_genus_melt$Genus, sep = ", ")

#Select a subset of columns
S_Ox_genus_agg <- S_Ox_genus_melt[,c("Seq_ID", "Spartina", "Compartment",
                                     "Location", "Genus",  "Abundance")] 

#Rearrange factor for plotting purposes
S_Ox_genus_agg$Compartment <- as.factor(S_Ox_genus_agg$Compartment)
S_Ox_genus_agg$Compartment <- factor(S_Ox_genus_agg$Compartment,
                                     levels(S_Ox_genus_agg$Compartment)[c(1,3,2)])

S_Ox_genus_agg$Spartina <- as.factor(S_Ox_genus_agg$Spartina)

#Plot!
fig_S_ox <- ggplot(S_Ox_genus_agg, aes (x = Compartment, y = Abundance*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 60), breaks = c(0, 2, 8, 20, 38, 55)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),  legend.position = "top") + labs(y = "Relative abundance (%)", x = "")

fig_S_ox

# Q15: What sulfur oxidizing Genus is concentrated in the root endosphere of Spartina alterniflora?
# What does the term Candidatus mean?  
# If you do a quick search in Google Scholar:
# As symbionts of what organism was this Genus previously described?


#Analyze now all putative S oxidizers
S_ox_ra <- Genus_ra_melt %>% filter(Genus %in% S_Ox)
S_ox_ra$Seq_ID <- as.factor(S_ox_ra$Seq_ID)

#Add the abundance of all putative S oxiding genus per sample
S_ox_ra_sum <- aggregate(S_ox_ra$Abundance, 
                         by=list(Seq_ID=S_ox_ra$Seq_ID, Spartina=S_ox_ra$Spartina,
                                 Compartment=S_ox_ra$Compartment, Location=S_ox_ra$Location),
                         FUN=sum)

names(S_ox_ra_sum)[5] <- "S_ox_ra_sum"

S_ox_ra_sum$S_ox_ra_sum <- S_ox_ra_sum$S_ox_ra_sum*100

##Plotting

S_ox_ra_sum$Compartment <- as.factor(S_ox_ra_sum$Compartment)
S_ox_ra_sum$Compartment <- factor(S_ox_ra_sum$Compartment,
                                  levels(S_ox_ra_sum$Compartment)[c(1,3,2)])

S_ox_ra_sum$Spartina <- as.factor(S_ox_ra_sum$Spartina)


g <- ggplot(data = S_ox_ra_sum, aes(x = Compartment, y = S_ox_ra_sum, fill = Spartina))
fig_S_ox_all <- g + geom_boxplot() + 
  expand_limits(y = 0) + theme_ra + labs(y = "Sulfur oxidizers (%)")  + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 65), breaks = c(0,2,8,20, 38, 60))

fig_S_ox_all

# Q16: In what microhabitat of the marsh (Tall, Medium, Short Spartina)
# you see a greater endospheric enrichment of S oxidizers? How do would you explain this?
# Based on the description of this system, where do expect you have greater porewater sulfide concentration? 


#Analysis of Sulfate Reducers, same as before with S oxidizers
S_Re_genus_phy <- subset_taxa(Genus_ra, Genus %in% S_Re)
S_Re_genus_phy_filt <- phyloseq_filter_prevalence(S_Re_genus_phy, prev.trh = 0, 
                                                  abund.trh = 1.2, abund.type = "total",
                                                  threshold_condition = "AND")

S_Re_genus_melt <- psmelt(S_Re_genus_phy_filt)
S_Re_genus_melt$Genus <- gsub("Sva0081_sediment_group", "Sva0081", S_Re_genus_melt$Genus)
S_Re_genus_melt$Genus <- paste(S_Re_genus_melt$Order, S_Re_genus_melt$Genus, sep = ", ")


S_Re_genus_agg <- S_Re_genus_melt[,c("Seq_ID", "Spartina", "Compartment",
                                     "Location", "Genus", "Abundance")]

S_Re_genus_agg$Compartment <- as.factor(S_Re_genus_agg$Compartment)
S_Re_genus_agg$Compartment <- factor(S_Re_genus_agg$Compartment,
                                     levels(S_Re_genus_agg$Compartment)[c(1,3,2)])

S_Re_genus_agg$Spartina <- as.factor(S_Re_genus_agg$Spartina)


fig_S_Re <- ggplot(S_Re_genus_agg, aes (x = Compartment, y = Abundance*100, fill = Spartina)) +
  geom_boxplot() + facet_wrap(~Genus, ncol = 1) + coord_flip() + theme_ra +
  scale_y_continuous(trans="mysqrt", limits=c(0, 43), breaks = c(0, 2.5, 9, 22, 40)) + 
  theme(axis.text.x = element_text(angle = 0), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18),  legend.position = "top") + labs(y = "Relative abundance (%)", x = "")

#Inspect the figure of putative sulfate reducers by compartment and marsh microhabitat
fig_S_Re

#Analyze all putative sulfate reducers 
S_re_ra <- Genus_ra_melt %>% filter(Genus %in% S_Re)
S_re_ra$Seq_ID <- as.factor(S_re_ra$Seq_ID )

S_re_ra_sum <- aggregate(S_re_ra$Abundance, by=list(Seq_ID=S_re_ra$Seq_ID, Spartina=S_re_ra$Spartina,
                                                    Compartment=S_re_ra$Compartment, Location=S_re_ra$Location,
                                                    Depth=S_re_ra$Depth), FUN=sum)

names(S_re_ra_sum)[6] <- "S_re_ra_sum"

S_re_ra_sum$S_re_ra_sum <- S_re_ra_sum$S_re_ra_sum*100

#Plot
S_re_ra_sum$Compartment <- as.factor(S_re_ra_sum$Compartment)
S_re_ra_sum$Compartment <- factor(S_re_ra_sum$Compartment,
                                  levels(S_re_ra_sum$Compartment)[c(1,3,2)])

S_re_ra_sum$Spartina <- as.factor(S_re_ra_sum$Spartina)

g <- ggplot(data = S_re_ra_sum, aes(x = Compartment, y = S_re_ra_sum, fill = Spartina))
fig_S_Re_all <- g + geom_boxplot() + 
  expand_limits(y = 0) + theme_ra + labs(y = "Sulfate reducers (%)")  + 
  scale_y_continuous(trans="mysqrt", limits=c(0, 55), breaks = c(0,3,12,28, 50))

fig_S_Re_all

# Q17: In what microhabitat of the marsh (Tall, Medium, Short Spartina) you see a
# greater endospheric enrichment of putative sulfate reducers? How do would you explain this?

fig_function <- ggpubr::ggarrange(fig_nitri_all, fig_S_ox_all, fig_S_Re_all,
                                  fig_nitri, fig_S_ox, fig_S_Re,
                                  labels = c("", "", "", "", "", ""),
                        nrow = 2, ncol = 3, heights = c(1,2), align = "v",
                        font.label = c(size = 28), common.legend = TRUE,
                        legend = "top")

fig_function

#Save jpeg figure summarizing all rel. abundance plots,
#you can use it for your report!

jpeg("../fig_function.jpeg", res = 300, width = 4800, height = 3000)
fig_function
dev.off()
