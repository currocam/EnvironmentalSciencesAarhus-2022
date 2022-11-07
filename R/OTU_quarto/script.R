#sessionInfo()
# args <- commandArgs(trailingOnly = TRUE)
# cat(args, sep = "\n")

if (!require("BiocManager")){install.packages("BiocManager")}
if (!require("phyloseq")){BiocManager::install("phyloseq")}
if (!require("tidyverse")){install.packages("tidyverse")}

otu <- readr::read_csv("OTU/data/02OTU.csv") %>%
  phyloseq::otu_table(taxa_are_rows = TRUE)
taxa <- readr::read_csv("OTU/data/02TAX.csv")%>%
  as.matrix() %>%
  phyloseq::tax_table()

#create metadata table from otu/tax table info
n_samples <- 6

meta <- tibble(
  rowname = colnames(otu),
  sample = rep(1:n_samples, ncol(otu)/n_samples),
  timepoint = c(rep("T0", n_samples), rep("T1", n_samples), rep("T2", n_samples), rep("T4", n_samples)),
  treatment = rep(c(rep("Control", 3), rep("Alfalfa", 3)), 4)
) %>%
  tidyr::unite("sample_time", c("timepoint", "sample"), remove = FALSE) %>%
  tidyr::unite("sample_type", c("sample", "treatment"), remove = FALSE) %>%
  tidyr::unite("treatment_timepoint", c("treatment", "timepoint"), remove = FALSE) %>%
  tidyr::unite("treatment_timepoint", c("treatment", "timepoint"), remove = FALSE) %>%
  tidyr::unite("sample_treatment_timepoint", c("sample", "treatment", "timepoint"), remove = FALSE) %>%
  tibble::column_to_rownames() %>%
  phyloseq::sample_data()
  
# Combine OTU table, taxonomy table and metadata into a phyloseq object.
(physeq <- phyloseq(otu, taxa, meta))

# Remove prefix
silva_prefix <- c(
  "k_", "p_", "c_", "o_", "f_", "g_", "s_"
)

tax_table(physeq) <- silva_prefix %>%
  paste(collapse = "|") %>%
  gsub("", tax_table(physeq))
# Remove anything labelled "unknown". This will not affect abundances,
# only the "name" of the species/genus/etc. level. 
tax_table(physeq) <- gsub("Unknown.*", "", tax_table(physeq))
# Replace spaces with underscores:
tax_table(physeq) <- gsub(" ", "_", tax_table(physeq))



###----Make unique OTU labels----
makeTaxLabel <- function(physeq){
  tax_table(physeq) %>%
    as.data.frame() %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        ~ na_if(.,"")
        )) %>%
    purrr::transpose() %>%
    purrr::map_chr(
      ~as.character(.x) %>%
        na.omit() %>%
        tail(1)
    ) %>%
    make.unique()
}
# Asign new names
taxa_names(physeq) <- makeTaxLabel(physeq)

head(otu_table(physeq))
tail(otu_table(physeq))

head(tax_table(physeq))
tail(tax_table(physeq))


head(tax_table(physeq), n = 30)

head(otu_table(physeq))

###----Ordination_all----

(physeq_nmds <- ordinate(physeq, method = "NMDS", distance = "bray"))

# Shepard plot - all
vegan::stressplot(physeq_nmds)

# NMDS plot - all
library('rwantshue')
scheme <- iwanthue(seed = 1, force_init = TRUE)
get_wants_hue <- function(n, singleton,palette = "colorblind_friendly"){
 singleton$hex(n,color_space = hcl_presets[[palette]])
}

(p1 <- physeq %>%
  plot_ordination(
    physeq_nmds,
    type ="samples",
    color = "timepoint",
    shape = "treatment",
    title = "NMDS, bray-curtis dissimilarity"
    )+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values= get_wants_hue(4, scheme)) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
)

#test adonis

library(vegan)
library(ggplot2)
library(grid)

metadata <- as(sample_data(physeq), "data.frame")

adonis(distance(physeq, method="bray") ~ treatment,
       data = metadata)

adonis(distance(physeq, method="bray") ~  treatment_timepoint,
       data = metadata)

###test pairwise adonis (fra Marie R.A.)


permanova <- vegan::adonis(t(otu)~ treatment,
                           data = metadata, permutations=999, method = "bray")

permanova$aov.tab
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

post_hoc_permanova <- pairwise.adonis(t(otu), metadata$treatment_timepoint, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova


###----Richness all----

p1 <- plot_richness(
  physeq, x="treatment_timepoint",
  color = "sample", shape ="treatment",
  measures=c("Shannon", "InvSimpson")
  )
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_all_calc_data_wo_plant.tsv', quote=FALSE, sep='\t')

###----Stacked Barplots all----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_all_mean <- merge_samples(physeq, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_all_mean <- transform_sample_counts(physeq_all_mean, function(x) x/ 3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_all_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_all_mean)

tax_table(physeq_all_mean)

unique(tax_table(physeq_all_mean)[,"Phylum"] )

unique(tax_table(physeq_all_mean)[,"Kingdom"] )

#agglomeration on the Phylum level
physeq_all_mean_phylum <- tax_glom(physeq_all_mean, taxrank = "Phylum")

#agglomeration on the Phylum level
physeq_all_mean_kingdom <- tax_glom(physeq_all_mean, taxrank = "Kingdom")


#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_all_mean
physeq_all_mean_phylum
physeq_all_mean_kingdom

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_all_mean_phylum)
sample_sums(physeq_all_mean_kingdom)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

##transform phyloseq object to a data frame (DF)
physeq_all_mean_phylumDF<- psmelt(physeq_all_mean_phylum)
physeq_all_mean_kingdomDF<- psmelt(physeq_all_mean_kingdom)

#inspect the dataframe
str(physeq_all_mean_phylumDF)

str(physeq_all_mean_kingdomDF)

#make the phyla characters, not factors
physeq_all_mean_phylumDF$Phylum <- as.character(physeq_all_mean_phylumDF$Phylum)

#make the kingdoms characters, not factors
physeq_all_mean_kingdomDF$Kingdom <- as.character(physeq_all_mean_kingdomDF$Kingdom)

#add new column with renamed low abundant taxa
physeq_all_mean_phylumDF <- physeq_all_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 1, "< 1%"))

#add new column with renamed low abundant taxa
physeq_all_mean_kingdomDF <- physeq_all_mean_kingdomDF %>% 
  mutate(Kingdom2 = replace(Kingdom, Abundance < 0.5, "< 0.5%"))


#check all phyla names
unique(physeq_all_mean_phylumDF$Phylum2)

#check all phyla names
unique(physeq_all_mean_kingdomDF$Kingdom2)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_all_mean_phylumDF <- physeq_all_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "", "unassigned Bacteria"))

#reorder the phyla so that they are stacked according to abundance
physeq_all_mean_phylumDF$Phylum2 <- reorder(physeq_all_mean_phylumDF$Phylum2,
                                            physeq_all_mean_phylumDF$Abundance)

ggplot(physeq_all_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = get_wants_hue(length(physeq_all_mean_phylumDF$Phylum2), scheme)) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#plot kingdoms
physeq_all_mean_kingdomDF$Kingdom <- reorder(physeq_all_mean_kingdomDF$Kingdom,
                                             physeq_all_mean_kingdomDF$Abundance)


ggplot(physeq_all_mean_kingdomDF, aes(Sample, Abundance, fill=Kingdom)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b25b3b",
                               "#9751e1",
                               "#59c13e",
                               "#583395",
                               "#b3b82f",
                               "#ce54b9",
                               "#6b9a3c",
                               "#7781d4",
                               "#db8b2a",
                               "#6991ba",
                               "#d7452e",
                               "#54ac6d",
                               "#cb4577",
                               "#4f9e8c",
                               "#a85857",
                               "#a98a3d",
                               "#8e5a83",
                               "#727346")) +
  labs(y= "Relative abundance [%]",
       fill= "Kingdom") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

###----ANOVA----

library(BiocManager)
BiocManager::install("microbiome")

library("microbiome"); packageVersion("microbiome")
library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq)
meta1 <- microbiome::meta(physeq)


permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_all.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453", "#800080"), cex=2, lwd=1) 
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()

###----Subsetting----

#Make a phyloseq object that only contains OTUs belonging to a specific group or clade.

physeq_bac <- subset_taxa(physeq, Kingdom=="Bacteria")
physeq_bac
head(otu_table(physeq_bac))
head(tax_table(physeq_bac))

physeq_arc <- subset_taxa(physeq, Kingdom=="Archaea")
physeq_arc
head(otu_table(physeq_arc))
head(tax_table(physeq_arc))

physeq_opi <- subset_taxa(physeq, Kingdom=="Opisthokonta")
physeq_opi
head(otu_table(physeq_opi))
head(tax_table(physeq_opi))

physeq_sar <- subset_taxa(physeq, Kingdom=="SAR")
physeq_sar
head(otu_table(physeq_sar))
head(tax_table(physeq_sar))

#physeq_pla <- subset_taxa(physeq, Kingdom=="Archaeplastida")
#physeq_pla
#head(otu_table(physeq_pla))
#head(tax_table(physeq_pla))

###----Richness bac----

p1 <- plot_richness(physeq_bac, x="treatment_timepoint", color = "sample", shape ="treatment", measures=c("Shannon", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_bac_calc_data.tsv', quote=FALSE, sep='\t')

###----Richness arc----

p1 <- plot_richness(physeq_arc, x="treatment_timepoint", color = "sample", shape ="treatment", measures=c("Shannon", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_arc_calc_data.tsv', quote=FALSE, sep='\t')

###----Richness sar----

p1 <- plot_richness(physeq_sar, x="treatment_timepoint", color = "sample", shape ="treatment", measures=c("Shannon", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_sar_calc_data.tsv', quote=FALSE, sep='\t')

###----Richness opi----

p1 <- plot_richness(physeq_opi, x="treatment_timepoint", color = "sample", shape ="treatment", measures=c("Shannon", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_opi_calc_data.tsv', quote=FALSE, sep='\t')

###----Richness pla----

p1 <- plot_richness(physeq_pla, x="treatment_timepoint", color = "sample", shape ="treatment", measures=c("Shannon", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill="white" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

# write the calculations in a file
write.table(p1$data, file='richness_pla_calc_data.tsv', quote=FALSE, sep='\t')

####----Ordination_bac----

physeq_nmds_bac <- ordinate(physeq_bac, method = "NMDS", distance = "bray")

# Shepard plot - bac
vegan::stressplot(physeq_nmds_bac)

# NMDS plot - bac

p1 = plot_ordination(physeq_bac, physeq_nmds_bac, type ="samples", color = "timepoint", shape = "treatment", title = "NMDS, bray-curtis dissimilarity")+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c("#2f85fe", "#e05436", "#009453", "#800080")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

####----Ordination_arc----

physeq_nmds_arc <- ordinate(physeq_arc, method = "NMDS", distance = "bray")

# Shepard plot - arc
vegan::stressplot(physeq_nmds_arc)

# NMDS plot - arc

p1 = plot_ordination(physeq_arc, physeq_nmds_arc, type ="samples", color = "timepoint", shape = "treatment", title = "NMDS, bray-curtis dissimilarity")+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c("#2f85fe", "#e05436", "#009453", "#800080")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

####----Ordination_sar----

physeq_nmds_sar <- ordinate(physeq_sar, method = "NMDS", distance = "bray")

# Shepard plot - sar
vegan::stressplot(physeq_nmds_sar)

# NMDS plot - sar

p1 = plot_ordination(physeq_sar, physeq_nmds_sar, type ="samples", color = "timepoint", shape = "treatment", title = "NMDS, bray-curtis dissimilarity")+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c("#2f85fe", "#e05436", "#009453", "#800080")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

####----Ordination_opi----

physeq_nmds_opi <- ordinate(physeq_opi, method = "NMDS", distance = "bray")

# Shepard plot - opi
vegan::stressplot(physeq_nmds_opi)

# NMDS plot - opi

p1 = plot_ordination(physeq_opi, physeq_nmds_opi, type ="samples", color = "timepoint", shape = "treatment", title = "NMDS, bray-curtis dissimilarity")+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c("#2f85fe", "#e05436", "#009453", "#800080")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

####----Ordination_pla----

physeq_nmds_pla <- ordinate(physeq_pla, method = "NMDS", distance = "bray")

# Shepard plot - pla
vegan::stressplot(physeq_nmds_pla)

# NMDS plot - pla

p1 = plot_ordination(physeq_pla, physeq_nmds_pla, type ="samples", color = "timepoint", shape = "treatment", title = "NMDS, bray-curtis dissimilarity")+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c("#2f85fe", "#e05436", "#009453", "#800080")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

####----ANOVA bac----

# library("microbiome"); packageVersion("microbiome")
# library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_bac)
meta1 <- microbiome::meta(physeq_bac)

permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_bac.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()

####----ANOVA arc----

# library("microbiome"); packageVersion("microbiome")
# library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_arc)
meta1 <- microbiome::meta(physeq_arc)

permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_arc.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()

####----ANOVA sar----

# library("microbiome"); packageVersion("microbiome")
# library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_sar)
meta1 <- microbiome::meta(physeq_sar)

permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_sar.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()


####----ANOVA opi----

# library("microbiome"); packageVersion("microbiome")
# library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_opi)
meta1 <- microbiome::meta(physeq_opi)

permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_opi.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()

####----ANOVA pla----
#NOT DONE HERE
# library("microbiome"); packageVersion("microbiome")
# library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_pla)
meta1 <- microbiome::meta(physeq_pla)

permanova <- vegan::adonis(t(otu1)~treatment_timepoint,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_pla.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)
#NOT DONE HERE
dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$treatment_timepoint, type="centroid")
mod
mod$treatment_timepoint <- meta$treatment_timepoint
TukeyHSD(mod)

plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

dev.off()


###----Stacked Barplots bac----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_bac_mean <- merge_samples(physeq_bac, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_bac_mean <- transform_sample_counts(physeq_bac_mean, function(x) x/3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_bac_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_bac_mean)

#agglomeration on the Phylum level
physeq_bac_mean_phylum <- tax_glom(physeq_bac_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_bac_mean
physeq_bac_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_bac_mean_phylum)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame (DF)
physeq_bac_mean_phylumDF<- psmelt(physeq_bac_mean_phylum)

#inspect the dataframe
str(physeq_bac_mean_phylumDF)

#make the phyla characters, not factors
physeq_bac_mean_phylumDF$Phylum <- as.character(physeq_bac_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_bac_mean_phylumDF <- physeq_bac_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 1, "< 1%"))

#check all phyla names
unique(physeq_bac_mean_phylumDF$Phylum2)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_bac_mean_phylumDF <- physeq_bac_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned Bacteria"))

#reorder the phyla so that they are stacked according to abundance
physeq_bac_mean_phylumDF$Phylum2 <- reorder(physeq_bac_mean_phylumDF$Phylum2,
                                            physeq_bac_mean_phylumDF$Abundance)

ggplot(physeq_bac_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#6b6bd6",
                               "#928e22",
                               "#ab70d2",
                               "#5659a8",
                               "#bf7324",
                               "#5488e3",
                               "#6a7f2f",
                               "#9363a9",
                               "#449c67",
                               "#d24952",
                               "#378f9a",
                               "#b75236",
                               "#5d87cb",
                               "#976d2f",
                               "#3b3a72",
                               "#506c30",
                               "#b34f72",
                               "#346237",
                               "#c25064",
                               "#26513e",
                               "#827db8",
                               "#453f1b",
                               "#498ab7",
                               "#712925",
                               "#518167",
                               "#4c3c62",
                               "#818d60",
                               "#9a80a3",
                               "#967044",
                               "#507491")) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

###----Stacked Barplots arc----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_arc_mean <- merge_samples(physeq_arc, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_arc_mean <- transform_sample_counts(physeq_arc_mean, function(x) x/3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_arc_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_arc_mean)

#agglomeration on the Phylum level
physeq_arc_mean_phylum <- tax_glom(physeq_arc_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_arc_mean
physeq_arc_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_arc_mean_phylum)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame (DF)
physeq_arc_mean_phylumDF<- psmelt(physeq_arc_mean_phylum)

#inspect the dataframe
str(physeq_arc_mean_phylumDF)

#make the phyla characters, not factors
physeq_arc_mean_phylumDF$Phylum <- as.character(physeq_arc_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_arc_mean_phylumDF <- physeq_arc_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 0.001, "< 0.001%"))

#check all phyla names
unique(physeq_arc_mean_phylumDF$Phylum)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_arc_mean_phylumDF <- physeq_arc_mean_phylumDF %>% 
  mutate(Phylum = replace(Phylum, Phylum == "NA", "unassigned arcteria"))

#reorder the phyla so that they are stacked according to abundance
physeq_arc_mean_phylumDF$Phylum <- reorder(physeq_arc_mean_phylumDF$Phylum,
                                           physeq_arc_mean_phylumDF$Abundance)

ggplot(physeq_arc_mean_phylumDF, aes(Sample, Abundance, fill=Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#928e22",
                               "#ab70d2",
                               "#5659a8",
                               "#bf7324",
                               "#5488e3",
                               "#6a7f2f",
                               "#9363a9",
                               "#449c67",
                               "#d24952",
                               "#378f9a",
                               "#b75236",
                               "#5d87cb",
                               "#976d2f",
                               "#3b3a72",
                               "#506c30",
                               "#b34f72",
                               "#346237",
                               "#c25064",
                               "#26513e",
                               "#827db8",
                               "#453f1b",
                               "#498ab7",
                               "#712925",
                               "#518167",
                               "#4c3c62",
                               "#818d60",
                               "#9a80a3",
                               "#967044")) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

###----Stacked Barplots sar----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_sar_mean <- merge_samples(physeq_sar, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_sar_mean <- transform_sample_counts(physeq_sar_mean, function(x) x/3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_sar_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_sar_mean)

#agglomeration on the Phylum level
physeq_sar_mean_phylum <- tax_glom(physeq_sar_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_sar_mean
physeq_sar_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_sar_mean_phylum)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame (DF)
physeq_sar_mean_phylumDF<- psmelt(physeq_sar_mean_phylum)

#inspect the dataframe
str(physeq_sar_mean_phylumDF)

#make the phyla characters, not factors
physeq_sar_mean_phylumDF$Phylum <- as.character(physeq_sar_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_sar_mean_phylumDF <- physeq_sar_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 0.01, "< 1%"))

#check all phyla names
unique(physeq_sar_mean_phylumDF$Phylum2)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_sar_mean_phylumDF <- physeq_sar_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned sarteria"))

#reorder the phyla so that they are stacked according to abundance
physeq_sar_mean_phylumDF$Phylum2 <- reorder(physeq_sar_mean_phylumDF$Phylum2,
                                            physeq_sar_mean_phylumDF$Abundance)

ggplot(physeq_sar_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#a74d57",
                               "#5a8a37",
                               "#697ce0",
                               "#a88229",
                               "#66479c",
                               "#61693c",
                               "#6a83bf",
                               "#b05233",
                               "#3f8777",
                               "#556086")) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

###----Stacked Barplots opi----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_opi_mean <- merge_samples(physeq_opi, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_opi_mean <- transform_sample_counts(physeq_opi_mean, function(x) x/3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_opi_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_opi_mean)

#agglomeration on the Phylum level
physeq_opi_mean_phylum <- tax_glom(physeq_opi_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_opi_mean
physeq_opi_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_opi_mean_phylum)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame (DF)
physeq_opi_mean_phylumDF<- psmelt(physeq_opi_mean_phylum)

#inspect the dataframe
str(physeq_opi_mean_phylumDF)

#make the phyla characters, not factors
physeq_opi_mean_phylumDF$Phylum <- as.character(physeq_opi_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_opi_mean_phylumDF <- physeq_opi_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 0.01, "< 1%"))

#check all phyla names
unique(physeq_opi_mean_phylumDF$Phylum2)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_opi_mean_phylumDF <- physeq_opi_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned opiteria"))

#reorder the phyla so that they are stacked according to abundance
physeq_opi_mean_phylumDF$Phylum2 <- reorder(physeq_opi_mean_phylumDF$Phylum2,
                                            physeq_opi_mean_phylumDF$Abundance)

ggplot(physeq_opi_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#a74d57",
                               "#5a8a37",
                               "#697ce0",
                               "#a88229",
                               "#66479c",
                               "#61693c",
                               "#6a83bf",
                               "#b05233",
                               "#3f8777",
                               "#556086")) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

###----Stacked Barplots pla----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_pla_mean <- merge_samples(physeq_pla, "treatment_timepoint")

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_pla_mean <- transform_sample_counts(physeq_pla_mean, function(x) x/3)

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_pla_mean)

#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level. 
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_pla_mean)

#agglomeration on the Phylum level
physeq_pla_mean_phylum <- tax_glom(physeq_pla_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_pla_mean
physeq_pla_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_pla_mean_phylum)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame (DF)
physeq_pla_mean_phylumDF<- psmelt(physeq_pla_mean_phylum)

#inspect the dataframe
str(physeq_pla_mean_phylumDF)

#make the phyla characters, not factors
physeq_pla_mean_phylumDF$Phylum <- as.character(physeq_pla_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_pla_mean_phylumDF <- physeq_pla_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 0.01, "< 1%"))

#check all phyla names
unique(physeq_pla_mean_phylumDF$Phylum2)

#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_pla_mean_phylumDF <- physeq_pla_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned plateria"))

#reorder the phyla so that they are stacked according to abundance
physeq_pla_mean_phylumDF$Phylum2 <- reorder(physeq_pla_mean_phylumDF$Phylum2,
                                            physeq_pla_mean_phylumDF$Abundance)

ggplot(physeq_pla_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ad4248",
                               "#5e87d3",
                               "#d9624b",
                               "#c5a73f",
                               "#d77bb7",
                               "#cd772c",
                               "#e1556e",
                               "#38dbda",
                               "#7d59bf",
                               "#cf3f80",
                               "#ac4258",
                               "#adae5f",
                               "#a27934",
                               "#7d245b",
                               "#c486da",
                               "#533583",
                               "#9b4429",
                               "#ba4aa3",
                               "#538133",
                               "#64c36f",
                               "#c45a83",
                               "#45ba8a",
                               "#91b23e",
                               "#6a7de6")) +
  labs(y= "Relative abundance [%]",
       fill= "Phlya") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/

