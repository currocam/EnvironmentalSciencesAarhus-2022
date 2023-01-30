###########################
# New version
# Code written by Thanassis Zervas 
# az@envs.au.dk
# Creative commons, 2021. 

###----Load packages----

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("tidyverse"); packageVersion("tidyverse")

###----Import data----
setwd("maria_plots/")
set.seed(42)

#BEWARE! Are your numbers using "." or "," as the decimal counter? If "," - use read.csv2() instead!
# otu <- read.csv(file = 'RNA_OTU_rel_tripl_new.csv', sep = ';', header = TRUE, fileEncoding="UTF-8-BOM")
otu <- read.csv(file = 'RNA_OTU_tripl_wo1.csv', sep = ';', header = TRUE, fileEncoding="UTF-8-BOM")
tax <- as.matrix(read.csv2(file = 'RNA_TAX_onlyanno.csv', sep = ';', header = TRUE, fileEncoding="UTF-8-BOM"))

#create metadata table from otu/tax table info
meta <- as.data.frame(colnames(otu))
meta$layer = c("AL","AL","AL", "TZ1", "TZ1", "TZ1","TZ2", "TZ2", "PF","AL","AL","AL","AL", "TZ1", "TZ1","TZ1","TZ2", "TZ2", "PF","AL","AL","AL","AL", "TZ1", "TZ1", "TZ1","TZ2", "TZ2", "PF")
meta$age = c("AY","AY","AY", "AM", "AM", "AO","AO", "AO", "AO","AY","AY","AY","AY", "AM", "AM","AO","AO", "AO", "AO","AY","AY","AY","AY", "AM", "AM", "AO","AO", "AO", "AO")
meta$depth = c("20","30","40", "50", "60", "70","80", "90", "100","10","20","30","40", "50", "60","70","80", "90", "100","10","20","30","40", "50", "60", "70","80", "90", "100")
meta$SOM = c("5.39","10.19","13.68","2.59","2.93","1.52","1.1","1","1.04","8.73","5.39","10.19","13.68","2.59","2.93","1.52","1.1","1","1.04","8.73","5.39","10.19","13.68","2.59","2.93","1.52","1.1","1","1.04")
meta$H2O = c("22.52","22.05","26.57","7.73","15.61","8.72","6.35","7.96","7.8","28.8","22.52","22.05","26.57","7.73","15.61","8.72","6.35","7.96","7.8","28.8","22.52","22.05","26.57","7.73","15.61","8.72","6.35","7.96","7.8")
meta$age_num = c("1.13","1.16","1.2","2635","3770","26500","22100","26200","26200","1.04","1.13","1.16","1.2","2635","3770","26500","22100","26200","26200","1.04","1.13","1.16","1.2","2635","3770","26500","22100","26200","26200")
meta$pH = c("4.02","4.29","4.25","4.64","4.13","4.86","4.48","4.57","4.91","4.22","4.02","4.29","4.25","4.64","4.13","4.86","4.48","4.57","4.91","4.22","4.02","4.29","4.25","4.64","4.13","4.86","4.48","4.57","4.91")
meta$protozoa = c("1.71","2.18","3.72","8.39","3.26","6.17","5.66","4.49","2.42","2.89","1.79","2.10","2.59","5.02","6.98","5.72","5.67","5.34","4.91","3.89","1.41","2.16","3.28","4.58","6.43","9.29","4.48","4.85","3.41")
meta$prokarya = c("94.26","94.46","90.70","83.02","88.80","85.87","83.51","80.88","89.45","81.77","93.40","91.89","93.30","86.66","81.25","86.12","84.13","83.86","82.96","77.93","95.55","94.20","90.18","83.84","82.92","76.37","80.65","71.65","86.27")
meta$bacpred = c("5.35","7.20","13.43","13.43","3.45","4.05","1.97","1.02","0.56","6.44","5.28","7.22","13.78","13.39","4.92","2.74","1.72","1.21","4.78","7.71","4.98","7.16","13.92","7.57","15.55","6.25","2.74","0.96","1.20")

meta
rownames(meta) = meta$`colnames(otu)`


###----Make Phyloseq object----
OTU = otu_table(otu, taxa_are_rows = TRUE); head(OTU)
TAX = tax_table(tax); head(TAX)
META <- sample_data(meta)
head(META)

# Combine OTU table, taxonomy table and metadata into a phyloseq object.
physeq <- phyloseq(OTU, TAX, META)
physeq
taxa_names(physeq) 
# OTUnames <- read_file("OTUs.txt", sep = ",")

backup_physeq_all <- physeq

# Check if everything looks good
sample_sums(physeq)
sample_names(physeq)
rank_names(physeq)
sample_variables(physeq)
otu_table(physeq)[1:3, 1:2]
taxa_names(physeq)[1:5] 

###----Tidy up names----
colnames(tax_table(physeq)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                 o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(physeq)

#Remove anything labelled "unknown". This will not affect abundances, only the "name" of the species/genus/etc. level. 
tax_table(physeq) <- gsub("Unknown.*", "", tax_table(physeq))     #this doesn't seem to work...
tax_table(physeq) <- gsub("incer*ae sedis*", "", tax_table(physeq))

###----Make unique OTU labels----
makeTaxLabel <- function(OTU, mydata){
  
  # Makes a string label using the lowest informative tax level
  #
  # Args:
  #   OTU: OTU number
  #   mydata: the phyloseq object with the tax table
  #
  # Returns:
  #   a tax name
  OTU <- as.character(OTU)  # the OTU numbers are stored as character not integer!
  taxstrings <- as.character(tax_table(mydata)[OTU])
  #taxstrings <- as.character(taxTab(mydata))
  empty_strings <- c("k_", "p_", "c_", "o_", "f_", "g_", "s_", " ", "", NA)
  tax_name <- NA
  tax_level <- length(taxstrings)  # start at lowest tax level
  
  
  while(is.na(tax_name) |
        (tax_name %in% empty_strings)){
    tax_name  <- taxstrings[tax_level]
    tax_level <- tax_level -1
  }
  tax_name
}

# remove chloroplast and undefined annotations on kingdom level
physeq <-
  subset_taxa(physeq, Kingdom != "Apusozoa kingdom incertae sedis" | is.na(Kingdom)& 
                Kingdom != "Calypogeia muelleriana (Chloroplast)"| is.na(Kingdom)& 
                Kingdom != "Elliptochloris bilobata (Chloroplast)"| is.na(Kingdom)& 
                Kingdom != "Carex siderosticta (Chloroplast)"| is.na(Kingdom)& 
                Kingdom != "Cryptomycota Holomycota kingdom incertae sedis"| is.na(Kingdom)& 
                Kingdom != "Breviatea kingdom incertae sedis"| is.na(Kingdom))


mynames = NULL
for (i in 1:length(taxa_names(physeq))){
  mynames <- rbind(mynames, c(makeTaxLabel(taxa_names(physeq)[i],physeq)))
}

mynames <- make.unique(mynames)

head(taxa_names(physeq))
taxa_names(physeq) = mynames

head(otu_table(physeq))
tail(otu_table(physeq))

head(tax_table(physeq))
tail(tax_table(physeq))

# Remove prefixes from tax table
tax_table(physeq) <- gsub("s__", "", tax_table(physeq))
tax_table(physeq) <- gsub("g__", "", tax_table(physeq))
tax_table(physeq) <- gsub("f__", "", tax_table(physeq))
tax_table(physeq) <- gsub("o__", "", tax_table(physeq))
tax_table(physeq) <- gsub("c__", "", tax_table(physeq))
tax_table(physeq) <- gsub("p__", "", tax_table(physeq))
tax_table(physeq) <- gsub("k__", "", tax_table(physeq))

# Replace spaces with underscores:
# tax_table(physeq) <- gsub(" ", "_", tax_table(physeq))

head(tax_table(physeq), n = 30)

head(otu_table(physeq))

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_rel <- transform_sample_counts(physeq, function(x) x/sum(x)*100)
sample_sums(physeq_rel)

#reating rel counts per depth
physeq_mean_rel <- merge_samples(physeq_rel, "depth")             # this is what will be depicted in the graph later

#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_mean_rel <- transform_sample_counts(physeq_mean_rel, function(x) x/sum(x)*100)       #careful when duplicates // /sum(x)*100
sample_sums(physeq_mean_rel)
head(tax_table(physeq_mean_rel))

# # Agglomerate at certain level.
# # This maneuver combines all duplicates at species level, OTU abundances are summed.
# physeq_p <- tax_glom(physeq_rel, taxrank = "Phylum", NArm = FALSE)
# physeq_p
# head(tax_table(physeq_p))
# sample_sums(physeq_p)
# 
# physeq_c <- tax_glom(physeq_rel, taxrank = "Class", NArm = FALSE)
# physeq_c
# head(tax_table(physeq_c))
# sample_sums(physeq_c)
# 
# physeq_g <- tax_glom(physeq_rel, taxrank = "Genus", NArm = FALSE)
# physeq_g
# head(tax_table(physeq_g))
# sample_sums(physeq_g)

###----Ordination_all----

physeq_nmds <- ordinate(physeq_rel, method = "NMDS", distance = "bray")

# Shepard plot - all
vegan::stressplot(physeq_nmds)

# NMDS plot - all
p1 = plot_ordination(physeq_rel, physeq_nmds, type =" samples", color = "layer", shape = "age", title = "NMDS, bray-curtis dissimilarity total rRNA OTUs")+
  geom_point(size=3)+
  theme_bw()+
  # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))+
  stat_ellipse()
print(p1)
#dev.print(pdf, 'NMDS_all.pdf')

# # NMDS plot - phyla
# p1 = plot_ordination(physeq_p, physeq_nmds, type =" samples", color = "layer", shape = "age", title = "NMDS, bray-curtis dissimilarity total rRNA phyla")+
#   geom_point(size=3)+
#   theme_bw()+
#   # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
#   theme(legend.title = element_blank())+
#   theme(strip.background = element_rect(fill="white" ))
# print(p1)
# dev.print(pdf, 'NMDS_all_phyla.pdf')
# # NMDS plot - class
# p1 = plot_ordination(physeq_c, physeq_nmds, type ="taxa", color = "Phylum", title = "NMDS, bray-curtis dissimilarity total rRNA")+
#   geom_point(size=3)+
#   theme_bw()+
#   # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
#   theme(legend.title = element_blank())+
#   theme(strip.background = element_rect(fill="white" ))
# print(p1)
# dev.print(pdf, 'NMDS_all_classes.pdf')


###----Richness all----

p1 <- plot_richness(physeq, x="layer", color =  "depth", measures=c("Observed","Shannon","Simpson","InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.05)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "right") +
  theme(strip.background = element_rect(fill="white" ))
p1
dev.print(pdf, 'richness_all.pdf')

# write the calculations in a file
write.table(p1$data, file='diversity_ALL.tsv', quote=FALSE, sep='\t')

###----Stacked Barplots all----

#check taxonomy level names
rank_names(physeq_mean_rel)

#agglomeration on the Phylum level
physeq_mean_rel_phylum <- tax_glom(physeq_mean_rel, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_mean_rel
physeq_mean_rel_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_mean_rel_phylum)
#After merging, we need to dived the numbers to the number of replicates, in this case "3"
physeq_mean_rel_phylum <- transform_sample_counts(physeq_mean_rel_phylum, function(x) x/sum(x)*100)       #careful when duplicates // /sum(x)*100

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to "< 1%".

#transform phyloseq object to a data frame (DF)
physeq_mean_rel_phylumDF<- psmelt(physeq_mean_rel_phylum)

#inspect the dataframe
str(physeq_mean_rel_phylumDF)

#make the phyla characters, not factors
physeq_mean_rel_phylumDF$Phylum <- as.character(physeq_mean_rel_phylumDF$Phylum)


#there are some reads that were assigned only to the kingdom level, 
# i.e. NA on the phylum level, so we will rename them
physeq_mean_rel_phylumDF <- physeq_mean_rel_phylumDF %>% 
  mutate(Phylum = replace(Phylum, Phylum == "NA", "unassigned"))


#check all phyla names
unique(physeq_mean_rel_phylumDF$Phylum)

physeq_mean_rel <- transform_sample_counts(physeq_mean_rel, function(x) x/sum(x)*100)       #careful when duplicates // /sum(x)*100

# add new column with renamed low abundant taxa
physeq_mean_rel_phylumDF <- physeq_mean_rel_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 3, "< 3%"))

#Checking the final results, the sum should be 1 for each sample
sample_sums(physeq_mean_rel)


#reorder the phyla so that they are stacked according to abundance
physeq_mean_rel_phylumDF$Phylum <- reorder(physeq_mean_rel_phylumDF$Phylum2,
                                           physeq_mean_rel_phylumDF$Abundance)

ggplot(physeq_mean_rel_phylumDF, aes(Sample, Abundance, fill=Phylum)) +
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = c( "#e9b5eb",
  #                               "#f6be0e",
  #                               "#6b00ac",
  #                               "#7adb86",
  #                               "#fe00b5",
  #                               "#52dcb3",
  #                               "#8262ff",
  #                               "#829000",
  #                               "#001586",
  #                               "#eb5f00",
  #                               "#5b7cff",
#                               "#a00013",
#                               "#00b4df",
#                               "#ff6f6e",
#                               "#0273c9",
#                               "#6e001a",
#                               "#a698ff",
#                               # "#8f0049",
#                               # "#005975",
#                               # "#aba9ff",
#                               "#ffa8ba")) +
labs(x= "depth [cm]",y= "Relative abundance [%]",
     fill= "Phyla") +
  theme_bw()
dev.print(pdf, 'barplot_all.pdf')
#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/


#### ANOVA ####

library("microbiome"); packageVersion("microbiome")
library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq)
meta1 <- microbiome::meta(physeq)

otu1
permanova_depth <- vegan::adonis(t(otu1)~depth, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_depth$aov.tab))
permanova_age <- vegan::adonis(t(otu1)~age, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_age$aov.tab))
permanova_layer <- vegan::adonis(t(otu1)~layer, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_layer$aov.tab))
permanova_SOM <- vegan::adonis(t(otu1)~SOM, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_SOM$aov.tab))
permanova_H2O <- vegan::adonis(t(otu1)~H2O, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_H2O$aov.tab))
permanova_pH <- vegan::adonis(t(otu1)~pH, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_pH$aov.tab))
permanova_age_num <- vegan::adonis(t(otu1)~age_num, data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova_age_num$aov.tab))
# permanova_myxococ <- vegan::adonis(t(otu1)~bacpred, data = meta1, permutations=999, method = "bray")
# print(as.data.frame(permanova_myxococ$aov.tab))
# permanova_protozoa <- vegan::adonis(t(otu1)~protozoa, data = meta1, permutations=999, method = "bray")
# print(as.data.frame(permanova_protozoa$aov.tab))

#same problem as with vegan: multicollinearity produces unrealistic significances!!!

write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_all.tsv', quote=FALSE, sep='\t')

#Calculate beta dispersion (in this case - the dispersion between the kits)

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$age, type="centroid")
mod
mod$age <- meta$age
TukeyHSD(mod)
write.table(print(results), file='tukey_age_all.tsv', quote=FALSE, sep='\t')

mod <- vegan::betadisper(dist, meta1$layer, type="centroid")
mod
mod$layer <- meta$layer
TukeyHSD(mod)
write.table(TukeyHSD(mod), file='tukey_layer_all.tsv', quote=FALSE, sep='\t')


par(mfrow = c(1,2))
plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA, layer", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="layer", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))
dev.print(pdf, 'PCOA_boxplot_all.pdf')
par(mfrow = c(1,1))
# dev.off()
mod <- vegan::betadisper(dist, meta1$depth, type="centroid")
mod
mod$depth <- meta$depth
TukeyHSD(mod)

write.table(TukeyHSD(mod), file='tukey_depth_all.tsv', quote=FALSE, sep='\t')

# #### COmmon and specialist taxa ####
# data(GlobalPatterns)
# 
# df = psmelt(physeq_rel)
# # needs refining, as currently only takes OTUs that occurr 10x
# otus_nsites_observed =df %>% group_by(OTU, depth) %>% filter(Abundance >0.001) %>% tally() %>% filter(n>0) %>% group_by(OTU) %>% tally()
# 
# common_ASVS = otus_nsites_observed %>% filter(n==10)
# common_ASVS #1630 OTUs occurr everywhere
# unique_ASVs = otus_nsites_observed %>% filter(n==1)
# unique_ASVs # 1 OTU occurs only in one sample
# head(common_ASVS)
# head(unique_ASVs)
# 
# physeq_unique <- subset_taxa(physeq_rel, OTU == "s__.723")
#### methanogens and trophs ####
# ## from Liu 2015 taxonomy methanogens
# methanogenic_orders =c("Methanopyrales", "Methanococcales", "Methanobacteriales", "Methanomicrobiales", "Methanosarcinales", "Methanocellales"   )
# methanogenic_genus =c("Methanobacterium","Methanobrevibacter","Methanosphaera","Methanothermobacter","Methanothermus","Methanococcus","Methanothermococcus","Methanocaldococcus","Methanotorris","Methanoculleus","Methanofollis","Methanogenium","Methanolacinia","Methanomicrobium","Methanoplanus","Methanospirillum","Methanocorpusculum","Methanocalculus","Methanolinea","Candidatus Methanoregula","Candidatus Methanosphaerula","Methanosarcina","Methanococcoides","Methanohalobium","Methanohalophilus","Methanolobus","Methanomethylovorans","Methanimicrococcus","Methanosalsum","Methanosaeta","Methermicoccus","Methanopyrus","Methanocella")
# #methanogenic_species 
# 
# arch_methanotroph_names= c( "ANME-2a-2b","ANME-2a", "ANME-2b" , "ANME-2c", "Candidatus Methanoperedens", "ANME-3", "ANME-1a", "ANME-1b")
# 
# 
# methanotroph_families = c("Methylococcaceae", "Methylocystaceae", "Beijerinckiaceae", "Methylacidiphilaceae", "Methylomirabilaceae")
# 
# #dedysh and knief 2018 Diversity and Phylogeny of Described Aerobic Methanotrophs
# methanotroph_genus = c("Methylococcus", "Methylomonas", "Methylobacter", "Methylomicrobium", "Methylosarcina", "Methylocaldum", "Methylogaea", "Methylosoma", "Methyloparacoccus", "Methyloglobulus", "Methyloprofundus", "Methylomagnum", "Methylosphaera", "Methylothermus", "Methylohalobius", "Methylomarinovum", "Methylosinus", "Methylocystis", "Methylocella", "Methylocapsa", "Methyloferula", "Methylacidiphilum", "Methylacidimicrobium", "Crenothrix", "Methylobacterium-Methylorubrum","Methylobacterium", "Methylorubrum","Candidatus Methylomirabilis")
# 
# 
# # r methanogen barplot on Orderlevel}
# methanogendf = psmelt(physeq) %>% filter(Order %in% methanogenic_orders
# ) %>% mutate(Genus = ifelse(is.na(Genus)| Genus=="uncultured",paste0("Unknown ",Family),Genus))
# 
# methanogendf2 = psmelt(physeq) %>% mutate(Order= ifelse(Order %in% methanogenic_orders, Order, "Other"),
#                                             Genus= ifelse(Order %in% methanogenic_orders, Genus, "Other"),
#                                             Genus = ifelse(is.na(Genus)| Genus=="uncultured",paste0("Unknown ",Family),Genus))
# methanogen_order= rev(c(unique(filter(methanogendf2, !(Genus %in% c("Other", "Candidatus Methanoperedens")) )$Genus), "Other","Candidatus Methanoperedens" ))
# methanogendf2$Genus <- factor(methanogendf2$Genus, levels = methanogen_order )
# 
# barplotter_col_man = function(x, tax_level, label_tax, col_codes){
#   ggplot(aes(x=depth,y=Abundance,fill = tax_level), data=x) + 
#     geom_bar(aes( fill=tax_level), stat="identity",width = 3, position = "stack") +
#     facet_grid(cols=vars(layer) )+
#     labs(y="Relative abundance %", x="Depth [cm]", fill=label_tax)+
#     theme_classic()+  coord_flip()+ 
#     scale_x_reverse(breaks=seq(0,110, by=10), minor_breaks=seq(0,110, by=5))+
#     theme_classic() +theme( axis.text.x=element_text(size=8),legend.position = "bottom")+
#     scale_fill_manual(values = col_codes, limits= rev(levels(tax_level)))+
#     guides(fill = guide_legend(ncol = 3)) +
#     ylim(0,100)
# }
# 
# barplotter = function(x, tax_level, label_tax){
#   ggplot(aes(x=depth,y=Abundance,fill = tax_level), data=x) + 
#     geom_bar(aes( fill=tax_level), stat="identity",width = 3, position = "stack") +
#     facet_grid(cols=vars(newcore) )+
#     labs(y="Relative abundance %", x="Depth [cm]", fill=label_tax)+
#     theme_classic()+  coord_flip()+ 
#     scale_x_reverse(breaks=seq(0,110, by=10), minor_breaks=seq(0,110, by=5))+
#     theme_classic() +theme( axis.text.x=element_text(size=8),legend.position = "bottom")+
#     scale_fill_manual(values = methanogen_cols, limits= rev(levels(tax_level)))+
#     guides(fill = guide_legend(ncol = 3)) +
#     ylim(0,100)
# }
# 
# barplot_methanogens = barplotter(methanogendf, methanogendf$Genus, "Genus") # +ylim(0,80)
# 
# 
# methanogen_cols= c("#ff0080",  "grey90",
#                    "#cc6190","#62be56","#cf5746","#b5b545","#688bd1","#dc9143","#52ba9f","#a966cc","#5c8b43","#997c3f")
# names(methanogen_cols) <- c(methanogen_order)
# 
# 
# barplot_methanogens2 = barplotter_col_man(methanogendf2, methanogendf2$Genus, "Genus", col_codes = methanogen_cols)
# 
# 
# methanogen_sums_sample = methanogendf %>% group_by(Sample, layer, depth, age) %>% summarise(sum=sum(Abundance)) %>% ungroup()
# methanogen_sums_cores= methanogendf %>%  group_by(layer, Genus) %>% summarise(sum=sum(Abundance)) %>% mutate(perc= sum/sum(sum)*100) %>% arrange(layer, desc(perc))
# 
# methanoperendens_sums = methanogendf %>% filter(Genus=="Candidatus Methanoperedens") %>% group_by(layer, age, depth) %>% summarise(sum=sum(Abundance)) %>% ungroup()
# barplot_methanogens2
# 
# #barplotter_col_man(methanogendf2, methanogendf2$Family, "Family", col_codes = colcodes25)



# #### Network ####
# 
# # # plot_net(physeq_p, maxdist = 0.2, type="taxa", color = "Kingdom", point_label = "Phylum")
# # 
# # ig <- make_network(physeq_mean_rel_genus, max.dist=0.2, type="taxa", distance="bray")
# # plot_network(ig, physeq_mean_rel_phylum, color = "Kingdom", type="taxa", label = "Phylum")
# # 
# # 
# # ig <- make_network(physeq_mean_rel_phylum, max.dist=0.2, type="samples", distance="bray")
# # plot_network(ig, physeq_mean_rel_phylum, color = "age", type="samples", label = "depth")
# 
# # only thaw layers
# subset_10 <-subset_samples(physeq_p, depth=="10")
# ig_10 <- make_network(subset_10, max.dist=0.1, type="taxa", distance="bray")
# subset_20 <-subset_samples(physeq_p, depth=="20")
# ig_20 <- make_network(subset_20, max.dist=0.1, type="taxa", distance="bray")
# subset_30 <-subset_samples(physeq_p, depth=="30")
# ig_30 <- make_network(subset_30, max.dist=0.1, type="taxa", distance="bray")
# subset_40 <-subset_samples(physeq_p, depth=="40")
# ig_40 <- make_network(subset_40, max.dist=0.1, type="taxa", distance="bray")
# subset_50 <-subset_samples(physeq_p, depth=="50")
# ig_50 <- make_network(subset_50, max.dist=0.1, type="taxa", distance="bray")
# subset_60 <-subset_samples(physeq_p, depth=="60")
# ig_60 <- make_network(subset_60, max.dist=0.1, type="taxa", distance="bray")
# subset_70 <-subset_samples(physeq_p, depth=="70")
# ig_70 <- make_network(subset_70, max.dist=0.1, type="taxa", distance="bray")
# subset_80 <-subset_samples(physeq_p, depth=="80")
# ig_80 <- make_network(subset_80, max.dist=0.1, type="taxa", distance="bray")
# subset_90 <-subset_samples(physeq_p, depth=="90")
# ig_90 <- make_network(subset_90, max.dist=0.1, type="taxa", distance="bray")
# subset_100 <-subset_samples(physeq_p, depth=="100")
# ig_100 <- make_network(subset_100, max.dist=0.1, type="taxa", distance="bray")
# 
# plot_network(ig_10, subset_10, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network10.pdf')
# plot_network(ig_20, subset_20, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network20.pdf')
# plot_network(ig_30, subset_30, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network30.pdf')
# plot_network(ig_40, subset_40, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network40.pdf')
# plot_network(ig_50, subset_50, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network50.pdf')
# plot_network(ig_60, subset_60, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network60.pdf')
# plot_network(ig_70, subset_70, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network70.pdf')
# plot_network(ig_80, subset_80, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network80.pdf')
# plot_network(ig_90, subset_90, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network90.pdf')
# plot_network(ig_100, subset_100, color = "Kingdom", type="taxa", label = "Phylum")
# dev.print(pdf, 'network100.pdf')



### ----Subsetting by Kingdoms----


#Make a phyloseq object that only contains OTUs belonging to a specific group or clade.

physeq_euk <- subset_taxa(physeq, Kingdom =="Opisthokonta" | Kingdom =="Excavata" | Kingdom =="Archaeplastida" | Kingdom =="SAR" | Kingdom =="Amoebozoa" | Kingdom =="Hacrobia"| is.na(Kingdom))
physeq_euk
head(otu_table(physeq_euk))
head(tax_table(physeq_euk))

physeq_bac <- subset_taxa(physeq, Kingdom=="Bacteria")
physeq_bac
head(otu_table(physeq_bac))
head(tax_table(physeq_bac))

physeq_arc <- subset_taxa(physeq, Kingdom=="Archaea")
physeq_arc
head(otu_table(physeq_arc))
head(tax_table(physeq_arc))

###-----Subsetting by protozoa/prokarya/myxcoccales-----
#### trying to subset predators ####
# based on https://david-barnett.github.io/microViz/reference/tax_select.html
#eukaryotic predators
protozoa <- c("Amoebozoa","Heteronematina","Jakobida","Kinetoplastea","Telonemia","Choanoflagellida","Cercozoa","Ciliophora")
#prokaryotic predators
bactpred <- c("Bdellovibrionales","Myxococcales", "Vampirovibrionales")
#all predators
pred <- c("Amoebozoa","Heteronematina","Jakobida","Kinetoplastea","Telonemia","Choanoflagellida","Cercozoa","Ciliophora","Bdellovibrionales","Myxococcales", "Vampirovibrionales")


ps<- physeq

tax_select <- function(ps,
                       tax_list,
                       ranks_searched = "all",
                       strict_matches = FALSE,
                       n_typos = 1,
                       deselect = FALSE) {
  if (inherits(ps, "ps_extra")) ps <- ps_get(ps)
  
  # collapse tax_list to a string of regex OR patterns
  taxaString <- paste(tax_list, collapse = "|")
  
  # get tax table to search
  Taxa <- phyloseq::tax_table(ps)
  
  if (!identical(ranks_searched, "all")) {
    # check valid taxonomic ranks given for searching
    if (any(!ranks_searched %in% phyloseq::rank_names(ps))) {
      stop(
        "Invalid rank names given: ", paste(ranks_searched, collapse = " "),
        "\n- Should be any/some of: ", paste(phyloseq::rank_names(ps), collapse = "/")
      )
    }
    
    Taxa <- Taxa[, ranks_searched, drop = FALSE]
  }
  
  if (isTRUE(strict_matches)) {
    selectionVec <- apply(Taxa, MARGIN = 1, FUN = function(r) any(r %in% tax_list))
  } else if (n_typos > 0) {
    selectionVec <- apply(Taxa, MARGIN = 1, FUN = function(r) {
      any(sapply(tax_list, FUN = function(x) {
        agrepl(x, r, max.distance = n_typos)
      }))
    })
  } else {
    selectionVec <- apply(Taxa, MARGIN = 1, FUN = function(r) any(grepl(taxaString, r)))
  }
  
  # include exact rownames/taxa names matches
  ROWNAMES <- rownames(Taxa)
  selectionVec <- selectionVec | ROWNAMES %in% tax_list
  
  # stop with error if no taxa matched! (but not "in deselect mode")
  if (isFALSE(deselect) && !any(selectionVec)) {
    stop("No taxa matched.")
  }
  
  # if taxa DEselection is requested, invert the logical vector, to cause removal of matching taxa
  if (isTRUE(deselect)) {
    selectionVec <- !selectionVec
  }
  
  ps <- phyloseq::prune_taxa(x = ps, taxa = selectionVec)
  
  return(ps)
}

physeq_protozoa <- physeq %>% tax_select(tax_list = protozoa)

physeq_bactpred <- physeq %>% tax_select(tax_list = bactpred) %>% tax_select(tax_list = "Sorangium", strict_matches = FALSE, deselect = TRUE)

physeq_pred <- physeq %>% tax_select(tax_list = pred) %>% tax_select(tax_list = "Sorangium", strict_matches = FALSE, deselect = TRUE)


#PRED#
sample_sums(physeq_pred)
physeq_pred_rel <- transform_sample_counts(physeq_pred, function(x) x/sum(x)*100)
sample_sums(physeq_pred_rel)

#reating rel counts per depth
physeq_pred_mean <- merge_samples(physeq_pred, "depth")             # this is what will be depicted in the graph later
sample_sums(physeq_pred_mean)
physeq_pred_mean_rel <- transform_sample_counts(physeq_pred_mean, function(x) x/sum(x)*100)
sample_sums(physeq_pred_mean_rel)

#BACT PRED
sample_sums(physeq_bactpred)
physeq_bactpred_rel <- transform_sample_counts(physeq_bactpred, function(x) x/sum(x)*100)
sample_sums(physeq_bactpred_rel)

#reating rel counts per depth
physeq_bactpred_mean <- merge_samples(physeq_bactpred, "depth")
sample_sums(physeq_bactpred_mean)
physeq_bactpred_mean_rel <- transform_sample_counts(physeq_bactpred_mean, function(x) x/sum(x)*100)
sample_sums(physeq_bactpred_mean_rel)

#PROTOZOA
sample_sums(physeq_protozoa)
physeq_protozoa_rel <- transform_sample_counts(physeq_protozoa, function(x) x/sum(x)*100)
sample_sums(physeq_protozoa_rel)

#reating rel counts per depth
physeq_protozoa_mean <- merge_samples(physeq_protozoa, "depth")
sample_sums(physeq_protozoa_mean)
physeq_protozoa_mean_rel <- transform_sample_counts(physeq_protozoa_mean, function(x) x/sum(x)*100)
sample_sums(physeq_protozoa_mean_rel)


###----Richness bac----
p1 <- plot_richness(physeq_bac, x="depth", color = "layer", shape = "age", measures=c("Observed","Shannon"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "right") +
  theme(strip.background = element_rect(fill="white" ))

p1

# write the calculations in a file
write.table(p1$data, file='diversity_bac.tsv', quote=FALSE, sep='\t')
dev.print(pdf, 'richness_bac.pdf')

###----Richness euk----
p1 <- plot_richness(physeq_euk, x="depth", color = "layer", shape = "age",measures=c("Observed","Shannon"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "right") +
  theme(strip.background = element_rect(fill="white" ))

p1

# write the calculations in a file
write.table(p1$data, file='diversity_euk.tsv', quote=FALSE, sep='\t')
dev.print(pdf, 'richness_euk.pdf')

###----Richness euk Predators----
p1 <- plot_richness(physeq_protozoa, x="depth", color = "layer", shape = "age", measures=c("Observed","Shannon"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "right") +
  theme(strip.background = element_rect(fill="white" ))

p1
dev.print(pdf, 'richness_prot.pdf')
write.table(p1$data, file='diversity_protozoa.tsv', quote=FALSE, sep='\t')

###----Richness bact pred----

p1 <- plot_richness(physeq_bactpred, x="depth", color = "layer", shape = "age", measures=c("Observed","Shannon","Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "right") +
  theme(strip.background = element_rect(fill="white" ))

p1
dev.print(pdf, 'richness_bacpred.pdf')


# write the calculations in a file
write.table(p1$data, file='diversity_bactPRED.tsv', quote=FALSE, sep='\t')



##########ORDINATION
###----Ordination_bac----

physeq_nmds_bac <- ordinate(physeq_bac, method = "NMDS", distance = "bray")

# Shepard plot - all
vegan::stressplot(physeq_nmds_bac)

# NMDS plot - all
p1 = plot_ordination(physeq_bac, physeq_nmds_bac, type =" samples", color = "layer", shape = "age", title = "NMDS Bateria total counts")+
  geom_point(size=3)+
  theme_bw()+
  # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)
dev.print(pdf, 'NMDS_bac.pdf')



###----Ordination_euk----

physeq_nmds_euk <- ordinate(physeq_euk, method = "NMDS", distance = "bray")

# Shepard plot - all
vegan::stressplot(physeq_nmds_euk)

# NMDS plot - all
p1 = plot_ordination(physeq_euk, physeq_nmds_euk, type =" samples", color = "layer", shape = "age", title = "NMDS EUK total counts")+
  geom_point(size=3)+
  theme_bw()+
  # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)
dev.print(pdf, 'NMDS_euk.pdf')


###----Ordination_protozoa----

physeq_nmds_protozoa <- ordinate(physeq_protozoa, method = "NMDS", distance = "bray")

# Shepard plot - all
vegan::stressplot(physeq_nmds_protozoa)

# NMDS plot - all
p1 = plot_ordination(physeq_protozoa, physeq_nmds_protozoa, type ="taxa", color = "Kingdom", title = "NMDS Protozoa total counts")+
  geom_point(size=3)+
  theme_bw()+
  # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)
dev.print(pdf, 'NMDS_prot.pdf')

###----Ordination_bactpred----

physeq_nmds_bactpred <- ordinate(physeq_bactpred, method = "NMDS", distance = "bray")

# Shepard plot - all
vegan::stressplot(physeq_nmds_bactpred)

# NMDS plot - all
p1 = plot_ordination(physeq_bactpred, physeq_nmds_bactpred, type ="taxa", color = "Order", title = "NMDS Bacterial Predators total counts")+
  geom_point(size=3)+
  theme_bw()+
  # scale_color_manual(values=c("#2f85fe", "#e05436", "#009453")) +
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill="white" ))
print(p1)

dev.print(pdf, 'NMDS_bacpred.pdf')




###----ANOVA bac----

library("microbiome"); packageVersion("microbiome")
library("vegan"); packageVersion("vegan")

otu1 <- microbiome::abundances(physeq_bac)
meta1 <- microbiome::meta(physeq_bac)


permanova <- vegan::adonis(t(otu1)~layer,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

# write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_bac.tsv', quote=FALSE, sep='\t')

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$age, type="centroid")
mod
mod$age <- meta$age
TukeyHSD(mod)


mod <- vegan::betadisper(dist, meta1$layer, type="centroid")
mod
mod$layer <- meta$layer
TukeyHSD(mod)


mod <- vegan::betadisper(dist, meta1$depth, type="centroid")
mod
mod$depth <- meta$depth
TukeyHSD(mod)


# plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA BAC", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

# dev.off()


###----ANOVA euk----
otu1 <- microbiome::abundances(physeq_euk)
meta1 <- microbiome::meta(physeq_euk)


permanova <- vegan::adonis(t(otu1)~layer,
                           data = meta1, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab))

# write.table(as.data.frame(permanova$aov.tab), file='permanova_bray_overall_sar.tsv', quote=FALSE, sep='\t')

dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$age, type="centroid")
mod
mod$age <- meta$age
TukeyHSD(mod)


mod <- vegan::betadisper(dist, meta1$layer, type="centroid")
mod
mod$layer <- meta$layer
TukeyHSD(mod)


mod <- vegan::betadisper(dist, meta1$depth, type="centroid")
mod
mod$depth <- meta$depth
TukeyHSD(mod)


# plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA SAR", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

# dev.off()



###----ANOVA protozoa----
otu1 <- microbiome::abundances(physeq_protozoa)
meta1 <- microbiome::meta(physeq_protozoa)


dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$age, type="centroid")
mod
mod$age <- meta$age
TukeyHSD(mod)

mod <- vegan::betadisper(dist, meta1$layer, type="centroid")
mod
mod$layer <- meta$layer
TukeyHSD(mod)

mod <- vegan::betadisper(dist, meta1$depth, type="centroid")
mod
mod$depth <- meta$depth
TukeyHSD(mod)


# plot(mod,  hull=FALSE, ellipse=TRUE, main = "PCoA", sub=NULL, col=c("#2f85fe", "#e05436", "#009453"), cex=2, lwd=1) #+
#boxplot(mod$distances ~ mod$group, main= "Distance to Centroid", xlab="kit", ylab="Distance", col= c("#2f85fe", "#e05436", "#009453"))

# dev.off()


###----ANOVA bactpred----
otu1 <- microbiome::abundances(physeq_bactpred)
meta1 <- microbiome::meta(physeq_bactpred)


dist <- vegan::vegdist(t(otu1), method="bray")
mod <- vegan::betadisper(dist, meta1$layer, type="centroid")
mod
mod$layer <- meta$layer
TukeyHSD(mod)

mod <- vegan::betadisper(dist, meta1$age, type="centroid")
mod
mod$age <- meta$age
TukeyHSD(mod)

mod <- vegan::betadisper(dist, meta1$depth, type="centroid")
mod
mod$depth <- meta$depth
TukeyHSD(mod)

# dev.off()




###----Stacked Barplots bac----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_bac_mean <- merge_samples(physeq_bac, "depth")

#After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_bac_mean <- transform_sample_counts(physeq_bac_mean, function(x) x/sum(x)*100)

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

#After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_bac_mean_phylum <- transform_sample_counts(physeq_bac_mean_phylum, function(x) x/sum(x)*100)
#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to "< 1%".

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
  # scale_fill_manual(values = c("#bd6b8f",
  #                              "#6db543",
  #                              "#7661cd",
  #                              "#c0b047",
  #                              "#c35abc",
  #                              "#60bf8b",
  #                              "#d13f73",
  #                              "#3e8149",
  #                              "#ca5340",
  #                              "#45b0cf",
  #                              "#cc8444",
#                              "#7882c9",
#                              "#7a7732",
#                              "#6db543")) +
labs(y= "Relative abundance [%]",
     fill= "Phyla") +
  theme_bw()
dev.print(pdf, 'barplot_bac.pdf')

###----Stacked Barplots euk----

#For stacked barplots, we will work with means of replicates using the merge_samples function. The default of this function is merge_samples(x, group, fun= mean).

physeq_bac_mean <- merge_samples(physeq_euk, "depth")

#After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_bac_mean <- transform_sample_counts(physeq_bac_mean, function(x) x/sum(x)*100)

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

#After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_bac_mean_phylum <- transform_sample_counts(physeq_bac_mean_phylum, function(x) x/sum(x)*100)
#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to "< 1%".

#transform phyloseq object to a data frame (DF)
physeq_bac_mean_phylumDF<- psmelt(physeq_bac_mean_phylum)

#inspect the dataframe
str(physeq_bac_mean_phylumDF)

#make the phyla characters, not factors
physeq_bac_mean_phylumDF$Phylum <- as.character(physeq_bac_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_bac_mean_phylumDF <- physeq_bac_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 5, "< 5%"))

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
  # scale_fill_manual(values = c("#bd6b8f",
  #                              "#6db543",
  #                              "#7661cd",
  #                              "#c0b047",
  #                              "#c35abc",
  #                              "#60bf8b",
  #                              "#d13f73",
  #                              "#3e8149",
  #                              "#ca5340",
  #                              "#45b0cf",
  #                              "#cc8444",
#                              "#7882c9",
#                              "#7a7732",
#                              "#6db543")) +
labs(y= "Relative abundance [%]",
     fill= "Phyla") +
  theme_bw()
dev.print(pdf, 'barplot_euk.pdf')

#more colours and combinations can be autogenerated here: https://medialab.github.io/iwanthue/
#this part was adapted from: https://mvuko.github.io/meta_phyloseq/


###----Stacked Barplots Protozoa----


#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level.
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_protozoa_mean)

#agglomeration on the Phylum level
physeq_protozoa_mean_phylum <- tax_glom(physeq_protozoa_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_protozoa_mean
physeq_protozoa_mean_phylum

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_protozoa_mean_phylum)

# #After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_protozoa_mean_phylum_rel <- transform_sample_counts(physeq_protozoa_mean_phylum, function(x) x/sum(x)*100)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to "< 1%".

#transform phyloseq object to a data frame (DF)
physeq_protozoa_mean_phylumDF<- psmelt(physeq_protozoa_mean_phylum_rel)

#inspect the dataframe
str(physeq_protozoa_mean_phylumDF)

#make the phyla characters, not factors
physeq_protozoa_mean_phylumDF$Phylum <- as.character(physeq_protozoa_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
physeq_protozoa_mean_phylumDF <- physeq_protozoa_mean_phylumDF %>%
  mutate(Phylum2 = replace(Phylum, Abundance < 1, "<1%"))

#check all phyla names
unique(physeq_protozoa_mean_phylumDF$Phylum2)

#there are some reads that were assigned only to the kingdom level,
# i.e. NA on the phylum level, so we will rename them
physeq_protozoa_mean_phylumDF <- physeq_protozoa_mean_phylumDF %>%
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned Protozoa"))

#reorder the phyla so that they are stacked according to abundance
physeq_protozoa_mean_phylumDF$Phylum2 <- reorder(physeq_protozoa_mean_phylumDF$Phylum2,
                                                 physeq_protozoa_mean_phylumDF$Abundance)

ggplot(physeq_protozoa_mean_phylumDF, aes(Sample, Abundance, fill=Phylum2)) +
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = c("#bd6b8f",
  #                              "#6db543",
  #                              "#7661cd",
  #                              "#c0b047",
  #                              "#c35abc",
  #                              "#60bf8b",
  #                              "#d13f73",
  #                              "#3e8149",
  #                              "#ca5340",
  #                              "#45b0cf",
  #                              "#cc8444",
#                              "#7882c9",
#                              "#7a7732")) +
labs(y= "Relative abundance [%]",
     fill= "Phyla") +
  theme_bw()
dev.print(pdf, 'barplot_prot.pdf')


###----Stacked Barplots bactpred----


#Stacked barplots can show taxonomic levels of choice. We will create a plot that shows the phylum level.
#First, we will agglomerate data to the phylum level using the tax_glom function, and then plot the agglomerated data.
#We can also agglomerate the data to other levels, depending on the object we are analyzing

#check taxonomy level names
rank_names(physeq_bactpred_mean)

#agglomeration on the Phylum level
physeq_bactpred_mean_order <- tax_glom(physeq_bactpred_mean, taxrank = "Order")

#check the number of taxa in the whole phyloseq object and then how many phyla they were assigned to
physeq_bactpred_mean
physeq_bactpred_mean_order

#check sample sums to make sure nothing was deleted due to NAs (should sum up to 1)
sample_sums(physeq_bactpred_mean_order)

# #After merging, we need to dived the numbers to the number of replicates, in this case "5"
physeq_bactpred_mean_order_rel <- transform_sample_counts(physeq_bactpred_mean_order, function(x) x/sum(x)*100)

#Since we get too many phyla to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. To do this, we will use the power of tidyverse again. First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to "< 1%".

#transform phyloseq object to a data frame (DF)
physeq_bactpred_mean_orderDF<- psmelt(physeq_bactpred_mean_order_rel)

#inspect the dataframe
str(physeq_bactpred_mean_orderDF)

#make the phyla characters, not factors
physeq_bactpred_mean_orderDF$Order <- as.character(physeq_bactpred_mean_orderDF$Order)

#add new column with renamed low abundant taxa
physeq_bactpred_mean_orderDF <- physeq_bactpred_mean_orderDF %>%
  mutate(Order2 = replace(Order, Abundance < 1, "<1%"))

#check all phyla names
unique(physeq_bactpred_mean_orderDF$Order2)

#there are some reads that were assigned only to the kingdom level,
# i.e. NA on the family level, so we will rename them
physeq_bactpred_mean_orderDF <- physeq_bactpred_mean_orderDF %>%
  mutate(Order2 = replace(Order2, Order2 == "NA", "unassigned bactpred"))

#reorder the phyla so that they are stacked according to abundance
physeq_bactpred_mean_orderDF$Order2 <- reorder(physeq_bactpred_mean_orderDF$Order2,
                                               physeq_bactpred_mean_orderDF$Abundance)

ggplot(physeq_bactpred_mean_orderDF, aes(Sample, Abundance, fill=Order2)) +
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = c("#bd6b8f",
  #                              # "#6db543",
  #                              # "#7661cd",
  #                              # "#c0b047",
  #                              # "#c35abc",
  #                              # "#60bf8b",
  #                              # "#d13f73",
  #                              # "#3e8149",
  #                              # "#ca5340",
  #                              # "#45b0cf",
  #                              # "#cc8444",
#                              # "#7882c9",
#                              "#7a7732")) +
labs(y= "Relative abundance [%]",
     fill= "Order") +
  theme_bw()
dev.print(pdf, 'barplot_bacpred.pdf')