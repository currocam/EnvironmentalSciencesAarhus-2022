###########################
# New version
# Code written by Thanassis Zervas 
# az@envs.au.dk
# Creative commons, 2021. 

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

# Agglomerate at certain level.
# This maneuver combines all duplicates at species level, OTU abundances are summed.
physeq_p <- tax_glom(physeq_rel, taxrank = "Phylum", NArm = FALSE)
physeq_p
head(tax_table(physeq_p))
sample_sums(physeq_p)

physeq_c <- tax_glom(physeq_rel, taxrank = "Class", NArm = FALSE)
physeq_c
head(tax_table(physeq_c))
sample_sums(physeq_c)

physeq_g <- tax_glom(physeq_rel, taxrank = "Genus", NArm = FALSE)
physeq_g
head(tax_table(physeq_g))
sample_sums(physeq_g)