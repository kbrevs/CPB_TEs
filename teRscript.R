### transposable elements ###
setwd("~/Documents/CPB")
library(tidyverse)

#### genes AND flanking regions ####
#open the intersection file
flanks <- read_delim("~/Documents/CPB/TEs/Feb2016/RM_custom/alpaths_genes_plus_flanks_intersects.gff", "\t", col_names = F)
#View(head(flanks))

#open the gene names file
gene_gff <- as.data.frame(read_csv("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/OGSv1.0_IDmapping.csv"))
#View(head(gene_gff))
#subset to get rid of stuff I think we dont need
gene_gff_rel <- dplyr::select(gene_gff, OGSv1.0_gene_ID, OGSv1.0_gene_name,OGSv1.0_transcript_name)
#View(head(gene_gff_rel))


#merge them somehow
flanks$OGSv1.0_gene_ID <- substr(flanks$X9,4,13)

intersections <- full_join(flanks, gene_gff_rel, by="OGSv1.0_gene_ID")
#View(head(intersections))

#now we need the actual TE names..... 
library(stringr)
te_names <- read_delim("~/Documents/CPB/TEs/TEnames.txt", "#", col_names=F)
names(te_names) <- c("TE_ID", "TE_type")
#View(te_names)

#OKAY, we need unique identifiers for each TE name, so - Dnd MOnsters I suppose
#http://miroz.com.hr/random/monsters.csv
dnd <- read_csv("/home/beetle/Downloads/srd35-db-v1.3/monsternames.csv")
#removing the dragons
dnd <- filter(dnd, !grepl("Dragon",Unique_Name))
dnd <- filter(dnd, !grepl("Elemental",Unique_Name))
dnd <- filter(dnd, !grepl("Were",Unique_Name))

# and give each TE its own name!
te_names <- bind_cols(te_names, dnd[1:nrow(te_names),])

#format TE names for merginf
intersections$TE_ID <- gsub('"', "", gsub(":","",str_extract(intersections$X18, ':(.*)"')))
#View(intersections)

intersections <- full_join(intersections, te_names, by="TE_ID")


#just genes

#ok, that looks squared away

#filter out all of the ones where the gene name is just the LDEC000000 (this is not necessary anymore now that we have the GO terms)
#str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE

#got rid of anything without a good name
#intersections_filtered <- filter(intersections, str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE)

#still want this tho
intersections_filtered <- filter(intersections, X3=="gene")

# #### STOP HERE PROBABLY ####
# 
# intersections_filtered <- filter(intersections_filtered, OGSv1.0_gene_name)
# 
# length(unique(intersections_filtered$OGSv1.0_gene_name))
# #24708 uniques.... pre filtered
# intersections_filtered <- filter(intersections_filtered, !grepl("receptor",OGSv1.0_gene_name))
# #number of TE intersections per gene
# genes_by_TE_count <- intersections_filtered %>% group_by(OGSv1.0_gene_name) %>% tally()
# View(genes_by_TE_count)
# 
# Gr48 <- filter(intersections_filtered, OGSv1.0_gene_name=="Gustatory receptor 48")
# TEs_in_Gr48 <- Gr48 %>% group_by(TE_name) %>% tally()
# 
# CH23 <- filter(intersections_filtered, OGSv1.0_gene_name=="Cadherin 23")
# TEs_in_CH23 <- CH23 %>% group_by(TE_name) %>% tally()
# 
# ACTH <- filter(intersections_filtered, OGSv1.0_gene_name=="Acetylcholine esterase 1")
# TEs_in_ACTH <- ACTH %>% group_by(TE_name) %>% tally()
# 
# 
# 
# #can we get GO terms at all?
# library(biomaRt)
# ensembl <- useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
# goterms <- getBM(attributes=c('external_gene_name','description','go_id', 'name_1006'), values= c('ACTB', 'TNF'), mart= ensembl)
# View(goterms)
# filter(goterms, grepl("Acetylcholinesterase", description))
# grep(goterms, "Acetylcholinesterase")
# 
# 
# 
# 
# 
# 
# #### JUST flanking regions ####
# #open the intersection file
# flanks <- read_delim("~/Documents/CPB/TEs/Feb2016/RM_custom/alpaths_flank_intersections.gff", "\t", col_names = F)
# View(head(flanks))
# 
# #open the gene names file
# gene_gff <- read_csv("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/OGSv1.0_IDmapping.csv")
# View(head(gene_gff))
# #subset to get rid of stuff I think we dont need
# gene_gff_rel <- select(gene_gff, OGSv1.0_gene_ID, OGSv1.0_gene_name,OGSv1.0_transcript_name)
# View(head(gene_gff_rel))
# 
# 
# #merge them somehow
# flanks$OGSv1.0_gene_ID <- substr(flanks$X9,4,13)
# 
# intersections <- full_join(flanks, gene_gff_rel, by="OGSv1.0_gene_ID")
# View(head(intersections))
# 
# #now we need the actual TE names..... 
# library(stringr)
# te_names <- read_delim("~/Documents/CPB/TEs/TEnames.txt", "#", col_names=F)
# names(te_names) <- c("TE_ID", "TE_name")
# View(te_names)
# intersections$TE_ID <- gsub('"', "", gsub(":","",str_extract(intersections$X18, ':(.*)"')))
# View(intersections)
# 
# intersections <- full_join(intersections, te_names, by="TE_ID")
# 
# 
# #just genes
# 
# #ok, that looks squared away
# 
# #filter out all of the ones where the gene name is just the LDEC000000
# filter(intersections)
# str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE
# 
# #got rid of anything without a good name
# intersections_filtered <- filter(intersections, str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE)
# 
# 
# intersections_filtered <- filter(intersections_filtered, X3=="gene")
# 
# 
# 
# length(unique(intersections_filtered$OGSv1.0_gene_name))
# #24708 uniques.... pre filtered
# 
# #number of TE intersections per gene
# genes_by_TE_count <- intersections_filtered %>% group_by(OGSv1.0_gene_name) %>% tally()
# View(genes_by_TE_count)
# 
# 
# 
# 
# #### and just inside genes ####
# #open the intersection file
# flanks <- read_delim("~/Documents/CPB/TEs/Feb2016/RM_custom/alpaths_annotations_intersections.gff", "\t", col_names = F)
# View(head(flanks))
# 
# #open the gene names file
# gene_gff <- read_csv("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/OGSv1.0_IDmapping.csv")
# View(head(gene_gff))
# #subset to get rid of stuff I think we dont need
# gene_gff_rel <- select(gene_gff, OGSv1.0_gene_ID, OGSv1.0_gene_name,OGSv1.0_transcript_name)
# View(head(gene_gff_rel))
# 
# 
# #merge them somehow
# flanks$OGSv1.0_gene_ID <- substr(flanks$X9,4,13)
# 
# intersections <- full_join(flanks, gene_gff_rel, by="OGSv1.0_gene_ID")
# View(head(intersections))
# 
# #now we need the actual TE names..... 
# library(stringr)
# te_names <- read_delim("~/Documents/CPB/TEs/TEnames.txt", "#", col_names=F)
# names(te_names) <- c("TE_ID", "TE_name")
# View(te_names)
# intersections$TE_ID <- gsub('"', "", gsub(":","",str_extract(intersections$X18, ':(.*)"')))
# View(intersections)
# 
# intersections <- full_join(intersections, te_names, by="TE_ID")
# 
# 
# #just genes
# 
# #ok, that looks squared away
# 
# #filter out all of the ones where the gene name is just the LDEC000000
# filter(intersections)
# str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE
# 
# #got rid of anything without a good name
# intersections_filtered <- filter(intersections, str_detect(intersections$OGSv1.0_gene_name, "[:digit:]{6}") == FALSE)
# 
# 
# intersections_filtered <- filter(intersections_filtered, X3=="gene")
# 
# 
# 
# length(unique(intersections_filtered$OGSv1.0_gene_name))
# #24708 uniques.... pre filtered
# 
# #number of TE intersections per gene
# genes_by_TE_count <- intersections_filtered %>% group_by(OGSv1.0_gene_name) %>% tally()
# View(genes_by_TE_count)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### try to draw gene pictures ####
# source("https://bioconductor.org/biocLite.R")
# biocLite("gviz")
# browseVignettes("rtracklayer")
# 
# #### open TEs ####
# 
# 
# # Gviz 
# library(rtracklayer)
# library(GenomicFeatures)
# library(Gviz)
# 
# gff<-readGFF("gffplottest.gff")
# gff$TE_ID <- gsub('"', "", gsub(":","",str_extract(gff$Target, ':(.*)"')))
# gff <- left_join(gff, te_names, by="TE_ID")
# colnames(gff)[11] <- "symbol"
# 
# te_atrack <- AnnotationTrack(gff, name="TEs")
# plotTracks(te_atrack)
# 
# te_gtrack <- GenomeAxisTrack()
# plotTracks(list(te_gtrack, te_atrack))
# 
# te_grtrack <- GeneRegionTrack(gff, name="TEs", 
#                               transcriptAnnotation = "symbol",
#                               background.title="brown", 
#                               background.panel="#FFFEDB",)
# 
# plotTracks(list(te_gtrack,te_grtrack))
# 
# #try a gene
# 
# # scaffold1 ###
# 
# LDEC000001intersections <- filter(intersections, OGSv1.0_gene_ID=="LDEC000001")
# genes_GFF <- as.data.frame(readGFF("lepdec_OGSv1.0.gff"))
# TE_GFF <- readGFF("Ldec.genome.10062013.fa.out.gff")
# #names
# TE_GFF$TE_ID <- gsub('"', "", gsub(":","",str_extract(TE_GFF$Target, ':(.*)"')))
# TE_GFF <- left_join(TE_GFF, te_names, by="TE_ID")
# colnames(TE_GFF)[11] <- "symbol"
# 
# 
# scaffold1_genes <- filter(genes_GFF, seqid=="Scaffold1")
# scaffold1_TEs <- filter(TE_GFF, seqid=="Scaffold1")
# 
# atrack_scaffold1 <- GeneRegionTrack(scaffold1_genes, name="Gene Model",
#                                     transcriptAnnotation="ID")
# gtrack_scaffold1 <- GenomeAxisTrack()
# grtrack_scaffold1 <- GeneRegionTrack(scaffold1_TEs, name="TEs", 
#                               transcriptAnnotation = "symbol",
#                               background.title="brown", 
#                               background.panel="#FFFEDB")
# 
# 
# plotTracks(list(gtrack_scaffold1,atrack_scaffold1,grtrack_scaffold1), from = 10000, to=20000)
# 
# 
# #### Scaffold 27 ###
# 
# scaffold27_genes <- filter(genes_GFF, seqid=="Scaffold27")
# colnames(scaffold27_genes)[10] <- "symbol"
# colnames(scaffold27_genes)[19] <- "orcs"
# scaffold27_exons <- filter(scaffold27_genes, type=="exon")
# scaffold27_TEs <- filter(TE_GFF, seqid=="Scaffold27")
# 
# #ld_annotations_GFF <- as.data.frame(readGFF("Ldec.genome.10062013.annotations.gff"))
# #old_scaffold27_genes <- filter(old_annotations_GFF, seqid=="Scaffold27")
# #colnames(old_scaffold27_genes)[9]<-"symbol"
# 
# atrack_scaffold27 <- GeneRegionTrack(scaffold27_exons, name="Exons",
#                                     transcriptAnnotation="symbol")
# atrack2_scaffold27 <- GeneRegionTrack(filter(filter(scaffold27_genes, symbol=="Acetylcholine esterase 1"), type=="gene"), name="Gene",
#                                       transcriptAnnotation="symbol",fill="lightblue")
# #oldtrack <- GeneRegionTrack(old_scaffold27_genes, name="Old", transcriptAnnotation="symbol")
# gtrack_scaffold27 <- GenomeAxisTrack()
# grtrack_scaffold27 <- GeneRegionTrack(scaffold27_TEs, name="TEs", 
#                                      transcriptAnnotation = "symbol",
#                                      background.title="brown", 
#                                      background.panel="#FFFEDB")
# 
# 
# 
# plotTracks(list(gtrack_scaffold27,atrack_scaffold27,atrack2_scaffold27,grtrack_scaffold27), 
#            from = 2157368, to=2174350,
#            main="title")
# displayPars(atrack2_scaffold27)
# 
# 
# 
# 
# 
# #Gustatory receptor 48
# 
# 
# GR48 <- filter(genes_GFF, Name=="Gustatory receptor 48")
# 
# scaffold122_genes <- filter(genes_GFF, seqid=="Scaffold122")
# colnames(scaffold122_genes)[10] <- "symbol"
# colnames(scaffold122_genes)[19] <- "orcs"
# scaffold122_exons <- filter(scaffold122_genes, type=="exon")
# scaffold122_TEs <- filter(TE_GFF, seqid=="Scaffold122")
# 
# #ld_annotations_GFF <- as.data.frame(readGFF("Ldec.genome.10062013.annotations.gff"))
# #old_scaffold122_genes <- filter(old_annotations_GFF, seqid=="Scaffold122")
# #colnames(old_scaffold122_genes)[9]<-"symbol"
# 
# #atrack_scaffold122 <- GeneRegionTrack(scaffold122_exons, name="Exons",
#                                       transcriptAnnotation="symbol")
# atrack2_scaffold122 <- GeneRegionTrack(scaffold122_genes, name="Gene",
#                                        transcriptAnnotation="symbol")
# #oldtrack <- GeneRegionTrack(old_scaffold122_genes, name="Old", transcriptAnnotation="symbol")
# gtrack_scaffold122 <- GenomeAxisTrack()
# grtrack_scaffold122 <- GeneRegionTrack(scaffold122_TEs, name="TEs", 
#                                        transcriptAnnotation = "symbol",
#                                        background.title="brown", 
#                                        background.panel="#FFFEDB")
# 
# 
# 
# plotTracks(list(gtrack_scaffold122,atrack2_scaffold122,grtrack_scaffold122), 
#            from = 722231, to=813264,
#            main="title")
# displayPars(atrack2_scaffold122)
# 
# 
# #### ABC transporter subfamily G member ###
# filter(genes_GFF, Name=="ABC transporter subfamily G member")
# scaffold1041_genes <- filter(genes_GFF, seqid=="Scaffold1041")
# colnames(scaffold1041_genes)[10] <- "symbol"
# colnames(scaffold1041_genes)[19] <- "orcs"
# scaffold1041_exons <- filter(scaffold1041_genes, type=="exon")
# scaffold1041_TEs <- filter(TE_GFF, seqid=="Scaffold1041")
# 
# #ld_annotations_GFF <- as.data.frame(readGFF("Ldec.genome.10062013.annotations.gff"))
# #old_scaffold1041_genes <- filter(old_annotations_GFF, seqid=="Scaffold1041")
# #colnames(old_scaffold1041_genes)[9]<-"symbol"
# 
# atrack_scaffold1041 <- GeneRegionTrack(scaffold1041_exons, name="Exons",
#                                        transcriptAnnotation="symbol")
# atrack2_scaffold1041 <- GeneRegionTrack(scaffold1041_genes, name="Gene",
#                                         transcriptAnnotation="symbol")
# #oldtrack <- GeneRegionTrack(old_scaffold1041_genes, name="Old", transcriptAnnotation="symbol")
# gtrack_scaffold1041 <- GenomeAxisTrack()
# grtrack_scaffold1041 <- GeneRegionTrack(scaffold1041_TEs, name="TEs", 
#                                         transcriptAnnotation = "symbol",
#                                         background.title="brown", 
#                                         background.panel="#FFFEDB")
# 
# 
# 
# plotTracks(list(gtrack_scaffold1041,atrack_scaffold1041,atrack2_scaffold1041,grtrack_scaffold1041), 
#            from = 100775, to=140851,
#            main="title")
# displayPars(atrack2_scaffold1041)
# 
# 
# #### is the bioconductor GO database useful? ####
# source("https://bioconductor.org/biocLite.R")
# biocLite("GO.db")
# library(GO.db)
# columns(GO.db)
# keys(GO.db) #this is all the IDs
# vals <- select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
# head(vals)
# vals2 <- select(GO.db, keys("GO:0000347"))
# head(vals2)
# get("GO:0007005", GOBPANCESTOR)
# #ok, but i want to use select...

#um how about just a csv of go terms. csv ####
goterms <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/go.csv", "\t",
                      col_names = c("term_ID",
                                    "term_name",
                                    "namespace",
                                    "alt_id",
                                    "definition",
                                    "comment",
                                    "synonym(s)",
                                    "xref(s)",
                                    "is_a(parents)",
                                    "disjoint_from",
                                    "relationship(parents)"))
#View(goterms)
# Remove obsolete terms
goterms <- filter(goterms, !grepl("obsolete",term_name))


#um how about just a csv of go terms - GO slim ####
goslim <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/goslim.csv", "\t",
                     col_names = c("term_ID",
                                   "term_name",
                                   "namespace",
                                   "alt_id",
                                   "definition",
                                   "comment",
                                   "synonym(s)",
                                   "xref(s)",
                                   "is_a(parents)",
                                   "disjoint_from",
                                   "relationship(parents)"))
View(goslim)

# now we need the correspondence matching file to convert go to goslim ####
goslim.gaf <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/subset.annotations.mapped.gaf", "\t",
                         col_names = as.character(seq(1:17)))
View(goslim.gaf)                    
length(unique(goslim.gaf$"5")) #141, there were 150 to begin with
#seem to be some missing, but maybe for stuff CPB doesnt do, lets check
#open the full list
slimterms <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/slimterms.csv", "\n", col_names="term")
slimterms_included <- as.data.frame(unique(goslim.gaf$"5"))
colnames(slimterms_included) <-"term"
#who is missing
anti_join(slimterms, slimterms_included, "term")
# 1 GO:0140014
# 2 GO:0043473
# 3 GO:0034330
# 4 GO:0030555
# 5 GO:0030533
# 6 GO:0021700
# 7 GO:0007568
# 8 GO:0006412
# 9 GO:0000278
#do not know why those are not showing up 

#### make the two column thing: GO term and GO slim
full.gaf <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/genes_orfs_Pfam-A_GO_GOterm_mapping.gaf", "\t",
                       col_names = as.character(seq(1:15)))
#full.gaf$"16" <- seq(1:nrow(full.gaf))
full.gaf$"2" <- paste(as.character(full.gaf$"2"),as.character(full.gaf$"5"), sep="-")
#I think I need to make a unique identifier, merge, then remove
write_delim(full.gaf, "~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/genes_orfs_Pfam-A_GO_GOterm_mapping_UNIQUE.gaf", "\t")

#import the new files, merge and purge!
# fullterms <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/genes_orfs_Pfam-A_GO_GOterm_mapping_UNIQUE.gaf", "\t",
#                        col_names = as.character(seq(1:17)))
# slimeterms <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/subset.annotations.mapped_UNIQUE.gaf", "\t",
#                                        col_names = as.character(seq(1:17)))
# length(unique(fullterms$"5"))
# length(unique(slimeterms$"5"))
# fma <- dplyr::select(fullterms, 2, 5)
# sma <- dplyr::select(slimeterms, 2, 5)
# colnames(fma) <- c("name", "GOterm")
# colnames(sma) <- c("name", "slimGOterm")
# matches <- full_join(fma, sma, by="name")
# View(matches)

# new plan, time to chop
slimeterms <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/annotations.mapped_UNIQUE.gaf", "\t",
                         col_names = as.character(seq(1:17)))
slimeterms <- dplyr::select(slimeterms, 2, 5)
#View(slimeterms)
slimeterms$old <- substr(slimeterms$"2", 31, 40)
slimeterms <- dplyr::select(slimeterms, 5, "old")
colnames(slimeterms) <- c("GOaway", "GOSLIM", "GOfine")
slimeterms <- dplyr::select(slimeterms, 2, 3)
View(slimeterms)
unique(slimeterms$GOfine)
#fix missing Gs
slimeterms$GOfine <- gsub("GG", "G", gsub("O:", "GO:", slimeterms$GOfine))
slimeterms <- unique(slimeterms)

#### hmm2go ####
hmm2go <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/genes_orfs_Pfam-A_GO.tsv", delim="\t",
                     col_names = c("Gene_ID", "pfam", "donno", "GO_activity", "GO", "desc"))
#View(head(hmm2go))

hmm2go2col <- read_delim("~/Documents/CPB/reference_genomes/lepdec_OGSv1.0/hmmer2go/genes_orfs_Pfam-A_GO_GOterm_mapping.tsv", delim="\t",
                         col_names = c("OGSv1.0_gene_ID", "GOs"))
#View(head(hmm2go2col))
#truncate names (what are those extra numbers?)
hmm2go2col$OGSv1.0_gene_ID <- substr(hmm2go2col$OGSv1.0_gene_ID, 0, 10)

View(head(intersections_filtered))

hmm2tes <- left_join(intersections_filtered, hmm2go2col, by="OGSv1.0_gene_ID")
# we need to merge the GO terms with the TE intersections - to see which GO terms are associated with TEs the most.
#View(hmm2tes)

#how many TEs are in genes with no GO term?
# what is the significance? not a gene... or just no hits?
sum(is.na(hmm2tes$GOs)) 
nrow(hmm2tes)

#there are many duplicated GO terms in the GOs column...
SplitFunction <- function(x) {
  b <- unlist(strsplit(x, '[,]'))
  c <- b[!duplicated(b)]
  return(paste(c, collapse=","))
}

hmm2tes$GOsDeDup <- sapply(hmm2tes$GOs, SplitFunction)

#look at gos!  I think this is outdated ####
# GOs <- dplyr::select(hmm2tes,GOsDeDup)
# GOs <- as.data.frame(GOs) %>% separate(GOsDeDup, into = paste("V", 1:10, sep = "_"), sep=",")
# GOs <- as.data.frame(stack(GOs)$values)
# head(GOs, n=30)
# colnames(GOs) <- "GO"
# GOs <- na.omit(GOs)

#goplot <- ggplot(GOs, aes(GO))
#goplot + geom_bar()

# gocount <- GOs %>% group_by(GO) %>% tally()

#merging time
#hmm2gosub <- unique(dplyr::select(hmm2go, GO, GO_activity, desc))
#View(hmm2gosub)
#gocountWnames <- left_join(gocount, hmm2gosub, by="GO")
#View(gocountWnames)

# goplot <- ggplot(head(nodescs), aes(GO_activity,n))
# goplot + geom_bar(stat="identity")
# 
# 
# nodescs <- unique(dplyr::select(gocountWnames, GO, n, GO_activity))
# View(nodescs)
# 
# # and by TE ####
# tecount <- hmm2tes %>% group_by(Unique_Name) %>% tally()
# tecount <- na.omit(tecount)
# teplot <- ggplot(tecount, aes(Unique_Name, n))
# teplot + geom_bar(stat="identity")


# #Okay, what GO terms is that DestructionBeetle close to? ####
# devastationBeetle <- filter(hmm2tes, Unique_Name=="DevastationBeetle")
# View(devastationBeetle)
# GOs <- dplyr::select(devastationBeetle,GOsDeDup)
# GOs <- as.data.frame(GOs) %>% separate(GOsDeDup, into = paste("V", 1:10, sep = "_"), sep=",")
# GOs <- as.data.frame(stack(GOs)$values)
# colnames(GOs) <- "GO"
# GOs <- na.omit(GOs)
# gocount <- GOs %>% group_by(GO) %>% tally()
# View(gocount)
# 
# 
# dreamlarva <- filter(hmm2tes, Unique_Name=="DreamLarva")
# View(dreamlarva)
# 



#### Experimental Tidying ####
# Ok, The files so far:
# hmm2tes: is our main dataframe here with output - with one row PER TE INSTANCE
# goterms: this is all the terms
# slimterms: this is the slim terms (not sure if we need, since just subset)
# slimeterms: the matching of fine-grained go terms to slim terms 

# I guess I'll make a... mega dataframe, splitting observations on "GODeDup"

megaframe <- hmm2tes %>% mutate(GOsDeDup = strsplit(as.character(GOsDeDup), ",")) %>% unnest(GOsDeDup)
megaframe <- rename(megaframe, term_ID=GOsDeDup)
slimeterms <- rename(slimeterms, term_ID=GOfine)

#slim terms
megaframe <- left_join(megaframe, slimeterms, by="term_ID")
#all the othe GO stuff
megaframe <- left_join(megaframe, goterms, by="term_ID")

View(megaframe)
object.size(megaframe)
write.csv(megaframe, "megaframe.csv")


megaframe<- read_csv("~/Documents/CPB/CPB_TEs/megaframe.csv")



#### count tables for TEs, GO terms, TE type, etc - SLIM ####
setwd("~/Documents/CPB/CPB_TEs/")
library("RCurl")
library("grid")
library("scales")
library("gridExtra")

te_type_table <- select(megaframe, TE_type) %>% group_by(TE_type) %>% tally() 
te_type_table <- na.omit(arrange(te_type_table, desc(n)))
colnames(te_type_table) <- c("TE Type", "Count")
pdf("TE_type_table.pdf", width=4, height=10)
grid.draw(tableGrob(te_type_table))
dev.off()



te_ind_table <- select(megaframe, Unique_Name, TE_type) %>% group_by(Unique_Name) %>% tally() 
te_ind_table <- na.omit(arrange(te_ind_table, desc(n)))
te_ind_table <- full_join(te_ind_table, te_names, by="Unique_Name")
te_ind_table <- select(te_ind_table, Unique_Name, n, TE_type)
colnames(te_ind_table) <- c("TE Unique ID", "Count", "TE Type")

pdf("TE_ind_table.pdf", width=10, height=60)
grid.draw(tableGrob(te_ind_table))
dev.off()







#### Make a graph of TE <-> GO_terms ####

library(igraph)
TE_GO_Edgelist <- na.omit(dplyr::select(megaframe, Unique_Name, GOSLIM, TE_type))
# TE_GO_Graph <- graph_from_data_frame(TE_GO_Edgelist)
# #owlbear = filter(dplyr::select(megaframe, Unique_Name, term_ID), Unique_Name=="Owlbear")
# #TE_GO_Graph <- graph_from_data_frame(owlbear)
# TE_GO_Graph <- graph_from_data_frame(head(TE_GO_Edgelist), n=400)
# 
# plot.igraph(TE_GO_Graph)
# is.bipartite(TE_GO_Graph)

# hg <- head(TE_GO_Edgelist, n=400)
# hg
# headgraph <- graph.data.frame(hg, directed = F)
# V(headgraph)$type <- c(rep(TRUE, nrow(unique(hg[,1]))), rep(FALSE, nrow(unique(hg[,2]))))#V(headgraph)$name %in% hg[,1]
# bipartite.projection(headgraph)
# is.bipartite(headgraph)
# plot.igraph(headgraph, layout=layout.bipartite,
#             vertex.color=c("red","green")[V(headgraph)$type+1])
# layout_as_bipartite(headgraph)
# 
# col <- c("steelblue", "orange")
# shape <- c("circle", "square")
# 
# plot(headgraph,
#      vertex.color = col[as.numeric(V(headgraph)$type)+1],
#      vertex.shape = shape[as.numeric(V(headgraph)$type)+1]
# )







hg <- TE_GO_Edgelist
hg <- hg %>% group_by_(.dots=c("Unique_Name", "GOSLIM")) %>% tally()
hg <- na.omit(hg) 

View(hg)

headgraph <- graph.data.frame(hg, directed = F)
V(headgraph)$type <- c(rep(TRUE, nrow(unique(hg[,1]))), rep(FALSE, nrow(unique(hg[,2]))))

#bipartite.projection(headgraph)
is.bipartite(headgraph) #just a check
is.directed(headgraph)

#### commmunity detection ####

# eb <- cluster_edge_betweenness(headgraph, weights = E(headgraph)$n)
# membership(eb)
# print(eb)

lc <- cluster_louvain(headgraph, weights = E(headgraph)$n)
membership(lc)
print(lc)
sizes(lc)

plot(lc, headgraph, col = membership(lc),
     mark.groups = communities(lc))
             

#### graph with full go terms ####
hg <- na.omit(dplyr::select(megaframe, Unique_Name, term_ID))

hg <- hg %>% group_by_(.dots=c("Unique_Name", "term_ID")) %>% tally()
hg <- filter(hg, term_ID!="NA")
hg <- na.omit(hg) 

View(hg)

headgraph <- graph.data.frame(hg, directed = F)
V(headgraph)$type <- c(rep(TRUE, nrow(unique(hg[,1]))), rep(FALSE, nrow(unique(hg[,2]))))

#bipartite.projection(headgraph)
is.bipartite(headgraph) #just a check
is.directed(headgraph)

#### degree distribution ####
degree_distribution()

#### commmunity detection ####

# eb <- cluster_edge_betweenness(headgraph, weights = E(headgraph)$n)
# membership(eb)
# print(eb)

lc <- cluster_louvain(headgraph, weights = E(headgraph)$n)
membership(lc)
print(lc)
sizes(lc)

plot(lc, headgraph, col = membership(lc),
     mark.groups = communities(lc))



#plot.igraph(headgraph, layout=layout.bipartite,
#            vertex.color=c("red","green")[V(headgraph)$type+1])
#layout_as_bipartite(headgraph) #this is just a layout

col <- c("steelblue", "orange")
shape <- c("circle", "square")

pdf("graphplot.pdf",width=60,height=40,paper='special') 

par(mar=c(0,0,0,0))

plot(headgraph,
     vertex.color = col[as.numeric(V(headgraph)$type)+1],
     vertex.shape = shape[as.numeric(V(headgraph)$type)+1],
     vertex.size=8,
     #E(headgraph)$label.cex <- .3,
     vertex.label=NA
     # layout = layout_with_fr(headgraph, weights=E(headgraph)$n)
)


dev.off()


lay <- layout_with_graphopt(headgraph)
plot(headgraph,
     vertex.color = col[as.numeric(V(eagle_G)$type)+1],
     vertex.shape = shape[as.numeric(V(eagle_G)$type)+1],
     vertex.size=8,
     layout=lay
     #E(headgraph)$label.cex <- .3,
     # layout = layout_with_fr(headgraph, weights=E(headgraph)$n)
)



#specific TE or GO plots

eagle <- filter(hg, Unique_Name %in% c("BehemothEagle","DevastationBeetle"))
eagle_G <- graph.data.frame(eagle, directed = F)
V(eagle_G)$type <- c(rep(TRUE, nrow(unique(eagle[,1]))), rep(FALSE, nrow(unique(eagle[,2]))))#V(headgraph)$name %in% hg[,1]
is.bipartite(eagle_G) #just a check

lay <- layout_with_graphopt(eagle_G)

plot(eagle_G,
     vertex.color = col[as.numeric(V(eagle_G)$type)+1],
     vertex.shape = shape[as.numeric(V(eagle_G)$type)+1],
     vertex.size=8,
     layout=lay
     #E(headgraph)$label.cex <- .3,
     # layout = layout_with_fr(headgraph, weights=E(headgraph)$n)
)


#### GO terms with TE frequency ####
slimcounts <- select(TE_GO_Edgelist,GOSLIM) %>% group_by(GOSLIM) %>% tally()
colnames(slimcounts) <- c("term_ID", "count")
ggplot(slimcounts) + geom_bar(stat="identity")

slimcounts <- left_join(slimcounts, goslim, by = "term_ID" ) 
View(slimcounts)
#slimcounts4table <- arrange(select(slimcounts, term_ID, count, term_name, namespace, definition), desc(count))
#write_csv(slimcounts4table, "slimcounts4tableTEST.csv")

#### standardize by... basepair? ####
# we want the 
# 1) number of basepairs per GO term
# 2) the number of genes per GO term
View(megaframe)

#slim
slim_bp <- select(megaframe, X4, X5, GOSLIM)
slim_bp <- unique(slim_bp)
slim_bp <- mutate(slim_bp, bp=X5-X4)
#count them? 
total_slim_bp <- slim_bp %>% group_by(GOSLIM) %>% summarise(total_bp = sum(bp))
colnames(total_slim_bp) <- c("term_ID", "total_bp")
#merge them
slimcounts <- left_join(slimcounts, total_slim_bp, by = "term_ID" ) 
#divide
slimcounts <- mutate(slimcounts, standardized_by_BP_counts = count/total_bp)

slimcounts4table <- arrange(select(slimcounts, term_ID, standardized_by_BP_counts, count, term_name, namespace, definition), desc(standardized_by_BP_counts))
write_csv(slimcounts4table, "slimcounts4tableTEST.csv")
slimcounts4table <- read_csv("slimcounts4tableTEST.csv")


#full go terms
gotermcounts <- na.omit(dplyr::select(megaframe, Unique_Name, term_ID)) %>% group_by(term_ID) %>% tally()
colnames(gotermcounts) <- c("term_ID", "count")

goterm_bp <- select(megaframe, X4, X5, term_ID)
goterm_bp <- unique(goterm_bp)
goterm_bp <- mutate(goterm_bp, bp=X5-X4)
#count them? 
total_goterm_bp <- goterm_bp %>% group_by(term_ID) %>% summarise(total_bp = sum(bp))
colnames(total_goterm_bp) <- c("term_ID", "total_bp")
#merge them
gotermcounts <- left_join(gotermcounts, total_goterm_bp, by = "term_ID" ) 
gotermcounts <- left_join(gotermcounts, goterms, by = "term_ID" ) 
#divide
gotermcounts <- mutate(gotermcounts, standardized_by_BP_counts = count/total_bp)

gotermcounts4table <- arrange(select(gotermcounts, term_ID, standardized_by_BP_counts, count, term_name,term_ID, namespace, definition), desc(standardized_by_BP_counts))
write_csv(gotermcounts4table, "gotermcounts4tableTEST.csv")

#### 2) genes per GO term #### 
slim_genes <- select(megaframe, OGSv1.0_gene_ID, GOSLIM)
slim_genes <- na.omit(unique(slim_genes))
slim_genes <- select(slim_genes,GOSLIM) %>% group_by(GOSLIM) %>% tally()
colnames(slim_genes) <- c("term_ID", "total_slim_genes")
slimcounts <- left_join(slimcounts, slim_genes, by = "term_ID" ) 

slimcounts <- mutate(slimcounts, standardized_by_gene_counts = count/total_slim_genes)
slimcounts4table <- arrange(select(slimcounts, term_ID, standardized_by_BP_counts, standardized_by_gene_counts,count, term_name, namespace, definition), desc(standardized_by_BP_counts))
#write_csv(slimcounts4table, "slimcounts4tableTEST.csv")

#full go terms
go_genes <- select(megaframe, OGSv1.0_gene_ID, term_ID)
go_genes <- na.omit(unique(go_genes))
go_genes <- select(go_genes,term_ID) %>% group_by(term_ID) %>% tally()
colnames(go_genes) <- c("term_ID", "total_GO_genes")
go_counts <- left_join(gotermcounts, go_genes, by = "term_ID" ) 

go_counts <- mutate(go_counts, standardized_by_gene_counts = count/total_GO_genes)
gocounts4table <- arrange(select(go_counts, term_ID, standardized_by_BP_counts, standardized_by_gene_counts,count, term_name, namespace, definition), desc(standardized_by_BP_counts))
View(gocounts4table)

#### bipartite package - ecology stuff ####
library(bipartite)

hg$genome = "CPBref"
colnames(hg) <- c("higher","lower","freq","webID")
hg <- data.frame(hg[,c(1,2,4,3)])
#bipgraph <- frame2webs(hg, varnames=c("GOSLIM", "Unique_Name", "genome", "n"))
bipgraph <- frame2webs(hg, type.out="array")

#neither of these work
visweb(bipgraph)
plotweb(bipgraph)

#scale between 0 and 1 
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
hg$freq01 <- range01(hg$freq)

#this is sorta not great
ggplot(data = hg, aes(x=higher, y=lower, fill=freq01)) + 
  geom_tile()

bipgraph <- frame2webs(hg,type.out = "list")
networklevel(bipgraph)

#### network statistics ####

hg2 <- acast(hg, GOSLIM ~ Unique_Name, sum)
View(hg2)
plotweb(hg2)
nlhg2 <- networklevel(hg2,index="ALLBUTDD")


hldata(Safariland)
View(Safariland)
typeof(Safariland)
typeof(as.integer(bipgraph))
nls <- networklevel(Safariland)
data(bipgraph)
bipgraph <- as.matrix(bipgraph)


degreedistr(hg2)

?heatmap
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "", Rowv=-rowSums(x), Colv=-rowMeans(x))




