library(phyloseq)
library(phylosmith)
library(dada2)
library(plyr)
library(ggplot2)
library(vegan)
library(reshape2)


setwd("~/Downloads/CABBI_RNA_to_DNA/")
#Generate phyloseq objects from the OTU table and taxonomy table generated from DADA2, 16s rRNA ata.
otu_rrna_16s <- readRDS("seq_table.RDS")
tax_rrna_16s <- readRDS("tax_table.RDS")
metadata <- read.csv(file="metadata_cabbi_rrna_16s.csv")
colnames(otu_rrna_16s) <- metadata$Samples
phylo.rrna <- phyloseq(otu_table(otu_rrna_16s, taxa_are_rows=TRUE), tax_table(tax_rrna_16s), metadata)

rownames(metadata) <- metadata$Sample
sample_data(phylo.rrna) <- metadata
phylo.rrna

otu_cabbi <- readRDS("seq_table_16srDNA.RDS")
tax_cabbi <- readRDS("tax_table_16srDNA.RDS")
ps <- phyloseq(otu_table(otu_cabbi, taxa_are_rows=TRUE), tax_table(tax_cabbi))
metadata.rdna<-read.csv(file="meatadata_cabbi_rDNA_481.csv")
rownames(metadata.rdna)<- metadata.rdna$Samples
sample_data(ps)<- metadata.rdna

#As you could see for the name of samples in DNA and RNA, the N#_P# are swapted. so i will make the order of these names are the same
#To swapt the name order in 16s rRNA sequencing
sample_names_s <- vector()
for(name in rownames(metadata)){
  begin <- strsplit(name, '_')[[1]][1:3]
  end <- strsplit(name, '_')[[1]][6:9]
  check <- strsplit(name, '_')[[1]][4:5]
  if(strsplit(check[1], '')[[1]][1] == 'N'){
    check <- check[c(2,1)]
  }
  sample <- paste(c(begin, check, end), sep = '_', collapse = '_')
  sample_names_s <- rbind(sample_names_s, sample)
sample_names(phylo.rrna) <- sample_names_s
head(sample_data(phylo.rrna))

#heatmap.rrna <- phylosmith::abundance_heatmap(phylo.rrna, classification = 'Phylum', treatment = c('Year'),  transformation = 'log10', colors = 'default')


## Investigating the differences between whole community and active community


merged.phylo <- merge_phyloseq(phylo.rrna, ps)
merged.phylo.1 <- prune_species(speciesSums(merged.phylo)>1, merged.phylo)
saveRDS(merged.phylo.1, file="data/merged_DNA_RNA_phyloseq.rds")
merged.phylo.1 <-readRDS(file="data/merged_DNA_RNA_phyloseq.rds")
merged.phylo.1
```

## Alpha diversity
```{r}
richness <- estimate_richness(merged.phylo.1, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
sam <- sample_data(merged.phylo.1)
sa <- data.frame(rownames(sam), sam$Year, sam$Nitrogen, sam$Days, sam$Nucelic)

rownames(richness) = gsub(pattern="X*", replacement="", x=rownames(richness))
merged <- merge(richness, sa, by.x= "row.names", by.y="rownames.sam.")
merged.melt <- melt(merged, id.vars=c("Row.names","sam.Nucelic"), measure.vars=c("Chao1","Shannon","InvSimpson"))
alpha.merged <- ggplot(merged.melt, aes(x=sam.Nucelic,y=value))+geom_boxplot()+facet_grid(~variable)+
  theme(axis.text.x = element_text(angle = 90, size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(angle=90),
        panel.background=element_rect(fill="white", colour="white"),
        axis.title = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title=element_text(size=16),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        strip.background = element_rect(colour = "black", fill = "white",size=2, linetype="solid"),
        strip.text=element_text(size=14))+
  xlab("")
compa <- list(c("RNA","DNA"))
alpha.merged.pvalue<- alpha.merged +stat_compare_means(comparisons = compa, label="p.signif", method="t.test")#+ #add pairwise comparisons p-value
  #stat_compare_means(label.y=4, method = "anova") 
ggsave(alpha.merged.pvalue, file="figures/alpha_diversity_merged_two_phylo.pdf", units="in", width=6, height=6)
alpha.merged.pvalue

```
The alpha diversity is significanlty higher associated with 16s rDNA sequencing, which makes sense, indicating the whole community has more species than the active communities.

So, what are different?
Let me take a look of the first 10 most abundant phyla associated with rDNA /rRNA
Rarefying
```{r}
set.seed(100)
merged.phylo.1.rare <- rarefy_even_depth(merged.phylo.1, sample.size=5000, rngseed=TRUE)
# removed 72 samples.
```

```{r}
merged.relative <- transform_sample_counts(merged.phylo.1.rare, function(x) x/sum(x))
merged.relative.glom <- tax_glom(merged.relative, taxrank="Phylum") 
Top10phylum.merged.glom = names(sort(taxa_sums(merged.relative.glom), decreasing=T))[1:10]
Top10phylum.merged.phylo <- prune_taxa(Top10phylum.merged.glom,merged.relative.glom)
```




