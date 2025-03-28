# Code by: Isabelle Du Plessis, 2023

require(ggplot2) # version 3.4.2
require(readr) # version 2.1.4
require(gridExtra) # version 2.3
require(ggpubr) # version 0.6.0
require(ggvenn) # version 0.1.10

# R version 4.0.2 (2020-06-22)

#### Load Data ####

# Load sample metadata
metadata = read.csv("../intermediate_files/metadata.csv")

# Load vOTU data
votus_cov75thres = read.delim("../intermediate_files/votus_cov75thres_1224.txt", header=TRUE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")
votus_cov75thres$sample = parse_number(votus_cov75thres$sample)

# Load host data
HostPrediction = read.delim("../intermediate_files/HostPrediction_0225.tsv")
hostdf = read.delim("../intermediate_files/hostProperties_0225.tsv")
colnames(hostdf) = c("hostorder", "sulfuroxidizing", "sulfatereducing", "ironoxidizing", "nitrifying")
#colnames(hostdf) = c("virus", "hostorder", "sulfuroxidizing", "sulfatereducing", "ironoxidizing", "nitrifying")


#### Create main dataframe ####

# Create dataframe with all vOTUs and all info
bigdf = data.frame(votus_cov75thres[, c(1,2)])
colnames(bigdf) = c("sample","virus")

# Initialize host vectors with domain, phylum, and order values
d <- rep("Unknown", nrow(bigdf))
p <- rep("Unknown", nrow(bigdf))
o <- rep("Unknown", nrow(bigdf))

# For each vOTU, find corresponding compartment, phenotype, and host prediction
for(i in 1:nrow(bigdf)){
  bigdf$Compartment[i]=metadata$Compartment[bigdf$sample[i]==metadata$Sample_ID] # Add compartment to bigdf
  bigdf$Phenotype[i]=metadata$Spartina[bigdf$sample[i]==metadata$Sample_ID] # Add phenotype to bigdf
  if(bigdf$virus[i] %in% HostPrediction$vOTUs){ # For vOTUs that have hosts predicted, add host taxonomy to host vectors
    index=which(bigdf$virus[i]==HostPrediction$vOTUs)
    d[i]=HostPrediction$Domain[index]
    p[i]=HostPrediction$Phylum[index]
    o[i]=HostPrediction$Order[index]
  }
}

# Set new bigdf columns equal to host vectors
bigdf$Domain=d 
bigdf$Phylum=p
bigdf$Order=o

# Factor compartments for plotting
bigdf$Compartment = factor(bigdf$Compartment,levels = c("Bulk sediment","Rhizosphere","Endosphere"))

# Main dataframe is complete

#### Phylum Stacked Bar Plot - Figure 5a ####

# Get phyla of predicted hosts, removing replicate viruses in the same compartment/phenotype
phylumdf = bigdf[,c(2,3,4,6)] # Get relevant columns: Compartment, Phenotype, and Phylum
phylumdf = unique(phylumdf[which(phylumdf$Phylum!="Unknown"),]) # Get unique rows to remove replicates, get only known hosts
phylumdf$value = 1 # For plotting

# Get proportion of vOTUs with hosts predicted out of all vOTUs per compartment/phenotype
phylumdf$hostcount = 0 # Initialize column for number of vOTUs with hosts predicted

# Get number of vOTUs with each phylum host prediction, per compartment/phenotype
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),])

phylumdf$totalcount = 0 # Initialize column for number of vOTUs total
tmp = unique(bigdf[,c(2,3,4)]) # Get unique viruses, compartment, and phenotype

# Get total number of vOTUs in each compartment/phenotype
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Bulk sediment" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Bulk sediment" & tmp$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Rhizosphere" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Rhizosphere" & tmp$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Endosphere" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Endosphere" & tmp$Phenotype=="Short"),])

# Plot phylum vs phenotype/compartment
phyla = ggplot(phylumdf, aes(interaction(Phenotype, Compartment), value, fill=Phylum)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_test() +
  ylab("Proportion of Phylum-Level Hosts") +
  scale_fill_manual(values = c("#89C5DA", "#DA5724", "#CE50CA", "#8569D5", "#C0717C", 
                               "#689030", "#673770", "#5F7FC7", "#38333E", "#508578", "#C84248", "#CD9BCD", "#7FDCC0",
                               "#D1A33D")) +
  scale_x_discrete("Sample Compartment and Phenotype", labels=c("Bulk Sediment\nShort", "Bulk Sediment\nTall", "Rhizosphere\nShort", "Rhizosphere\nTall", "Root\nShort", "Root\nTall")) +
  theme(text = element_text(size = 15)) + theme(axis.text.x = element_text(size = 10)) +
  ggtitle("Predicted Target Host Phyla") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5)) +
  geom_text(aes(label = paste(round((hostcount/totalcount)*100), "%")), vjust = -0.4)

pdf("../Figures/Figure4.pdf",width=10)
print(phyla)
dev.off()
