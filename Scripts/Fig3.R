# Plot vOTU network
# Only includes VC subclusters and edges to the vOTUs found in this study
# Includes vOTUs classified as clustered and Overlap
# Code by Hannah Snyder and Marian Dominguez-Mirazo

#Import libraries
library(ggplot2) #version ‘3.4.3’
library(GGally) #version ‘2.1.2’
library(tidyverse) #version ‘2.0.0’

theclusters = read.csv("../intermediate_files/relevant_clusters.list",header = F, col.names = c("votus","cluster"))
clusters = unique(theclusters$cluster)
la = read.csv("../intermediate_files/c1.ntw",header=F)

# Get cluster data into a format that includes: 
# our votus, other votus in the same cluster, cluster
save_here = matrix(nrow=0,ncol=3)
for(i in 1:length(clusters)){
  sublist = theclusters[theclusters$cluster==clusters[i],]
  this_ours = sublist$votus[grep("SMV25_",sublist$votus)]
  comb_votus = rep(sublist$votus,length(this_ours))
  #comb_votus = rep(sublist$votus,length(sublist$votus)) # include non ours links
  comb_ours = sort(rep(this_ours,length(sublist$votus))) 
  #comb_ours = sort(rep(sublist$votus,length(sublist$votus))) # include non ours links
  comb_cluster = rep(clusters[i],length(comb_ours))
  save_here = rbind(save_here,cbind(comb_ours,comb_votus,comb_cluster))
}
# Remove same votu in both sides 
save_here2 = save_here[save_here[,1]!=save_here[,2],]
# Invert both sides
save_here2 = unique(rbind(save_here2,cbind(save_here[,2],save_here[,1],save_here[,3])))
#Match to data in c1.ntw
new_la = la[match(paste(save_here[,1],save_here[,2]),paste(la$V1,la$V2)),]
# Remove nas
newnew_la = new_la[!is.na(new_la$V3),]
# Create a matrix with edges cause dataframes dont work for some reason
mmm = matrix(c(newnew_la$V1,newnew_la$V2),byrow=F,ncol=2)
# Create dataframe for edges
ntwk_our= newnew_la
colnames(ntwk_our) = c("OTU1","OTU2","Score")

# Create nodes
nodes = ggnet2(mmm[,1:2], mode = "fruchtermanreingold", layout.par = list(list=(niter=2000)))$data 
# Create edges
edges <- ntwk_our %>% 
  mutate(Pair = paste(OTU1, OTU2, sep = ".")) %>% 
  gather(key = "Member", value = "label", -Pair, -Score) %>% 
  inner_join(nodes, by = "label")


# include habitat information 
info_nodes = nodes
info_nodes$source = "RefSeq"
info_nodes$source[grep("EarthsVirome",info_nodes$label)] = "earths virome"
info_nodes$source[grep("gary_all",info_nodes$label)] = "gary all"
info_nodes$source[grep("alaska_puertorico",info_nodes$label)] = "alaska puertorico"
info_nodes$source[grep("borton",info_nodes$label)] = "borton"
info_nodes$source[grep("GOV",info_nodes$label)] = "GOV"
info_nodes$source[grep("biochar",info_nodes$label)] = "biochar"
info_nodes$source[grep("SPRUCE_viral_seq",info_nodes$label)] = "SPRUCE viral seq"
info_nodes$source[grep("merson",info_nodes$label)] = "merson"
info_nodes$source[grep("SMV25_",info_nodes$label)] = "our"

habitats = read.csv("../intermediate_files/PIGEON_v1_habitat_origin.csv")
habitats_needed = habitats[,c(-2,-3)]

#Remove the refseq viruses and our viruses as they will not appear in the habitat dataframe
our = info_nodes[info_nodes$source == "our",]
our$contig_origin = "This study"

#Create a df that has the habitat for known votu collection sites
known_nodes = info_nodes[(info_nodes$source != "our" & info_nodes$source != "Ref-seq"),]
known_habitats = merge(known_nodes, habitats_needed,  by.x = "label", by.y = "pigeon_header", all = TRUE)

#combine all of the votu data back into one dataframe
complete = our
complete = rbind(complete, known_habitats)

# alter certain habitat names so we can join them as one later
complete$contig_origin[grep("Marine",complete$contig_origin)] = "Aquatic"
complete$contig_origin[grep("Host",complete$contig_origin)] = "Plant-associated"
complete$contig_origin[grep("Plants",complete$contig_origin)] = "Plant-associated"
complete$contig_origin[grep("Peat",complete$contig_origin)] = "Wetland"
complete$contig_origin[grep("Wetland Soil",complete$contig_origin)] = "Wetland"
complete$contig_origin[grep("SPRUCE",complete$source)] = "Wetland"
complete$contig_origin[grep("Terrestrial",complete$contig_origin)] = "Terrestrial"
complete$contig_origin[grep("soil",complete$contig_origin)] = "Other soils"
complete$contig_origin[grep("Soil",complete$contig_origin)] = "Other soils"
#getting rid of extra rows not related to our samples
complete = na.omit(complete)
# Change order in legend
complete$contig_origin = factor(complete$contig_origin,c("This study", "RefSeq","Aquatic","Plant-associated","Wetland","Terrestrial","Other soils"))


# Plot
thisplot = ggplot() +
  geom_line(data = edges, aes(x,y,group = Pair), color = "black", linewidth = 0.5, alpha = 0.1) +
  geom_point(data = complete, aes(x, y,color=contig_origin),alpha = 0.8, shape = 16, size = 2) +
  theme_minimal() + 
  scale_color_manual(values=c("RefSeq" = "gray","This study" = "red", "Aquatic" = "deepskyblue4", "Plant-associated" = "seagreen4", "Other soils" = "navajowhite3", "Terrestrial" = "burlywood4", "Wetland" = "darkorchid4" ),
                     name = "Sequence origin") +
  theme_minimal() +
  theme(text = element_text(size = 15), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(2, 'line'),
        aspect.ratio = 1)

pdf(file = ("../Figures/Fig3.pdf"), width=8, height = 5)
print(thisplot)
dev.off()

