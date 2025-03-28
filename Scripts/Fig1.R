# Code for Figure 1 b,c,d
# by Isabelle Du Plessis, Hannah Snyder, and Marian Dominguez-Mirazo, 2025

require(ggplot2) # version 3.4.2
require(reshape2) # version 1.4.4
require(vegan) # version 2.5-6
require(grid) # version 4.0.2
require(gridExtra) # version 2.3

#### Load Data ####

path = getwd() # Set path

# load sample metadata
metadata <- read.csv("../intermediate_files/metadata.csv")
#change sample id so that we can compare
metadata$Sample_ID = paste("sample_",metadata$Sample_ID,sep="")
#change order so that we can compare
metadata = metadata[order(metadata$Sample_ID),]

# load votu data
votus_cov75thres <- read.delim("../intermediate_files/votus_cov75thres_1224.txt", header=TRUE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")

### Get TPM from RPKM calculations
votus_cov75thres$tpm = 0 # Initialize tpm column
# Formula: TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6
for (i in 1:nrow(votus_cov75thres)){
  rpkm = votus_cov75thres[i,5]
  sum_rpkm = sum(votus_cov75thres$rpkm)
  tpm = rpkm/sum_rpkm*10^6
  votus_cov75thres[i,6] = tpm
}

##############################################
## Figure 1b: Cumulative vOTUs / Saturation Plot
##############################################

# Permutate the sampling order to generate a saturation plot
perms = 100
# Create a presence-absence table for vOTU per sample
le_tmp = table(votus_cov75thres[,1:2])
# initialize storage
main = matrix(0,nrow=perms,ncol=24)

# For each permutation
for(perm in 1:perms){
  # Record of which vOTUs have been counted
  tots = rep(0,769)
  # Permutate sample order 
  tmp = le_tmp[sample(1:24),]
  # for each sample
  for(i in 1:24){
    # Get vOTUs found in this sample and remove the ones we have already counted before
    tmp2= tmp[i,] - tots
    # Prevent negative numbers (for vOTUs found in previous samples but not in this one)
    tmp2[tmp2<0]=0
    # Count how many new vOTUs
    main[perm,i] = sum(tmp2) #change the 1
    # Update the record of observed vOTUs
    tots = tots + tmp2
    # Prevent double counting vOTU apperance
    tots[tots>1]=1
  }
  # Cumsum of observed vOTUs 
  main[perm,] = cumsum(main[perm,])
}
# Average saturation line
lemean = colSums(main)/perms
# Make into ggplot friendly format
df_mean = data.frame("sample"=1:24,"means"=lemean)
df = melt(main)

# Plot
saturation = ggplot(df_mean,aes(sample, means)) +
  geom_point(data=df, aes(Var2, value, group='Var1'),color = "darkgrey") +
  geom_line(color="dodgerblue2", linewidth = 1) +
  xlab("Number of Samples") +
  scale_y_continuous("Cumulative Number of vOTUs") +
  theme_test() + ggtitle("") + 
  theme(text = element_text(size = 20), plot.title = element_text(face = "bold", size = 15, hjust = .5),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5),"pt"))

# Save graph 
pdf("../Figures/Figure1b_saturation.pdf",height=5, width = 5)
print(saturation)
dev.off()

##############################################
### Figure 1c: vOTU richness
##############################################

# Count how many vOTUs appear in each sample
tmp = melt(table(votus_cov75thres$sample))
colnames(tmp) = c("sample","votus")
# Add metadata 
tmp = cbind(tmp, metadata[,2:5])
tmp$Compartment = factor(tmp$Compartment,levels = c("Bulk sediment","Rhizosphere","Endosphere"))

#Plot 
votu_graph = ggplot(data=tmp,aes(Compartment, votus, color=Compartment, group=interaction(Compartment,Spartina))) +
  geom_point(aes(shape = Spartina), size=2.5, position = position_dodge(width = .9)) +
  theme_test() + 
  stat_summary(fun = mean, geom = "tile", color = "black", height = ((max(tmp$votus)-min(tmp$votus))*.003), width = .8, position = position_dodge(width = .9)) +
  scale_color_manual(values = c("Bulk sediment" = "#853512",
                                "Endosphere" = "#558A78",
                                "Rhizosphere" = "#EEAA23"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  scale_x_discrete(NULL, labels = c("Bulk Sediment", "Rhizosphere", "Root")) + 
  ylab("Number of vOTUs") +
  theme(text = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5.5, b = 0, l = 0))) +
  theme(legend.position = "none") +
  labs(color = NULL, shape = NULL) +
  theme(plot.margin = unit(c(15, 15, 30, 5.5),"pt")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Save plot 
pdf("../Figures/Fig1c_richness.pdf",height = 4,width = 5)
print(votu_graph)
dev.off()

##############################################
### Figure 1c: Statistical Tests 
##############################################

# Get the number of vOTUs in samples per compartment and phenotype combination
BT = tmp$votus[tmp$Compartment =="Bulk sediment"& tmp$Spartina=="Tall"]
BS = tmp$votus[tmp$Compartment =="Bulk sediment"& tmp$Spartina=="Short"]
RhT = tmp$votus[tmp$Compartment =="Rhizosphere"& tmp$Spartina=="Tall"]
RhS = tmp$votus[tmp$Compartment =="Rhizosphere"& tmp$Spartina=="Short"]
RoT = tmp$votus[tmp$Compartment =="Endosphere"& tmp$Spartina=="Tall"]
RoS = tmp$votus[tmp$Compartment =="Endosphere"& tmp$Spartina=="Short"]

# Wilcoxon test
# Short vs tall
wilcox.test(BT,BS) #NS
wilcox.test(RhT,RhS) # p-value = 0.02857
wilcox.test(RoT,RoS) # p-value = 0.02857
# Tall vs tall
wilcox.test(BT,RhT) #NS
wilcox.test(BT,RoT) #NS
wilcox.test(RhT,RoT) #NS
# Short vs short
wilcox.test(BS,RhS) # p-value = 0.02857
wilcox.test(BS,RoS) #NS
wilcox.test(RhS,RoS) #NS but almost lol

##############################################
### Figure 1d: NMDS
##############################################

sample = votus_cov75thres$sample
votu = votus_cov75thres$votus
tpm = votus_cov75thres$tpm
# Create a new data frame
newdf2 = data.frame(sample,votu,tpm)

votu <- unique(newdf2$votu)
samples <- sort(unique(newdf2$sample)) 
# Initialize a matrix with n rows for sample number
# And n cols for votu number
df1 <- matrix(nrow = length(samples), ncol = length(votu), dimnames = list(samples,votu))

# Make dataframe into matrix format to run NMDS
for(r in 1:nrow(newdf2)){
  samp <- newdf2[r, 1]
  tax <- newdf2[r, 2]
  df1[samp,tax] <- newdf2[r, 3]
} # 1, 2, 3 here relate the the column number in the raw data in which the sample name, species name and data are in

df1[is.na(df1)] <- 0   #convert NA's to 0

# Run NMDS
nmds <- metaMDS(df1, distance = "bray", autotransform = FALSE)
data.scores = as.data.frame((scores(nmds)$sites))

# Add information on compartment and phenotype
data.scores$Sample = samples
height_list = c()
compartment_list = c()
for (i in (data.scores$Sample)){
  if (i %in% (c("sample_1","sample_7","sample_13","sample_19"))){
    compartment_list = append(compartment_list,"Bulk")
    height_list = append(height_list,"Tall")
  }
  if (i %in% (c("sample_2","sample_8","sample_14","sample_20"))){
    compartment_list = append(compartment_list,"Rhizosphere")
    height_list = append(height_list,"Tall")
  }
  if (i %in%(c("sample_3","sample_9","sample_15","sample_21"))){
    compartment_list = append(compartment_list, "Endosphere")
    height_list = append(height_list, "Tall")
  }
  if (i %in%c("sample_4","sample_10","sample_16","sample_22")){
    compartment_list = append(compartment_list, "Bulk")
    height_list = append(height_list, "Short")
  }
  if (i %in%(c("sample_5","sample_11","sample_17","sample_23"))){
    compartment_list = append(compartment_list, "Rhizosphere")
    height_list = append(height_list, "Short")
  }
  if (i %in% c("sample_6","sample_12","sample_18","sample_24")){
    compartment_list = append(compartment_list, "Endosphere")
    height_list = append(height_list, "Short")
  }
}
data.scores$Compartment = compartment_list
data.scores$Height = height_list

nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2.5, aes( shape = Height, colour = Compartment)) + theme_test() +
  scale_colour_manual(values = c("Bulk" = "#853512", "Endosphere" = "#558A78","Rhizosphere" = "#EEAA23"), limits = c("Bulk", "Rhizosphere", "Endosphere"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  ggtitle("") + labs(color = "", shape = "") +
  theme(text = element_text(size = 15), plot.title = element_text(face = "bold", size = 15, hjust = .5))

pdf("../Figures/Figure1d_NMDS.pdf",height=4, width = 5.5)
print(nmds_plot)
dev.off()

