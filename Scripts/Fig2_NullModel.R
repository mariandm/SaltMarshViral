### Code for Null model of Figure 2
### By Marian Dominguez-Mirazo, Isabelle Du Plessis, Hannah Snyder

require(reshape2)
require(pracma)

######################################
#### Load and manipulate observations

# Load data
metadata = read.csv("../intermediate_files/metadata.csv")
metadata$Sample_ID =paste("sample",metadata$Sample_ID,sep = "_")
votus_cov75thres = read.table("../intermediate_files//votus_cov75thres_1224.txt")
votus_cov75thres$Compartment = metadata$Compartment[match(votus_cov75thres[,1],metadata$Sample_ID)]
votus_cov75thres$Phenotype = metadata$Spartina[match(votus_cov75thres[,1],metadata$Sample_ID)]

# As in Figure2 script, only use the vOTUs that appear in more than one sample
this_votus_cov75thres = votus_cov75thres[votus_cov75thres[,2]%in%rownames(table(votus_cov75thres[,2]))[table(votus_cov75thres[,2])>1],]
# Keep only relevant information
tmp = this_votus_cov75thres[,c(2,6,7)]
# How many vOTUs appear in n samples
cnt_table = table(table(tmp$V2))

# Get a probability table per phenotype/compartment combination (it's counts really)
prob_table = melt(table(tmp[,2:3]))

# Create an observed table of compartment and phenotype multicombination
tmp2 = tmp
colnames(tmp2) = c("vOTU","variable","variable")
new_compress = t(2^(4:0) * t(table(unique(rbind(tmp2[,1:2],tmp2[,c(1,3)])))))
match_compartment = data.frame("Compartment" = c("bulk","endo","rhizo","bulk_endo","bulk_rhizo","endo_rhizo","bulk_endo_rhizo"),
                               "Code"=c(16,8,4,24,20,12,28))
match_phenotype = data.frame("Phenotype"=c("short","tall","both"),
                             "Code" = c(2,1,3))
this_table = data.frame("Compartment" = match_compartment$Compartment[match(rowSums(new_compress[,1:3]),match_compartment$Code)],
                        "Phenotype" = match_phenotype$Phenotype[match(rowSums(new_compress[,4:5]),match_phenotype$Code)])

observed = table(this_table)


######################################
#### Run the null model

#### Get per compartment / phenotype combination probability
ps = prob_table$value/sum(prob_table$value)
# Corresponding comaprtment / phenotype label 
p_labels = paste(prob_table[,1],prob_table[,2],sep="_")

# Allocate space for per combination probability
# At the end, this should look like prob_table
expected_prob = data.frame("Compartment"=prob_table[,1],
                           "Phenotype"=prob_table[,2],
                           "count"=0)
# Allocate space, cols multi phenotypes, rows multi compartments
mean_fill = matrix(0,ncol=3,nrow=7)

# Count how many times each votus appears
votu_cnt = table(tmp[,1])
#corresponding votu names 
votu_names = rownames(votu_cnt)

nsims = 1e4
for(k in 1:nsims){
  # initialize per simulation matrix
  fill = matrix(0,ncol=length(p_labels),nrow=length(votu_names))
  colnames(fill) = p_labels
  rownames(fill) = votu_names
  
  # loop across each votu
  for(j in 1:length(votu_names)){
    this_votu = votu_names[j]
    this_n = votu_cnt[j]

    # Loop through votu number of occurrences
    for(i in 1:this_n){
      match_idx = which(cumsum(ps)>runif(1))[1]
      fill[this_votu,p_labels[match_idx]] = fill[this_votu,p_labels[match_idx]] + 1
    }
  }
  
  # Count combination appearances
  expected_prob$count = expected_prob$count + colSums(fill)
  
  # Get overlap distribution using binary numbers
  compress = matrix(0,ncol=5,nrow=length(votu_names))
  colnames(compress) = c("bulk","endo","rhizo","tall","short")
  rownames(compress) = votu_names
  # Is it in bulk, any phenotype
  compress[,"bulk"] = (rowSums(fill[,c(1,4)])>0)*1
  # Is it in endo, any phenotype
  compress[,"endo"] = (rowSums(fill[,c(2,5)])>0)*1
  # Is it in rhizo, any phenotype
  compress[,"rhizo"] = (rowSums(fill[,c(3,6)])>0)*1
  # Is it in tall, any compartment
  compress[,"tall"] = (rowSums(fill[,4:6])>0)*1
  # Is it in short, any compartment
  compress[,"short"] = (rowSums(fill[,1:3])>0)*1
  binary = 2^(4:0)
  new_compress = t(t(compress)*binary)
  # Match compartment/phenotype presence with binary code
  match_compartment = data.frame("Compartment" = c("bulk","endo","rhizo","bulk_endo","bulk_rhizo","endo_rhizo","bulk_endo_rhizo"),
                                 "Code"=c(16,8,4,24,20,12,28))
  match_phenotype = data.frame("Phenotype"=c("tall","short","both"),
                               "Code" = c(2,1,3))
  # Count each type
  this_table = data.frame("Compartment" = match_compartment$Compartment[match(rowSums(new_compress[,1:3]),match_compartment$Code)],
                          "Phenotype" = match_phenotype$Phenotype[match(rowSums(new_compress[,4:5]),match_phenotype$Code)])
  mean_fill = mean_fill + table(this_table)
}

# Normalize by number of simulations
expected_prob_normalized = expected_prob
expected_prob_normalized$count = expected_prob_normalized$count/nsims
expected = mean_fill/nsims

# Check that the combination probability tables look the same
# Then the null model works
prob_table
expected_prob_normalized
# Beautiful
# Save null model 
write.table(expected,"../intermediate_files/nullmodel.txt",sep="\t",quote = F)

