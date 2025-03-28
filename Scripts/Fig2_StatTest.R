### Code for Statistical tests of Figure 2
### By Marian Dominguez-Mirazo, Isabelle Du Plessis, Hannah Snyder

# Load data
metadata = read.csv("../intermediate_files/metadata.csv")
metadata$Sample_ID =paste("sample",metadata$Sample_ID,sep = "_")
votus_cov75thres = read.table("../intermediate_files//votus_cov75thres_1224.txt")
votus_cov75thres$Compartment = metadata$Compartment[match(votus_cov75thres[,1],metadata$Sample_ID)]
votus_cov75thres$Phenotype = metadata$Spartina[match(votus_cov75thres[,1],metadata$Sample_ID)]

# As in Figure2 script, only use the vOTUs that appear in more than one sample
this_votus_cov75thres = votus_cov75thres[votus_cov75thres[,2]%in%rownames(table(votus_cov75thres[,2]))[table(votus_cov75thres[,2])>1],]
# Keep only relevant information
tmp2 = this_votus_cov75thres[,c(2,6,7)]

# Create an observed table of compartment and phenotype multicombination
colnames(tmp2) = c("vOTU","variable","variable")
new_compress = t(2^(4:0) * t(table(unique(rbind(tmp2[,1:2],tmp2[,c(1,3)])))))
match_compartment = data.frame("Compartment" = c("bulk","endo","rhizo","bulk_endo","bulk_rhizo","endo_rhizo","bulk_endo_rhizo"),
                               "Code"=c(16,8,4,24,20,12,28))
match_phenotype = data.frame("Phenotype"=c("short","tall","both"),
                             "Code" = c(2,1,3))
this_table = data.frame("Compartment" = match_compartment$Compartment[match(rowSums(new_compress[,1:3]),match_compartment$Code)],
                        "Phenotype" = match_phenotype$Phenotype[match(rowSums(new_compress[,4:5]),match_phenotype$Code)])

observed = table(this_table)

# Load null model
expected = read.table("../intermediate_files/nullmodel.txt",row.names = 1)

###############################################
######### Statistical tests 

# For phenotype sharedness (one phenotype, two phenotypes)
# uy si very significative
expected_phenotype_both = sum(expected[,1])
expected_phenotype_one = sum(expected[,2:3])
observed_phenotype_both = sum(observed[,1])
observed_phenotype_one = sum(observed[,2:3])
chisq = (observed_phenotype_both - expected_phenotype_both)^2/expected_phenotype_both + (observed_phenotype_one - expected_phenotype_one)^2/expected_phenotype_one
pchisq(chisq, 1, lower.tail = FALSE)

# For distribution of short and tall appearing in a single phenotype
# Significance here too
observed_short_one = sum(observed[,2])
observed_tall_one = sum(observed[,3])
expected_short_one = sum(expected[,2])/sum(expected[,2:3]) * observed_phenotype_one
expected_tall_one = sum(expected[,3])/sum(expected[,2:3]) * observed_phenotype_one
chisq = (observed_short_one - expected_short_one)^2/expected_short_one + (observed_tall_one - expected_tall_one)^2/expected_tall_one
pchisq(chisq, 1, lower.tail = FALSE)

# For Compartment sharedness (one compartment, two compartments, three compartments)
# Significance here too
observed_compartment_one = sum(observed[c(1,4,6),])
expected_compartment_one = sum(expected[c(1,5,7),])
observed_compartment_two = sum(observed[c(3,5),])
expected_compartment_two = sum(expected[c(2,4,6),])
observed_compartment_three = sum(observed[2,])
expected_compartment_three = sum(expected[3,])
chisq = (observed_compartment_one - expected_compartment_one)^2/expected_compartment_one + (observed_compartment_two - expected_compartment_two)^2/expected_compartment_two + (observed_compartment_three - expected_compartment_three)^2/expected_compartment_three
pchisq(chisq, 2, lower.tail = FALSE)

# For distribution of compartments appearing in a single compartment
# Significance here too
observed_bulk = sum(observed[1,])
observed_endo = sum(observed[4,])
observed_rhizo = sum(observed[6,])
expected_bulk = sum(expected[1,])/sum(expected[c(1,5,7),]) * observed_compartment_one
expected_endo = sum(expected[5,])/sum(expected[c(1,5,7),]) * observed_compartment_one
expected_rhizo = sum(expected[7,])/sum(expected[c(1,5,7),]) * observed_compartment_one

chisq = (observed_bulk - expected_bulk)^2/expected_bulk + (observed_rhizo - expected_rhizo)^2/expected_rhizo + (observed_endo - expected_endo)^2/expected_endo
pchisq(chisq, 2, lower.tail = FALSE)

# For distribution of compartments appearing in two compartments
# Significance
expected_bulk_endo = sum(expected[2,]) / sum(expected[c(2,4,6),]) * expected_compartment_two
expected_bulk_rhizo = sum(expected[4,]) / sum(expected[c(2,4,6),]) * expected_compartment_two
expected_endo_rhizo = sum(expected[6,]) / sum(expected[c(2,4,6),]) * expected_compartment_two
observed_bulk_endo = 0
observed_bulk_rhizo = sum(observed[3,])
observed_endo_rhizo = sum(observed[5,])

chisq = (observed_bulk_endo - expected_bulk_endo)^2/expected_bulk_endo + (observed_bulk_rhizo - expected_bulk_rhizo)^2/expected_bulk_rhizo+
  (observed_endo_rhizo - expected_endo_rhizo)^2/expected_endo_rhizo
pchisq(chisq, 2, lower.tail = FALSE)

