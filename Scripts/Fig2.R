### Code for Figure 2
### By Isabelle Du Plessis, Hannah Snyder, and Marian Dominguez-Mirazo

#load packages
library(ggplot2) # version 3.4.2
library(ggpattern) # version 1.0.1
library(ggpubr) # version 0.6.0
library(grid) # version 4.2.1
library(gridExtra) # version 2.3
library(dplyr) # version 1.1.3

# R version 4.0.2 (2020-06-22)


#Figure 2
#load in votu metadata
votus_cov75thres = read.table("../intermediate_files/votus_cov75thres_1224.txt",header = T)

######################################
########FIGURE 2 a ###################

#calculate the number of times a votu is repeated in the dataframe, save in occurrence column
votus <-votus_cov75thres %>%
  group_by(votus) %>%
  dplyr::mutate(Occurrence = n())
votus = votus[,c("sample","votus","Occurrence")] # Keep only necessary data
colnames(votus) = c("sampleID","vOTU_ID","Occurrence")

#use the distinct combinations of vOTU_ID and occurrence to make a barplot using ggplot2
distinct_votus = votus %>% distinct(vOTU_ID,Occurrence, .keep_all = FALSE)

#here we make a list of the votu_ids that occur more than one time
multiple_occurrence = (distinct_votus %>% filter(Occurrence > 1))$vOTU_ID

#here we start making the bargraph to display the number of times a votu occurs
distinct_votus$ColorID = "Multiple"
distinct_votus$ColorID[distinct_votus$Occurrence == 1] <- "One"
distinct_votus = distinct_votus[,-1]
distinct_votus2 <-distinct_votus %>%
  group_by(Occurrence) %>%
  dplyr::mutate(Count = n())
distinct_votus2 <- distinct_votus2%>% distinct()

Figure2a <- ggplot(distinct_votus2, aes(x=Occurrence, y=Count, fill=ColorID)) +
  geom_bar(position='stack', stat='identity') + theme_test() + 
  geom_text(aes(y= Count+12,label=Count),size=5) +
  scale_fill_manual(values = c("One" = "gray", "Multiple" = "black")) +
  scale_x_continuous(limits = c(0.3,24.7),breaks = c(1,5,10,15,20,24),expand = c(0,0)) + 
  scale_y_continuous(limits =c(0,350),expand = c(0,0))+
  theme(legend.position = "none")  + 
  xlab("Sample Occurrence") + 
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank()
        ) +
  ylab("vOTU count")


##########################################
####FIGURE 2 b/c  Organization ###########

#we next made a dataframe including just the vOTU_IDs that occur more than once,
#and put this into donut plots separated by phenotype and compartment
#next, make a list of the sample_ids in order that correspond to every votu
Compartment = c()
Spartina = c()
for (i in votus$sampleID){
  if (i %in% (c("sample_1","sample_7","sample_13","sample_19"))){
    Spartina = append(Spartina, "Tall Height")
    Compartment = append(Compartment, "Bulk Sediment")
  }
  if (i %in%(c("sample_2","sample_8","sample_14","sample_20"))){
    Spartina = append(Spartina, "Tall Height")
    Compartment = append(Compartment, "Rhizosphere")  
  }
  if (i %in%(c("sample_3","sample_9","sample_15","sample_21"))){
    Spartina = append(Spartina, "Tall Height")
    Compartment = append(Compartment, "Root")
  }
  if (i %in%c("sample_4","sample_10","sample_16","sample_22")){
    Spartina = append(Spartina, "Short Height")
    Compartment = append(Compartment, "Bulk Sediment")
  }
  if (i %in%(c("sample_5","sample_11","sample_17","sample_23"))){
    Spartina = append(Spartina, "Short Height")
    Compartment = append(Compartment, "Rhizosphere")
  }
  if (i %in% c("sample_6","sample_12","sample_18","sample_24")){
    Spartina = append(Spartina, "Short Height")
    Compartment = append(Compartment, "Root")
  }
}

# Load the null model output
null_model = read.table("../intermediate_files/nullmodel.txt",row.names = 1)

####################################################
##################### FIGURE 2 b ###################

#make a separate dataframe for just the votuID and compartment data
#then only keep the votus that had an occurrence value > 1
all_phenotype_df = data.frame(votus$vOTU_ID, Spartina)
tall = all_phenotype_df %>% filter(Spartina == "Tall Height")
short = all_phenotype_df %>% filter(Spartina == "Short Height")
phenotype_df = all_phenotype_df %>% filter(votus.vOTU_ID %in% multiple_occurrence)
unique_df = (phenotype_df %>% distinct())

#make a column that includes the number of times each is repeated
df2 <- unique_df %>%
  group_by(votus.vOTU_ID) %>%
  dplyr::mutate(Occurrence = n())
df2 <- df2[,-2]
df2 = (df2 %>% distinct())

#Save votus names with two phenotypes
twopheno_names = df2$votus.vOTU_ID[df2$Occurrence==2]
onepheno_df = unique_df[!(unique_df$votus.vOTU_ID%in%twopheno_names),]

#count up the number of votus that occur in only one phenotype and only two phenotypes
#make a dataframe with the data for these, and create a donut plot for the phenotypic data
tall_phenotype = sum(onepheno_df$Spartina=="Tall Height")
short_phenotype = sum(onepheno_df$Spartina=="Short Height")
two_phenotype= length(twopheno_names)

data_observed_pheno = data.frame(
  "Phenotype" = c("Tall", "Short","Both"), 
  "Count"=c(tall_phenotype,short_phenotype,two_phenotype),
  "Model"="Observed"
) 

# Copy the structure for the null model
data_expected_pheno = data.frame("Phenotype" = c("Tall","Short","Both"),
                           "Count" = c(sum(null_model[,3]),sum(null_model[,2]),sum(null_model[,1])),
                           "Model" = "Expected")


data_pheno = rbind(data_expected_pheno,data_observed_pheno)
data_pheno$Phenotype = factor(data_pheno$Phenotype,levels=c("Tall","Short","Both"))
data_pheno$Model = factor(data_pheno$Model,levels=c("Observed","Expected")) 


Figure2b_pheno = ggplot() +
  geom_col(data=data_pheno,aes(y=Count,x=Model,fill=Phenotype),color="black") +
  scale_fill_manual(values = c("Both"="#808080", "Short"="darkseagreen3","Tall"="darkseagreen4")) + 
  theme_minimal() +
  theme(axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15)) +
  scale_y_continuous(breaks=c(seq(0,400,100),444)) +
  ylab("vOTU count") +
  xlab("")

pheno_legend = get_legend(Figure2b_pheno)
Figure2b_pheno_nolegend = Figure2b_pheno + theme(legend.position = "none")

####################################################
##################### FIGURE 2 c ###################

# Expected data
data_expected_comp = data.frame("Compartment" = c("Bulk sediment","Rhizosphere","Root","BRh","BRo","RhRo","All compartments"),
                           "Count" = c(sum(null_model[1,]),
                                       sum(null_model[7,]),
                                       sum(null_model[5,]),
                                       sum(null_model[4,]),
                                       sum(null_model[2,]),
                                       sum(null_model[6,]),
                                       sum(null_model[3,])),
                           "Model"="Expected")


######## Data compartment counting
# starting on the compartment data here
all_compartent_df = data.frame( votus$vOTU_ID, Compartment)
compartent_df = all_compartent_df %>% filter(votus.vOTU_ID %in% multiple_occurrence)
#make a separate dataframe for just the votuID and compartment data
#make a dataframe of just the votuID and Compartment combinations that are unique
unique_df = (compartent_df %>% distinct())

#make a column that includes the number of times each is repeated:
df2 <-unique_df %>%
  group_by(votus.vOTU_ID) %>%
  dplyr::mutate(Occurrence = n())
df2 <- df2[,-2]
df2 = (df2 %>% distinct())

#count up the number of votus that occur in one, two, and three compartments
#make a dataframe with the data for these, and create a donut plot for the compartment data
one_compartment = sum(df2$Occurrence == "1") 
two_compartment = sum(df2$Occurrence == "2") 
three_compartment = sum(df2$Occurrence == "3") 

## vOTUs that occur in 1 Compartment
votu_one_comp = df2 %>% filter(Occurrence == 1)
comp_one = all_compartent_df %>% filter(votus.vOTU_ID %in% votu_one_comp$votus.vOTU_ID)
comp_one = (comp_one %>% distinct())
rhi = comp_one %>% filter(Compartment == "Rhizosphere")
bulk = comp_one %>% filter(Compartment == "Bulk Sediment")
end = comp_one %>% filter(Compartment == "Root")
rhi_number = length(unique(rhi$votus.vOTU_ID))
end_number = length(unique(end$votus.vOTU_ID))
bulk_number= length(unique(bulk$votus.vOTU_ID))


## vOTUs that occur in 2 Compartments
votu_two_comp = df2 %>%  filter(Occurrence == 2)
two_list = votu_two_comp$votus.vOTU_ID
df_two = all_compartent_df %>%  filter(votus.vOTU_ID %in% two_list)

bulk_rhi = df_two %>% filter(Compartment == "Bulk Sediment" | Compartment == "Rhizosphere")
bulk_rhi = (bulk_rhi %>% distinct())
bulk_rhi = bulk_rhi %>%
  group_by(votus.vOTU_ID) %>%
  dplyr::mutate(Occurrence = n())
bulk_rhi_two = bulk_rhi %>% filter(Occurrence == 2)
bulk_rhi_two_number = length(unique(bulk_rhi_two$votus.vOTU_ID))

end_rhi = df_two %>% filter(Compartment == "Root" | Compartment == "Rhizosphere")
end_rhi = (end_rhi %>% distinct())
end_rhi = end_rhi %>%
  group_by(votus.vOTU_ID) %>%
  dplyr::mutate(Occurrence = n())
end_rhi_two = end_rhi %>% filter(Occurrence == 2)
end_rhi_two_number = length(unique(end_rhi_two$votus.vOTU_ID))

end_bulk = df_two %>% filter(Compartment == "Root" | Compartment == "Bulk Sediment")
end_bulk = (end_bulk %>% distinct())
end_bulk = end_bulk %>%
  group_by(votus.vOTU_ID) %>%
  dplyr::mutate(Occurrence = n())
end_bulk_two = end_bulk %>% filter(Occurrence == 2)
end_bulk_two_number = length(unique(end_bulk_two$votus.vOTU_ID))

# observed data frame structure for donut plot
data_observed_comp = data.frame("Compartment" = c("Bulk sediment","Rhizosphere","Root","BRh","BRo","RhRo","All compartments"),
                                "Count" = c(bulk_number,
                                            rhi_number,
                                            end_number,
                                            bulk_rhi_two_number,
                                            end_bulk_two_number,
                                            end_rhi_two_number,
                                            three_compartment),
                                "Model" = "Observed")

data_comp = rbind(data_expected_comp,data_observed_comp)
data_comp$Compartment = factor(data_comp$Compartment,levels = c("Bulk sediment","Rhizosphere","Root","BRh","RhRo","BRo","All compartments"))
data_comp$Model = factor(data_comp$Model, levels=c("Observed","Expected"))

Figure2c_comp = ggplot() +
  geom_col_pattern(data=data_comp,aes(y=Count,x=Model,
                                      fill=Compartment, 
                                      pattern_fill = Compartment,
                                      pattern_color=Compartment),
                    pattern = "stripe", 
                    colour = 'black',
                    pattern_angle = 45,
                    pattern_density = 0.5,
                    pattern_size = 0) +
  scale_fill_manual(values = c("Bulk sediment"="#853512", "Rhizosphere"="#EEAA23","Root"="#558A78","BRh"="#853512","BRo"="#853512","RhRo"="#EEAA23","All compartments"="#404040")) +
  scale_pattern_fill_manual(values = c("Bulk sediment"="#853512", "Rhizosphere"="#EEAA23","Root"="#558A78","BRh"="#EEAA23","BRo"="#558A78","RhRo"="#558A78","All compartments"="#404040")) +
  scale_pattern_color_manual(values = c("Bulk sediment"="#853512", "Rhizosphere"="#EEAA23","Root"="#558A78","BRh"="#853512","BRo"="#853512","RhRo"="#EEAA23","All compartments"="#404040")) +
  theme_minimal() +
  theme(legend.position="none",
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        panel.grid.major.x = element_blank()) +
  scale_y_continuous(breaks=c(seq(0,400,100),444)) +
  ylab("vOTU count") +
  xlab("")
  
Figure2b_pheno_nolegend = Figure2b_pheno_nolegend + theme(plot.margin=unit(c(2,0,0,0.5), "cm"))
Figure2c_comp = Figure2c_comp + theme(plot.margin=unit(c(2,0,0,0.5), "cm"))

for_legend_df = data_comp[data_comp$Compartment%in%c("Bulk sediment","Rhizosphere","Root","All compartments"),]
for_legend_df$Compartment = factor(for_legend_df$Compartment,levels=c("Bulk sediment","Rhizosphere","Root","All compartments"))

for_legend = ggplot() +
  geom_col(data=for_legend_df,aes(y=Count,x=Model,fill=Compartment),color="#404040") +
  scale_fill_manual(values = c("Bulk sediment"="#853512", "Rhizosphere"="#EEAA23","Root"="#558A78","All compartments"="#404040")) +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=12))

comp_legend = get_legend(for_legend)


pdf("../Figures/Figure2.pdf",height=11,width=12)
grid.arrange(Figure2a,Figure2b_pheno_nolegend,pheno_legend,Figure2c_comp,comp_legend,
             layout_matrix=matrix(c(1,1,1,1,2,3,4,5),ncol=4,byrow=T),
             widths=c(0.35,0.15,0.35,0.15),heights=c(0.35,0.65))
dev.off()

