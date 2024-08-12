library(tidyverse)
library(purrr) #map function
library(tools) #needed for md5sum function
library(lme4) #mixed linear models
library(lmerTest) #to get p-values from summary() and anova()
library(corrr) #calculate correlation
library(ggraph) #used for graph visualization
library(igraph) #more graph tools
library(ggcorrplot) #map of correlation between variables
library(dplyr)
# Used previously to help determine best model
library(car) 
library(MASS)
library(lmtest)
library(performance)



InputDataFileName <- "TILs_6wk_TN_RB_SBetc_Data&ANOVA_UPDATED.csv"
InputDataFileName.md5 <- md5sum(InputDataFileName)
print(InputDataFileName.md5) #put in README metadata file

RootData <- read.csv(InputDataFileName, na.strings=c("-")) %>%
  mutate_all(trimws) %>% #strip whitespace from everything
  mutate(Rep_GH2019_1to3=as.factor(Rep_GH2019_1to3)) %>% #Year is a factor
  mutate(Ranking_from_cooling_pad=as.numeric(Ranking_from_cooling_pad)) %>% #Rep is a factor
  mutate(Genotype=as.factor(Genotype)) %>% #Genotype is a factor
  mutate(across(TN_sixwk:LfLngth_TipHt_Minus_CollarHt, as.numeric)) %>% #response variables set to numeric
  mutate(across(TN_sixwk:LfLngth_TipHt_Minus_CollarHt, ~scale(., center = TRUE, scale = TRUE) %>% as.numeric)) %>% #scale to get Z-scores
  filter(!(Genotype %in% c("Katy", "FRNC", "PI312777", "Rondo", "TQNG"))) # Drop checks


traitlist <- RootData %>% dplyr::select(TN_sixwk:LfLngth_TipHt_Minus_CollarHt) %>% colnames()


# identify outliers of traits using previous grouping technique
Rootdata.outliers <- RootData %>%
  arrange(Genotype, Rep_GH2019_1to3, Ranking_from_cooling_pad) 
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(Rootdata.outliers, file=paste0(outdir,"/Rootdata_arrangedbyRep_GH2019_1to3_CoolingPad.csv"))

# By trait, if value is x>3 st. deviations then put TRUE to indicate outlier
for(trait in traitlist){
  newcol <- paste0(trait, " outliers")
  
  RootData <- RootData %>%
    arrange(Genotype, Ranking_from_cooling_pad) %>%
    mutate(!!newcol := ifelse(abs(!!sym(trait)) > 3, TRUE, FALSE))
}
# Create and save data frame of logical outlier list
outlier_col <- paste0(traitlist, " outliers")
Rootdata.outliers_df <- dplyr::select(RootData, Genotype, Ranking_from_cooling_pad, all_of(outlier_col))
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(Rootdata.outliers_df, file=paste0(outdir,"Root_outliers_ByCoolPad.csv"))



lm1.calc <- function(lmdataset, lmtrait){
  lm.formula <- 
    paste0(lmtrait," ~ Ranking_from_cooling_pad + Rep_GH2019_1to3 + (1|Genotype)")
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset)
  return(lmer.result)
}

lm2.calc <- function(lmdataset, lmtrait){
  lm.formula <- 
    paste0(lmtrait," ~ Ranking_from_cooling_pad + Rep_GH2019_1to3 + (Rep_GH2019_1to3 * Ranking_from_cooling_pad) + (1|Genotype)")
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset)
  return(lmer.result)
}



#traitlist <- RootData %>% dplyr::select(TN_sixwk:LfLngth_TipHt_Minus_CollarHt) %>% colnames()
#print(traitlist)

lm1.traits.ROOT <- traitlist %>% 
  set_names() %>% 
  map(~lm1.calc(RootData, .x))

lm2.traits.ROOT <- traitlist %>% 
  set_names() %>% 
  map(~lm2.calc(RootData, .x))

blups.calc <- function(traitlms, traitname){
  blups <- ranef(traitlms[[traitname]])$Genotype 
  colnames(blups) <- c(traitname)  #rename column from "(intercept)" to traitname
  blups <- blups %>% mutate(Genotype = rownames(blups))
  return(blups )
}


BLUPS.traits.ROOTS.lm1 <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm1.traits.ROOT, .x))

BLUPS.traits.ROOTS.lm2 <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm2.traits.ROOT, .x))

#For the first model: create and save file
BLUPS.ROOTS.df.lm1 <- as.list(BLUPS.traits.ROOTS.lm1) %>%
  reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by Genotype
  dplyr::select(all_of(c("Genotype", traitlist))) #sets order of columns
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.ROOTS.df.lm1, file=paste0(outdir,"/ROOTBLUPS.lm1.csv"))

#For the second model: create and save file
BLUPS.ROOTS.df.lm2 <- as.list(BLUPS.traits.ROOTS.lm2) %>%
  reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by Genotype
  dplyr::select(all_of(c("Genotype", traitlist))) #sets order of columns
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.ROOTS.df.lm2, file=paste0(outdir,"/ROOTBLUPS.lm2.csv"))


performance_Lfnu <- model_performance(lm1.traits.ROOT$LfNu)
print(performance_Lfnu)

calculate_performance <- function(trait_models) {
  map(trait_models, ~model_performance(.x))
}
performance.lm1.traits.ROOT <- calculate_performance(lm1.traits.ROOT)
create_summary <- function(performance_list, trait_names) {
  performance_df <- map2_dfr(performance_list, trait_names, ~as.data.frame(.x) %>% mutate(Trait = .y))
  return(performance_df)
}
performance_summary <- create_summary(performance.lm1.traits.ROOT, traitlist)
write.csv(performance_summary, file = paste0(outdir, "/Model_performance_summary_lm1.csv"))


# Drop checks to produce correlation graphs
BLUPS.Root.lm1.df.NOchecks <- BLUPS.ROOTS.df.lm1 %>%
  filter(!(Genotype %in% c("FRNC", "Katy", "PI312777", "LMNT","Rondo")))

BLUPS.Root.lm2.df.NOchecks <- BLUPS.ROOTS.df.lm2 %>%
  filter(!(Genotype %in% c("FRNC", "Katy", "PI312777", "LMNT","Rondo")))

#plot correlations between element BLUPs among TILs
#Start of calculation and graphs used for flooded data set
#
#
BLUPcorr_Root_lm1 <- BLUPS.Root.lm1.df.NOchecks %>% 
  correlate() %>% 
  stretch()

graph_BLUPcorr_Root.lm1 <- BLUPcorr_Root_lm1 %>%
  filter(abs(r) > .3) %>%
  graph_from_data_frame(directed = FALSE)

# plots simple correlation graph - only names and connection
ggraph(graph_BLUPcorr_Root.lm1) +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name))

# graph with more detail - includes color, nodes, and edge properties
corrNodegraph_BLUPS.Root.lm1 <- ggraph(graph_BLUPcorr_Root.lm1) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = sqrt(abs(r)), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between root traits (mixed model 1)")
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_BLUPS.Root.lm1, width = 10, height = 5, file=paste0(outdir,"/corrNodegraph_BLUPS.Root.lm1.png"))

# correlation calculation for ggcorrplot
BLUPS.Root.lm1.df.NumONLY <- subset(BLUPS.Root.lm1.df.NOchecks, select = -c(Genotype) )
Root_lm1_corr <- round(cor(BLUPS.Root.lm1.df.NumONLY, use = "pairwise.complete.obs"), 2)
Root_lm1BLUPp.mat <- cor_pmat(BLUPS.Root.lm1.df.NumONLY, use = "pairwise.complete.obs")
# graph of correlation 
corrmap_Root.lm1 <-ggcorrplot(Root_lm1_corr, method = "square",
                             title = "Correlation - BLUPS of elemental variables (F)",
                             legend.title = "Pearson Corr", lab = TRUE,
                             lab_col = "black",
                             lab_size = 1.75,
                             ggtheme = theme_minimal(),
                             colors = c("red2", "white", "dodgerblue")
)
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrmap_Root.lm1, width = 6, height = 6, file=paste0(outdir,"/CorrMap_Root_lm1.png"))
corrmap_Root.lm1

#Save the correlation plots as a table of pairwise correlations
indicies <- which(upper.tri(Root_lm1_corr), arr.ind = TRUE)
Root_lm1_corrpairs <- data.frame(
  row = rownames(Root_lm1_corr)[indicies[, 1]],
  col = rownames(Root_lm1_corr)[indicies[, 2]],
  cor = Root_lm1_corr[indicies]
)
names(Root_lm1_corrpairs) <- c("Variable1", "Variable2", "Correlation")
Root_lm1_corrpairs <- arrange(Root_lm1_corrpairs, Variable1, Variable2)
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(Root_lm1_corrpairs, file=paste0(outdir,"/Root_lm1_corrpairs.csv"))

# correllogram for p-val. Crosses denote insignificance at default level .05
Root_lm1hypothesis_Graph <- ggcorrplot(Root_lm1_corr,
                                    title = "Hypothesis Matrix for Flooded (p=.05)",
                                    method = "square",
                                    outline.color = "black",
                                    colors = c("red2", "white", "dodgerblue"),
                                    p.mat = Root_lm1BLUPp.mat)
outdir <- "Outputs_09-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(Root_lm1hypothesis_Graph, width = 6, height = 6, file=paste0(outdir,"/Root_lm1hypothesis_Graph.png"))

