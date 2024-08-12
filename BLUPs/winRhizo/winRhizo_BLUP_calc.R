library(tidyverse)
library(lme4)
library(tools)
library(purrr)
library(corrr)
library(dplyr)
library(ggplot2)
library(visNetwork)
library(igraph)
library(ggraph)
library(ggcorrplot)


InputDataFileName <- "winRhizo_data.csv"
InputDataFileName.md5 <- md5sum(InputDataFileName)
print(InputDataFileName.md5) #put in README metadata file

winRhizo.data <- read.csv(InputDataFileName, na.strings=c("-")) %>%
  #mutate(sample..tub.plate..21.plates.per.slant. = sub("^[^_]*_", "", sample..tub.plate..21.plates.per.slant.)) %>%
  mutate(Genotype = sub("TIL:", "", Genotype))  %>%
  mutate(Genotype=as.factor(Genotype)) %>% #GENOTYPE is a factor
  mutate(across(LatRtAvLngth_RtSystm_mm:C_thickRtSA, as.numeric)) %>%
  mutate(across(LatRtAvLngth_RtSystm_mm:C_thickRtSA,  ~scale(., center = TRUE, scale = TRUE) %>% as.numeric)) %>%
  mutate(Check_TIL_OtherCV = as.factor(Check_TIL_OtherCV)) %>%
  mutate(X2.tub_SET = as.factor(X2.tub_SET)) %>%
  #mutate(across(where(is.numeric), ~ifelse(. < 0, NA, .))) %>%
  select(-LatRt_SizeClasses, -LatRt_MaxDiam) #Dropping these seems to let the hypothesis matrices to be saved

# identify outliers of traits using previous grouping technique
wizRhizo.outliers <- winRhizo.data %>%
  arrange(Genotype, X2.tub_SET) 
outdir <- "Outputs_08_11"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(wizRhizo.outliers, file=paste0(outdir,"/winRhizo_arrangedbyX2tub_set.csv"))


library(optimx) #additional optimizers for lmer/glmer functions
library(nloptr) #more optimizers
# ~  Check_TIL_OtherCV + (1|X2.tub_SET) + (1|Genotype)
lm.calc1 <- function(lmdataset, lmtrait){
  lm.formula <- 
    paste0(lmtrait,"~ Check_TIL_OtherCV + (1|X2.tub_SET) + (1|Check_TIL_OtherCV/Genotype) ") #Ridding of nest removes all errors for all traits + ABC only traits
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset, REML = TRUE, control = lmerControl(optCtrl = list(optimizer = "BFGS"))) #, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb')
  #lmer.result <- nlme::lme(as.formula(lm.formula), data=lmdataset)
  return(lmer.result)
}


# glm.calc1 <- function(glmdataset, glmtrait){
#   glm.formula <-
#     paste0(glmtrait," ~  Check_TIL_OtherCV + (1|X2.tub_SET) + (1|Check_TIL_OtherCV/Genotype)")
#   glmer.result <- glmer(as.formula(glm.formula), data = glmdataset, family = gaussian, control=glmerControl(optimizer="Nelder_Mead", optCtrl = list(maxfun = 2e5))) #,optCtrl=list(maxfun=2e5)
#   return(glmer.result)
# }


lm.calc2 <- function(lmdataset, lmtrait){
  lm.formula <-
    paste0(lmtrait," ~ Check_TIL_OtherCV + (1|X2.tub_SET) + (1|Check_TIL_OtherCV/Genotype)")
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset, REML = TRUE, control = lmerControl(optCtrl = list("L-BFGS-B", maxfun = 10000, ftol_abs = 1e-5))) #L-BFGS-B
  return(lmer.result)
}

#Select traits to perform analysis
traitlist <- winRhizo.data %>% dplyr::select(LatRtAvLngth_RtSystm_mm:C_thickRtSA) %>% colnames()
print(traitlist)

ABCtraitlist <- winRhizo.data %>% dplyr::select(B_LatLngthSum_cm:C_thickRtSA) %>% colnames()
print(ABCtraitlist)
#ABCtraitlist <- ABCtraitlist[-c(4, 5, 6)] #Drop TotRt_SA_cm2, TotRt_AvgDiam_cm, TotRt_Vol_cm3

lm1.traits.winRhizo <- traitlist %>% 
  set_names() %>% 
  map(~lm.calc1(winRhizo.data, .x))
lm1.ABCtraits.winRhizo <- ABCtraitlist %>% 
  set_names() %>% 
  map(~lm.calc1(winRhizo.data, .x))
#For the model comparison.
lm2.traits.winRhizo <- traitlist %>% 
  set_names() %>% 
  map(~lm.calc2(winRhizo.data, .x))

# glm1.traits.winRhizo <- traitlist %>% 
#   set_names() %>% 
#   map(~glm.calc1(winRhizo.data, .x))

blups.calc <- function(traitlms, traitname){
  blups <- ranef(traitlms[[traitname]])$Genotype 
  colnames(blups) <- c(traitname)  #rename column from "(intercept)" to traitname
  blups <- blups %>% mutate(Genotype = rownames(blups))
  return(blups)
}


BLUPS.traits.winRhizo.lm1 <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm1.traits.winRhizo, .x))
BLUPS.ABCtraits.winRhizo.lm1 <- ABCtraitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm1.ABCtraits.winRhizo, .x))

BLUPS.traits.winRhizo.lm2 <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm2.traits.winRhizo, .x))
# BLUPS.traits.winRhizo.glm1 <- traitlist %>% 
#   set_names() %>% 
#   map(~blups.calc(glm1.traits.winRhizo, .x))


#Create and save data frame for BLUPs
BLUPS.traits.winRhizo.df_lm1 <- as.list(BLUPS.traits.winRhizo.lm1) %>%
  reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by GENOTYPE
  dplyr::select(all_of(c("Genotype", traitlist))) #sets order of columns
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.traits.winRhizo.df_lm1, file=paste0(outdir,"/winRhizo.BLUPS_lm1.csv"))
#For the ABC traits
BLUPS.ABCtraits.winRhizo.df_lm1 <- as.list(BLUPS.ABCtraits.winRhizo.lm1) %>%
  reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by GENOTYPE
  dplyr::select(all_of(c("Genotype", ABCtraitlist))) #sets order of columns
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.ABCtraits.winRhizo.df_lm1, file=paste0(outdir,"/winRhizo.ABCBLUPS_lm1.csv"))

BLUPS.traits.winRhizo.df_lm2 <- as.list(BLUPS.traits.winRhizo.lm2) %>%
  reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by GENOTYPE
  dplyr::select(all_of(c("Genotype", traitlist))) #sets order of columns
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.traits.winRhizo.df_lm2, file=paste0(outdir,"/winRhizo.BLUPS_lm2.csv"))

# BLUPS.traits.winRhizo.df_glm1 <- as.list(BLUPS.traits.winRhizo.glm1) %>%
#   reduce(full_join, by="Genotype") %>% #join BLUPs of each trait by GENOTYPE
#   dplyr::select(all_of(c("Genotype", traitlist))) #sets order of columns
# outdir <- "Outputs_08_24"
# if (!dir.exists(outdir)){
#   dir.create(outdir)
# }
# write_csv(BLUPS.traits.winRhizo.df_glm1, file=paste0(outdir,"/winRhizo.BLUPS_glm.csv"))



#Remove checks
BLUPS.winRhizo.noCHECKS <- BLUPS.traits.winRhizo.df_lm1 %>%
  filter(!(Genotype %in% c("francis", "katy", "LMNT","TQNG")))
BLUPS.winRhizo.noCHECKS.lm2 <- BLUPS.traits.winRhizo.df_lm2 %>%
  filter(!(Genotype %in% c("francis", "katy", "LMNT","TQNG")))
ABCBLUPS.winRhizo.noCHECKS <- BLUPS.ABCtraits.winRhizo.df_lm1 %>%
  filter(!(Genotype %in% c("francis", "katy", "LMNT","TQNG")))

BLUPS.winRhizo.numONLY <- BLUPS.winRhizo.noCHECKS %>%
  select(-c(Genotype))
BLUPS.winRhizo.numONLY.lm2 <- BLUPS.winRhizo.noCHECKS.lm2 %>%
  select(-c(Genotype))
ABCBLUPS.winRhizo.numONLY <- ABCBLUPS.winRhizo.noCHECKS %>%
  select(-c(Genotype))

#Correlation calculation
BLUPcorr_winRhizo <- BLUPS.winRhizo.numONLY %>% 
  correlate() %>% 
  stretch()
BLUPcorr_winRhizo.lm2 <- BLUPS.winRhizo.numONLY.lm2 %>% 
  correlate() %>% 
  stretch()
ABCBLUPcorr_winRhizo <- ABCBLUPS.winRhizo.numONLY %>% 
  correlate() %>% 
  stretch()

#Filter edges based on correlation value and create graph object
graph_BLUPcorr_winRhizo <- BLUPcorr_winRhizo %>%
  filter(abs(r) > .9) %>%
  graph_from_data_frame(directed = FALSE)
graph_BLUPcorr_winRhizo.lm2 <- BLUPcorr_winRhizo.lm2 %>%
  filter(abs(r) > .9) %>%
  graph_from_data_frame(directed = FALSE)
graph_ABCBLUPcorr_winRhizo <- ABCBLUPcorr_winRhizo %>%
  filter(abs(r) > .9) %>%
  graph_from_data_frame(directed = FALSE)

#Simple example of ggraph function. Not necessary for analysis
ggraph(graph_BLUPcorr_winRhizo, layout = "kk") +
  geom_edge_link() +
  geom_node_point(repel = TRUE) +
  geom_node_text(aes(label = name))


#Make detailed node graph. This is for the full trait list.
corrNodegraph_BLUPS.winRhizo <- ggraph(graph_BLUPcorr_winRhizo) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = sqrt(abs(r)), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between root traits from winRhizo")

outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_BLUPS.winRhizo, width = 24, height = 20, file=paste0(outdir,"/correlationNodes_winRhizo.png"))

corrNodegraph_BLUPS.winRhizo.lm2 <- ggraph(graph_BLUPcorr_winRhizo.lm2) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = sqrt(abs(r)), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between root traits from winRhizo")

outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_BLUPS.winRhizo.lm2, width = 24, height = 20, file=paste0(outdir,"/correlationNodes_winRhizo_lm2.png"))

#Node graph for ABC
corrNodegraph_ABCBLUPS.winRhizo <- ggraph(graph_ABCBLUPcorr_winRhizo) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = sqrt(abs(r)), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between root traits from winRhizo")

outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_ABCBLUPS.winRhizo, width = 24, height = 20, file=paste0(outdir,"/ABCcorrelationNodes_winRhizo.png"))


#correlation calculation for ggcorrplot
BLUPS.winRhizo.df.NumONLY <- subset(BLUPS.winRhizo.noCHECKS, select = -c(Genotype) )#, LatRt_MaxDiam, LatRt_SizeClassesMay need to drop additional variables. 0's produce bad graphs!
winRhizo_corr <- round(cor(BLUPS.winRhizo.df.NumONLY, use = "pairwise.complete.obs"), 2)
wRBLUPp.mat <- cor_pmat(BLUPS.winRhizo.df.NumONLY, use = "pairwise.complete.obs")

ABCBLUPS.winRhizo.df.NumONLY <- subset(ABCBLUPS.winRhizo.noCHECKS, select = -c(Genotype) ) #May need to drop additional variables. 0's produce bad graphs!
winRhizo_ABCcorr <- round(cor(ABCBLUPS.winRhizo.df.NumONLY, use = "pairwise.complete.obs"), 2)
ABCwRBLUPp.mat <- cor_pmat(ABCBLUPS.winRhizo.df.NumONLY, use = "pairwise.complete.obs")

# graph of correlation
corrmap_winRhizo <-ggcorrplot(winRhizo_corr, method = "square",
                              title = "Correlation - BLUPS of elemental variables (F)",
                              legend.title = "Pearson Corr", lab = TRUE,
                              lab_col = "black",
                              lab_size = 3,
                              ggtheme = theme_minimal(),
                              colors = c("red2", "white", "dodgerblue"),
                              sig.level = .01
)
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrmap_winRhizo, width = 20, height = 20, file=paste0(outdir,"/winRhizocorr_map.png"))

#Save the correlation plots as a table of pairwise correlations
indicies <- which(upper.tri(winRhizo_corr), arr.ind = TRUE)
winRhizo_corrpairs <- data.frame(
  row = rownames(winRhizo_corr)[indicies[, 1]],
  col = rownames(winRhizo_corr)[indicies[, 2]],
  cor = winRhizo_corr[indicies]
)
names(winRhizo_corrpairs) <- c("Variable1", "Variable2", "Correlation")
winRhizo_corrpairs <- arrange(winRhizo_corrpairs, Variable1, Variable2)
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(winRhizo_corrpairs, file=paste0(outdir,"/winRhizo_corrpairs.csv"))




#For ABC traits. The correlation matrix
ABCcorrmap_winRhizo <-ggcorrplot(winRhizo_ABCcorr, method = "square",
                              title = "Correlation - BLUPS of elemental variables (F)",
                              legend.title = "Pearson Corr", lab = TRUE,
                              lab_col = "black",
                              lab_size = 3,
                              ggtheme = theme_minimal(),
                              colors = c("red2", "white", "dodgerblue"),
                              sig.level = .01
)
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(ABCcorrmap_winRhizo, width = 20, height = 20, file=paste0(outdir,"/winRhizoABCcorr_map.png"))


# correllogram for p-val. Crosses denote insignificance at default level .05
winRhizohypothesis_Graph <- ggcorrplot(winRhizo_corr,
                                       title = "Hypothesis Matrix for Flooded (p=.05)",
                                       method = "square",
                                       outline.color = "black",
                                       colors = c("red2", "white", "dodgerblue"),
                                       p.mat = wRBLUPp.mat)
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(winRhizohypothesis_Graph, width = 20, height = 20, file=paste0(outdir,"/winRhizohypothesis_Graph.png"))
#Hypotheisis matrix for ABC traits
winRhizoABChypothesis_Graph <- ggcorrplot(winRhizo_ABCcorr,
                                       title = "Hypothesis Matrix for Flooded (p=.05)",
                                       method = "square",
                                       outline.color = "black",
                                       colors = c("red2", "white", "dodgerblue"),
                                       p.mat = ABCwRBLUPp.mat)
outdir <- "Outputs_08_24"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(winRhizoABChypothesis_Graph, width = 20, height = 20, file=paste0(outdir,"/ABCwinRhizohypothesis_Graph.png"))


# identify outliers of traits using previous grouping technique
# winRhizo.outliers <- winRhizo.data %>%
#   arrange(Genotype, Check_TIL_OtherCV) %>%
#   mutate(is_outlier = ifelse(abs(RtDepth_cm) > 3, TRUE, FALSE))
# 
# # apply to entire traitlist
# for(trait in traitlist){
#   newcol <- paste0(trait, " outliers")
#   
#   winRhizo.data <- winRhizo.data %>%
#     arrange(X2.tub_SET, Genotype) %>%
#     mutate(!!newcol := ifelse(abs(!!sym(trait)) > 3, TRUE, FALSE))
# }
# # Create data frame of the logical vectors
# outlier_col <- paste0(traitlist, " outliers")
# winRhizo.outlier_df <- select(winRhizo.data, Genotype, Check_TIL_OtherCV, X2.tub_SET, all_of(outlier_col))
# outdir <- "Outputs_08_24"
# if (!dir.exists(outdir)){
#   dir.create(outdir)
# }
# write_csv(winRhizo.outlier_df, file=paste0(outdir,"/winRhizo_outliers.df.csv"))


# Scaled values >3 st. deviations away are classified as TRUE. Done for all traits
for(trait in traitlist){
  newcol <- paste0(trait, " outliers")
  
  winRhizo.data <- winRhizo.data %>%
    arrange(Genotype, X2.tub_SET)  %>%
    mutate(!!newcol := ifelse(abs(!!sym(trait)) > 3, TRUE, FALSE))
}
# Create data and save data frame consisting of only the logical outlier columns
outlier_col <- paste0(traitlist, " outliers")
winRhizo.outliers_df <- dplyr::select(winRhizo.data, Genotype, X2.tub_SET, tub, Check_TIL_OtherCV, all_of(outlier_col))
outdir <- "Outputs_08_11"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(winRhizo.outliers_df, file=paste0(outdir,"/winRhizo.outliers_df.csv"))

