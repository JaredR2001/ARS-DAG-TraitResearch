library(tidyverse)
library(purrr) #map function
library(tools) #needed for md5sum function
library(lme4) #mixed linear models
library(lmerTest) #to get p-values from summary() and anova()
library(corrr) #calculate correlation
library(ggraph) #used for graph visualization
library(igraph) #more graph tools
library(ggcorrplot) #map of correlation between variables
library(scales)


InputDataFileName <- "TILIonomicsData.csv"
InputDataFileName.md5 <- md5sum(InputDataFileName)
print(InputDataFileName.md5) #put in README metadata file

Ionomicsdata <- read.csv(InputDataFileName, na.strings=c("-")) %>%
  select(-X) %>% #drop extra empty column
  mutate_all(trimws) %>% #strip white space from everything
  mutate(YEAR=as.factor(YEAR)) %>% #Year is a factor
  mutate(REP_Irig_nu=as.factor(REP_Irig_nu)) %>% #Rep is a factor
  mutate(GENOTYPE=as.factor(GENOTYPE)) %>% #GENOTYPE is a factor
  mutate(Irrigation_Treatment = str_replace(Irrigation_Treatment, "flooded", "Flooded")) %>%
  #mutate(Ni = ifelse(Ni > 1, NA, Ni)) %>% # remove high Ni values only if including checks
  mutate(across(DHD:Zn, as.numeric)) %>% #response variables set to numeric
  filter(!(GENOTYPE %in% c("SABR", "DG100", "CPRS", "TQNG"))) #Drop checks for BLUP calculations

# list genotypes
Ionomicsdata$GENOTYPE %>% unique %>% sort

# split by flood and unflooded
Unflood_data <- filter(Ionomicsdata, Irrigation_Treatment == "UnFlooded") 
Unflood_data <- Unflood_data %>%
  mutate(across(DHD:Zn, ~scale(., center = TRUE, scale = TRUE) %>% as.numeric)) %>%
  mutate(Rb = ifelse(Rb < -2, NA, Rb)) %>% # remove high Rb values, if not including checks...adding leads to missing BLUPs for unflooded
  mutate(Rb = ifelse(Rb > 2, NA, Rb)) %>%
  mutate(Fe = ifelse(Fe < -4, NA, Fe))
#boxplot(Unflood_data$Rb)

Flood_data <- filter(Ionomicsdata, Irrigation_Treatment == "Flooded") 
Flood_data <- Flood_data %>%
  mutate(across(DHD:Zn, ~scale(., center = TRUE, scale = TRUE) %>% as.numeric)) #scale to get Z-scores


lm.calc <- function(lmdataset, lmtrait){
  lm.formula <- 
    paste0(lmtrait," ~ YEAR + REP_Irig_nu + (1|GENOTYPE)")
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset)
  return(lmer.result)
}

# function to run the mixed linear model: specify data set and trait name.
lm.calc.noyear <- function(lmdataset, lmtrait){
  lm.formula <- 
    paste0(lmtrait," ~ REP_Irig_nu + (1|GENOTYPE)")
  lmer.result <- lmer(as.formula(lm.formula), data=lmdataset)
  return(lmer.result)
}
# Example of how to use the above function individually
example.lm.calc <- lm.calc.noyear(Flood_data, "As")

# get list of traits to pipe into furrr map functions
traitlist <- Ionomicsdata %>% dplyr::select(DHD:Zn) %>% colnames()
print(traitlist)

# get a dataframe with lmer results for each trait piped into the map function
lm.traits.flood <- traitlist %>% 
  set_names() %>% 
  map(~lm.calc.noyear(Flood_data, .x)) #.x a the trait name from piped list

lm.traits.unflood <- traitlist %>% 
  set_names() %>% 
  map(~lm.calc.noyear(Unflood_data, .x))

# how to access lmer results of a trait
summary(lm.traits.flood$As)
anova(lm.traits.flood$As)
# alternative way to access by trait
lm.traits.flood[["As"]]

# run BLUPs for an individual trait 
BLUPS.As <- ranef(lm.traits.flood$As)$GENOTYPE

# function to calculate BLUPs:specify lmer results dataframe and a trait name
blups.calc <- function(traitlms, traitname){
  blups <- ranef(traitlms[[traitname]])$GENOTYPE 
  colnames(blups) <- c(traitname)  #rename column from "(intercept)" to traitname
  blups <- blups %>% mutate(GENOTYPE = rownames(blups))
  return(blups )
}
# Example of how to use the above function individually
example.blups <- blups.calc(lm.traits.flood, "As")


# get a dataframe with BLUP results for each trait piped into the map function
BLUPS.traits.flood <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm.traits.flood, .x))
# also run for unflooded
BLUPS.traits.unflood <- traitlist %>% 
  set_names() %>% 
  map(~blups.calc(lm.traits.unflood, .x))

# Example of accessing BLUPs from dataframe
#BLUPS.traits.flood$As

# put BLUPs from all traits together as columns
BLUPS.flood.df <- as.list(BLUPS.traits.flood) %>%
  reduce(full_join, by="GENOTYPE") %>% #join BLUPs of each trait by GENOTYPE
  dplyr::select(all_of(c("GENOTYPE", traitlist))) #sets order of columns

outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.flood.df, file=paste0(outdir,"/FloodBLUPS.csv"))

# also run for unflooded
BLUPS.unflooded.df <- as.list(BLUPS.traits.unflood) %>%
  reduce(full_join, by="GENOTYPE") %>% #join BLUPs of each trait by GENOTYPE
  dplyr::select(all_of(c("GENOTYPE", traitlist))) #sets order of columns

outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.unflooded.df, file=paste0(outdir,"/UnfloodedBLUPS.csv"))

# Imputing step for Strontium values
imputed.BLUPS.Sr <- BLUPS.unflooded.df # Save BLUPs in new .df so we can manipulate it

# Model creation using only complete values in BLUPS.unflooded.df
Sr.BLUP.MODEL <- lm(Sr ~ 1 + Ca + Mg, data = na.omit(BLUPS.unflooded.df))
# Predicts the missing values in Sr
imputed.BLUPS.Sr$Sr[is.na(imputed.BLUPS.Sr$Sr)] <- predict(Sr.BLUP.MODEL, newdata = imputed.BLUPS.Sr)[is.na(imputed.BLUPS.Sr$Sr)]

# Logical vector of non-empty rows 
nonempty_rows <- !is.na(imputed.BLUPS.Sr$Sr)
# Copy over non-empty entries to data frame - only for unflooded (flooded field had all Sr BLUPs)
imputed.BLUPS.predicted <- imputed.BLUPS.Sr[nonempty_rows, ]
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write.csv(imputed.BLUPS.predicted, file=paste0(outdir,"/imputedSr.unflooded.BLUPS.csv"))

# Step to impute DHD values. Sr completed and saved
# Create temporary df for manipulation of flooded and unflooded BLUPS
impute.BLUPS.temp.u <- imputed.BLUPS.predicted # Most recent unflood BLUPs 
impute.BLUPS.temp.f <- BLUPS.flood.df 

#Creating model for DHD and adding names to columns
DHD.BLUP.MODEL.j <- lmer(DHD ~  Irrigation_Treatment + REP_Irig_nu + (1|GENOTYPE), data = Ionomicsdata)


# Get the random effects from model 
all.DHD.BLUPS <- ranef(DHD.BLUP.MODEL.j)$GENOTYPE # BLUPs calculated using the entire ionomics data set
colnames(all.DHD.BLUPS)[1] <- "DHD"
all.DHD.BLUPS <- all.DHD.BLUPS %>%
  mutate(GENOTYPE = rownames(.))


# Unflooded
# Apply new DHD BLUPS to DHD column: Unflooded. Added imputed BLUPS file.
impute.BLUPS.temp.u$DHD <- all.DHD.BLUPS$DHD[match(impute.BLUPS.temp.u$GENOTYPE, all.DHD.BLUPS$GENOTYPE)]

# Change name of DHD column
colnames(impute.BLUPS.temp.u)[colnames(impute.BLUPS.temp.u) == "DHD"] <- "DHD.j"
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(impute.BLUPS.temp.u, file=paste0(outdir,"/UnfloodedBLUPS.imputed.DHDj.csv"))

#Apply new BLUPs to NAs in DHD column
imputed.BLUPS.Sr$DHD <- ifelse(is.na(imputed.BLUPS.Sr$DHD), all.DHD.BLUPS$DHD[match(imputed.BLUPS.Sr$GENOTYPE, all.DHD.BLUPS$GENOTYPE)], imputed.BLUPS.Sr$DHD)
colnames(imputed.BLUPS.Sr)[colnames(imputed.BLUPS.Sr) == "DHD"] <- "DHD.u"
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(imputed.BLUPS.Sr, file=paste0(outdir,"/UnfloodedBLUPS.imputed.DHDu.csv"))



# Apply new BLUPs to entire DHD column. DHD was calculated using entire Ionomics sheet.
impute.BLUPS.temp.f$DHD <- all.DHD.BLUPS$DHD[match(impute.BLUPS.temp.f$GENOTYPE, all.DHD.BLUPS$GENOTYPE)]

# Rename DHD column
colnames(impute.BLUPS.temp.f)[colnames(impute.BLUPS.temp.f) == "DHD"] <- "DHD.j"
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(impute.BLUPS.temp.f, file=paste0(outdir,"/FloodedBLUPS.imputed.DHDj.csv"))

# Impute the DHD BLUPs but only to the NA entries
BLUPS.flood.df$DHD <- ifelse(is.na(BLUPS.flood.df$DHD), all.DHD.BLUPS$DHD[match(BLUPS.flood.df$GENOTYPE, all.DHD.BLUPS$GENOTYPE)], BLUPS.flood.df$DHD)
colnames(BLUPS.flood.df)[colnames(BLUPS.flood.df) == "DHD"] <- "DHD.f"
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(BLUPS.flood.df, file=paste0(outdir,"/FloodedBLUPS.imputed.DHDf.csv"))


# remove check varieties. Not needed after removing checks from BLUP calc. 
BLUPS.flood.df.NOchecks <- BLUPS.flood.df %>%
  filter(!(GENOTYPE %in% c("SABR", "DG100", "CPRS", "LMNT","TQNG")))

BLUPS.unflooded.df.NOchecks <- imputed.BLUPS.Sr %>%
  filter(!(GENOTYPE %in% c("SABR", "DG100", "CPRS", "LMNT","TQNG")))



# plot correlations between element BLUPs among TILs
# Start of calculation and graphs used for flooded data set
#
#
BLUPcorr_Flood <- BLUPS.flood.df.NOchecks %>% 
  correlate() %>% 
  stretch()

graph_BLUPcorr_Flood <- BLUPcorr_Flood %>%
  filter(abs(r) > .3) %>%
  graph_from_data_frame(directed = FALSE)

# plots simple correlation graph - only names and connection
ggraph(graph_BLUPcorr_Flood) +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name))

# graph with more detail - includes color, nodes, and edge properties
corrNodegraph_BLUPS.flood <- ggraph(graph_BLUPcorr_Flood) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = sqrt(abs(r)), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between elements (flooded field)")
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_BLUPS.flood, width = 10, height = 5, file=paste0(outdir,"/Nodegraph_Flooded_08-12.png"))

# correlation calculation for ggcorrplot
BLUPS.flood.df.NumONLY <- subset(BLUPS.flood.df.NOchecks, select = -c(GENOTYPE) )
Flooded_corr <- round(cor(BLUPS.flood.df.NumONLY, use = "pairwise.complete.obs"), 2)
FBLUPp.mat <- cor_pmat(BLUPS.flood.df.NumONLY, use = "pairwise.complete.obs")
# graph of correlation. the correlation matrix
corrmap_Flooded <-ggcorrplot(Flooded_corr, method = "square",
                             title = "Correlation - BLUPS of elemental variables (F)",
                             legend.title = "Pearson Corr", lab = TRUE,
                             lab_col = "black",
                             lab_size = 1.75,
                             ggtheme = theme_minimal(),
                             colors = c("red2", "white", "dodgerblue")
)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrmap_Flooded, width = 6, height = 6, file=paste0(outdir,"/CorrMap_Flooded_08-12.png"))
corrmap_Flooded

#Save the correlation plots as a table of pairwise correlations
indicies <- which(upper.tri(Flooded_corr), arr.ind = TRUE)
Flooded_corrpairs <- data.frame(
  row = rownames(Flooded_corr)[indicies[, 1]],
  col = rownames(Flooded_corr)[indicies[, 2]],
  cor = Flooded_corr[indicies]
)
names(Flooded_corrpairs) <- c("Variable1", "Variable2", "Correlation")
Flooded_corrpairs <- arrange(Flooded_corrpairs, Variable1, Variable2)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(Flooded_corrpairs, file=paste0(outdir,"/Flooded_corrpairs_08-12.csv"))

# correllogram for p-val. = .05. Crosses denote not signficant (in this case not different from 0)
Floodhypothesis_Graph <- ggcorrplot(Flooded_corr,
                                    title = "Hypothesis Matrix for Flooded (p=.05)",
                                    method = "square",
                                    outline.color = "black",
                                    colors = c("red2", "white", "dodgerblue"),
                                    p.mat = FBLUPp.mat)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(Floodhypothesis_Graph, width = 6, height = 6, file=paste0(outdir,"/Floodhypothesis_Graph_08-12.png"))



#Start of graphs for unflooded data set
#
#
BLUPcorr_Unflood <- BLUPS.unflooded.df.NOchecks %>% 
  correlate() %>% 
  stretch()

graph_BLUPcorr_Unflood <- BLUPcorr_Unflood %>%
  filter(abs(r) > .3) %>%
  graph_from_data_frame(directed = FALSE)

# plots simple correlation graph - only names and connection
ggraph(graph_BLUPcorr_Unflood) +
  geom_edge_link() +
  geom_node_point() +
  geom_node_text(aes(label = name))

# graph with more detail - includes color, nodes, and edge properties
corrNodegraph_BLUPS.unflooded <- ggraph(graph_BLUPcorr_Unflood) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("red2", "white", "dodgerblue")) +
  geom_node_point(color = "gray40", size = 2) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between elements (unflooded field)")
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrNodegraph_BLUPS.unflooded, width = 10, height = 5, file=paste0(outdir,"/Nodegraph_Unflooded_08-12.png"))

# correlation calculation for ggcorrplot
BLUPS.unflood.df.NumONLY <- subset(BLUPS.unflooded.df.NOchecks, select = -c(GENOTYPE) )
Unflooded_corr <- round(cor(BLUPS.unflood.df.NumONLY, use = "pairwise.complete.obs"), 2) # Change to complete.obs if theres replacement error
UfBLUPp.mat <- cor_pmat(BLUPS.unflood.df.NumONLY, use = "pairwise.complete.obs")
# actual graph of correlation map
corrmap_Unflooded <- ggcorrplot(Unflooded_corr, method = "square",
                                title = "Correlation - BLUPS of elemental variables (U)",
                                legend.title = "Pearson Corr", lab = TRUE,
                                lab_col = "black",
                                lab_size = 1.75,
                                ggtheme = theme_minimal(),
                                colors = c("red2", "white", "dodgerblue")
)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(corrmap_Unflooded, width = 6, height = 6, file=paste0(outdir,"/CorrMap_Unflooded_08-12.png"))
corrmap_Unflooded

#Save the correlation plots as a table of pairwise correlations
indicies <- which(upper.tri(Unflooded_corr), arr.ind = TRUE)
Unflooded_corrpairs <- data.frame(
  row = rownames(Flooded_corr)[indicies[, 1]],
  col = rownames(Flooded_corr)[indicies[, 2]],
  cor = Unflooded_corr[indicies]
)
names(Unflooded_corrpairs) <- c("Variable1", "Variable2", "Correlation")
Unflooded_corrpairs <- arrange(Unflooded_corrpairs, Variable1, Variable2)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
write_csv(Unflooded_corrpairs, file=paste0(outdir,"/Unflooded_corrpairs_08-12.csv"))

# correllogram for p-val. Crosses denote insignificance at default level .05
Unfloodhypothesis_Graph <-ggcorrplot(Unflooded_corr,
                                     title = "Hypothesis Matrix for Unflooded (p=.05)",
                                     method = "square",
                                     outline.color = "black",
                                     colors = c("red2", "white", "dodgerblue"),
                                     p.mat = UfBLUPp.mat)
outdir <- "Outputs_08-12"
if (!dir.exists(outdir)){
  dir.create(outdir)
}
ggsave(Unfloodhypothesis_Graph, width = 6, height = 6, file=paste0(outdir,"/Unfloodhypothesis_Graph_08-12.png"))
#try ggcorrplot, corrr, ggraph
#https://drsimonj.svbtle.com/how-to-create-correlation-network-plots-with-corrr-and-ggraph




###############################################################################################
# Uncomment the lines below to obtain boxplots with overlayed LMNT check distributions
# Start of distribution graphs across TIL and checks
#
# library(ggplot2)
# traitlist[1] <- "DHD" # Not necessary, however, useful when looking at different imputed DHD
# 
# # Create new column for checks/TIL classification. 
# # To create the boxplots change "Flood_data" to desired trait file, e.g. "Unflood_data"
# Flood_data$CheckType <- sapply(Flood_data$GENOTYPE, function(x){
#   if (x == "LMNT") {
#     return("LMNT")
#   } else if (x == "TQNG") {
#     return("TQNG")
#   } else {
#     return("TIL")
#   }
# })
# 
# 
# # Try and place box-plots across the distribution of variables. First attempt
# library(ggstance)
# lapply(traitlist, function(trait) {
# 
# checks_data <- Flood_data[Flood_data$CheckType %in% c("LMNT", "TQNG"), ]
# 
# p <- ggplot(Flood_data, aes_string(x = trait)) +
#   geom_histogram(aes(y = ..density..), fill = 'blue', color = 'white', alpha = .5) +
#   geom_boxploth(data = checks_data, aes_string(x = trait, y = trait, fill = "CheckType"),
#                position = position_dodge(.8), width = .1) +
#   scale_fill_manual(values = c("LMNT" = "red", "TQNG" = "green"),
#                     name = "Observation Type") +
#   labs(title = paste("Distribution of ", trait),
#        x = trait,
#        y = "Density") +
#   theme_minimal()
# print(p)
# 
# # Uncomment to save the overlayed figures
# #file_name <- paste0(folder_name, "/", trait, "_overlay.png")
# #ggsave(filename = file_name, plot = p)
# })
# 
# 
# # Create both boxplot and histogram . Line up values and match scale for both
# library(gridExtra)
# library(scales)
# library(patchwork)
# 
# for (trait in traitlist) {
#   min_val <- min(Unflood_data[[trait]], na.rm = TRUE) # Find min value for a variable
#   max_val <- max(Unflood_data[[trait]], na.rm = TRUE) # Find max value for a variable
#   breaks_for_graph <- pretty(c(min_val, max_val)) # Calculate the breaks to keep consistent for traits
#   checktype_colors <- c(LMNT = "blue", TQNG = "red") 
#   
#   # Plots the histogram
#   p1 <- ggplot(Unflood_data[Unflood_data$CheckType == 'TIL',], aes_string(x = trait)) +
#     geom_histogram(fill = 'dodgerblue', color = 'white') +
#     labs(title = "Histogram for TIL") +
#     coord_cartesian(xlim = c(min_val, max_val)) +
#     scale_x_continuous(breaks = breaks_for_graph, limits = c(min_val, max_val))
#   
#   # Plots the box plot for the specified checks. Can add and remove with "c('TQNG', 'LMNT')"
#   p2 <-
#     ggplot(Unflood_data[Unflood_data$CheckType %in% c('TQNG', 'LMNT'),], aes_string(x = "CheckType", y = trait)) +
#     geom_boxplot(aes(fill = CheckType)) +
#     labs(title = "Boxplot for LMNT and TQNG") +
#     coord_cartesian(ylim = c(min_val, max_val)) +
#     coord_flip() +
#     scale_y_continuous(breaks = breaks_for_graph, limits = c(min_val, max_val)) +
#     scale_fill_manual(values = checktype_colors)
#   
#   # Could be removed? Might be unecessary if using patchwork
#   p1_dim <- get_dim(p1)
#   set_dim(p2, p1_dim)
#   aligned_plots <- align_patches(p1, p2)
#   
#   #grid.arrange(p1, p2, nrow = 2)
#   
#   #library(patchwork)
#   combined_plot <- p1 / p2
#   print(combined_plot)
#   ggsave(paste0(trait, "_plot.png"), combined_plot, width = 1280, height = 720, units = c("px"))
#   
#   
# }
# 
# plot(residuals(lm.traits.unflood$As) ~ predict(lm.traits.unflood$As))
