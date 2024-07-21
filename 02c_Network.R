### SLEEP_PTSD_ED ###

# load required packages
if (!require("bootnet")) {install.packages("bootnet"); require("bootnet")}
if (!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if (!require("forcats")) {install.packages("forcats"); require("forcats")}
if (!require("magrittr")) {install.packages("magrittr"); require("magrittr")}
if (!require("NetworkComparisonTest")) {install.packages("NetworkComparisonTest"); require("NetworkComparisonTest")}
if (!require("networktools")) {install.packages("networktools"); require("networktools")}
if (!require("qgraph")) {install.packages("qgraph"); require("qgraph")}

# Centrality function 
# https://reisrgabriel.com/blog/2021-10-11-compare-centrality/
compareCentrality <- function(net1, net2,
                              include = c("Strength",
                                          "Closeness",
                                          "Betweenness",
                                          "ExpectedInfluence",
                                          "all",
                                          "All"),
                              orderBy = c("Strength",
                                          "Closeness",
                                          "Betweenness",
                                          "ExpectedInfluence"),
                              decreasing = T,
                              legendName = '',
                              net1Name = 'Network 1',
                              net2Name = 'Network 2'){
  
  if(include == "All" | include == "all"){
    include = c("Strength",
                "Closeness",
                "Betweenness",
                "ExpectedInfluence")
  }
  
  df <- centralityTable(net1, net2) %>% filter(measure %in% include)
  
  df %>% 
    mutate(graph = case_when(graph == 'graph 1' ~ net1Name,
                             graph == 'graph 2' ~ net2Name),
           graph = as.factor(graph),
           node = as.factor(node)) %>% 
    
    mutate(node = fct_reorder(node, value)) %>% 
    
    ggplot(aes(x = node, y = value, group = graph)) +
    
    geom_line(aes(linetype = graph), size = 1) +
    
    labs(x = '', y = '') +
    
    scale_linetype_discrete(name = legendName) +
    
    coord_flip() +
    
    facet_grid(~measure) +
    
    theme_bw()
  
}

difference_value <- function(NCT, alpha = 0.05){
  
  diff_edges <- NCT$einv.pvals %>% dplyr::filter(`p-value` <= alpha)
  
  for (i in 1:nrow(diff_edges)) {
    var_1 <- as.character(diff_edges[i, 1])
    var_2 <- as.character(diff_edges[i, 2])
    
    value_net_1 <- NCT$nw1[var_1, var_2]
    value_net_2 <- NCT$nw2[var_1, var_2]
    
    abs_difference <- abs(value_net_1 - value_net_2)
    p_value <- diff_edges$`p-value`[i]
    
    cat("Test Edge", i, "\n----\n")
    cat(var_1, "and", var_2)
    cat("\nNetwork 1:", value_net_1,
        "\nNetwork 2:", value_net_2)
    cat("\nAbsolute difference:", abs_difference,
        "with p-value =", p_value, "\n----\n")
  }
}

# Import data

data <- read.csv("Created_Data/ED_ptsd/SLEEPPTSDEDS_imputed_11-21-23.csv")
data <- data %>% dplyr::select(-id)
data <- na.omit(data)

# Check data 
summary(data)
glimpse(data)
colnames(data)
colnames(data)[13:14] = c("bdi_16", "bdi_20")

# Set seed 
set.seed(1210)


#### Network with combined sample ###
network_all <- estimateNetwork(data %>% dplyr::select(-PTSDyn), 
                               default = "EBICglasso", 
                               corMethod = "cor", 
                               corArgs = list(method = "spearman", use = "pairwise.complete.obs")) # Spearman correlations

summary(network_all)
network_all$graph
plot(network_all)

# Plot graph 
# Number of nodes
n = 16 

# Symptom clusters
groupsint = list("Eating Disorder" = c(1:6),
               "Sleep Disturbances" = c(13:16),
               "PTSD" = c(7:12))

# Define column names and item descriptions
names2 <- c("Binge", "Purge", "Fast", "Restrict", "Fear Weight", "Feel Fat", 
            "Intrusive", "Avoid Thinking", "Avoid Activities", "Unable Pos", 
            "Alert", "Startle", "Sleep Pattern", "Fatigue", "Nightmare", 
            "Sleep Trouble")

items <- c(
  "How many times did you have a sense of having lost control over your eating?", 
  "How many times have you made yourself sick (vomit) as a means of controlling your shape/weight?",
  "Have you gone for long periods of time (8 waking hours or more) without eating anything at all in order to influence your shape or weight?", 
  "Have you tried to exclude from your diet any foods that you like in order to influence your shape or weight?", 
  "Have you had a definite fear that you might gain weight?",
  "Have you felt fat?", 
  "Repeated, disturbing memories, thoughts, or images of a stressful experience from the past?",
  "Avoid thinking about or talking about a stressful experience from the past or avoid having feelings related to it", 
  "Avoid activities or situations because they remind you of a stressful experience from the past?", 
  "Feeling emotionally numb or unable to have loving feelings for those close to you?", 
  "Being super alert or watchful on guard?", 
  "Feeling jumpy or easily startled.", 
  "Changes in sleep pattern.", 
  "Tiredness or fatigue.", 
  "Repeated, disturbing dreams of a stressful experience from the past?", 
  "Trouble falling or staying asleep?"
)

# Plot 1
plot_all <- plot(network_all, 
                 labels = names, 
                 layout = "spring", 
                 vsize = 8, 
                 cut = 0, 
                 border.width = .7, 
                 border.color = "black", 
                 groups = groupsint, 
                 color = c("lightcoral", "slategray2", "darkseagreen"),
                 label.color = "black", 
                 legend.cex = .3,
                 negDashed = TRUE)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# Plot 2
plot_allv2 <- plot(network_all, 
                   groups = groupsint, 
                   labels = names,
                   layout = "spring", 
                   theme = "colorblind", 
                   palette = "colorblind", 
                   
                   # Output Arguments
                   width = 7 * 1.4, 
                   height = 7,
                   
                   # Graphical Arguments 
                   vsize = 7, # size of nodes 
                   border.width = 0.1, # controls width of border
                   
                   # Node Labels 
                   label.cex = 0.9, 
                   label.color = "black",
                   label.prop = 1, # proportion of the width of the node that the label scales
                   negDashed = TRUE, 
                   
                   # Legend 
                   legend = TRUE,
                   legend.cex = 0.5 # scalar of the legend 
)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# Plot 3
plot_allv3 <- plot(network_all, 
                   groups = groupsint, 
                   labels = names2,
                   layout = "spring", 
                   theme = "colorblind", 
                   palette = "colorblind", 
                   
                   # Output Arguments
                   width = 7 * 1.4, 
                   height = 7,
                   
                   # Graphical Arguments 
                   vsize = 7, # size of nodes 
                   border.width = 0.1, # controls width of border
                   
                   # Node Labels 
                   label.cex = 0.9, 
                   label.color = "black",
                   label.prop = 1, # proportion of the width of the node that the label scales
                   negDashed = TRUE, 
                   
                   # Legend 
                   legend = TRUE,
                   legend.cex = 0.5 # scalar of the legend 
)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# Plot 4
plot_allv4 <- plot(network_all, 
                   groups = groupsint, 
                   labels = names2,
                   layout = "spring", 
                   theme = "colorblind", 
                   palette = "colorblind", 
                   
                   # Output Arguments
                   width = 7 * 1.4, 
                   height = 7 * 2.4,
                   
                   # Graphical Arguments 
                   vsize = 7, # size of nodes 
                   border.width = 0.1, # controls width of border
                   
                   # Node Labels 
                   label.cex = 0.9, 
                   label.color = "black",
                   label.prop = 1, # proportion of the width of the node that the label scales
                   negDashed = TRUE, 
                   
                   # Legend 
                   legend = TRUE,
                   legend.cex = 0.4, # scalar of the legend 
                   legend.mode = "style1", 
                   nodeNames = items, 
                   GLratio = 2.5, # relative size of graph compared to the layout
                   layoutScale = c(0.7, 0.9), 
                   layoutOffset = c(-0.4, 0) 
)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# Centrality graph
centralityPlot(network_all, scale="raw")

centralityPlot(network_all, scale="raw", labels=names2)
centralityPlot(network_all, scale="raw", labels=names2, orderBy="Strength")

# Check that you get the same thing as above 
centralityPlot(plot_all, scale="raw", orderBy="Strength", labels=names2)
centralityPlot(plot_all, scale="z-scores", orderBy="Strength", labels=names2)

# Centrality Table
cen_tab = centralityTable(network_all, standardized=FALSE) # raw values
cen_tab_strength = subset(cen_tab, measure == "Strength")
cen_tab_strength

# Constructing a partial correlation matrix
All_pcormat <-getWmat(network_all)
write.csv(All_pcormat, "Output/All_pcormat.csv")

bridge(network_all$graph, 
       communities = c('1','1','1','1','1','1',
                       '2','2','2','2','2','2',
                       '3','3','3','3'), 
       useCommunities = "all", 
       directed = NULL, 
       nodes = NULL)

# Double check
bridge(plot_allv3, communities = c('1','1','1','1','1','1',
                                   '2','2','2', '2','2','2',
                                   '3','3','3','3'), 
       useCommunities = "all", 
       directed = NULL, 
       nodes = NULL)

# Create bridge graph 
intbridge <- bridge(network_all$graph, 
                    communities=c('1','1','1','1','1','1',
                                  '2','2','2','2','2','2',
                                  '3','3','3','3'), 
                    useCommunities = "all",
                    directed = NULL,
                    nodes = NULL)

# Double check
intbridge2 <- bridge(plot_allv3, 
                     communities=c('1','1','1','1','1','1',
                                   '2','2','2', '2','2','2',
                                   '3','3','3','3'), 
                     useCommunities = "all", 
                     directed = NULL, 
                     nodes = NULL)

plot(intbridge2, 
     order = "value", 
     standardized = TRUE, 
     include = "Bridge Expected Influence (1-step)")

# Save bridge graph as pdf
pdf("Output/All_bridge.pdf", width = 15)
plot(intbridge)
plot(intbridge2)
dev.off()

# Check network stability
All_boot <- bootnet(network_all,
                    boots = 1000, 
                    nCores = 4)

All_boot_case <- bootnet(network_all,
                         boots = 1000, 
                         nCores = 4, 
                         type = "case")

save(All_boot, file = "Output/All_boot.Rdata")
save(All_boot_case, file = "Output/All_boot_case.Rdata")

# Plot edge weight CI
pdf("Output/Stability_all.pdf")
plot(All_boot, labels = FALSE, order = "sample") 
plot(All_boot, labels = TRUE, order = "sample") 
dev.off()

# Edge weights diff test
pdf("Output/Edge_difftest_all.pdf")
plot(All_boot, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

# Plot centrality stability
pdf("Output/Stability_all.pdf") 
plot(All_boot_case)
dev.off()

# Centrality stability coefficient
corStability(All_boot_case)

# Centrality diff test
pdf("Output/Centrality_difftest_all.pdf")
plot(All_boot, "strength", order="sample", labels=TRUE, names=names) 
dev.off()

caseDroppingBoot <- bootnet(network_all, 
                            boots = 1000, 
                            type = "case", 
                            statistics = c("bridgeStrength", "bridgeCloseness", "bridgeBetweenness", "bridgeExpectedInfluence"), 
                            communities = groupsint)

corStability(caseDroppingBoot)

plot(caseDroppingBoot, statistics="bridgeExpectedInfluence")

plot(caseDroppingBoot, statistics="bridgeStrength")

(nonParametricBoot <- bootnet(network_all, 
                             boots = 1000, 
                             type = "nonparametric", 
                             statistics = c("bridgeStrength", "bridgeExpectedInfluence"), 
                             communities = groupsint))
table(nonParametricBoot$bootTable$type)

pdf("Output/BEIdifftest_all.pdf")
plot(nonParametricBoot, statistics = "bridgeExpectedInfluence", plot = "difference")
dev.off()

pdf("Output/BSdifftest_all.pdf")
plot(nonParametricBoot, statistics = "bridgeStrength", plot = "difference")
dev.off()


### PTSD AND ED NETWORK ###

# Select participants with a 1 for PTSDyn variable
network_ptsd <- estimateNetwork(data %>% filter(PTSDyn == 1) %>%dplyr::select(-PTSDyn), 
                                default = "EBICglasso", 
                                corMethod = "cor", 
                                corArgs = list(method = "spearman", use = "pairwise.complete.obs")) # Spearman correlation
n = 16 #number of nodes
colnames(data)

# Symptom clusters
groupsint = list("Eating Disorder" = c(1:6),
                 "Sleep Disturbances" = c(13:16),
                 "PTSD" = c(7:12))

# Define column names and item descriptions
names2 <- c("Binge", "Purge", "Fast", "Restrict", "Fear Weight", "Feel Fat", 
            "Intrusive", "Avoid Thinking", "Avoid Activities", "Unable Pos", 
            "Alert", "Startle", "Sleep Pattern", "Fatigue", "Nightmare", 
            "Sleep Trouble")

items <- c(
  "How many times did you have a sense of having lost control over your eating?", 
  "How many times have you made yourself sick (vomit) as a means of controlling your shape/weight?",
  "Have you gone for long periods of time (8 waking hours or more) without eating anything at all in order to influence your shape or weight?", 
  "Have you tried to exclude from your diet any foods that you like in order to influence your shape or weight?", 
  "Have you had a definite fear that you might gain weight?",
  "Have you felt fat?", 
  "Repeated, disturbing memories, thoughts, or images of a stressful experience from the past?",
  "Avoid thinking about or talking about a stressful experience from the past or avoid having feelings related to it", 
  "Avoid activities or situations because they remind you of a stressful experience from the past?", 
  "Feeling emotionally numb or unable to have loving feelings for those close to you?", 
  "Being super alert or watchful on guard?", 
  "Feeling jumpy or easily startled.", 
  "Changes in sleep pattern.", 
  "Tiredness or fatigue.", 
  "Repeated, disturbing dreams of a stressful experience from the past?", 
  "Trouble falling or staying asleep?"
)

summary(network_ptsd)
network_ptsd$graph
plot(network_ptsd)


# Plot 1 
plot_ptsd <- plot(network_ptsd, 
                  labels = names, 
                  layout = "spring", 
                  vsize = 8, 
                  cut = 0, 
                  border.width = .7, 
                  border.color = "black", 
                  groups = groupsint, 
                  color = c("lightcoral", "slategray2", "darkseagreen"),
                  label.color = "black", 
                  legend.cex = .3,
                  negDashed = TRUE)
title(paste0("ED and PTSD (n = ", nrow(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn)), ")"))

# Plot 2
plot_ptsdv2 <- plot(network_ptsd, 
                    groups = groupsint, 
                    labels = names2,
                    layout = "spring", 
                    theme = "colorblind", 
                    palette = "colorblind", 
                    
                    # Output Arguments
                    width = 7 * 1.4, 
                    height = 7,
                    
                    # Graphical Arguments 
                    vsize = 7, # size of nodes 
                    border.width = 0.1, # controls width of border
                    
                    # Node Labels 
                    label.cex = 0.9, 
                    label.color = "black",
                    label.prop = 1, # proportion of the width of the node that the label scales
                    negDashed = TRUE, 
                    
                    # Legend 
                    legend = TRUE,
                    legend.cex = 0.5 # scalar of the legend 
)
title(paste0("ED and PTSD (n = ", nrow(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn)), ")"))

# Plot 3
plot_ptsdv3 <- plot(network_ptsd, 
                    groups = groupsint, 
                    labels = names2,
                    layout = "spring", 
                    theme = "colorblind", 
                    palette = "colorblind", 
                    
                    # Output Arguments
                    width = 7 * 1.4, 
                    height = 7 * 2.4,
                    
                    # Graphical Arguments 
                    vsize = 7, # size of nodes 
                    border.width = 0.1, # controls width of border
                    
                    # Node Labels 
                    label.cex = 0.9, 
                    label.color = "black",
                    label.prop = 1, # proportion of the width of the node that the label scales
                    negDashed = TRUE, 
                    
                    # Legend 
                    legend = TRUE,
                    legend.cex = 0.35, # scalar of the legend 
                    legend.mode = "style1", 
                    nodeNames = items, 
                    GLratio = 2.5, # relative size of graph compared to the layout
                    layoutScale = c(0.8, 0.9), 
                    layoutOffset = c(-0.35, 0) 
)
title(paste0("ED and PTSD (n = ", nrow(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn)), ")"))


# Centrality graph
centralityPlot(network_ptsd, scale="raw")

centralityPlot(network_ptsd, scale="raw", labels=names2)
centralityPlot(network_ptsd, scale="raw", labels=names2, orderBy="Strength")

# Check that you get the same thing as above 
centralityPlot(plot_ptsd, scale="raw", orderBy="Strength", labels=names2)
centralityPlot(plot_ptsd, scale="z-scores", orderBy="Strength", labels=names2)


# Centrality Table
cen_tab_ptsd = centralityTable(network_ptsd, standardized=FALSE) # raw values
cen_tab_ptsd_strength = subset(cen_tab_ptsd, measure == "Strength")
cen_tab_ptsd_strength

# Constructing a partial correlation matrix
PTSD_pcormat <-getWmat(network_ptsd)
write.csv(PTSD_pcormat, "PTSD_pcormat.csv")

colnames(data)
bridge(network_ptsd$graph, 
       communities = c('1','1','1','1','1','1',
                       '2','2','2', '2','2','2',
                       '3','3','3','3'),
       useCommunities = "all", 
       directed = NULL, 
       nodes = NULL)

bridge(plot_ptsdv2, 
       communities = c('1','1','1','1','1','1',
                       '2','2','2', '2','2','2',
                       '3','3','3','3'),
       useCommunities = "all", 
       directed = NULL, 
       nodes = NULL)

# Create bridge graph
intbridge_ptsd <- bridge(network_ptsd$graph, 
                         communities = c('1','1','1','1','1','1',
                                         '2','2','2','2','2','2',
                                         '3','3','3','3'),
                         useCommunities = "all", 
                         directed = NULL, 
                         nodes = NULL)

(intbridge_ptsd2 <- bridge(plot_ptsdv2, 
                          communities = c('1','1','1','1','1','1',
                                          '2','2','2','2','2','2',
                                       '3','3','3','3'),
                          useCommunities = "all", 
                          directed = NULL,  
                          nodes = NULL))

plot(intbridge_ptsd2, 
     order = "value", 
     standardized = FALSE, 
     include = "Bridge Expected Influence (1-step)")

# Save bridge graph as pdf
pdf("Output/PTSD_bridge.pdf", width = 15)
plot(intbridge_ptsd)
plot(intbridge_ptsd2)
dev.off()

# Network stability
PTSD_boot <- bootnet(network_ptsd, 
                     boots = 1000,
                     nCores = 4)

PTSD_boot_case <- bootnet(network_ptsd, 
                          boots = 1000,
                          nCores = 4, 
                          type = "case")

save(PTSD_boot, file = "Output/PTSD_boot.Rdata")
save(PTSD_boot_case, file = "Output/PTSD_boot_case.Rdata")

# Plot edge weight CI
pdf("Output/Stability_ptsd.pdf")
plot(PTSD_boot, labels = FALSE, order = "sample") 
plot(PTSD_boot, labels = TRUE, order = "sample") 
dev.off()

# Edge weights diff test
pdf("Output/Edge_difftest_ptsd.pdf")
plot(PTSD_boot, 
     "edge", 
     plot = "difference", 
     onlyNonZero = TRUE, 
     order = "sample")
dev.off()

# Plot centrality stability
pdf("Output/Stability_ptsd.pdf") 
plot(PTSD_boot_case)
dev.off()

# Centrality stability coefficient
corStability(PTSD_boot_case)

# Centrality diff test
pdf("Output/Centrality_difftest_ptsd.pdf")
plot(PTSD_boot, 
     "strength",
     order = "sample", 
     labels = TRUE,
     names = names2) 
dev.off()

# Case dropping bootstrap for PTSD network
caseDroppingBoot_ptsd <- bootnet(network_ptsd, 
                                 boots = 1000, 
                                 type = "case", 
                                 statistics = c("bridgeStrength", "bridgeCloseness", "bridgeBetweenness", "bridgeExpectedInfluence"), 
                                 communities = groupsint)

# Correlation stability for case dropping bootstrap
corStability(caseDroppingBoot_ptsd)

# Plot bridgeExpectedInfluence for case dropping bootstrap
plot(caseDroppingBoot_ptsd, statistics = "bridgeExpectedInfluence")

# Plot bridgeStrength for case dropping bootstrap
plot(caseDroppingBoot_ptsd, statistics = "bridgeStrength")

# Non-parametric bootstrap for PTSD network
nonParametricBoot_ptsd <- bootnet(network_ptsd, 
                                  boots = 1000, 
                                  type = "nonparametric", 
                                  statistics = c("bridgeStrength", "bridgeExpectedInfluence"), 
                                  communities = groupsint)


# Save bridgeExpectedInfluence difference plot
pdf("Output/BEIdifftest_ptsd.pdf")
plot(nonParametricBoot_ptsd, statistics = "bridgeExpectedInfluence", plot = "difference")
dev.off()

# Save bridgeStrength difference plot
pdf("Output/BSdifftest_ptsd.pdf")
plot(nonParametricBoot_ptsd, statistics = "bridgeStrength", plot = "difference")
dev.off()

### ED Only network ###

# Estimate Network for ED group without PTSD
network_ed <- estimateNetwork(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn), 
                              default = "EBICglasso", 
                              corMethod = "cor", 
                              corArgs = list(method = "spearman", use = "pairwise.complete.obs")) # use Spearman correlation

# Summary 
summary(network_ed)
network_ed$graph

# Plot
plot(network_ed)

# Plot 1
plot_ed <- plot(network_ed, 
                labels = names, 
                layout = "spring", 
                vsize = 8, 
                cut = 0, 
                border.width = .7, 
                border.color = "black", 
                groups = groupsint, 
                color = c("lightcoral", "slategray2", "darkseagreen"),
                label.color = "black", 
                legend.cex = .5,
                negDashed = TRUE)
title(paste0("ED_Only (n = ", nrow(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn)), ")"))

# Plot 2
plot_edv2 <- plot(network_ed, 
                  groups = groupsint, 
                  labels = names2,
                  layout = "spring", 
                  theme = "colorblind", 
                  palette = "colorblind", 
                  
                  # Output Arguments
                  width = 7 * 1.4, 
                  height = 7,
                  
                  # Graphical Arguments 
                  vsize = 7, # size of nodes 
                  border.width = 0.1, # controls width of border
                  
                  # Node Labels 
                  label.cex = 0.9, 
                  label.color = "black",
                  label.prop = 1, # proportion of the width of the node that the label scales
                  negDashed = TRUE, 
                  
                  # Legend 
                  legend = TRUE,
                  legend.cex = 0.5 # scalar of the legend 
)
title(paste0("ED_Only (n = ", nrow(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn)), ")"))

# Plot 3
plot_edv3 <- plot(network_ed, 
                  groups = groupsint, 
                  labels = names2,
                  layout = "spring", 
                  theme = "colorblind", 
                  palette = "colorblind", 
                  
                  # Output Arguments
                  width = 7 * 1.4, 
                  height = 7 * 2.4,
                  
                  # Graphical Arguments 
                  vsize = 7, # size of nodes 
                  border.width = 0.1, # controls width of border
                  
                  # Node Labels 
                  label.cex = 0.9, 
                  label.color = "black",
                  label.prop = 1, # proportion of the width of the node that the label scales
                  negDashed = TRUE, 
                  
                  # Legend 
                  legend = TRUE,
                  legend.cex = 0.38, # scalar of the legend 
                  legend.mode = "style1", 
                  nodeNames = items, 
                  GLratio = 2.5, # relative size of graph compared to the layout
                  layoutScale = c(0.75, 0.9), 
                  layoutOffset = c(-0.4, 0) 
)
title(paste0("ED_Only (n = ", nrow(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn)), ")"))

# Centrality plots
centralityPlot(network_ed, scale = "raw")
centralityPlot(network_ed, scale = "raw", labels = names2)
centralityPlot(network_ed, scale = "raw", labels = names2, orderBy = "Strength")

# Check centrality plots with plot_ed
centralityPlot(plot_ed, scale = "raw", orderBy = "Strength", labels = names2)
centralityPlot(plot_ed, scale = "z-scores", orderBy = "Strength", labels = names2)

# Centrality table
cen_tab_ed <- centralityTable(network_ed, standardized = FALSE) # raw values
cen_tab_ed_strength <- subset(cen_tab_ed, measure == "Strength")
print(cen_tab_ed_strength)

# Constructing a partial correlation matrix
ED_pcormat <- getWmat(network_ed)
write.csv(ED_pcormat, "ED_pcormat.csv")

# Calculate bridge centrality measures
bridge(network_ed$graph,
       communities = c('1','1','1','1','1','1',
                       '2','2','2','2','2','2',
                       '3','3','3','3'), 
       useCommunities = "all", 
       directed = NULL, 
       nodes = NULL)

bridge(plot_edv2, 
       communities = c('1','1','1','1','1','1',
                       '2','2','2','2','2','2',
                       '3','3','3','3'), 
       useCommunities = "all",
       directed = NULL,
       nodes = NULL)

# Create bridge graph
intbridge_ed <- bridge(network_ed$graph, communities = c('1','1','1','1','1','1',
                                                         '2','2','2', '2','2','2',
                                                         '3','3','3','3'), 
                       useCommunities = "all",
                       directed = NULL, 
                       nodes = NULL)

intbridge_ed2 <- bridge(plot_edv2, communities = c('1','1','1','1','1','1',
                                                   '2','2','2','2','2','2',
                                                   '3','3','3','3'), 
                        useCommunities = "all", 
                        directed = NULL, 
                        nodes = NULL)

# Compare bridge centrality results
print(intbridge_ed)
print(intbridge_ed2)

# Save bridge graph 
pdf("Ouput/ED_bridge.pdf", width = 15)
plot(intbridge_ed)
plot(intbridge_ed2)
dev.off()

# Plot bridge expected influence for intbridge_ed2
plot(intbridge_ed2, 
     order = "value", 
     standardized = FALSE, 
     include = "Bridge Expected Influence (1-step)")

# Stability analysis
ED_boot <- bootnet(network_ed, boots=1000,nCores=4)
ED_boot_case <- bootnet(network_ed, boots=1000,nCores=4, type="case")

save(ED_boot, file = "Output/ED_boot.Rdata")
save(ED_boot_case, file = "Output/ED_boot_case.Rdata")

# Plot edge weight confidence intervals
pdf("Stability_ed.pdf")
plot(ED_boot, labels = FALSE, order = "sample") 
plot(ED_boot, labels = TRUE, order = "sample") 
dev.off()

# Edge weights difference test
pdf("Edge_difftest_ed.pdf")
plot(ED_boot, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

# Plot centrality stability
pdf("Centrality_stability_ed.pdf") 
plot(ED_boot_case)
dev.off()

# Centrality stability coefficient
corStability(ED_boot_case)

## Centrality difference test
pdf("Output/Centrality_difftest_ed.pdf")
plot(ED_boot, "strength", order = "sample", labels = TRUE, names = names2) 
dev.off()

# Case dropping bootstrap for ED network
caseDroppingBoot_ed <- bootnet(network_ed, 
                               boots = 1000, 
                               type = "case", 
                               statistics = c("bridgeStrength", "bridgeCloseness", "bridgeBetweenness", "bridgeExpectedInfluence"), 
                               communities = groupsint)
corStability(caseDroppingBoot_ed)

# Plot case dropping bootstrap results
pdf("Output/Case_Dropping_Bootstrap_ED.pdf")
plot(caseDroppingBoot_ed, statistics = "bridgeExpectedInfluence")
plot(caseDroppingBoot_ed, statistics = "bridgeStrength")
dev.off()

# Non-parametric bootstrap for ED network
nonParametricBoot_ed <- bootnet(network_ed, 
                                boots = 1000,
                                type = "nonparametric", 
                                statistics = c("bridgeStrength", "bridgeExpectedInfluence"), 
                                communities = groupsint)

# Plot non-parametric bootstrap difference tests
pdf("Output/BEIdifftest_ED.pdf")
plot(nonParametricBoot_ed, statistics = "bridgeExpectedInfluence", plot = "difference")
dev.off()

pdf("Output/BSdifftest_ED.pdf")
plot(nonParametricBoot_ed, statistics = "bridgeStrength", plot = "difference")
dev.off()

### NETWORK COMPARISON TEST ###

## Compare Centrality
# From the function - https://reisrgabriel.com/blog/2021-10-11-compare-centrality/
#   
# Compare centrality measures

compareCentrality(network_ed, network_ptsd,
                  include = "Strength",
                  legendName = "ED Only vs PTSD and ED",
                  net1Name = "ED Only",
                  net2Name = "PTSD and ED")

##Checking invariant network structure
M <-
  max(
    abs(c(network_ed$graph) - c(network_ptsd$graph))
  )

cat("The biggest edge difference is:", M)

# Checking invariant global strength
S <-
  abs(
    sum(
      abs(c(network_ed$graph)) -
        abs(c(network_ptsd$graph))
    )
  )/2

cat("Strength difference between the two networks is:", S)

# Running and reporting NetworkComparisonTest
set.seed(123) # random permutation seed
nct_results <- NCT(network0, network1,
                   it = 1000,
                   progressbar = T)

nct_results

# Checking invariant edge strength
(nct_test_edges <- NCT(network0, network1, 
                      it = 1000, test.edges = T,
                      p.adjust.methods = "BH",
                      progressbar = F))

# Checking which edges are different 
difference_value(nct_test_edges)

# Image of all three graphs with nodes in same locations
layout1 = averageLayout(network_all, network_ptsd, network_ed)

plot(network_all, layout=layout1, cut = 0, labels = names2) 
plot(network_ptsd, layout=layout1, cut = 0, labels = names2)
plot(network_ed, layout=layout1, cut = 0, labels = names2)

### Plot Network Graphs

# Full Sample
plot(network_all, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.5, # scalar of the legend 
)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# PTSD and ED
plot(network_ptsd, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.5, # scalar of the legend 
)
title(paste0("ED and PTSD (n = ", nrow(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn)), ")"))

# ED Only
plot(network_ed, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.5, # scalar of the legend 
)
title(paste0("ED Only (n = ", nrow(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn)), ")"))

### Large Network Graphs

# Full Sample
plot(network_all, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7 * 2.4,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.35, # scalar of the legend 
     legend.mode = "style1", 
     nodeNames = items, 
     GLratio = 2.5, # relative size of graph compared to the layout
     layoutScale = c(0.8, 0.9), 
     layoutOffset = c(-0.3, 0) 
)
title(paste0("Full Sample (n = ", nrow(data %>% dplyr::select(-PTSDyn)), ")"))

# PTSD and ED
plot(network_ptsd, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7 * 2.4,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.35, # scalar of the legend 
     legend.mode = "style1", 
     nodeNames = items, 
     GLratio = 2.5, # relative size of graph compared to the layout
     layoutScale = c(0.8, 0.9), 
     layoutOffset = c(-0.3, 0) 
)
title(paste0("ED and PTSD (n = ", nrow(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn)), ")"))

# ED Only
plot(network_ed, 
     groups = groupsint, 
     labels = names2,
     layout = layout1, 
     theme = "colorblind", 
     palette = "colorblind", 
     
     # Output Arguments
     width = 7 * 1.4, 
     height = 7 * 2.4,
     
     # Graphical Arguments 
     vsize = 7, # size of nodes 
     border.width = 0.1, # controls width of border
     
     # Node Labels 
     label.cex = 0.9, 
     label.color = "black",
     label.prop = 1, # proportion of the width of the node that the label scales
     negDashed = TRUE, 
     
     # Legend 
     legend = TRUE,
     legend.cex = 0.35, # scalar of the legend 
     legend.mode = "style1", 
     nodeNames = items, 
     GLratio = 2.5, # relative size of graph compared to the layout
     layoutScale = c(0.8, 0.9), 
     layoutOffset = c(-0.3, 0) 
)
title(paste0("ED Only (n = ", nrow(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn)), ")"))

### Centrality Plots

centralityPlot(list("PTSD and ED" = plot_ptsdv2, "ED" = plot_edv2), 
               labels = names2, 
               scale = "z-scores")

centralityPlot(list("Total Sample" = plot_allv3, "PTSD and ED" = plot_ptsdv2, "ED" = plot_edv2), 
               labels = names2, 
               scale = "z-scores")

centralityPlot(list("Total Sample" = plot_allv3, "PTSD and ED" = plot_ptsdv2, "ED" = plot_edv2), 
               labels = names2, 
               scale = "z-scores", 
               include = "All")

### Network Simulation

network_all2 <- estimateNetwork(data %>% dplyr::select(-PTSDyn), 
                                default = "EBICglasso", 
                                corMethod = "cor", 
                                corArgs = list(method = "spearman", use = "pairwise.complete.obs"), 
                                refit = TRUE) # Use Spearman correlation

simRes_all = netSimulator(network_all2$graph, 
                          dataGenerator = ggmGenerator(ordinal = FALSE), 
                          default = "EBICglasso", 
                          nCases = c(100, 250, 500, 1000, 2500), 
                          tuning = 0.5, 
                          nReps = 100, 
                          nCores = 4)
simRes_all
plot(simRes_all) 
# N=100 achieves a correlation between the "true" and estimated networks above 
plot(simRes_all, yvar = c("strength", "closeness", "betweenness"))

### Network Comparison Tests

network_ptsd <- estimateNetwork(data %>% filter(PTSDyn == 1) %>% dplyr::select(-PTSDyn), 
                                default = "EBICglasso", 
                                corMethod = "cor", 
                                corArgs = list(method = "spearman", use = "pairwise.complete.obs")) # Use Spearman correlation

network_ed <- estimateNetwork(data %>% filter(PTSDyn == 0) %>% dplyr::select(-PTSDyn), 
                              default = "EBICglasso", 
                              corMethod = "cor", 
                              corArgs = list(method = "spearman", use = "pairwise.complete.obs")) # Use Spearman correlation

res_nct = NCT(network_ptsd, 
              network_ed, 
              it = 2000, 
              test.edges = TRUE, 
              edges = "all", 
              p.adjust.methods = "BH")

# Global strength results - differences in the level of connectivity 
res_nct$glstrinv.pval
res_nct$glstrinv.real 

# Omnibus test results 
res_nct$nwinv.pval
res_nct$nwinv.real 

# Edge results 
res_nct$einv.pvals
res_nct$einv.pvals[which(res_nct$einv.pvals[, 3] <= 0.05),]

# Plot Network Graphs with Maximum Edge Value

maxi = max(network_ptsd$graph, network_ed$graph)

qgraph(network_ptsd$graph, maximum = maxi, theme = "colorblind", labels = names2)
title("PTSD and ED", line = -2)

qgraph(network_ed$graph, maximum = maxi, theme = "colorblind", labels = names2)
title("ED", line = -2)

### bridge - PTSD: 
intbridge_ptsd2
intbridge_ptsd

### bridge - ED: 
intbridge_ed2
intbridge_ed

### bridge: all 
intbridge
intbridge2