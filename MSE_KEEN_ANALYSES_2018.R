###################################################################################
#                                                                                ##
# Title of Script                                                                ##
# Data are current as of 2018-07-10                                              ##
# Data source: Marine Subtidal Ecology Course - Kelp Ecosystem Ecology Network - ##
# Friday Harbor Labs - University of Oregon                                      ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2018-07-10                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:

# Data analyses from the 2018 Subtidal Marine Ecology Course taught at Friday 
# Harbor Labs by Aaron Galloway, Alex Lowe, Pema Kitaeff, and Mo Turner. Data
# quantifies kelp and kelp-associated communitites at the South West corner of 
# Shaw Island, WA in June 2018. Data were collected according to the Kelp Ecosystem
# Ecology Network experimental protocol, but using 3 removal plots, 3 control plots
# and 4 quadrats of each kind (UPC, Open Quadrat) per plot. Definitive data was 
# collected by Galloway and Lowe. Redundant data was collected by students for 
# comparison. 

# Required Files (check that script is loading latest version):
# QUAD-Data.csv
# UPC_SUBS_CODES.csv
# WA_SP_CODES.csv

# Associated Scripts:
# None

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# DATA MANIPULATION                                                               #
# GENERATE FIGURES                                                                #
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2018-07-09 Script created by Ross Whippo

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(vegan)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

quaddata <- read.csv("QUAD_Data.csv")
codes <- read.csv("WA_SP_CODES.csv")
subs <- read.csv("UPC_SUBS_CODES.csv")

# create column designating 'INSTRUCTOR' vs 'student' data
quaddata <- quaddata %>%
  mutate(COLLECTION_TYPE = ifelse(OBSERVER == "Lowe" | OBSERVER == "Galloway", 
                                  "INSTRUCTOR", "STUDENT"))
quaddata$COLLECTION_TYPE <- as.factor(quaddata$COLLECTION_TYPE)

# create column of full species names
names(codes)[names(codes)=="SP_CODE"] <- "MO_SP_CODE"
quaddata <- left_join(quaddata, codes, by="MO_SP_CODE")

quaddata$GroupingCode <- as.character(quaddata$GroupingCode)
quaddata$GroupingCode[quaddata$GroupingCode == ""] <- "Abiotic substrate"
quaddata$GroupingCode[quaddata$GroupingCode == "Mollusks (general)"] <- "Molluscs (general)"
quaddata$GroupingCode <- as.factor(quaddata$GroupingCode)

# open quadrat data
QUADS <- quaddata %>%
  subset(DATA_TYPE == "QUAD")

# UPC quadrat data
UPC <- quaddata %>%
  subset(DATA_TYPE == "UPC")

###################################################################################
# DATA MANIPULATION                                                               #
###################################################################################

# summarize counts of organisms for quads by group
quad_count_groups <- QUADS %>%
  group_by(COLLECTION_TYPE, GroupingCode) %>%
  summarise(sum(COUNT))
names(quad_count_groups)[names(quad_count_groups)=="sum(COUNT)"] <- "TOT_COUNT"

# summarize counts of organisms and substrate for UPC by group
UPC_count_groups <- UPC %>%
  subset(DATA_TYPE == "UPC") %>%
  group_by(COLLECTION_TYPE, GroupingCode) %>%
  summarise(sum(COUNT)) 
names(UPC_count_groups)[names(UPC_count_groups)=="sum(COUNT)"] <- "TOT_COUNT"


# Create matrix for community composition analysis
quad_matrix <- QUADS %>%
  group_by(COLLECTION_TYPE, PLOT, MO_SP_CODE) %>%
  summarise(sum(COUNT))
names(quad_matrix)[names(quad_matrix)=="sum(COUNT)"] <- "TOT_COUNT"

KEEN_mat_quad_full <- quad_matrix %>%
  spread(MO_SP_CODE, TOT_COUNT, fill = 0)
KEEN_quad_mat_run <- KEEN_mat_quad_full[,c(3:61)]

# MDS of community
KEEN_quad_mds <- metaMDS(KEEN_quad_mat_run)
KEEN_quad_mds_points <- KEEN_quad_mds$points
KEEN_quad_mds_points <- data.frame(KEEN_quad_mds_points)
plot_data_mds <- data.frame(KEEN_quad_mds_points, KEEN_mat_quad_full[,1:2])
vec_quad <- envfit(KEEN_quad_mds$points, KEEN_quad_mat_run, perm = 1000)
vec_quad_df <- as.data.frame(vec_quad$vectors$arrows*sqrt(vec_quad$vectors$r))
vec_quad_df$species <- rownames(vec_quad_df)
library(plyr)
chulls_quad <- ddply(plot_data_mds, .(COLLECTION_TYPE), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)
#plot_data_mds$PLOT <- as.factor(plot_data_mds$PLOT)

###################################################################################
# ANALYSES                                                                        #
###################################################################################

# All Richness
All_rich <- quaddata %>%
  select(PLOT, QUAD, DATA_TYPE, COUNT, COLLECTION_TYPE, MO_SP_CODE, GroupingCode) %>%
  group_by(PLOT, QUAD, DATA_TYPE, COLLECTION_TYPE, MO_SP_CODE, GroupingCode) %>%
  summarise(sum(COUNT))
names(All_rich)[names(All_rich)=="sum(COUNT)"] <- "TOT_COUNT"

# Extract Quadrats
All_quad <- All_rich %>%
  filter(DATA_TYPE == "QUAD")
# Extract Animals
algae <- c("Brown algae", "Green algae", "Red algae fleshy (fuct groups)
", "Red algae fleshy", "Red algae coraline")
Anim_quad <- All_quad %>%
  filter(!GroupingCode %in% algae)
# total redundant grouping codes
Anim_quad <- Anim_quad %>%
  ungroup() %>%
  select(-c(GroupingCode, DATA_TYPE)) %>%
  group_by(PLOT, QUAD, COLLECTION_TYPE, MO_SP_CODE, TOT_COUNT) %>%
  summarise(sum(TOT_COUNT))
names(Anim_quad)[names(Anim_quad)=="sum(TOT_COUNT)"] <- "NEW_COUNT"
Anim_quad <- Anim_quad %>%
  select(-TOT_COUNT)

# spread quad animals
Anim_spread_quad <- Anim_quad %>%
  group_by(COLLECTION_TYPE) %>%
  spread(MO_SP_CODE, NEW_COUNT, fill = 0)

Anim_pool_quad <- Anim_spread_quad[,4:53]
Anim_quad_plots <-Anim_spread_quad[,c(1:3)]
Anim_quad_plots <-Anim_quad_plots %>%
  unite_("PLOTQUAD_TYPE", c("PLOT", "QUAD", "COLLECTION_TYPE"), sep = "_", remove = FALSE)
attach(Anim_quad_plots)
Anim_quad_rich <- Anim_pool_quad %>%
  specpool(PLOTQUAD_TYPE)
Anim_quad_rich <- Anim_quad_rich %>%
  bind_cols(Anim_quad_plots)
Anim_quad_rich$PLOT <- as.factor(Anim_quad_rich$PLOT)


# extract quad algae
Algae_quad <- All_quad %>%
  filter(GroupingCode %in% algae)

# spread quad animals
Algae_spread_quad <- Algae_quad %>%
  group_by(COLLECTION_TYPE) %>%
  spread(MO_SP_CODE, TOT_COUNT, fill = 0)

Algae_pool_quad <- Algae_spread_quad[,6:14]
Algae_quad_plots <-Algae_spread_quad[,c(1:2,4)]
Algae_quad_plots <-Algae_quad_plots %>%
  unite_("PLOTQUAD_TYPE", c("PLOT", "QUAD", "COLLECTION_TYPE"), sep = "_", remove = FALSE)
attach(Algae_quad_plots)
Algae_quad_rich <- Algae_pool_quad %>%
  specpool(PLOTQUAD_TYPE)
Algae_quad_rich <- Algae_quad_rich %>%
  bind_cols(Algae_quad_plots)
Algae_quad_rich$PLOT <- as.factor(Algae_quad_rich$PLOT)

###################################################################################
# GENERATE FIGURES                                                                #
###################################################################################

############### FIGURE 1
# Total contribtion of Group Codes for Quadrats and UPC by Students and
# Instructors

# Figure of counts of observed groups by collection type for open quads
COUNT_GROUPS_TYPE_QUAD <- ggplot(quad_count_groups, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = GroupingCode)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE) +
  labs(x="", y="Count") 
#COUNT_GROUPS_TYPE_QUAD

# Figure of counts of observed of groups by collection type for UPC
COUNT_GROUPS_TYPE_UPC <- ggplot(UPC_count_groups, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = GroupingCode)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE, option = "plasma") +
  labs(x="", y="Count") +
  ylim(0, 1250)
#COUNT_GROUPS_TYPE_UPC


Figure1 <- ggarrange(COUNT_GROUPS_TYPE_QUAD, COUNT_GROUPS_TYPE_UPC, 
                     labels = c("A", "B"),
                     ncol = 1, nrow = 2,
                     common.legend = FALSE, 
                     legend = "right",
                     align = "v")
annotate_figure(Figure1, bottom = text_grob("Figure 1: Total contribution of taxonomic groupings to A) open quadrat and B) UPC \n surveys done by students in comparison to the definitive dataset.", size = 10))

# best size 600x850


############### FIGURE 2
# Visualization of community composition by Students and Instructors from open
# quadrats

# MDS
QUAD_MDS <- ggplot(plot_data_mds, aes(x=MDS1, y=MDS2, pch = COLLECTION_TYPE, color = factor(PLOT))) + 
  labs(color = "PLOT") +
  theme_minimal() +
  geom_polygon(data=chulls_quad, aes(x=MDS1, y=MDS2, group=COLLECTION_TYPE), fill=NA, color = "grey") +
  geom_point(size = 4) +
  scale_color_viridis_d()

Figure2 <- ggarrange(QUAD_MDS, 
                     ncol = 1, nrow = 1,
                     common.legend = FALSE, 
                     legend = "right")
annotate_figure(Figure2, bottom = text_grob("Figure 2: MDS of open quadrat communities surveyed by students comparted to surveys \n conducted by instructors.", size = 10))

# best size 700x500

############### FIGURE 3
# Richness of animals and algae
# animal quadrats

QUAD_RICH_ANIM <- ggplot(Anim_quad_rich, aes(x = PLOT, y = Species, fill=COLLECTION_TYPE)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  ylab("Observed Animal Species Richness") +
  xlab("Plot Number") 
  
QUAD_RICH_ALGAE <- ggplot(Algae_quad_rich, aes(x = PLOT, y = Species, fill=COLLECTION_TYPE)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  ylab("Observed Algal Species Richness") +
  xlab("Plot Number") 

Figure3 <- ggarrange(QUAD_RICH_ANIM, QUAD_RICH_ALGAE, 
                     ncol = 1, nrow = 2,
                     labels = c("A", "B"),
                     common.legend = TRUE, 
                     legend = "right")
annotate_figure(Figure3, bottom = text_grob("Figure 3: Observed measures of species richness for A) animals and B) algae \n in open quadrats measured by instructors and students", size = 10))

# best size: 500x700

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#