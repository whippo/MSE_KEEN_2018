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
library(stringr)

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

# summarize counts of organisms for quads by group without CALI
nocali_QUADS <- QUADS %>%
  filter(!MO_SP_CODE == "CALI")

quad_count_nocali <- nocali_QUADS %>%
  group_by(COLLECTION_TYPE, GroupingCode) %>%
  summarise(sum(COUNT))
names(quad_count_nocali)[names(quad_count_nocali)=="sum(COUNT)"] <- "TOT_COUNT"

# summarize counts of organisms for quads by group for brown algae only
brown_QUADS <- QUADS %>%
  filter(GroupingCode == "Brown algae")
brown_QUADS$MO_SP_CODE <- as.factor(brown_QUADS$MO_SP_CODE)
# rename and regroup brown algae
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="SAJU"] <- "Juv. Kelp"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="UKKJ"] <- "Juv. Kelp"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="SALANI"] <- "Saccharina spp."
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="SALA"] <- "Saccharina spp."
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="SACO"] <- "Saccharina spp."
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="COCO"] <- "C. costata"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="AGFI"] <- "A. fimbriatum"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="AGCL"] <- "A. clathratum"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="DELI"] <- "D. ligulatum"
levels(brown_QUADS$MO_SP_CODE)[levels(brown_QUADS$MO_SP_CODE)=="SAMU"] <- "S. muticum"

quad_count_browns <- brown_QUADS %>%
  group_by(COLLECTION_TYPE, MO_SP_CODE) %>%
  summarise(sum(COUNT))
names(quad_count_browns)[names(quad_count_browns)=="sum(COUNT)"] <- "TOT_COUNT"
names(quad_count_browns)[names(quad_count_browns)=="MO_SP_CODE"] <- "Species"

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



# SPECIES ACCUMULATION

# Instructors
Anim_quad_instructors <- Anim_spread_quad %>%
  filter(COLLECTION_TYPE == "INSTRUCTOR")
Anim_quad_inst_pool <- Anim_quad_instructors[,4:53]
inst_accum <- specaccum(Anim_quad_inst_pool)

# Students
Anim_quad_students <- Anim_spread_quad %>%
  filter(COLLECTION_TYPE == "STUDENT")
Anim_quad_stud_pool <- Anim_quad_students[,4:53]
stud_accum <- specaccum(Anim_quad_stud_pool)



# UNDER/OVER ESTIMATION
# Use All_quad for as starting point, use Anim_quad for animal quadrats and 
# Algae_quad for algae quadrats

All_UPC <- All_rich %>%
  filter(DATA_TYPE == "UPC")
# Extract Algae
Algae_UPC <- All_UPC %>%
  filter(GroupingCode %in% algae)

# Extract Substrate
Subs_UPC <- All_UPC %>%
  filter(GroupingCode == "Abiotic substrate")


# sum animal quads
Anim_sums <- Anim_quad %>%
  group_by(MO_SP_CODE, COLLECTION_TYPE) %>%
  summarise(sum(NEW_COUNT))
names(Anim_sums)[names(Anim_sums) == "sum(NEW_COUNT)"] <- "NEW_COUNT"

# Spread out animal quads
Anim_spread_overcount <- Anim_sums %>%
  spread(COLLECTION_TYPE, NEW_COUNT, fill = 0)
# Add column of differences
Anim_spread_overcount$diff <- Anim_spread_overcount$STUDENT - Anim_spread_overcount$INSTRUCTOR
# identify percentiles
quantile(Anim_spread_overcount$diff, prob = seq(0, 1, length = length(Anim_spread_overcount$diff)), type = 5)
# top 5% = <24, bottom 5% = >-19

Anim_quad_topbot <- Anim_spread_overcount %>%
  filter(diff > 23 | diff < -18)
Anim_quad_topbot <- left_join(Anim_quad_topbot, codes, by="MO_SP_CODE")


# sum algae quads
Algae_sums <- Algae_quad %>%
  group_by(MO_SP_CODE, COLLECTION_TYPE) %>%
  summarise(sum(TOT_COUNT))
names(Algae_sums)[names(Algae_sums) == "sum(TOT_COUNT)"] <- "NEW_COUNT"

# Spread out algae quads
Algae_spread_overcount <- Algae_sums %>%
  spread(COLLECTION_TYPE, NEW_COUNT, fill = 0)
# Add column of differences
Algae_spread_overcount$diff <- Algae_spread_overcount$STUDENT - Algae_spread_overcount$INSTRUCTOR
# identify percentiles
quantile(Algae_spread_overcount$diff, prob = seq(0, 1, length = length(Algae_spread_overcount$diff)), type = 5)
# top 12.5% = <5, bottom 12.5% = >-65

Algae_quad_topbot <- Algae_spread_overcount %>%
  filter(diff > 6 | diff < -66)
Algae_quad_topbot <- left_join(Algae_quad_topbot, codes, by="MO_SP_CODE")





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

# Figure of counts of observed groups by collection type for open quads NO CALI
COUNT_GROUPS_TYPE_NOCALI <- ggplot(quad_count_nocali, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = GroupingCode)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE) +
  labs(x="", y="Count") 
#COUNT_GROUPS_TYPE_NOCALI

# Figure of counts of observed of groups by collection type for UPC
COUNT_GROUPS_TYPE_UPC <- ggplot(UPC_count_groups, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = GroupingCode)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE, option = "plasma") +
  labs(x="", y="Count") +
  ylim(0, 1250)
#COUNT_GROUPS_TYPE_UPC

# Figure of counts of observed groups by collection type for open quads browns only
COUNT_GROUPS_TYPE_BROWNS <- ggplot(quad_count_browns, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = Species)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE, option = 3, begin = 0.2) +
  labs(x="", y="Count") 
#COUNT_GROUPS_TYPE_NOCALI


Figure1 <- ggarrange(ggarrange(COUNT_GROUPS_TYPE_QUAD, COUNT_GROUPS_TYPE_NOCALI,
                               labels = c("A", "B"),
                               ncol = 1, nrow = 2,
                               common.legend = TRUE,
                               legend = "left"), ggarrange(COUNT_GROUPS_TYPE_UPC, COUNT_GROUPS_TYPE_BROWNS,
                     labels = c("C", "D"),
                     ncol = 1, nrow = 2,
                     common.legend = FALSE, 
                     legend = "right",
                     align = "v"), 
                     ncol = 2, nrow = 1)
annotate_figure(Figure1, bottom = text_grob("Figure 1: Total contribution of taxonomic groupings to open quadrat surveys A) with Calliostoma and B) without Callistoma, in C) \n UPC quadrats, and D) of brown algal species in open quadrat surveys done by students in comparison to an instructor dataset.", size = 10))

# best size 1000x800


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

# algal quadrats
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

############### FIGURE 4
# Species accumulation curves for animals in open quads
# animal quadrats

# function to create alpha channel in plot() function
col2alpha <- function(col, alpha) {
  col_rgb <- col2rgb(col)/255
  rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
}

par(oma = c(2,2,2,2))
plot(inst_accum, ci.type="poly", col="#FDE725FF", lwd=2, ci.lty=0, ci.col=col2alpha("#2D708EFF", "1"), ylim = c(0,45), xlab = "Quadrats Sampled", ylab = "Number of Species")
plot(stud_accum, ci.type="poly", col="#FDE725FF", lwd=2, ci.lty=0, ci.col=col2alpha("#73D055FF", "0.5"), add = TRUE) +
  mtext("Figure 4: New species as a function of number of quadrats sampled shown with \n 95% confidence intervals for instructors (blue) and students (green).", side = 1, line = 5, cex = 0.8)



# best size: 550x550

############### FIGURE 5
# Over and underestimation of groups comparing student to instructor data
# animal quadrats

Anim_quad_topbot_fig <- ggplot(Anim_quad_topbot, aes(x = Genus, y = diff)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylim(-100, 350) +
  ylab("Difference in Counts") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  geom_hline(yintercept = 0, color = "#FDE725FF", size = 2)

Algae_quad_topbot_fig <- ggplot(Algae_quad_topbot, aes(x = Genus, y = diff)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylim(-120, 50) +
  ylab("Difference in Counts") +
  geom_hline(yintercept = 0, color = "#73D055FF", size = 2)

Figure5 <- ggarrange(Anim_quad_topbot_fig, Algae_quad_topbot_fig, 
                     ncol = 1, nrow = 2,
                     labels = c("A", "B"))
annotate_figure(Figure5, bottom = text_grob("Figure 5: Absolute difference in total count differences between instructor and student open quadrats for \n A) animals and B) algae. Only differences magnitude for the top ~5% for animals and ~10% for algae are shown. \n Instructor values are considered baseline (0), positive values indicate student overcounting, negative values \n indicate undercounting.", size = 10))

# best size: 660x800




#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

### SCRATCH PAD

# DRILL IN ON BROWN ALGAE

smallspp <- c("Molluscs (general)", "Urchins", "Brown algae", "Stars")

subset_stack_data <- quad_count_groups %>%
  filter(GroupingCode %in% smallspp)

# Figure of counts of observed groups by collection type for open quads
ggplot(subset_stack_data, aes(x = COLLECTION_TYPE, y = TOT_COUNT)) + geom_bar(stat = "identity", aes(fill = GroupingCode)) + 
  theme_minimal() + 
  scale_fill_viridis(discrete=TRUE) +
  labs(x="", y="Count") 
#COUNT_GROUPS_TYPE_QUAD


# NO CALLIOSTOMA

# summarize counts of organisms for quads by group
nocali_QUADS <- QUADS %>%
  filter(!MO_SP_CODE == "CALI")

quad_count_groups <- nocali_QUADS %>%
  group_by(COLLECTION_TYPE, GroupingCode) %>%
  summarise(sum(COUNT))
names(quad_count_groups)[names(quad_count_groups)=="sum(COUNT)"] <- "TOT_COUNT"
