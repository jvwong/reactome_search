library(dplyr)
library(stringr)
rm(list=ls(all=TRUE)) 

#############################################
### ============ Data Processing ============
#############################################

### ============ Declare directories and files =========
base_dir <- "~/Projects/pc/reactome_search"
data_dir <- file.path(base_dir, "data")
reactome_query_path <- file.path(data_dir, "reactome_empty_queries.csv")
pcapps_query_path <- file.path(data_dir, "pcapps_queries_with_reactome.csv")

### ============ Load files =========

reactome_query <- read.csv(reactome_query_path, 
                           header = TRUE, 
                           check.names = FALSE,
                           na.strings = c("", "NA")) 
pcapps_query <- read.csv(pcapps_query_path,
                         header = TRUE, 
                         check.names = FALSE,
                         colClasses=c("TERM"="character", 
                                      "PATHWAYS"="integer",
                                      "INTERACTIONS"="integer"),
                         na.strings = c("", "NA"))

### ============ Filter rows and columns =========

reactome_query <- select(reactome_query, Term, Hits) %>% na.omit()
reactome_query <- reactome_query[ with( reactome_query,  grepl("^[0-9]+$", Hits) ), ]
reactome_query <- droplevels.data.frame(reactome_query)
pcapps_query <- pcapps_query %>% na.omit()

### ============ Merge data =========
reactome.w.pcapps.queries <- merge(reactome_query, pcapps_query, by.x = c("Term"), by.y = c("TERM") )
numeric_reactome_query_hits <- as.numeric(levels(reactome.w.pcapps.queries$Hits))[reactome.w.pcapps.queries$Hits]
reactome.w.pcapps.queries$Hits <- factor(reactome.w.pcapps.queries$Hits, levels = as.character( sort( unique( numeric_reactome_query_hits ) ) ) ) 
num_reactome_queries <- dim(reactome.w.pcapps.queries)[1]
levels_Hits <- levels(reactome.w.pcapps.queries$Hits)

########################################
### ============ Analysis ==============
########################################
library(ggplot2)
library(gridExtra)
library(grid)
### ============ Overall hits for pathways/interactions? ==============
reactome.w.pcapps.queries.pathways <- reactome.w.pcapps.queries %>% filter(PATHWAYS > 0)
num_with_pathways <- dim(reactome.w.pcapps.queries.pathways)[1]
reactome.w.pcapps.queries.interactions <- reactome.w.pcapps.queries %>% filter(INTERACTIONS > 0)
num_with_interactions <- dim(reactome.w.pcapps.queries.interactions)[1]

reactome.w.pcapps.queries.both <- reactome.w.pcapps.queries %>% filter(PATHWAYS > 0, INTERACTIONS > 0)
num_with_both <- dim(reactome.w.pcapps.queries.both)[1]
reactome.w.pcapps.queries.either <- reactome.w.pcapps.queries %>% filter(PATHWAYS > 0 | INTERACTIONS > 0)
num_with_either <- dim(reactome.w.pcapps.queries.either)[1]

fraction_of_queries <- c(num_with_pathways, num_with_interactions, num_with_both)/num_reactome_queries
fractionages_df <- data.frame( fraction = fraction_of_queries, type = c("Pathway","Interaction", "Pathway AND Interaction"))

fraction_df_plot <- ggplot(fractionages_df, aes(x=factor(type), y=fraction)) + 
  geom_bar(stat="identity", width = 0.7) + 
  labs(x = "PC Result Type", y = "Fraction w/ PC Result") + 
  theme(text = element_text(size=16),
       axis.text.x = element_text(face = "bold", size = 14),
       axis.text.y = element_text(face = "bold", size = 14))

grid.arrange(fraction_df_plot, nrow=1)
### ============ Highest occurring queries  ==============


### Distribution of query 'Hits'
hit_distribution_df_xticks <- c(1, 10, 20, 30, 40, 50, 75, 101, 253)
hit_distribution_df <- ggplot(reactome.w.pcapps.queries, aes(x=Hits)) + 
  geom_bar() +
  labs(x = "Query occurence", y = "Counts") + 
  theme(text = element_text(size=16), 
        axis.text.x = element_text(face = "bold", size = 14), 
        axis.text.y = element_text(face = "bold", size = 14)) +
  scale_x_discrete(breaks=hit_distribution_df_xticks)

grid.arrange(hit_distribution_df, nrow=1)

### Most frequent hits
n_top_levels <- 10
top_Hits_df <- reactome.w.pcapps.queries[ which(reactome.w.pcapps.queries$Hits %in% tail(levels_Hits, n=n_top_levels ) ), ]
grid.table(top_Hits_df, rows = NULL)


### ============ Drill-down into PC search results  ==============
pathways_df_plot <- ggplot(reactome.w.pcapps.queries.pathways, aes(x=PATHWAYS)) + 
  geom_bar() + 
  labs(x = "No. pathways in search result", y = "Counts") + 
  theme(text = element_text(size=16),
      axis.text.x = element_text(face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 14))

interactions_df_plot <- ggplot(reactome.w.pcapps.queries.interactions, aes(x=INTERACTIONS)) + 
  geom_bar() + 
  labs(x = "No. interactions in search result", y = "Counts") + 
  theme(text = element_text(size=16),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14))

grid.arrange(pathways_df_plot, interactions_df_plot, nrow=2)

### Most frequent hits
n_top_interactions <- 100
top_INTERACTIONS_df <- reactome.w.pcapps.queries.interactions[ which(reactome.w.pcapps.queries.interactions$INTERACTIONS > n_top_interactions ), ]
grid.table(top_INTERACTIONS_df, rows = NULL)



