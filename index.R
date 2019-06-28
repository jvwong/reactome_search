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

########################################
### ============ Analysis ==============
########################################
library(ggplot2)

### ============ Overall hits for pathways/interactions? ==============
num_with_pathways <- sum(reactome.w.pcapps.queries$PATHWAYS > 0)
num_with_interactions <- sum(reactome.w.pcapps.queries$INTERACTIONS > 0)

reactome.w.pcapps.queries.both <- reactome.w.pcapps.queries %>% filter(PATHWAYS > 0, INTERACTIONS > 0)
num_with_both <- dim(reactome.w.pcapps.queries.both)[1]
reactome.w.pcapps.queries.either <- reactome.w.pcapps.queries %>% filter(PATHWAYS > 0 | INTERACTIONS > 0)
num_with_either <- dim(reactome.w.pcapps.queries.either)[1]

fraction_of_queries <- c(num_with_pathways, num_with_interactions, num_with_both)/num_reactome_queries
fractionages_df <- data.frame( fraction = fraction_of_queries, type = c("Pathway","Interaction", "Pathway AND Interaction"))

fraction_df_plot <- ggplot(fractionages_df, aes(x=factor(type), y=fraction)) + geom_bar(stat="identity", width = 0.7)
fraction_df_plot + labs(x = "PC Result Type", y = "Fraction w/ PC Result") + theme(text = element_text(size=16),
                         axis.text.x = element_text(face = "bold", size = 12),
                         axis.text.y = element_text(face = "bold", size = 12))


### ============ Popular queries ==============

# Distribution of query 'Hits'
hit_distribution_df <- ggplot(reactome.w.pcapps.queries, aes(x=Hits)) + geom_bar()
hit_distribution_df + labs(x = "X", y = "Hits") 



