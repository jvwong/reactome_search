library(dplyr)
rm(list=ls(all=TRUE)) 

### ============ Declare directories and files =========
base_dir <- "~/Projects/pc/reactome_search"
data_dir <- file.path(base_dir, "data")
reactome_query_path <- file.path(data_dir, "reactome_empty_queries.csv")
pcapps_query_path <- file.path(data_dir, "pcapps_queries_with_reactome.csv")

### ============ Load files =========

reactome_query <- read.csv(reactome_query_path, na.strings = c("", "NA")) 
pcapps_query <- read.csv(pcapps_query_path, na.strings = c("", "NA"))

### ============ Filter rows and columns =========

reactome_query <- select(reactome_query, Term, Hits) %>% na.omit()
pcapps_query <- pcapps_query %>% na.omit()

### ============ Merge data =========
reactome.w.pcapps.queries <- merge(reactome_query, pcapps_query, by.x = c("Term"), by.y = c("TERM") )

