# conversions to add ENTREZ to the test data

setwd("~/tmp/shiny/scripts/")

library("tidyverse")
library("dplyr")
library("limma")
library("GO.db")
library("org.Hs.eg.db")

# dds_df_1$Gene <- as.character(dds_df_1$Gene)
dds_df_1 <- dds_df_1 %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname())

dds_df_2 <- dds_df_2 %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname())

dds_df_3 <- dds_df_3 %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname())

dds_df_4 <- dds_df_4 %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname())

## we need to save the results with ENTREZ IDs from now on

## get all GO IDs, and all ENTREZ IDs associated with each ID
go_list <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
                  "ENTREZID", "GO", multiVals = "list")
go_vector <- lapply(go_list, as.vector)
ezs <- sapply(go_vector, paste0, collapse = ";")
go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

### then map GO terms to genes between tables via ENTREZID
go_1 <- go_1 %>% dplyr::left_join(go_df, by ="GOID")
go_2 <- go_2 %>% dplyr::left_join(go_df, by ="GOID")
go_3 <- go_3 %>% dplyr::left_join(go_df, by ="GOID")
go_4 <- go_4 %>% dplyr::left_join(go_df, by ="GOID")

##### to do later ############################################################

# kegg_1
# kegg_2
# kegg_3
# kegg_4

saveRDS(dds_df_1, file = "../results/dds_df_1.rds")
saveRDS(dds_df_2, file = "../results/dds_df_2.rds")
saveRDS(dds_df_3, file = "../results/dds_df_3.rds")
saveRDS(dds_df_4, file = "../results/dds_df_4.rds")

saveRDS(go_1, file = "../results/go_1.rds")
saveRDS(go_2, file = "../results/go_2.rds")
saveRDS(go_3, file = "../results/go_3.rds")
saveRDS(go_4, file = "../results/go_4.rds")

# saveRDS(kegg_1, file = "./app/data/kegg_1.rds")
# saveRDS(kegg_2, file = "./app/data/kegg_2.rds")
# saveRDS(kegg_3, file = "./app/data/kegg_3.rds")
# saveRDS(kegg_4, file = "./app/data/kegg_4.rds")
