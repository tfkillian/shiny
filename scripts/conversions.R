# conversions to add ENTREZ to the test data

## things to do for tomorrow

## ENTREZID not unique

## in the table with the GO results

## use str_split



# 0) link list of results files for app template
# 1) each row selected in GO table makes a dynamic table of all ENTREZ results
# 2) download button for dynamic table?
# 3) volcano plot with genes

## convert kegg tables next

## Dev questions?
## do you want me to be able to select as many rows or just one
## color of highlighted genes = green

## is there a particular format that of .csv to be saved as?

## how are the lists linked? 

setwd("~/tmp/shiny/scripts/")

library("tidyverse")
library("dplyr")
library("limma")
library("GO.db")
library("org.Hs.eg.db")

# dds_df_1$Gene <- as.character(dds_df_1$Gene)
dds_df_1 <- dds_df_1 %>% dplyr::select(-neg_log10_padj) %>% 
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>%
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  F651, F652, F653, F741, F742, F743)

dds_df_2 <- dds_df_2 %>% dplyr::select(-neg_log10_padj) %>% 
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>%
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  H651, H652, H653, H741, H742, H743)

dds_df_3 <- dds_df_3 %>% dplyr::select(-neg_log10_padj) %>% 
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>% 
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  S651, S652, S653, S741, S742, S743)

dds_df_4 <- dds_df_4 %>% dplyr::select(-neg_log10_padj) %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>% 
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  F651, F652, F653, F741, F742, F743,
                  H651, H652, H653, H741, H742, H743,
                  S651, S652, S653, S741, S742, S743)

## we need to save the results with ENTREZ IDs from now on

## get all GO IDs, and all ENTREZ IDs associated with each ID
go_list <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
                  "ENTREZID", "GO", multiVals = "list")
go_vector <- lapply(go_list, as.vector)
ezs <- sapply(go_vector, paste0, collapse = ";")
go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

### then map GO terms to genes between tables via ENTREZID
go_1 <- go_1 %>% dplyr::select(-ENTREZID)
go_2 <- go_2 %>% dplyr::select(-ENTREZID)
go_3 <- go_3 %>% dplyr::select(-ENTREZID)
go_4 <- go_4 %>% dplyr::select(-ENTREZID)

go_1 <- go_1 %>% dplyr::left_join(go_df, by ="GOID")
go_2 <- go_2 %>% dplyr::left_join(go_df, by ="GOID")
go_3 <- go_3 %>% dplyr::left_join(go_df, by ="GOID")
go_4 <- go_4 %>% dplyr::left_join(go_df, by ="GOID")

go_1$ENTREZID_in_term <- as.character(go_1$ENTREZID)
go_2$ENTREZID_in_term <- as.character(go_2$ENTREZID)
go_3$ENTREZID_in_term <- as.character(go_3$ENTREZID)
go_4$ENTREZID_in_term <- as.character(go_4$ENTREZID)

go_1 <- go_1 %>% dplyr::select(-ENTREZID)
go_2 <- go_2 %>% dplyr::select(-ENTREZID)
go_3 <- go_3 %>% dplyr::select(-ENTREZID)
go_4 <- go_4 %>% dplyr::select(-ENTREZID)

##### to do later ############################################################

# kegg_1
# kegg_2
# kegg_3
# kegg_4

saveRDS(dds_df_1, file = "../results/res1.rds")
saveRDS(dds_df_2, file = "../results/res2.rds")
saveRDS(dds_df_3, file = "../results/res3.rds")
saveRDS(dds_df_4, file = "../results/res4.rds")

saveRDS(go_1, file = "../results/go_1.rds")
saveRDS(go_2, file = "../results/go_2.rds")
saveRDS(go_3, file = "../results/go_3.rds")
saveRDS(go_4, file = "../results/go_4.rds")

# saveRDS(kegg_1, file = "./app/data/kegg_1.rds")
# saveRDS(kegg_2, file = "./app/data/kegg_2.rds")
# saveRDS(kegg_3, file = "./app/data/kegg_3.rds")
# saveRDS(kegg_4, file = "./app/data/kegg_4.rds")
