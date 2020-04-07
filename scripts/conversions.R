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
# go_list <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
#                   "ENTREZID", "GO", multiVals = "list")
# go_vector <- lapply(go_list, as.vector)
# ezs <- sapply(go_vector, paste0, collapse = ";")
# go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

## this also misses some genes
# go_list <- as.list(org.Hs.egGO2ALLEGS)
# go_vector <- lapply(go_list, as.vector)
# ezs <- sapply(go_vector, paste0, collapse = ";")
# go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

## okay, so apparently you should use this method from goana
obj <- paste0("org.Hs.egGO2ALLEGS")
orgPkg <- paste0("org.Hs.eg.db")
egGO2ALLEGS <- getFromNamespace(obj,orgPkg)
GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)
go_set <- GeneID.PathID %>%
          dplyr::select(go_id, gene_id) %>%
          dplyr::rename(GOID = go_id)

# go_32501 <- go_set %>% dplyr::filter(GOID == "GO:0032501") 
# 
# go_ez <- go_32501 %>%
#          group_by(GOID) %>%
#          dplyr::summarise(ENTREZIDs = paste0(gene_id, collapse = ";"))
# 
# go_list <- go_set %>%
#            group_by(GOID) %>%
#            summarise(ENTREZIDs = paste(gene_id, collapse = ";"))

go_list <- go_set %>%
    group_by(GOID) %>%
    unique() %>%
    summarise(ENTREZID_in_term = paste(gene_id, collapse = ";"))
    #summarise(ENTREZIDs = paste(gene_id, collapse = ";"), n = dplyr::n())

go_df <- go_list

### then map GO terms to genes between tables via ENTREZID
go_1 <- go_1 %>% dplyr::select(-ENTREZID_in_term)
go_2 <- go_2 %>% dplyr::select(-ENTREZID_in_term)
go_3 <- go_3 %>% dplyr::select(-ENTREZID_in_term)
go_4 <- go_4 %>% dplyr::select(-ENTREZID_in_term)

go_1 <- go_1 %>% dplyr::left_join(go_df, by ="GOID")
go_2 <- go_2 %>% dplyr::left_join(go_df, by ="GOID")
go_3 <- go_3 %>% dplyr::left_join(go_df, by ="GOID")
go_4 <- go_4 %>% dplyr::left_join(go_df, by ="GOID")

# go_1$ENTREZID_in_term <- as.character(go_1$ENTREZID)
# go_2$ENTREZID_in_term <- as.character(go_2$ENTREZID)
# go_3$ENTREZID_in_term <- as.character(go_3$ENTREZID)
# go_4$ENTREZID_in_term <- as.character(go_4$ENTREZID)

# go_1 <- go_1 %>% dplyr::select(-ENTREZID)
# go_2 <- go_2 %>% dplyr::select(-ENTREZID)
# go_3 <- go_3 %>% dplyr::select(-ENTREZID)
# go_4 <- go_4 %>% dplyr::select(-ENTREZID)

###################### KEGG conversion ###################################
# https://www.researchgate.net/post/How_i_can_get_a_list_of_KEGG_pathways_and_its_list_of_genes
library("org.Hs.eg.db")
mapped <- mappedkeys(org.Hs.egPATH2EG)
L <- as.list(org.Hs.egPATH2EG[mapped])
Kegg_ID <- names(L)
Gene_IDs <- sapply(L, paste, collapse = ";")
kegg_genes <- cbind(Kegg_ID, Gene_IDs)
kegg_genes <- as.data.frame(kegg_genes)
names(kegg_genes) <- c("PathwayID", "ENTREZID_in_path")
rownames(kegg_genes) <- NULL
path_id_kegg <- kegg_genes
# write.table(cbind(Kegg_ID, Gene_IDs), file="KEGG to Genes.txt", sep="\t", row.names=FALSE, col.names=FALSE)

### this method seems to provide erroneous ENRTEZIDs 
# kegg_organism <- "mmu"
# GK <- getGeneKEGGLinks(species.KEGG = kegg_organism)
# GK_list <- as.list(GK)
# ezs <- sapply(GK_list, paste0, collapse = ";")

## get all KEGG PathwayIDs, and all ENTREZ IDs associated with each PathwayID
# GK_desc <- getKEGGPathwayNames(species.KEGG = kegg_organism, remove = TRUE)
# GK_col <- GK %>% group_by(PathwayID) %>% summarise(GeneIDs = paste(GeneID, collapse = ";"))
# path_id_kegg <- dplyr::full_join(GK_col, GK_desc, by = "PathwayID")
# names(path_id_kegg) <- c("PathwayID", "ENTREZ_GeneIDs", "Pathway")
# saveRDS(as.data.frame(path_id_kegg), file = "./app/data/kegg_path_2_id_mouse.rds")

# kegg_1$PathwayID <- gsub("path:mmu", "", kegg_1$PathwayID)
# kegg_2$PathwayID <- gsub("path:mmu", "", kegg_2$PathwayID)
# kegg_3$PathwayID <- gsub("path:mmu", "", kegg_3$PathwayID)
# kegg_4$PathwayID <- gsub("path:mmu", "", kegg_4$PathwayID)

kegg_1 <- kegg_1 %>% dplyr::select(-ENTREZID_in_path) # %>% dplyr::select(-PathwayID)
kegg_2 <- kegg_2 %>% dplyr::select(-ENTREZID_in_path) # %>% dplyr::select(-PathwayID)
kegg_3 <- kegg_3 %>% dplyr::select(-ENTREZID_in_path) # %>% dplyr::select(-PathwayID)
kegg_4 <- kegg_4 %>% dplyr::select(-ENTREZID_in_path) # %>% dplyr::select(-PathwayID)

kegg_1 <- kegg_1 %>% left_join(path_id_kegg, by = "PathwayID")
kegg_2 <- kegg_2 %>% left_join(path_id_kegg, by = "PathwayID")
kegg_3 <- kegg_3 %>% left_join(path_id_kegg, by = "PathwayID")
kegg_4 <- kegg_4 %>% left_join(path_id_kegg, by = "PathwayID")

kegg_1$ENTREZID_in_path <- as.character(kegg_1$ENTREZID_in_path)
kegg_2$ENTREZID_in_path <- as.character(kegg_2$ENTREZID_in_path)
kegg_3$ENTREZID_in_path <- as.character(kegg_3$ENTREZID_in_path)
kegg_4$ENTREZID_in_path <- as.character(kegg_4$ENTREZID_in_path)

# names(kegg_1) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_2) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_3) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_4) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")

############################# save ######################################

saveRDS(dds_df_1, file = "../results/res1.rds")
saveRDS(dds_df_2, file = "../results/res2.rds")
saveRDS(dds_df_3, file = "../results/res3.rds")
saveRDS(dds_df_4, file = "../results/res4.rds")

saveRDS(go_1, file = "../results/go_1.rds")
saveRDS(go_2, file = "../results/go_2.rds")
saveRDS(go_3, file = "../results/go_3.rds")
saveRDS(go_4, file = "../results/go_4.rds")

saveRDS(kegg_1, file = "../results/kegg_1.rds")
saveRDS(kegg_2, file = "../results/kegg_2.rds")
saveRDS(kegg_3, file = "../results/kegg_3.rds")
saveRDS(kegg_4, file = "../results/kegg_4.rds")
