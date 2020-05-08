##' This shiny app allows count matrix data analyzed for differential expression
##' to be visually explored. Different sample group comparisons resulting from
##' the DE analysis can be selected from the dropdown menu in each tab of the
##' user interface.
##' 
##' The first tab contains a description of the analysis and a table of
##' experimental groups used in the differential expression analysis
##' 
##' The second tab contains a PCA displaying the projection of the experimental
##' groups used in the differential expression analysis, projected onto the
##' first two principal components.
##' 
##' The third tab displays a *MA plot*, *Volcano plot*, a data table of results
##' from the DE analysis and a plot of normalized counts of a selected gene
##' determined by user click input. The plots and table shown in the user
##' interface are updated automatically with new user input.
##' 
##' The fourth tab displays results of the Gene Set Enrichment Analysis (GSEA)
##' in an html table.
##'
##' @title Interactive Shiny Omics Analysis Explorer
##' 
##' @author Theo Killian
##' 
##' The following variables *must be updated* when changing input files
##'
##' @param dds_df_1 a tibble of an RNAseq DESeq2 analysis with normalized counts
##' 
##' @param "dds_df_1" a string corresponding to the name of the .rds file above
##' 
##' @param col_1 a dataframe of the coldata an RNAseq DESeq2 analysis
##' 
##' @param "col_1" a string corresponding to the name of the .rds file above
##' 
##' @param pca_1 a .png file of the PCA taken from the  RNAseq DESeq2 analysis
##' 
##' @param "pca_1" a string corresponding to the name of the .png file above
##' 
##' @param go_1 a dataframe of the GO enrichment of an RNAseq DESeq2 analysis
##' 
##' @param "go_1" a string corresponding to the name of the .rds file above
##' 
##' @param serv_lib path to R libaries on the server
##' 
##' @param loc_lib path to R libaries on your local
##' 
##' @param serv_res this is the path to the DE results on the server
##' 
##' @param loc_res this is the path to the DE results on your local
##' 
##' The following variables *must not be updated* when changing input files
##' 
##' @param click_value reactive variable that matches click user input with the
##' row of the dataframe of the data point being clicked on
##' 
##' @param newValue variable that sets that gives reactive variable value() the
##' user click input
##' 
##' @param plot1_click reactive variable accepting click input on MA plot
##' 
##' @param plot2_click reactive variable accepting click input on volcano plot
##' 
##' @param plot3_click reactive variable accepting click input on data table
##' 
##' @param plot_graph_1 reactive output for MA plot that generates plot
##' 
##' @param plot_graph_2 reactive output for Volcano plot that generates plot
##' 
##' @param plot_graph_3 reactive output for count plot that generates plot
##' 
##' @param input_1 a reactive variable that selects the dropdown comparison user
##' input used to display the coldata in table_1 on tab 1
##' 
##' @param input_2 a reactive variable that selects the dropdown comparison user
##' input used to display the PCA on tab 2
##' 
##' @param input_3 a reactive variable that selects the dropdown comparison user
##' input used to display the DE results, Volcano Plot, MA Plot and results
##' table_2 on tab 3
##' 
##' @param input_4 a reactive variable that selects the dropdown comparison user
##' input used to display the GO enrichment in table_3 on tab 4
##' 
##' @param table_1 a reactive output that captures the value of `input_1` and
##' selects which is the coldata dataframe corresponding to the experimental
##' comparison found in the dropdown bar on tab 1
##'
##' @param table_2 a reactive output that captures the value of `input_3` and
##' selects which is the DE result dataframe corresponding to the experimental
##' comparison found in the dropdown bar on tab 3
##' 
##' @param table_3 a reactive output that captures the value of `input_4` to be
##' selects which is the GO enrichment dataframe corresponding to the
##' experimental comparison found in the dropdown bar on tab 4
##' 

## load libraries
library("DT")
library("shiny")
library("dplyr")
library("tidyr")
library("ggplot2")

####################### set current directory variables ########################
## The functions below streamline the variables that you must change between
## running the app on local, and running the app on the cluster. You must update
## the `my_directory` path on the server every time you start a new project, and
## the local path to library and results (if running on a local)
## NOTE: the server MUST USE absolute paths NEVER USE relative paths!!!


## project name
proj_name <- "202001_GECE_scRNAseq"

## define log2 fold change threshold
l2FC_thr <- 1

## define p-adjusted value threshold
padj_thr <- 0.05

## path to R libaries on the server
serv_lib <- "/data/cbio/rlib/3.6"

## path to R libaries on your local
loc_lib <- "~/R/x86_64-pc-linux-gnu-library/3.6"

## this is the path to the DE results on the server
serv_res <- "/srv/shiny-private-server/cbio/202001_GECE_scRNAseq/app/data/"

## this is the path to the DE results on your local
loc_res <- "~/tmp/202001_GECE_scRNAseq/app/data/"  

## This function switches R libraries 
setLibs <- function(x) {
    if (x == TRUE) {
    my_lib <- serv_lib
    } else {
    my_lib <- loc_lib
    }
    return(my_lib)
}

## This function switches directory paths for results
setPaths <- function(x) {
    if (x == TRUE) {
    my_directory <- serv_res
    } else {
    my_directory <- loc_res   
    }
    return(my_directory)
}

################################################################################
#### if this is running on the server, change these functions to TRUE !!!!! ####
### if this is running on a local machine, then set these functions to FALSE ###
myDirectory <- setPaths(FALSE)
.libPaths(setLibs(FALSE))
################################################################################

####### Load results of DEseq2 analysis to be displayed in Shiny app ###########
## The function read_count_matrix() returns a count matrix dataframe with a new
## added 'neg_log10_padj' column. This function expects a RNAseq count matrix
## data file in .tsv format derived from a DESeq2 analysis. The input file must
## have the following columns: `Gene`, `geneName`, `baseMean`, `log2FoldChange`,
## `pvalue`, `padj`, and at least 2 columns of normalized counts. The headers of
## the columns of the normalized counts must have some similar characters
## Example: df <- read_count_matrix(paste0(myDirectory, "results_file.rds"))

read_count_matrix <- function(dds_df) {
        dds_df <- readRDS(dds_df) %>% as_tibble() %>%
            dplyr::rename(Gene = Row.names) %>%
            dplyr::mutate(neg_log10_padj = -log10(padj)) %>% 
            dplyr::arrange(padj)
        return(dds_df)
}

## files to be read which need to correspond to your experimental comparisons ##
## big results file (res object + normalized counts)
dds_1 <- read_count_matrix(paste0(myDirectory, "dds_CD4_vs_Unknown.rds"))
dds_2 <- read_count_matrix(paste0(myDirectory, "dds_Treatment_alone.rds"))
dds_3 <- read_count_matrix(paste0(myDirectory, "dds_Treg_vs_CD4.rds"))
dds_4 <- read_count_matrix(paste0(myDirectory, "dds_Treg_vs_Unknown.rds"))
dds_5 <- read_count_matrix(paste0(myDirectory, "dds_Treatment_in_Treg.rds"))

## saved colData
# col_1 <- readRDS(paste0(myDirectory, "colData.rds"))

## saved PCA png files
## NOTE:: this MUST be the absolute path to work on the server
umap_1 <- (paste0(myDirectory, "UMAP_W_PredictedCellType2.png"))

## saved limma::goana GSEA results
# go_1 <- readRDS(paste0(myDirectory, "go_1.rds"))
# go_2 <- readRDS(paste0(myDirectory, "go_2.rds"))

############################ user interface ###################################
# this function defines the Shiny UI. new tabs can be added with tabPanel() and
# new rows can be added with flowRow() within the tabs and new columns can be
# added within the within rows with column(). For more info, see shiny tutorial
# https://shiny.rstudio.com/tutorial/
# The string that appearing in the drop down menu and can be changed to the name
# of whatever comparison is required, however "dds_df_1" should remain the name
# of whatever the .rds file is named

ui <- shinyUI(fluidPage(
              titlePanel("Explore Single-Cell RNAseq Analysis of 202001_GECE_scRNAseq
                         with Shiny"),
              mainPanel(
              tabsetPanel(type = "tabs",

              ########## second panel ##########
              tabPanel(title = "UMAP",
              fluidRow(column(width = 12,
              # h4("PCA is a method of visually identifying the similarity or difference between samples. PCA rotates the data cloud
              #    onto an orthagonal basis determined by the dimensions of maximal variance. The first two Principal Components (PCs)
              #    usually hold the majority of the variance of the data. The following plot shows the count matrix samples projected
              #    onto the two largest Principal Components (PC1 and PC2)."),
              # selectInput(inputId = "pca",
              #             label = "Please choose comparison",
              #             choices = c("pH 6.5 vs 7.4 for F cells" = "pca_1",
              #                         "pH 6.5 vs 7.4 for H cells" = "pca_2",
              #                         "pH 6.5 vs 7.4 for S cells" = "pca_3",
              #                         "pH 6.5 vs 7.4 for all cells" = "pca_4")),
              ## second panel row 1
              fluidRow(
              column(width = 12,
              h5("UMAP of experimental comparison"),
              plotOutput("image_1"),
              align = "left")
              )))),
              
              ########## third panel ##########              
              tabPanel("DE Analysis",
              h4("Interactive plots illustrating the results of the DESeq2 analysis "),
              selectInput(inputId = "comparison",
              label = "Please choose comparison",
              choices = c("CD4 vs Unknown" = "dds_1",
                          "Treatment alone" = "dds_2",
                          "Treg vs CD4" = "dds_3",
                          "Treg_vs_Unknown" = "dds_4",
                          "Treatment_in_Treg" = "dds_5")),
              ## third panel row 1
              fluidRow(
              # column(width = 4,
              # h4("MA Plot"),
              # h5("Threshold (red) set at p-adjusted value < 0.05"),
              # plotOutput("plot_graph_1",
              #            click = "plot1_click",
              #            height = 350)),
              column(width = 6,
              h4("Volcano Plot"),
              h5("Threshold (red) set at log2 fold change > 1 & p-adjusted value < 0.05"),
              plotOutput("plot_graph_2",
                         click = "plot2_click",
                         height = 350)),
              column(width = 6,
              h4("Count Plot"),
              h5("Counts for each sample"),
              plotOutput("plot_graph_3",
              height = 350)),
              ## third panel row 2
              fluidRow(
              column(width = 12,
              h4("DE results table"),
              h5("This table displays results of the statistic analysis. Note: specific
                 results can be queried in the search bar."),
              DT::dataTableOutput("table"),
              align = "left")
              )))
              # ,
              # 
              # ########## fourth panel ##########
              # tabPanel(title = "Enrichment Analysis",
              # fluidRow(
              # column(width = 12,
              # h4("Gene Set Enrichment Analysis (GSEA) was performed using limma::goana(), which tests for over-representation of gene
              # ontology (GO) terms or KEGG pathways in one or more sets of genes. The enriched gene set is here defined as the set of all
              # genes resulting from the DE analysis that possess log2 fold change > 1 & p-adjusted value < 0.05. This enriched gene set is
              # subjected to a hypergeometric test for differential enrichment (DE) against a gene 'universe', which is defined as all genes
              # from the original count matrix that possess ENTREZ identifiers. Note: the p-value for over-representation of the GO or KEGG
              # terms in the set (P.DE) is not p-adjusted."),
              # selectInput(inputId = "go",
              #             label = "Please choose comparison",
              #             choices = c("pH 6.5 vs 7.4 for F cells" = "go_1",
              #                         "pH 6.5 vs 7.4 for H cells" = "go_2",
              #                         "pH 6.5 vs 7.4 for S cells" = "go_3",
              #                         "pH 6.5 vs 7.4 for all cells" = "go_4")),
              # ## fourth panel row 1
              # fluidRow(
              # column(width = 12,
              # h4("limma::goana results table"),
              # h5("This table displaying results of the GO enrichment analysis. Note: specific GO terms can be queried in the search bar."),
              # h5("Note: Term = GO term, Ont = ontology that the GO term belongs to (e.g. 'BP', 'CC' and 'MF'), N = number of genes in the GO term,
              #    DE = number of genes in the DE set, P.DE = non-adjusted p-value for over-representation of the GO term in the set."),
              # DT::dataTableOutput("table_3")))
              #)))
              ))
))

################################## server ######################################
## this function defines the Shiny "server" all functions that the UI refers to
## that involve reactive variables go here
server <- function(input, output) {
    
##################### selecting dataset from dropdown menu #####################
## gets the dropdown menu input for the dataset plotted for the DE analysis tab
## NOTE: "dds_df_1" = dds_df_1  both should be the name of the data table to be
## selected correctly
    
    selected_df <- reactive({switch(input$comparison, 
                                    "dds_1" = dds_1,
                                    "dds_2" = dds_2,
                                    "dds_3" = dds_3,
                                    "dds_4" = dds_4,
                                    "dds_5" = dds_5
    )})

########### create reactive variable value to select clicked object ############
## This reactive variable captures the row of the gene clicked which matches
## the selected dataframe. The initial value is set to 1, (i.e. the gene with
## most sig. padj score). The reactive value will be used to plot the count plot
## and highlight the clicked gene in blue in the other 2 plots automatically.
## NOTE: nearPoints will select the nearest point to a click. threshold MUST BE
## set to 1000 so app will not crash when user clicks far from a selected point!
     
    ## reactive click variable used to    
    click_value <- reactiveVal(value = 1, label = 1)
    
    ## this reactive function records clicks on the MA plot 
    # observeEvent(input$plot1_click, ignoreNULL = TRUE, {
    #              point_gene <- nearPoints(selected_df(), input$plot1_click, threshold = 1000)
    #              ## if your data doesn't have "genes" this needs to be some sort
    #              ## of unique ID like UNIPROT or PEPTIDE ID
    #              newValue <- which(grepl(point_gene$Gene, selected_df()$Gene))
    #              if(is.null(newValue)) {
    #              click_value(1) 
    #              } else {
    #              click_value(newValue)     
    #              }})
    ## this reactive function records clicks on the Volcano plot 
    observeEvent(input$plot2_click, ignoreNULL = TRUE, {
                 point_gene <- nearPoints(selected_df(), input$plot2_click, threshold = 1000)
                 ## if your data doesn't have "genes" this needs to be some sort
                 ## of unique ID like UNIPROT or PEPTIDE ID
                 newValue <- which(grepl(point_gene$Gene, selected_df()$Gene))
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
    
    ## this reactive function records clicks a row of the results 
    observeEvent(input$table_rows_selected, ignoreNULL = TRUE, {
                 newValue <-  input$table_rows_selected
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
 
#################### First plot (MA plot) ######################################
## This MA plot is generated automatically from dataset selected from dropdown,
## with threshold at padj < 0.05. Points meeting this threshold are shown in
## red, all other points in black, and the selected gene (determined by user
## click input) is highlighted as a transparent blue circle
    
    # output$plot_graph_1 <- renderPlot({
    #     ggplot(selected_df(), aes(x = baseMean, y = log2FoldChange)) +
    #         geom_point(aes(colour = padj < padj_thr), size = 0.5) +
    #         geom_point(data = selected_df()[click_value(), ],
    #                    colour = "blue", alpha = 0.4, size = 5) +
    #         scale_colour_manual(name = 'padj < 0.05',
    #                             values = setNames(c('red', 'black'), c(TRUE, FALSE))) +
    #         scale_x_continuous(trans = "log10", limits = c(0.1, 300000)) +
    #         geom_smooth(colour = "red") +
    #         geom_abline(slope = 0,intercept = 0, colour = "blue") +
    #         theme(plot.title = element_text(hjust = 0.5))
    # })
    
##################### Second plot (Volcano plot) ###############################
## This Volcano plot is generated automatically from dataset selected from the
## dropdown menu, with threshold at log2FoldChange > 1 and padj < 0.05. Points
## meeting this threshold are shown in red, all other points in black, and the
## selected gene (determined by user click input) is highlighted as a
## transparent blue circle
    
    output$plot_graph_2 <- renderPlot({ 
        res_df <- as.data.frame(selected_df(), na.rm = TRUE)
        res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > l2FC_thr & res_df$padj < padj_thr)
        vol <- ggplot(data = res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = threshold)) +
            geom_point(size = 0.5) +
            geom_point(data = selected_df()[click_value(), ], colour = "blue", alpha = 0.4, size = 5) +
            xlab("log2 fold change") + ylab("-log10 p-value") +
            theme(plot.title = element_text(hjust = 0.5))
        vol + scale_color_manual(values = c("#000000", "#FF0000"))
    })
    
######################### Third plot (Count plot) ##############################
## This plot is generated from user clicks on the other two plots or data table
## row, showing the normalized counts of the selected gene for that comparison.
## NOTE: when updating with different dataframes, these names should be
## changed to some part of character string of the header of the comparison
## groups that are being compared. If pseudocounts are desired, change to log1p
    
        output$plot_graph_3 <- renderPlot({ 
        as_tibble(selected_df()[click_value(), ]) %>% 
        dplyr::select(Gene, # reads and groups the count columns by prefix
                      starts_with("Macr"), starts_with("Unkn"), starts_with("naiv"),
                      starts_with("Treg"), starts_with("Isot"), starts_with("58A2")) %>%
        head(1) %>% tidyr::gather(sample, counts, -c(Gene)) %>%
        mutate(group = substr(sample, 1, 4)) %>% ## group by prefix
        ggplot(aes(x = group, y = counts, colour = group)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.1, alpha = 0.3) +
        ggtitle(selected_df()[click_value(), ]$Gene) # display clicked gene as
    })

########################## Display saved PCA.png files #########################
## this function reads file paths of pre-generated PCA ggplots from the DE
## analysis. Note, that this function reads file PATHS, not files! These file
## paths MUST BE absolute paths
    
    # input_2 <- reactive(get(input$pca))
    # output$image_1 <- renderImage({
    #     filename <- file.path(input_2())
    #     list(src = filename, width = 800, height = 800)
    # }, deleteFile = FALSE)    
    
    output$image_1 <- renderImage({
        filename <- file.path(umap_1)
        list(src = filename, width = 800, height = 600
        )
    }, deleteFile = FALSE)
    
########################## Plot datatable DE results ###########################
## this table is generated from the dataset selected from dropdown, showing
## variables Gene, geneName, baseMean, log2FoldChange, pvalue, padj (arranged by
## highest padj), and only one row can be selected at a time
## NOTE: "selection = 'single'" means that only one row can be selected. this
## app does not currently support selecting multiple rows
    
    input_3 <- reactive(get(input$comparison))
    output$table <- DT::renderDataTable(input_3() %>%
                      dplyr::select(Gene:padj) %>%
                      mutate_at(2:4, round, 3),
                      selection = 'single')

########################## Plot datatable GO results ###########################
## this table is generated from a saved dataframes displaying the limma::goana
## GO enrichment results
    
    # input_4 <- reactive(get(input$go))
    # output$table_3 <- DT::renderDataTable(input_4())

} ### nothing EVER goes below this line because shiny will throw an error!!! ###
shinyApp(ui = ui, server = server)