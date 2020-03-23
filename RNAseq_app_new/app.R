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
##' @param l2FC log2 fold change threshold
##' 
##' @param thresh_padj p-adjusted value threshold
##' 
##' @param proj_name project name
##' 
##' @param thresh_ma MA threshold
##' 
##' @param thresh_vol Volcano/Gene Enrichment threshold
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

## load libraries
library("DT")
library("shiny")
library("dplyr")
library("tidyr")
library("ggplot2")
library("stringr") ###### make sure that this works on the server ###################### do we use this?? ######

####################### set current directory variables ########################
## The functions below streamline the variables that you must change between
## running the app on local, and running the app on the cluster. You must update
## the `my_directory` path on the server every time you start a new project, and
## the local path to library and results (if running on a local)
## NOTE: the server MUST USE absolute paths NEVER USE relative paths!!!

## define log2 fold change threshold
l2FC <- 1

## define p-adjusted value threshold
thresh_padj <- 0.05

## project name
proj_name <- "2020XX_XXXX_XXXX"

## MA threshold
thresh_ma <- paste0("Threshold set at p-adjusted value < ", thresh_padj)

## Volcano/Gene Enrichment threshold
thresh_vol <- paste0("Threshold set at log2 fold change > ", l2FC,
                     " & p-adjusted value < ", thresh_padj)

## path to R libaries on the server
serv_lib <- "/data/cbio/rlib/3.6"

## path to R libaries on your local
loc_lib <- "~/R/x86_64-pc-linux-gnu-library/3.6"

## this is the path to the DE results on the server
# serv_res <- paste0("/srv/shiny-private-server/cbio/", proj_name, "/app/data/")
serv_res <- "/srv/shiny-private-server/cbio/test/results/"
# serv_res <- "/srv/shiny-private-server/cbio/test/app/data/"

## this is the path to the DE results on your local
# loc_res <- paste0("~/tmp/", proj_name, "/app/data/")
loc_res <- "~/tmp/shiny/results/"  

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
        mutate(neg_log10_padj = -log10(padj))
        return(dds_df)
}

## files to be read which need to correspond to your experimental comparisons ##
## big results file (res object + normalized counts)
dds_df_1 <- read_count_matrix(paste0(myDirectory, "res1.rds"))
dds_df_2 <- read_count_matrix(paste0(myDirectory, "res2.rds"))
dds_df_3 <- read_count_matrix(paste0(myDirectory, "res3.rds"))
dds_df_4 <- read_count_matrix(paste0(myDirectory, "res4.rds"))

## saved colData
col_1 <- readRDS(paste0(myDirectory, "col_4.rds"))

## saved PCA png files
## NOTE:: this MUST be the absolute path to work on the server
pca_1 <- (paste0(myDirectory, "PCA4.png"))

############################################################################################## HERE! ############
## GO and KEGG need to be arranged by %>% arrange(P.ADJ.DE)
## need to calculate P.ADJ.DE for these files.

## NOTE: in order for the dynamic GO and KEGG tables to work, the results must
## have ENTREZ IDs as a column in each results dataframe!!!

## saved limma::goana GO results
go_1 <- readRDS(paste0(myDirectory, "go_1.rds")) # %>% arrange(P.DE)
go_2 <- readRDS(paste0(myDirectory, "go_2.rds")) # %>% arrange(P.DE)
go_3 <- readRDS(paste0(myDirectory, "go_3.rds")) # %>% arrange(P.DE)
go_4 <- readRDS(paste0(myDirectory, "go_4.rds")) # %>% arrange(P.DE)

## saved limma::kegga GO results
kegg_1 <- readRDS(paste0(myDirectory, "kegg_1.rds")) %>% arrange(P.DE)
kegg_2 <- readRDS(paste0(myDirectory, "kegg_2.rds")) %>% arrange(P.DE)
kegg_3 <- readRDS(paste0(myDirectory, "kegg_3.rds")) %>% arrange(P.DE)
kegg_4 <- readRDS(paste0(myDirectory, "kegg_4.rds")) %>% arrange(P.DE)

###################### link results files as lists #############################
## to "link" all results dataframes so that they can be selected and appear
## correctly in each dropdown menu, for each comparison, all results dataframes
## for that comparison are added as objects in a list and this list is what is
## selected by the dropdown menu and then the dataframe objects are what are
## passed to their respective plots and tables

list_1 <- list(dds_df_1, go_1, kegg_1)
list_2 <- list(dds_df_2, go_2, kegg_2)
list_3 <- list(dds_df_3, go_3, kegg_3)
list_4 <- list(dds_df_4, go_4, kegg_4)

############################ user interface ###################################
# This function defines the Shiny UI. New tabs can be added with tabPanel() and
# new rows can be added with flowRow() within the tabs and new columns can be
# added within the within rows with column(). The string that appearing in the
# drop down menu and can be changed to the name of whatever comparison is
# required, however "list_1" etc should not change

ui <- shinyUI(fluidPage(## title panel
              titlePanel(paste0("Explore RNAseq Analysis of ", proj_name, " data with Shiny")),
              mainPanel(tabsetPanel(type = "tabs",
                          
              ########## first panel ##########
              tabPanel(title = "Description",
              fluidRow(column(width = 12,
              ## first panel row 1
              fluidRow(column(width = 12,
              h4("Experimental groups of selected DE comparison"),
              DT::dataTableOutput("table_1"))
              )))),
              
              ########## second panel ##########
              tabPanel(title = "PCA",
              fluidRow(column(width = 12,
              ## second panel row 1
              fluidRow(column(width = 12,
              h5("PCA of experimental comparison"),
              plotOutput("image_1"),
              align = "left")
              )))),
              
              ########## third panel ##########              
              tabPanel("DE Analysis",
              h4(paste0("Interactive plots illustrating the results of the DESeq2 analysis of the ", proj_name, " data")),
              selectInput(inputId = "deseq",
              label = "Please choose comparison",
              choices = c("Comparison 1" = "list_1",
                          "Comparison 2" = "list_2",
                          "Comparison 3" = "list_3",
                          "Comparison 4" = "list_4")),
              ## third panel row 1
              fluidRow(
              column(width = 4,
              h4("MA Plot"),
              h5(paste0(thresh_ma)),
              plotOutput("plot_graph_1",
                         click = "plot1_click",
                         height = 350)),
              column(width = 4,
              h4("Volcano Plot"),
              h5(paste0(thresh_vol)),
              plotOutput("plot_graph_2",
                         click = "plot2_click",
                         height = 350)),
              column(width = 4,
              h4("Count Plot"),
              h5("Counts for each sample"),
              plotOutput("plot_graph_3",
              height = 350)),
              ## third panel row 2
              fluidRow(column(width = 12,
              h4("DE results table"),
              h5("This table displays results of the statistic analysis. Note: specific genes can be queried in the search bar."),
              DT::dataTableOutput("table_2"),
              align = "left")
              ))),
              
              ########## fourth panel ##########
              tabPanel(title = "GO Term Enrichment Analysis",
              fluidRow(column(width = 12,
              h4(paste0("GO Enrichment Analysis was performed using limma::goana(), which tests for over-representation of gene
                 ontology (GO) terms for specified enriched sets of genes. The Enriched Gene Set ", thresh_vol, " This enriched
                 gene set is subjected to a hypergeometric test for differential enrichment (DE) against a gene 'universe', which
                 is defined as all genes from the original count matrix that possess ENTREZ identifiers.")),
              selectInput(inputId = "go",
                          label = "Please choose comparison",
                          choices = c("Comparison 1" = "list_1",
                                      "Comparison 2" = "list_2",
                                      "Comparison 3" = "list_3",
                                      "Comparison 4" = "list_4")),
              # fourth panel row 1
              fluidRow(column(width = 12,
              h4("Volcano Plot"),
              h5(paste0(thresh_vol)),
              plotOutput("plot_graph_4", height = 350))),
              br(),
              # fourth panel row 2
              fluidRow(column(width = 12,
              h4("limma::goana results table"),
              h5("This table displaying results of the GO term enrichment analysis. Note: specific GO terms can be queried in the search bar."),
              h5("Note: Term = GO term, Ont = ontology that the GO term belongs to (e.g. 'BP', 'CC' and 'MF'), N = number of genes in the GO term,
                 DE = number of genes in the DE set, P.DE = non-adjusted p-value for over-representation of the GO term in the set."),
              # Button ########################## download button needs to be adjusted !! ############
              br(),
              br(),
              br(),
              DT::dataTableOutput("table_3")),
              br(),
              br(),
              br(),
              br(),
              fluidRow(column(width = 12,
              h4("Genes found in significant GO IDs"),
              h5("Download data of genes found in significant GO IDs"),
              # br(), br(), br(),
              downloadButton('downloadData1', 'Download data'),
              br(),
              br(),
              br(),
              DT::dataTableOutput("table_4")))
              )))),
              ########## fifth panel ##########
              tabPanel(title = "KEGG Pathway Enrichment Analysis",
              fluidRow(column(width = 12,
              h4(paste0("KEGG Pathway Enrichment Analysis was performed using limma::kegga(), which tests for over-representation of KEGG pathways in one or more sets of genes.
                 The Enriched Gene Set ", thresh_vol, " This enriched gene set is subjected to a hypergeometric test for differential enrichment (DE) against a gene 'universe',
                 which is defined as all genes from the original count matrix that possess ENTREZ identifiers.")),
              selectInput(inputId = "kegg",
                          label = "Please choose comparison",
                          choices = c("Comparison 1" = "list_1",
                                      "Comparison 2" = "list_2",
                                      "Comparison 3" = "list_3",
                                      "Comparison 4" = "list_4")),
              # fifth panel row 1
              fluidRow(column(width = 12,
              h4("Volcano Plot"),
              h5(paste0(thresh_vol)),
              plotOutput("plot_graph_5", height = 350))),
              br(),
              br(),
              br(),
              ## fifth panel row 2
              fluidRow(column(width = 12,
              h4("limma::kegga results table"),
              h5("This table displaying results of the KEGG pathway enrichment analysis. Note: specific KEGG pathways can be queried in the search bar."),
              # Button ############################### download button needs to be adjusted ! ############3
              br(),
              br(),
              br(),
              DT::dataTableOutput("table_5"),
              br(),
              br(),
              br(),
              br(),
              h4("Genes found in significant KEGG terms"),
              h5("Download data of genes found in significant GO IDs"),
              # br(), br(), br(),
              downloadButton('downloadData2', 'Download data'),
              br(),
              br(),
              br(),
              DT::dataTableOutput("table_6"))
              )))
              )))
))

################################## server ######################################
## this function defines the Shiny "server" all functions that the UI refers to
## that involve reactive variables go here
server <- function(input, output) {
    
################## selecting dataset from dropdown menu DESEQ ##################
## gets the dropdown menu input for the DE analysis tab
    
    selected_df1 <- reactive({switch(input$deseq, 
                                    "list_1" = list_1,
                                    "list_2" = list_2,
                                    "list_3" = list_3,
                                    "list_4" = list_4
    )})

################## selecting dataset from dropdown menu GO RES #################
## gets the dropdown menu input for the GO analysis tab
    
    selected_df2 <- reactive({switch(input$go,
                                    "list_1" = list_1,
                                    "list_2" = list_2,
                                    "list_3" = list_3,
                                    "list_4" = list_4
    )})

################## selecting dataset from dropdown menu KEGG RES ###############
## gets the dropdown menu input for the KEGG analysis tab
    
    selected_df3 <- reactive({switch(input$kegg,
                                     "list_1" = list_1,
                                     "list_2" = list_2,
                                     "list_3" = list_3,
                                     "list_4" = list_4
    )})

########### create reactive variable value to select clicked object ############
## This reactive variable captures the row of the gene clicked which matches
## the selected dataframe. The initial value is set to 1, (i.e. the gene with
## most sig. padj score). The reactive value will be used to plot the count plot
## and highlight the clicked gene in blue in the other 2 plots automatically.
## NOTE: nearPoints will select the nearest point to a click. threshold MUST BE
## set to 1000 so app will not crash when user clicks far from a selected point!
## if your data doesn't have "genes" for point_gene$Gene then this needs to be
## some sort of unique ID like UNIPROT or PEPTIDEID as every row must have a
## unique identifier for this reactive variable to function correctly
     
    ## reactive click variable used to    
    click_value <- reactiveVal(value = 1, label = 1)
    
    ## this reactive function records clicks on the DESeq2 tab MA plot 
    observeEvent(input$plot1_click, ignoreNULL = TRUE, {
                 point_gene <- nearPoints(selected_df1()[[1]], input$plot1_click, threshold = 1000)
                 newValue <- which(grepl(point_gene$Gene, selected_df1()[[1]]$Gene))
                 if(is.null(newValue)) {
                 click_value(1) 
                 } else {
                 click_value(newValue)     
                 }})
    
    ## this reactive function records clicks on the DESeq2 tab Volcano plot 
    observeEvent(input$plot2_click, ignoreNULL = TRUE, {
                 point_gene <- nearPoints(selected_df1()[[1]], input$plot2_click, threshold = 1000)
                 newValue <- which(grepl(point_gene$Gene, selected_df1()[[1]]$Gene))
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
    
    ## this reactive function records clicks a row of the DESeq2 table results 
    observeEvent(input$table_2_rows_selected, ignoreNULL = TRUE, {
                 newValue <-  input$table_2_rows_selected
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
    
########################## GO TERM Gene Table ##################################
## This reactive variable captures the row of the gene clicked which matches
## the selected dataframe for the GO results tab. The initial value is set to 1,
## (i.e. the GO TERM with most sig. padj score). The reactive value will be used
## to generate a table of the list of genes corresponding to the row of the
## clicked GO term, furthermore the genes corresponding to the clicked GO term
## will be highlighted in a volcano plot
    
    click_value2 <- reactiveVal(value = 1, label = 1) 
    
    ## this reactive function records clicks a row of the results in GO table
    observeEvent(input$table_3_rows_selected, ignoreNULL = TRUE, {
        newValue <-  input$table_3_rows_selected
        if(is.null(newValue)) {
            click_value2(1)
        } else {
            click_value2(newValue)
        }})

######################## KEGG Pathway Gene Table ###############################
## This reactive variable captures the row of the gene clicked which matches
## the selected dataframe for the KEGGresults tab. The initial value is set to
## 1, (i.e. the KEGG Pathway with most sig. padj score). The reactive value will
## be used to generate a table of the list of genes corresponding to the row of
## the clicked KEGG Pathway, furthermore the genes corresponding to the clicked
## KEGG Pathway will be highlighted in a volcano plot
    
    click_value3 <- reactiveVal(value = 1, label = 1) 
    
    ## this reactive function records clicks a row of the results in GO table
    observeEvent(input$table_5_rows_selected, ignoreNULL = TRUE, {
        newValue <-  input$table_5_rows_selected
        if(is.null(newValue)) {
            click_value3(1)
        } else {
            click_value3(newValue)
        }})
    
    ## if you want to know the value of what the cursor is clicking on, these
    ## reactive variables will display the output of that
    # output$click_value <- renderText({click_value()})
    # output$click_value2 <- renderText({click_value2()})

#################### First plot (MA plot) ######################################
## This MA plot is generated automatically from dataset selected from dropdown,
## with threshold at padj < 0.05. Points meeting this threshold are shown in
## red, all other points in black, and the selected gene (determined by user
## click input) is highlighted as a transparent blue circle
    
    output$plot_graph_1 <- renderPlot({
        ggplot(selected_df1()[[1]], aes(x = baseMean, y = log2FoldChange)) +
        geom_point(aes(colour = padj < thresh_padj), size = 0.5) +
        geom_point(data = selected_df1()[[1]][click_value(), ], colour = "blue", alpha = 0.4, size = 5) +
        scale_colour_manual(name = paste0('padj < ', thresh_padj), values = setNames(c('red', 'black'), c(TRUE, FALSE))) +
        scale_x_continuous(trans = "log10", limits = c(0.1, 300000)) +
        geom_smooth(colour = "red") +
        geom_abline(slope = 0,intercept = 0, colour = "blue") +
        theme(plot.title = element_text(hjust = 0.5))
    })
    
##################### Second plot (Volcano plot) ###############################
## This Volcano plot is generated automatically from dataset selected from the
## dropdown menu, with threshold at log2FoldChange > 1 and padj < 0.05. Points
## meeting this threshold are shown in red, all other points in black, and the
## selected gene (determined by user click input) is highlighted as a
## transparent blue circle
    
    output$plot_graph_2 <- renderPlot({ 
        res_df <- as.data.frame(selected_df1()[[1]], na.rm = TRUE)
        res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > l2FC & res_df$padj < thresh_padj)
        vol <- ggplot(data = res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = threshold)) +
               geom_point(size = 0.5) +
               geom_point(data = selected_df1()[[1]][click_value(), ], colour = "blue", alpha = 0.4, size = 5) +
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
        as_tibble(selected_df1()[[1]][click_value(), ]) %>% 
        dplyr::select(Gene, geneName,
                      starts_with("F65"), starts_with("H65"), starts_with("S65"),
                      starts_with("F74"), starts_with("H74"), starts_with("S74")) %>%
        head(1) %>% tidyr::gather(sample, counts, -c(Gene, geneName)) %>%
        mutate(group = substr(sample, 1, 4)) %>% 
        ggplot(aes(x = group, y = counts, color = group)) +
        geom_point(size = 3) +
        ggtitle(selected_df1()[[1]][click_value(), ]$geneName)
    })
    
########################## Plot datatable experimental coldata #################
## this table is generated from saving the experimental colData used in the dds
## object. the columns will change from experiment to experiment
    
    output$table_1 <- DT::renderDataTable(col_1)

########################## Display saved PCA.png files #########################
## this function reads file paths of pre-generated PCA ggplots from the DE
## analysis. Note, that this function reads file PATHS, not files! These file
## paths MUST BE absolute paths

    output$image_1 <- renderImage({
        filename <- file.path(pca_1)
        list(src = filename, width = 800, height = 800)
    }, deleteFile = FALSE) 
    
####################### Generate datatable DE results ##########################
## this table is generated from the dataset selected from dropdown, showing
## variables Gene, geneName, baseMean, log2FoldChange, pvalue, padj (arranged by
## highest padj), and only one row can be selected at a time
## NOTE: "selection = 'single'" means that only one row can be selected. this
## app does not currently support selecting multiple rows

    input_3 <- reactive(selected_df1()[[1]])
    output$table_2 <- DT::renderDataTable(input_3() %>%
                      dplyr::select(-neg_log10_padj) %>%
                      # mutate_at(3:5, round, 3) %>%
                      # mutate_at(8:last_col(), round, 3) %>%
                      arrange(padj), selection = 'single')

####################### Generate datatable GO results ##########################
## this table is generated from a saved dataframe displaying the limma::goana
## GO enrichment results
    
    input_4 <- reactive(selected_df2()[[2]])
    output$table_3 <- DT::renderDataTable(input_4(), selection = 'single')

######## Generate dynamic datatable of ENTREZ genes based on GO results ########
## this table is generated from a saved dataframe displaying the limma::goana
## GO enrichment results

    input_5 <- reactive(selected_df2()[[2]][click_value2(), ]$ENTREZID_in_term %>%
                        grepl(selected_df2()[[1]]$ENTREZID) %>% which())
    
    output$table_4 <- DT::renderDataTable((selected_df2()[[1]][selected_df2()[[1]]$ENTREZID %in% input_5(), ]))
    
    download_1 <- reactive((selected_df2()[[1]][selected_df2()[[1]]$ENTREZID %in% input_5(), ]))
    
############################### GO Volcano plot ################################
## This Volcano plot is generated automatically from dataset selected from the
## dropdown menu, with threshold at log2FoldChange > 1 and padj < 0.05. Points
## meeting this threshold are shown in red, all other points in black, and the
## selected genes corresponding to the GOID clicked (determined by user click
## input in the table below it) are highlighted in GREEN

    output$plot_graph_4 <- renderPlot({
        res_df <- as.data.frame(selected_df2()[[1]], na.rm = TRUE) 
        res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > l2FC & res_df$padj < thresh_padj)
        vol <- ggplot(data = res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = threshold)) +
            geom_point(size = 0.5) +
            geom_point(data = selected_df2()[[1]][selected_df2()[[1]]$ENTREZID %in% input_5(), ],
                       colour = "green", size = 0.7) +
            xlab("log2 fold change") + ylab("-log10 p-value") +
            theme(plot.title = element_text(hjust = 0.5))
        vol + scale_color_manual(values = c("#000000", "#FF0000"))
    })

######################### Generate datatable KEGG results ######################
## this table is generated from a saved dataframes displaying the limma::kegga
## KEGG enrichment results
    
    input_6 <- reactive(selected_df3()[[3]])
    output$table_5 <- DT::renderDataTable(input_6(), selection = 'single')
    
######## Generate dynamic datatable of ENTREZ genes based on KEGG results ######
## this table is generated from a saved dataframe displaying the limma::kegga
## KEGG term enrichment results
    
    input_7 <- reactive(selected_df3()[[3]][click_value3(), ]$ENTREZIDs_in_path %>%
                            grepl(selected_df3()[[1]]$ENTREZID) %>% which())
    
    output$table_6 <- DT::renderDataTable((selected_df3()[[1]][selected_df3()[[1]]$ENTREZID %in% input_7(), ]))
    
    download_2 <- reactive((selected_df3()[[1]][selected_df3()[[1]]$ENTREZID %in% input_7(), ]))
        
############################# KEGG Volcano plot ################################
## This Volcano plot is generated automatically from dataset selected from the
## dropdown menu, with threshold at log2FoldChange > 1 and padj < 0.05. Points
## meeting this threshold are shown in red, all other points in black, and the
## selected genes corresponding to the KEGG Pathway clicked (determined by user
## click input in the table below it) are highlighted in GREEN
    
    output$plot_graph_5 <- renderPlot({ 
        res_df <- as.data.frame(selected_df3()[[1]], na.rm = TRUE) 
        res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > l2FC & res_df$padj < thresh_padj)
        vol <- ggplot(data = res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = threshold)) +
            geom_point(size = 0.5) +
            geom_point(data = selected_df3()[[1]][selected_df3()[[1]]$ENTREZID %in% input_7(), ],
                       colour = "green", size = 0.7) +
            xlab("log2 fold change") + ylab("-log10 p-value") +
            theme(plot.title = element_text(hjust = 0.5))
        vol + scale_color_manual(values = c("#000000", "#FF0000"))
    })
    
######################### Download Buttons #####################################
## this reactive variable passes the selected table download to the download
## button. The downloadHandler takes a filename argument, which tells the web
## browser what filename to default to when saving.
## https://shiny.rstudio.com/reference/shiny/1.0.4/downloadButton.html
    
    ## DOWNLOAD BUTTON 1 Downloadable csv of selected dataset GO dynamic table
    output$downloadData1 <- downloadHandler(
        filename = function() {
            paste("GOID_genes", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(download_1(), file, sep = ",", row.names = FALSE)
        }
    )
    
    ## DOWNLOAD BUTTON 2 Downloadable csv of selected dataset KEGG dynamic table
    output$downloadData2 <- downloadHandler(
        filename = function() {
            paste("KEGG_path_genes", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(download_2(), file, sep = ",", row.names = FALSE)
        }
    )

} ### nothing EVER goes below this line because shiny will throw an error!!! ###
shinyApp(ui = ui, server = server)