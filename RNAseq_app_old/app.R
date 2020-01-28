##' This shiny app allows RNAseq count matrix data analyzed by *DESeq2* to be
##' visually explored. Different sample group comparisons resulting from the
##' *DESeq2* analysis can be selected from the dropdown menu in the left side
##' panel of the user interface. This app displays a *MA plot*, *Volcano plot*,
##' a data table of results from the *DESeq2* analysis and a plot of normalized
##' counts of a selected gene determined by user click input. The plots and
##' table shown in the user interface are updated automatically with new user
##' input.
##'
##' @title Interactive Shiny RNAseq DESeq2 Analysis Explorer
##' 
##' @author Theo Killian
##' 
##' The following variables *must be updated* when changing input files
##'
##' @param dds_df_1 a tibble of an RNAseq DESeq2 analysis with normalized counts
##' 
##' @param dds_df_2 a tibble of an RNAseq DESeq2 analysis with normalized counts
##' 
##' @param dds_df_3 a tibble of an RNAseq DESeq2 analysis with normalized counts
##' 
##' @param "dds_df_1" a string corresponding to the name of the .rds file
##' 
##' @param "dds_df_2" a string corresponding to the name of the .rds file
##' 
##' @param "dds_df_3" a string corresponding to the name of the .rds file
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
##' @param inputx a reactive variable that captures the dataframe selected from
##' user input
##' 
##' @param table a reactive output that captures the value of `inputx` to be
##' shown in the user interface

## set libpath
#.libPaths("/data/cbio/rlib/3.6") # comment out this line if not running on server
.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6") # local libpath

## load libraries
library("DT")
library("shiny")
library("dplyr")
library("tidyr")
library("ggplot2")

####### Load results of DEseq2 analysis to be displayed in Shiny app ###########
## The function read_count_matrix() returns a count matrix dataframe with a new
## added 'neg_log10_padj' column. This function expects a RNAseq count matrix
## data file in .tsv format derived from a DESeq2 analysis. The input file must
## have the following columns: `Gene`, `geneName`, `baseMean`, `log2FoldChange`,
## `pvalue`, `padj`, and at least 2 columns of normalized counts. The headers of
## the columns of the normalized counts must have some similar characters
##
## Example: df <- read_count_matrix("../results/name_of_count_matrix_file.tsv")

read_count_matrix <- function(dds_df) {
        dds_df <- readRDS(dds_df) %>% as_tibble() %>% 
        mutate(neg_log10_padj = -log10(padj)) %>% 
            dplyr::select(-stat)
        return(dds_df)
}

## relative path on local machine
dds_df_1 <- read_count_matrix("~/tmp/201910_FATH_RNA1/results/res1_F_only_201910_FATH_RNA1.rds")
dds_df_2 <- read_count_matrix("~/tmp/201910_FATH_RNA1/results/res1_H_only_201910_FATH_RNA1.rds")
dds_df_3 <- read_count_matrix("~/tmp/201910_FATH_RNA1/results/res1_S_only_201910_FATH_RNA1.rds")
dds_df_4 <- read_count_matrix("~/tmp/201910_FATH_RNA1/results/res1_pH_only_201910_FATH_RNA1.rds")

## absolute path on server
# dds_df_1 <- read_count_matrix("/srv/shiny-private-server/cbio/201910_FATH_RNA1/results/res1_F_only_201910_FATH_RNA1.rds")
# dds_df_2 <- read_count_matrix("/srv/shiny-private-server/cbio/201910_FATH_RNA1/results/res1_H_only_201910_FATH_RNA1.rds")
# dds_df_3 <- read_count_matrix("/srv/shiny-private-server/cbio/201910_FATH_RNA1/results/res1_S_only_201910_FATH_RNA1.rds")
# dds_df_4 <- read_count_matrix("/srv/shiny-private-server/cbio/201910_FATH_RNA1/results/res1_pH_only_201910_FATH_RNA1.rds")

############################ user interface ###################################
# The string that appearing in the drop down menu and can be changed to the name
# of whatever comparison is required, however "dds_df_1" should remain the name
# of whatever the .tsv file is named

ui <- shinyUI(fluidPage(titlePanel("Explore 201910_FATH_RNA1 RNAseq Analysis with Shiny"),
                        h5("DESeq2 analysis of RNAseq samples"),
              selectInput(inputId = "comparison",
                          label = "Please choose comparison",
                          choices = c("pH 6.5 vs 7.4 for F cells" = "dds_df_1",
                                      "pH 6.5 vs 7.4 for H cells" = "dds_df_2",
                                      "pH 6.5 vs 7.4 for S cells" = "dds_df_3",
                                      "pH 6.5 vs 7.4 for all cells" = "dds_df_4")),
              fluidRow(column(width = 4,
                              h4("MA Plot"),
                              h5("Threshold set at padj < 0.05"),
                              plotOutput("plot_graph_1",
                                         click = "plot1_click",
                                         height = 350)),
                       column(width = 4,
                              h4("Volcano Plot"),
                              h5("Threshold set at LFC > 1 & padj < 0.05"),
                              plotOutput("plot_graph_2",
                                         click = "plot2_click",
                                         height = 350)),
                       column(width = 4,
                              h4("Count Plot"),
                              h5("Counts for each sample"),
                              plotOutput("plot_graph_3",
                                         height = 350)),
              fluidRow(column(width = 4,
                              h4("Results table"),
                              h5("This table showing results of the statistic analysis. Note: specific genes can be queried in the search bar."),
                              DT::dataTableOutput("table"),
                              offset = 1,
                              align = "left")))))

################################## server ######################################

server <- function(input, output) {
    
##################### selecting dataset from dropdown menu #####################
# gets the dropdown menu input for the dataset. NOTE: "dds_df_1" = dds_df_1 both
# should be the name of the data table to be selected
    
    selected_df <- reactive({switch(input$comparison, 
                                    "dds_df_1" = dds_df_1,
                                    "dds_df_2" = dds_df_2,
                                    "dds_df_3" = dds_df_3,
                                    "dds_df_4" = dds_df_4
    )})
        
########## create reactive variable value() to select clicked object ###########
## This reactive variable captures the row of the gene clicked which matches
## the selected dataframe. The initial value is set to 1, (i.e. the gene with
## most sig. padj score). The reactive value will be used to plot the count plot
## and highlight the clicked gene in blue in the other 2 plots automatically.
## NOTE: nearPoints will select the nearest point to a click. threshold must be
## set to 1000 so app will not crash when user clicks far from a selected point!
        
    click_value <- reactiveVal(value = 1, label = 1)
 
    observeEvent(input$plot1_click, ignoreNULL = TRUE, {
                 point_gene <- nearPoints(selected_df(), input$plot1_click,
                                          threshold = 1000)
                 newValue <- which(grepl(point_gene$Gene, selected_df()$Gene))
                 if(is.null(newValue)) {
                 click_value(1) 
                 } else {
                 click_value(newValue)     
                 }})
    
    observeEvent(input$plot2_click, ignoreNULL = TRUE, {
                 point_gene <- nearPoints(selected_df(), input$plot2_click,
                                          threshold = 1000)
                 newValue <- which(grepl(point_gene$Gene, selected_df()$Gene))
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
    
    observeEvent(input$table_rows_selected, ignoreNULL = TRUE, {
                 newValue <-  input$table_rows_selected
                 if(is.null(newValue)) {
                 click_value(1)
                 } else {
                 click_value(newValue)
                 }})
    
    # output$click_value <- renderText({click_value()})
    
#################### First plot (MA plot) ######################################
## This MA plot is generated automatically from dataset selected from dropdown,
## with threshold at padj < 0.05. Points meeting this threshold are shown in
## red, all other points in black, and the selected gene (determined by user
## click input) is highlighted as a transparent blue circle
    
    output$plot_graph_1 <- renderPlot({
        ggplot(selected_df(), aes(x = baseMean, y = log2FoldChange)) +
            geom_point(aes(colour = padj < 0.05), size = 0.5) +
            geom_point(data = selected_df()[click_value(), ], colour = "blue",
                       alpha = 0.4, size = 5) +
            scale_colour_manual(name = 'padj < 0.05',
                                values = setNames(c('red', 'black'),
                                                  c(TRUE, FALSE))) +
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
        res_df <- as.data.frame(selected_df(), na.rm = TRUE)
        res_df$threshold <- as.factor(abs(res_df$log2FoldChange) >
                                          1 & res_df$padj < 0.05)
        vol <- ggplot(data = res_df, aes(x = log2FoldChange, y = neg_log10_padj,
                                          color = threshold)) +
            geom_point(size = 0.5) +
            geom_point(data = selected_df()[click_value(), ], colour = "blue",
                       alpha = 0.4, size = 5) +
            xlab("log2 fold change") + ylab("-log10 p-value") +
            theme(plot.title = element_text(hjust = 0.5))
        vol + scale_color_manual(values=c("#000000", "#FF0000"))
    })
    
######################### Third plot (Count plot) ##############################
## This plot is generated from user clicks on the other two plots or data table
## row, showing the normalized counts of the selected gene for that comparison.
## NOTE: when updating with different dataframes, "E16" and "E18" should be
## changed to some part of character string of the header of the comparison
## groups that are being compared.
    
        output$plot_graph_3 <- renderPlot({ 
        as_tibble(selected_df()[click_value(), ]) %>% 
        dplyr::select(Gene, geneName,
                      starts_with("F651"), starts_with("H651"), starts_with("S651"),
                      starts_with("F652"), starts_with("H652"), starts_with("S652"),
                      starts_with("F743"), starts_with("H743"), starts_with("S743")) %>%
        head(1) %>% tidyr::gather(sample, counts, -c(Gene, geneName)) %>%
        mutate(group = substr(sample, 1, 4)) %>% 
        ggplot(aes(x = group, y = log1p(counts), fill = group)) +
        geom_bar(stat = "identity") +
        ggtitle(selected_df()[click_value(), ]$geneName)
    })
    
########################## Plot datatable (DT) #################################
## this table is generated from the dataset selected from dropdown, showing
## variables Gene, geneName, baseMean, log2FoldChange, pvalue, padj (arranged by
## highest padj), and only one row can be selected at a time
    
    inputx <- reactive(get(input$comparison))
    output$table <- DT::renderDataTable(inputx() %>%
    dplyr::select(-neg_log10_padj) %>%
                  arrange(padj), selection = 'single')
}

shinyApp(ui = ui, server = server)