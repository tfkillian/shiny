##' This shiny app allowscount matrix data to be visually explored. Different
##' sample group comparisons resulting from the analysis can be selected from
##' the dropdown menu in the left side panel of the user interface. This app
##' displays a  *Volcano plot*, a data table of results from the analysis and
##' a plot of counts of a selected gene determined by user click input. The
##' plots and table shown in the user interface are updated automatically with
##' new user input.
##'
##' @title Interactive Shiny Proteomic Data Analysis Explorer
##' 
##' @author Theo Killian
##' 
##' The following variables *must be updated* when changing input files
##'
##' @param dds_df_1 a tibble of counts
##' 
##' @param dds_df_2 a tibble of counts
##' 
##' @param "dds_df_1" a string corresponding to the name of the .tsv file
##' 
##' @param "dds_df_2" a string corresponding to the name of the .tsv file
##' 
##' The following variables *must not be updated* when changing input files
##' 
##' @param click_value reactive variable that matches click user input with the
##' row of the dataframe of the data point being clicked on
##' 
##' @param newValue variable that sets that gives reactive variable value() the
##' user click input
##' 
##' @param plot2_click reactive variable accepting click input on volcano plot
##' 
##' @param plot3_click reactive variable accepting click input on data table
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

####### Load results of  analysis to be displayed in Shiny app ###########
read_count_matrix <- function(dds_df) {
    dds_df <- readRDS(dds_df) %>% as_tibble()
    return(dds_df)
}

# set current directory variable
# myDirectory <- "/srv/shiny-private-server/cbio/201910_PHOS_Shengju/"

## relative path on local machine
dds_df_1 <- read_count_matrix("/home/tfkillian/tmp/201910_PHOS_Shengju/results/res_201910_PHOS_Shengju.rds")
dds_df_1 <- dds_df_1 %>% dplyr::mutate(neg_log10_padj = -log10(padj_limma)) %>% arrange(padj_limma)

## absolute path on server
#dds_df_1 <- read_count_matrix("/srv/shiny-private-server/cbio/201910_PHOS_Shengju/results/res_201910_PHOS_Shengju.rds")
#dds_df_1 <- dds_df_1 %>% dplyr::mutate(neg_log10_padj = -log10(padj_limma)) %>% arrange(padj_limma)

############################ user interface ###################################
# NOTE: "Control vs 2hr", etc is the string that appearing in the drop down menu
# and can be changed to the name of whatever comparison is required, however
# "dds_df_1" should remain the name of whatever the .csv file is named

ui <- shinyUI(fluidPage(titlePanel("Explore 201910_PHOS_Shengju Human Proteomic Data with Shiny"),
                        h5("Statistical analysis of 3 controls vs. 3 exp. samples using limma"),
                        selectInput(inputId = "comparison",
                                    label = "Please choose comparison",
                                    choices = c(#"Raw Data" = "dds_df_1",
                                                "Quantile Normalized" = "dds_df_1")),
                        fluidRow(column(width = 4,
                                        h4("Volcano Plot"),
                                        h5("Threshold set at LFC > 1 & padj < 0.05"),
                                        plotOutput("plot_graph_2",
                                                   click = "plot2_click",
                                                   height = 350)),
                                 column(width = 4,
                                        h4("Count Plot"),
                                        h5("Quantile normalized counts for each sample"),
                                        plotOutput("plot_graph_3",
                                                   height = 350)),
                                 fluidRow(column(width = 4,
                                                 h4("Results table"),
                                                 h5("This table showing results of the statistic analysis, as well as the quantile normalized counts for each sample. Note: specific genes can be queried in the search bar."),
                                                 DT::dataTableOutput("table"),
                                                 offset = 1,
                                                 align = "left")))))

################################## server ######################################

server <- function(input, output) {
    
    ##################### selecting dataset from dropdown menu #####################
    # gets the dropdown menu input for the dataset. NOTE: "dds_df_1" = dds_df_1 both
    # should be the name of the data table to be selected
    
    selected_df <- reactive({switch(input$comparison, 
                                    "dds_df_1" = dds_df_1
    )})
    
    ########## create reactive variable value() to select clicked object ###########
    ## This reactive variable captures the row of the gene clicked which matches
    ## the selected dataframe. The initial value is set to 1, (i.e. the gene with
    ## most sig. padj score). The reactive value will be used to plot the count plot
    ## and highlight the clicked gene in blue in the other 2 plots automatically.
    ## NOTE: nearPoints will select the nearest point to a click. threshold must be
    ## set to 1000 so app will not crash when user clicks far from a selected point!
    
    click_value <- reactiveVal(value = 1, label = 1)
    
    observeEvent(input$plot2_click, ignoreNULL = TRUE, {
        point_gene <- nearPoints(selected_df(), input$plot2_click,
                                 threshold = 100000000)
        newValue <- which(grepl(point_gene$logFC, selected_df()$logFC))
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
    
    ##################### Second plot (Volcano plot) ###############################
    ## This Volcano plot is generated automatically from dataset selected from the
    ## dropdown menu, with threshold at log2FoldChange > 2 and padj < 0.05. Points
    ## meeting this threshold are shown in red, all other points in black, and the
    ## selected gene (determined by user click input) is highlighted as a
    ## transparent blue circle
    
    output$plot_graph_2 <- renderPlot({ 
        res_df <- as.data.frame(selected_df(), na.rm = TRUE)
        # res_df$threshold <- as.factor(abs(res_df$log_fold_change) >
        #                                   2 & res_df$padj < 0.05)
        res_df$threshold <- as.factor(abs(res_df$logFC) >
                                          2 & res_df$padj_limma < 0.05)
        # vol <- ggplot(data = res_df, aes(x = log_fold_change, y = neg_log10_padj,
        #                                  color = threshold)) +
        vol <- ggplot(data = res_df, aes(x = logFC, y = neg_log10_padj,
                                         color = threshold)) +
            geom_point(size = 0.5) +
            geom_point(data = selected_df()[click_value(), ], colour = "blue",
                       alpha = 0.4, size = 5) +
            xlab("log2 fold change") + ylab("-log10 padj-value") +
            theme(plot.title = element_text(hjust = 0.5))
        vol + scale_color_manual(values=c("#000000", "#FF0000"))
    })
    
    ######################### Third plot (Count plot) ##############################
    ## This plot is generated from user clicks on the other two plots or data table
    ## row, showing the normalized counts of the selected gene for that comparison.
    ## NOTE: when updating with different dataframes, "control" and "disease" should be
    ## changed to some part of character string of the header of the comparison
    ## groups that are being compared. Also, (sample, 1, 3) should be changed to
    ## position of the characters in the character string in the header of the
    ## comparison groups being compared
    
    output$plot_graph_3 <- renderPlot({ 
        as_tibble(selected_df()[click_value(), ]) %>% 
            dplyr::select(peptide_id, uniprot_id, starts_with("c"), starts_with("s"))%>%
            head(1) %>% tidyr::gather(sample, counts, -c(peptide_id, uniprot_id)) %>%
            dplyr::mutate(treatment = substr(sample, 1, 2)) %>%
            ggplot(aes(x = treatment, y = counts, color = treatment)) +
            geom_point(size = 1) +
            ggtitle(selected_df()[click_value(), ]$uniprot_id)
    })
    
    ########################## Plot datatable (DT) #################################
    ## this table is generated from the dataset selected from dropdown, showing
    ## variables Gene, geneName, baseMean, log2FoldChange, pvalue, padj (arranged by
    ## highest padj), and only one row can be selected at a time
    
    inputx <- reactive(get(input$comparison))
    output$table <- DT::renderDataTable(inputx() %>%
    dplyr::select(peptide_id, mod, psm, uniprot_id, pvalue_limma, padj_limma, logFC,
                  c1, c2, c3, s1, s2, s3) %>%
                  mutate_at(7:13, funs(round(., 3))) %>% 
                  arrange(padj_limma), selection = 'single')
}

shinyApp(ui = ui, server = server)