### UMAP Shiny Visualization of 202001_GECE_scRNAseq data

library("shiny")
library("SingleCellExperiment")
library("dplyr")
library("tidyr")
library("ggplot2")

## this is the path to the DE results on the server
serv_res <- "/srv/shiny-private-server/cbio/202001_GECE_scRNAseq/umap_vis_app/data/"

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
myDirectory <- setPaths(TRUE)
.libPaths(setLibs(TRUE))
################################################################################

load(paste0(myDirectory, "sceZinbCellAnnot.RData"))

# load("~/tmp/202001_GECE_scRNAseq/Analysis_Chris/out/sceZinbCellAnnot.RData")

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "feat",
                        label = "Select a feature",
                        choices = rownames(sceZinb)),
            selectInput(inputId = "factor",
                        label = "Factor",
                        selected = "Treatment",
                        choices = colnames(colData(sceZinb))),
            selectInput(inputId = "subset",
                        label = "Subset variable",
                        choices = colnames(colData(sceZinb))),
            uiOutput("subsetValue"),
            selectInput(inputId = "assay",
                        label = "Assay type ",
                        selected = "logcounts",
                        choices = assayNames(sceZinb)),
            selectInput(inputId = "dimred",
                        label = "Dimension reduction",
                        selected = "UMAP_W",
                        choices = reducedDimNames(sceZinb))
        ),
        mainPanel(
            h3("Interactive Dimensionality Reduction Visualization of 202001_GECE_scRNAseq data"),
            h4("Dimension reduction plot"),
            plotOutput("reducedDimPlot",
                       height = 350),
            h4("Count plot"),
            plotOutput("countPlot",
                       height = 350)
            
        )
    )
)

server <- function(input, output) {
    
    getSubset <- function() {
        if ("all" %in% input$subsetValue) {
            return(1:ncol(sceZinb))
        } else {
            return(which(colData(sceZinb)[, input$subset] %in% input$subsetValue))
        }
    }
    
    output$subsetValue <- renderUI({
        selectInput(inputId = "subsetValue",
                    label = "Value",
                    selected = "all",
                    choices = c("all",
                                unique(colData(sceZinb)[, input$subset])) %>%
                        setNames(nm = .))
    })
    
    output$reducedDimPlot <- renderPlot({
        feat <- input$feat
        dimred <- input$dimred
        col <- input$factor
        sel <- getSubset()
        df <- data.frame(DR1 = reducedDim(sceZinb, dimred)[sel, 1],
                         DR2 = reducedDim(sceZinb, dimred)[sel, 2], 
                         exprs = assay(sceZinb, input$assay)[feat, sel],
                         fact = colData(sceZinb)[sel, col])
        df %>%
            group_by(fact) %>%
            slice(chull(DR1, DR2)) -> hull
        ggplot(df, aes(x = DR1, y = DR2)) +
            geom_polygon(data = hull, aes(fill = fact), alpha = 0.25) +
            geom_point(aes(x = DR1, y = DR2, col = exprs))
    })
    
    output$countPlot <- renderPlot({
        sel <- getSubset()
        feat <- input$feat
        data.frame(y = assay(sceZinb, input$assay)[feat, sel],
                   weights = assay(sceZinb, "weights")[feat, sel],
                   x = colData(sceZinb)[sel, input$factor]) %>%
            ggplot(aes(x = x, y = y, fill = x, 
                       alpha = weights, weights = weights)) +
            geom_violin() +
            geom_jitter(width = 0.25) +
            theme(legend.position = "none") +
            ggtitle(feat)
    })
}
    
shinyApp(ui, server)