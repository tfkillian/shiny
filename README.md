# shiny

This repository contains Shiny web app templates for visualizaling bioinformatic
data analysis results. So far, there are two types of Shiny apps, which
are described below:

-   RNAseq viewer
    This type of app displays 3 interactive plots from an RNAseq analysis. These
    plots are, MA plot, Volcano plot and a count plot, which is dynamically
    plotted by user click input. There is also an html table displaying the
    results of the analysis which can be queried, and interacts with the display
    plots above. An advanced version of this app allows gene ontology terms to be
    selected, which genes relating to that GO terms and KEGG pathways appears in corresponding
    volcano plots.

-   Proteomic viewer
    This type of app displays 2 interactive plots from a proteomic mass spectrometry analysis. These
    plots are, Volcano plot and an expression plot, which is dynamically
    plotted by user click input. There is also an html table displaying the
    results of the analysis which can be queried, and interacts with the display
    plots above.
