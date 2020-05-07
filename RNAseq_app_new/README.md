# Abstract

This shiny app allows count matrix data analyzed for differential expression
to be visually explored. Different sample group comparisons resulting from
the DE analysis can be selected from the dropdown menu in each tab of the
user interface.
 
The first tab contains a description of the analysis and a table of
experimental groups used in the differential expression analysis
 
The second tab contains a PCA displaying the projection of the experimental
groups used in the differential expression analysis, projected onto the
first two principal components.
 
The third tab displays a *MA plot*, *Volcano plot*, a data table of results
from the DE analysis and a plot of normalized counts of a selected gene
determined by user click input. The plots and table shown in the user
interface are updated automatically with new user input.
 
The fourth tab displays results of the GO term Enrichment Analysis
in an html table along with a volcano plot that displays the selected GO terms
 
The fifth tab displays results of the KEGG pathway Enrichment Analysis
in an html table along with a volcano plot that displays the selected KEGG
pathways

NOTE: the app expects results in the following format:
de_res_*.rds
go_*.rds
kegg_*.rds

(e.g. DE results from the *first* comparison are to be called "de_res_1.rds",
DE results from the *second* comparison are to be called "de_res_2.rds")

The de_res_*.rds files must have a neg_log10_padj column, which is
calculated from -log10(padj)

The go_*.rds files are expected to have a "ENTREZID_in_term" column, with a
string of ";" separated ENTREZIDs for each GO term

The kegg_*.rds file expected to have a "ENTREZID_in_path" column, with a
string of ";" separated ENTREZIDs for each KEGG pathway

