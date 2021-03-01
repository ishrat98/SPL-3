#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(HDF5Array)
library(dplyr) 

library(scater)
library(plotly)
library(reshape2)
library(circlize)
library(sSeq)
library(shinycssloaders)
library(ComplexHeatmap)
library(tibble)
library(paletteer)
library(umap)
library(cluster)
library(dplyr)
library(tidyverse)
library(viridisLite)
library(HDF5Array)
library(formattable)
library(DT)
library(colourpicker)
library(shinyWidgets)
library(shinyjs)
library(monocle)
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)
library(TSCAN)
library(scales)
library(viridis)
library(Matrix)
# Define UI for application that draws a histogram
# shinyUI(fluidPage(
# 
#     # Application title
#     titlePanel("Old Faithful Geyser Data"),
# 
#     # Sidebar with a slider input for number of bins
#     sidebarLayout(
#         sidebarPanel(
#             sliderInput("bins",
#                         "Number of bins:",
#                         min = 1,
#                         max = 50,
#                         value = 30)
#         ),
# 
#         # Show a plot of the generated distribution
#         mainPanel(
#             plotOutput("distPlot")
#         )
#     )
# ))

cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="F:/SPL-3/updated app/cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")
ui <- dashboardPage(
  #skin = "purple",
  dashboardHeader(
    title = "UoM Single cell"
  ),
  
  # Sidebar #############################
  dashboardSidebar(
    width = 230,
    
    sidebarMenu(
      menuItem("Home", tabName = "overview", icon = icon("home")),
      menuItem("Projection", tabName = "projection", icon = icon("th")),
      menuItem("Summary", tabName = "Summary", icon = icon("align-justify")),
      menuItem("Gene expression", tabName = "Gene_expressionAll", icon = icon("dna"),
               menuSubItem('Single gene expression', tabName = "Gene_expression"),
               menuSubItem('Multiple gene expression', tabName = "Gene_expressionMultiple")),
      menuItem("Highly expressed genes", tabName = "HEG", icon = icon("filter")),
      menuItem("Marker genes", tabName = "MarkerGenes", icon = icon("hornbill")),
      menuItem("Enriched pathway", tabName = "Enriched_pathway", icon = icon("hubspot")),
      menuItem("Heatmap", tabName = "multiple_cluster_heatmap", icon = icon("columns"),
               menuSubItem('Cluster heatmap', tabName = 'clusterHeatmap'),
               menuSubItem('Sample heatmap', tabName = 'sampleHeatmap'),
               menuSubItem('Sample & cluster heatmap', tabName = 'sampleClusterHeatmap')),
      menuItem("Bubble plot", tabName = "bubblePlot", icon = icon("dot-circle"),
               menuSubItem('Cluster bubble plot', tabName = 'clusterBubblePlot'),
               menuSubItem('Sample bubble plot', tabName = 'sampleBubblePlot')),
      menuItem('Differential Expression', tabName = 'DE_menus', icon = icon('line-chart'), 
               menuSubItem('DE between clusters', tabName = 'DE_between_clusters'),
               menuSubItem('DE between samples', tabName = 'DE_between_samples'),
               menuSubItem('DE between sample & clusters', tabName = 'DE_between_sample_and_clusters'),
               menuSubItem('DE between manual selection', tabName = 'DE_between_manual_selection')),
      menuItem('Trajectory', tabName = 'trajectory', icon = icon('route'),
               menuSubItem('FirstLook', tabName = 'trajectory_FirstLook'),
               menuSubItem('Slingshot', tabName = 'trajectory_slingshot'),
               menuSubItem('Monocle', tabName = 'trajectory_monocle'),
               menuSubItem('Monocle3', tabName = 'trajectory_monocle3'),
               menuSubItem('TSCAN', tabName = 'trajectory_TSCAN'),
               menuSubItem('Slicer', tabName = 'trajectory_slicer')),
      menuItem('Analysis info', tabName = 'analysisInfo', icon = icon('info'))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
        /* logo */
                              .skin-blue .main-header .logo {
                              background-color: #800080;
                              }
                              
                              /* logo when hovered */
                              .skin-blue .main-header .logo:hover {
                              background-color: #800080;
                              }
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: #800080;
                              }        
                              
                              
                              
                              /* other links in the sidebarmenu when hovered */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #800080;
                              }
                              /* toggle button when hovered  */                    
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #800080;
                              }
                              .box.box-solid.box-primary>.box-header {
                               color:#fff;
                               background:#800080
                              }
                              
                              .box.box-solid.box-primary{
                              border-bottom-color:#800080;
                              border-left-color:#800080;
                              border-right-color:#800080;
                              border-top-color:#800080
                              
                              }

                              '))),
    
    tabItems(
      tabItem(tabName = 'overview',
              
              fluidRow(
                column(12,tags$h1('Load data')),
                column(12,
                       fileInput('fileInput', 'HDF5 file location', multiple = FALSE, accept = c(".rds",".crb",".cerebro","h5"),
                                 width = NULL, buttonLabel = "Browse...",
                                 placeholder = "No file selected")),
                
                br(),br(),br(),br(),br(),br(),br(),
                column(4, align="center",
                       box(
                         title = tags$p(style = "font-size: 200%;font-weight: 900;", dim(cdScFiltAnnot)[2]), width=20,  status = "info", "cells"
                       )),
                column(4,align="center",
                       box(
                         title = tags$p(style = "font-size: 200%;font-weight: 900;", nlevels(as.factor(cdScFiltAnnot$Sample))), width=20,  status = "info","samples"
                         
                       )),
                column(4, align="center",
                       box(
                         title = tags$p(style = "font-size: 200%;font-weight: 900;", nlevels(as.factor(cdScFiltAnnot$Clusters))), width=20,status = "info", "clusters"
                         
                       ))
                
                #downloadButton("exportTsne", label = "Download t-SNE"),
                #downloadButton("exportUmap", label = "Download UMAP")
                
              )
              
      ),
      tabItem(tabName = 'projection',
              fluidRow(
                box(
                  title = p(tags$span("Projection", style="padding-right:8px;"), 
                            actionButton("projectionInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  column(width=6, selectInput("projection", "Projection",
                                              choices = c('tSNE','UMAP'),
                                              multiple = FALSE,
                                              selectize = FALSE)),
                  column(width = 6, selectInput("colorCellsBy", "Color cells by",
                                                choices = c('Sample','Cluster','nUMI','nGene','percent_mt','CellType'),
                                                multiple = FALSE,
                                                selectize = FALSE)),
                  column(width = 6, sliderInput("plotOverviewDotSize", "Dot size:", 0, 10, 3.5, 0.5)),
                  column(width = 6, sliderInput("plotOverviewDotOpacity", "Dot opacity:", 0, 1, 1, 0.1)),
                  column(width=12,plotlyOutput("tsnePlotCluster", width = "100%") %>% withSpinner(type = getOption("spinner.type", default = 8))),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionSample", 
                    label = "Select/deselect samples", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionCluster", 
                    label = "Select/deselect clusters", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxProjectionCellType", 
                    label = "Select/deselect celltype", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ))
                )
              )
      ),
      tabItem(tabName = 'Summary',
              fluidRow(
                box(
                  title = p(tags$span("Samples by clusters", style="padding-right:8px;"), 
                            actionButton("sampleByclustersInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("tableSampleClust") %>% withSpinner(type = getOption("spinner.type", default = 8)) %>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = p(tags$span("Percent of clusters in sample", style="padding-right:8px;"), 
                            actionButton("percentClusterInSampleInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("sampleClusterProp", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = p(tags$span("UMI in sample", style="padding-right:8px;"), 
                            actionButton("umiinSampleInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("samplenUMI", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = p(tags$span("Number of genes expressed", style="padding-right:8px;"), 
                            actionButton("numberOfGenesExpressedInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ),status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("samplenGene", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = p(tags$span("UMI in clusters", style="padding-right:8px;"), 
                            actionButton("umiInClustersInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("clusternUMI", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = p(tags$span("Number of genes expressed in clusters", style="padding-right:8px;"), 
                            actionButton("numberOfGenesExpressedInClustersInfo", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                  ), status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("clusternGene", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
              
      ),
      tabItem(tabName = 'HEG',
              fluidRow(
                box(title='Highly expressed genes for samples', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    selectInput('hegChooseSample', 'Sample', 
                                choices = levels(as.factor(cdScFiltAnnot$Sample)), 
                                selected = '', multiple = FALSE,
                                selectize = FALSE),
                    DT::dataTableOutput("hegTableSample") %>% withSpinner(type = getOption("spinner.type", default = 8))%>% withSpinner(type = getOption("spinner.type", default = 8))
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Highly expressed genes for clusters", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  selectInput('hegChooseCluster', 'Cluster', 
                              choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                              selected = '', multiple = FALSE,
                              selectize = FALSE),
                  DT::dataTableOutput("hegTableCluster") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'MarkerGenes',
              fluidRow(
                box(title='Marker genes for samples', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    selectInput('markerChooseSample', 'Sample', 
                                choices = levels(as.factor(cdScFiltAnnot$Sample)), 
                                selected = '', multiple = FALSE,
                                selectize = FALSE),
                    DT::dataTableOutput("markerTableSample") %>% withSpinner(type = getOption("spinner.type", default = 8))
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Marker genes for clusters", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  selectInput('markerChooseCluster', 'Cluster', 
                              choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                              selected = '', multiple = FALSE,
                              selectize = FALSE),
                  DT::dataTableOutput("markerTableCluster") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'Enriched_pathway',
              fluidRow(
                box(title='Pathway enrichment for samples', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    column(width=6, selectInput('ENChooseSample', 'Sample', 
                                                choices = levels(as.factor(cdScFiltAnnot$Sample)), 
                                                selected = '', multiple = FALSE,
                                                selectize = FALSE)),
                    column(width=6,selectInput('ENChooseSampleGO', 'Sample', 
                                               choices = levels(as.factor(cdScFiltAnnot$Sample)), 
                                               selected = '', multiple = FALSE,
                                               selectize = FALSE)),
                    column(width=12,DT::dataTableOutput("ENTableSample") %>% withSpinner(type = getOption("spinner.type", default = 8)))
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Pathway enrichment for clusters", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  column(width=6,selectInput('ENChooseCluster', 'Cluster', 
                                             choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                                             multiple = FALSE,
                                             selectize = FALSE)),
                  # column(width=6, selectInput('ENChooseCluster', 'Cluster', 
                  #                             choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                  #                             selected = '', multiple = FALSE,
                  #                             selectize = FALSE)),
                  column(width = 6, selectInput('ENChooseClusterGO', 'Cluster', 
                                                choices = levels(as.factor(my.clusters[order(my.clusters)])), 
                                                selected = '', multiple = FALSE,
                                                selectize = FALSE)),
                  column(width=12,DT::dataTableOutput("ENTableCluster") %>% withSpinner(type = getOption("spinner.type", default = 8)))
                )
              )
      ),
      tabItem(tabName = 'Gene_expression',
              fluidRow(
                box(
                  title = "Gene expression", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  column(width=6,selectInput("geneName", "Gene Name", selected = GeneNameSorted[1], 
                                             choices = GeneNameSorted,
                                             multiple = FALSE,
                                             selectize = FALSE)),
                  column(width=6, selectInput("geneExprProjection", "Projection",
                                              choices = c('tSNE','UMAP'),
                                              multiple = FALSE,
                                              selectize = FALSE)),
                  column(width=6, sliderInput("geneExpressionplotOverviewDotSize", "Dot size:", 0, 10, 0.5, 0.5)),
                  column(width=6,sliderInput("geneExpressionplotOverviewDotOpacity", "Dot opacity:", 0, 1, 1, 0.1)),
                  column(width=6, colourpicker::colourInput("colmaxgeneExp", "Select colour for maximum value", "firebrick1")),
                  column(width=6, colourpicker::colourInput("colmingeneExp", "Select colour for minimum", "gray88")),
                  column(width=12,plotOutput("PlotGeneExpr", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionSample", 
                    label = "Select/deselect samples", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionCluster", 
                    label = "Select/deselect clusters", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  )),
                  column(width = 4, pickerInput(
                    inputId = "checkboxGeneExpressionCellType", 
                    label = "Select/deselect celltype", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ))
                ),
                box(
                  title = "Violin plot sample", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("violinPlot", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                box(
                  title = "Violin plot cluster", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("violinPlotClusterOrig", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'Gene_expressionMultiple',
              fluidRow(
                box(
                  title = "Multiple gene expression", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 4,
                  selectizeInput("geneNameMultiple", "List of genes", 
                                 selected = '', 
                                 choices = GeneNameSorted,
                                 multiple = TRUE),
                  selectInput("geneExprProjectionMultiple", "Projection",
                              choices = c('tSNE','UMAP'),
                              multiple = FALSE,
                              selectize = FALSE),
                  sliderInput("geneExpressionplotOverviewDotSizeMultiple", "Dot size:", 0, 10, 0.5, 0.5),
                  sliderInput("geneExpressionplotOverviewDotOpacityMultiple", "Dot opacity:", 0, 1, 1, 0.1),
                  colourpicker::colourInput("colmaxgeneExpMultiple", "Select colour for maximum value", "firebrick1"),
                  colourpicker::colourInput("colmingeneExpMultiple", "Select colour for minimum", "gray88"),
                  pickerInput(
                    inputId = "checkboxGeneExpressionSampleMultiple", 
                    label = "Select/deselect samples", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ),
                  pickerInput(
                    inputId = "checkboxGeneExpressionClusterMultiple", 
                    label = "Select/deselect clusters", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ),
                  pickerInput(
                    inputId = "checkboxGeneExpressionCellTypeMultiple", 
                    label = "Select/deselect celltype", 
                    choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                    selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                    options = list(
                      `actions-box` = TRUE, 
                      size = 10,
                      `selected-text-format` = "count > 20"
                    ), 
                    multiple = TRUE
                  ),
                  tags$form(
                    actionButton("buttonForMultipleGeneExpresion", "Visualize gene expression", styleclass = "primary")
                  )
                ),
                box(
                  title = "Multiple gene panel", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("PlotGeneExprMultiple", width = "100%") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'clusterHeatmap',
              fluidRow(
                box(title='Input for cluster heatmap', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseClustersHeatmap', 'Select your clusters', 
                                   choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('ClusterGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForClusterHeatmap", "Generate Heatmap", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Heatmap", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotClusterHeatmap", height = "600px"),
                  column(width=4, colourpicker::colourInput("colmaxClustHeatmap", "Maximum value", "red")),
                  column(width=4, colourpicker::colourInput("colmidClustHeatmap", "Middle value", "white")),
                  column(width=4, colourpicker::colourInput("colminClustHeatmap", "Minimum value", "blue"))
                )
              )
      ),
      tabItem(tabName = 'sampleHeatmap',
              fluidRow(
                box(title='Input for sample heatmap', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseSampleHeatmap', 'Select your samples', 
                                   choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('SampleGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForSampleHeatmap", "Generate Heatmap", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Heatmap", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotSampleHeatmap", height = "600px"),
                  column(width=4, colourpicker::colourInput("colmaxSampleHeatmap", "Maximum value", "red")),
                  column(width=4, colourpicker::colourInput("colmidSampleHeatmap", "Middle value", "white")),
                  column(width=4, colourpicker::colourInput("colminSampleHeatmap", "Minimum value", "blue"))
                )
              )
      ),
      tabItem(tabName = 'sampleClusterHeatmap',
              fluidRow(
                box(title='Input for sample & cluster heatmap', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseSampleClusterHeatmapSample', 'Select your samples', 
                                   choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    selectizeInput('ChooseSampleClusterHeatmapCluster', 'Select your clusters', 
                                   choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('SampleClusterGeneList', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForSampleClusterHeatmap", "Generate Heatmap", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Heatmap", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotSampleClusterHeatmap", height = "600px"),
                  column(width=4, colourpicker::colourInput("colmaxSampleClusterHeatmap", "Maximum value", "red")),
                  column(width=4, colourpicker::colourInput("colmidSampleClusterHeatmap", "Middle value", "white")),
                  column(width=4, colourpicker::colourInput("colminSampleClusterHeatmap", "Minimum value", "blue"))
                )
              )
      ),
      tabItem(tabName = 'clusterBubblePlot',
              fluidRow(
                box(title='Input for cluster bubbleplot', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseClustersBubblePlot', 'Select your clusters', 
                                   choices = as.factor(cdScFiltAnnot$Clusters[order(cdScFiltAnnot$Clusters)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('ClusterGeneListClusterBubblePlot', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForClusterBubbleplot", "Generate Bubbleplot", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Bubble plot", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotClusterBubbleplot", height = "600px")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'sampleBubblePlot',
              fluidRow(
                box(title='Input for sample bubbleplot', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 4,
                    selectizeInput('ChooseSampleBubblePlot', 'Select your samples', 
                                   choices = as.factor(cdScFiltAnnot$Sample[order(cdScFiltAnnot$Sample)]), 
                                   selected = '', 
                                   multiple = TRUE),
                    textAreaInput('sampleGeneListSampleBubblePlot', 'Paste your Gene list', value = '', width = "180px", height = "250px"),
                    tags$form(
                      actionButton("buttonForSampleBubbleplot", "Generate Bubbleplot", styleclass = "primary")
                    )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "Bubble plot", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 8,
                  plotOutput("plotSampleBubbleplot", height = "600px")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_menus',
              h2('Selected Sub-Item Four')
      ),
      tabItem(tabName = 'DE_between_clusters',
              fluidRow(
                box(title='Input for DE between clusters', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    column(width = 6,
                           selectizeInput("clust1", "Select Cluster/Clusters for Group1",
                                          choices = as.factor(my.clusters[order(my.clusters)]),
                                          multiple = TRUE)),
                    column(width = 6,
                           selectizeInput("clust2", "Select Cluster/Clusters for Group2",
                                          choices = as.factor(my.clusters[order(my.clusters)]),
                                          multiple = TRUE)),
                    column(width = 4,
                           tags$form(
                             actionButton("buttonForTwoClustDE", "Run DE", styleclass = "primary")
                           )),
                    column(width = 12,
                           plotlyOutput("AllClustTwoClustComp") %>% withSpinner(type = getOption("spinner.type", default = 8))),
                    column(width = 12,
                           plotOutput("plotSelectClust") %>% withSpinner(type = getOption("spinner.type", default = 8)) )
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "DE results", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("mytableTwoClust") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_between_samples',
              fluidRow(
                box(title='Input for DE between samples', status = "primary", solidHeader = TRUE,
                    collapsible = TRUE, width = 12,
                    column(width = 6,
                           selectInput("sample1", "Select Sample-1",
                                       choices = as.factor(colData(cdScFiltAnnot)[,'Sample']), multiple = TRUE)),
                    column(width = 6,
                           selectInput("sample2", "Select Sample-2",
                                       choices = as.factor(colData(cdScFiltAnnot)[,'Sample']), multiple = TRUE)),
                    column(width = 4,
                           tags$form(
                             actionButton("buttonForTwoSampleDE", "Run DE for two Samples", styleclass = "primary")
                           )),
                    column(width = 12,
                           plotlyOutput("AllClustTwoSampleComp") %>% withSpinner(type = getOption("spinner.type", default = 8)) ),
                    column(width = 12,
                           plotOutput("plotSelectSample") %>% withSpinner(type = getOption("spinner.type", default = 8)))
                    
                    #downloadButton("exportTsne", label = "Download t-SNE"),
                    #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "DE results", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("mytableTwoSample") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_between_sample_and_clusters',
              fluidRow(
                box(title = p(tags$span("Input selection", style="padding-right:8px;"), 
                              actionButton("DE_between_sample_and_clustersSelectionInfo", "help", 
                                           class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 4,
                pickerInput(
                  inputId = "checkboxMultiSelectionSampleSel1", 
                  label = "Selection-1 samples", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  #selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  selected = '', 
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                ),
                pickerInput(
                  inputId = "checkboxMultiSelectionClusterSel1", 
                  label = "Selection1 clusters", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                  #selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                  selected = '',
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                ),
                pickerInput(
                  inputId = "checkboxMultiSelectionSampleSel2", 
                  label = "Selection-2 samples", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  #selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  selected = '', 
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                ),
                pickerInput(
                  inputId = "checkboxMultiSelectionClusterSel2", 
                  label = "Selection-2 clusters", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                  #selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                  selected = '',
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                ),
                selectInput("selectProjectionMixtureSelection", "Projection by",
                            choices = c('tSNE','UMAP'),
                            multiple = FALSE,
                            selectize = FALSE),
                selectInput("colorCellsByMixtureSelection", "Color cells by",
                            choices = c('Sample','Cluster'),
                            multiple = FALSE,
                            selectize = FALSE),
                
                tags$form(
                  actionButton("buttonForDEMixedSelection", "Run DE", styleclass = "primary")
                )
                ),
                #downloadButton("exportTsne", label = "Download t-SNE"),
                #downloadButton("exportUmap", label = "Download UMAP")
                box(title = p(tags$span("Cell projection", style="padding-right:8px;"), 
                              actionButton("DE_between_sample_and_clustersProjectionInfo", "help", 
                                           class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE, width = 8,
                plotOutput("AllClustMixedSelection") %>% withSpinner(type = getOption("spinner.type", default = 8)) 
                ),
                box(title = p(tags$span("Selected cells", style="padding-right:8px;"), 
                              actionButton("DE_between_sample_and_clustersSelectinCellsInfo", "help", 
                                           class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE, width = 8,
                plotOutput("selectedCellsClustMixedSelection") %>% withSpinner(type = getOption("spinner.type", default = 8)) 
                ),
                box(title = "DE results", status = "primary", solidHeader = TRUE, 
                    collapsible = TRUE, width = 12,
                    DT::dataTableOutput("mytableMixedSelection") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      tabItem(tabName = 'DE_between_manual_selection',
              fluidRow(
                box(title=p(tags$span("Manual selection of cells", style="padding-right:8px;"), 
                            actionButton("DE_in_manual_selection_input", "help", 
                                         class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                column(width = 12, h5('Manually select the two groups of cells by dragging your mouse over the cells. 
                                         Then click the button')), 
                column(width = 6,       
                       tags$form(
                         actionButton("buttonForDEManualSelection", "Run DE between selected cells", styleclass = "primary")
                       )),
                column(width=6, useShinyjs(),
                       actionButton("reset", "Reset all selection")),
                column(width = 4, pickerInput(
                  inputId = "checkboxMultiSelectionSample", 
                  label = "Select/deselect samples", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                )),
                column(width = 4, pickerInput(
                  inputId = "checkboxMultiSelectionCluster", 
                  label = "Select/deselect clusters", 
                  choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))), 
                  selected = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                  options = list(
                    `actions-box` = TRUE, 
                    size = 10,
                    `selected-text-format` = "count > 20"
                  ), 
                  multiple = TRUE
                )),
                column(width = 4, selectInput("colorCellsByManualSelection", "Color cells by",
                                              choices = c('Sample','Cluster'),
                                              multiple = FALSE,
                                              selectize = FALSE)),
                column(width = 12,
                       plotlyOutput("AllClustManualSelection") %>% withSpinner(type = getOption("spinner.type", default = 8)) )
                
                #downloadButton("exportTsne", label = "Download t-SNE"),
                #downloadButton("exportUmap", label = "Download UMAP")
                ),
                
                box(
                  title = "DE results", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  DT::dataTableOutput("mytableManualSelection") %>% withSpinner(type = getOption("spinner.type", default = 8))
                )
              )
      ),
      
      tabItem(tabName = 'trajectory_FirstLook',
              box(
                title = "PC1", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_FirstLookOT1", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "First Principal Component", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("FirstPrincipalComponent", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      tabItem(tabName = 'trajectory_slingshot',
              box(
                title = "Slingshot", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_slingshotOT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "First Slingshot Psedutime", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_FirstSlingshot", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "Second Slingshot Psedutime", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_SecondSlingshot", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      
      tabItem(tabName = 'trajectory_monocle',
              fluidRow(
                box(title = p(tags$span("Input parameters", style="padding-right:8px;"), 
                              actionButton("trajectory_input", "info",
                                           class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 4,
                selectInput(
                  inputId = "trajectory_projection_info",
                  label = "Trajectory", 
                  #choices = names(cdScFiltAnnot$trajectory$monocle2),
                  # choices = unique(levels(as.factor(cdScFiltAnnot$trajectory$monocle2))),
                  # choices = names(as.factor(cdScFiltAnnot$trajectory$monocle2)
                  # choices = unique(levels(as.factor(cdScFiltAnnot$samples))), 
                  # selected = unique(levels(as.factor(cdScFiltAnnot$samples))),
                  choices = c('tSNE','UMAP'),
                  multiple = FALSE,
                  selectize = FALSE
                ),
                pickerInput(
                  inputId = "trajectory_samples_to_display", 
                  label = "Samples to display",
                  # choices = cdScFiltAnnot$sample_names,
                  # selected = cdScFiltAnnot$sample_names,
                  choices = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  #selected = unique(levels(as.factor(cdScFiltAnnot$Sample))), 
                  selected = '',
                  options = list("actions-box" = TRUE),
                  multiple = TRUE
                  
                ),
                pickerInput(
                  inputId = "trajectory_clusters_to_display", 
                  label = "Clusters to display",
                  # choices = cdScFiltAnnot$cluster_names,
                  # selected = cdScFiltAnnot$cluster_names,
                  choices = unique(levels(as.factor(cdScFiltAnnot$Clusters))),
                  selected = '',
                  options = list("actions-box" = TRUE),
                  multiple = TRUE
                  
                ),
                sliderInput(
                  "trajectory_percentage_cells_to_show",
                  label = "Show % of cells",
                  min = scatter_plot_percentage_cells_to_show[["min"]],
                  max = scatter_plot_percentage_cells_to_show[["max"]],
                  step = scatter_plot_percentage_cells_to_show[["step"]],
                  value = scatter_plot_percentage_cells_to_show[["default"]]
                ),
                
                
                selectInput(
                  "trajectory_dot_color",
                  label = "Color cells by",
                  # choices = c("state","pseudotime",names(cdScFiltAnnot$cells)[! names(cdScFiltAnnot$cells) %in% c("cell_barcode")])
                  # choices = unique(levels(as.factor(cdScFiltAnnot$cellType))), 
                  # selected = unique(levels(as.factor(cdScFiltAnnot$cellType))),
                  choices = c('state','pseudotime','Sample','Cluster','nUMI','nGene','percent_mt','CellType'),
                  multiple = FALSE,
                  selectize = FALSE
                ),
                
                sliderInput(
                  "trajectory_dot_size",
                  label = "Dot size",
                  min = scatter_plot_dot_size[["min"]],
                  max = scatter_plot_dot_size[["max"]],
                  step = scatter_plot_dot_size[["step"]],
                  value = scatter_plot_dot_size[["default"]]
                ),
                sliderInput(
                  "trajectory_dot_opacity",
                  label = "Dot opacity",
                  min = scatter_plot_dot_opacity[["min"]],
                  max = scatter_plot_dot_opacity[["max"]],
                  step = scatter_plot_dot_opacity[["step"]],
                  value = scatter_plot_dot_opacity[["default"]]
                )
                
                
                ),
                
                box(title = p(tags$span("Trajectory", style="padding-right:8px;"),
                              actionButton("trajectory_projection_info", "info",
                                           class = "btn-xs", title = "Additional information for this panel")
                ), status = "primary", solidHeader = TRUE, width = 8,
                plotOutput("trajectory_projection") %>% withSpinner(type = getOption("spinner.type", default = 8))
                ),
                
                box(
                  title = "Distribution along pseudotime", status = "primary", solidHeader = TRUE,
                  collapsible = TRUE, width = 12,
                  plotOutput("trajectory_density_plott", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
                )
                
              )
      ),
      
      tabItem(tabName = 'trajectory_monocle3',
              box(
                title = "Monocle2", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_monocle2OT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "Component", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_monocle2Component", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "Monocle Psedutime", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_monocle2Psedutime", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      
      tabItem(tabName = 'trajectory_TSCAN',
              box(
                title = "PCA Dimension 1", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_TSCAN_1", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              ),
              box(
                title = "TSCAN Psedutime", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_TSCAN_Pseudotime", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      tabItem(tabName = 'trajectory_slicer',
              box(
                title = "Slicer", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotlyOutput("trajectory_slicerOT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      tabItem(tabName = 'analysisInfo',
              h2('Analysis details')
              #includeMarkdown("../snRNA_seq_dataset_Hiseq127_JenniferScott/snRNA_seq_dataset_Hiseq127_JenniferScott.md")
      )
      
    )
  )
)

