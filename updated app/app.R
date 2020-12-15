library(shiny)
library(shinydashboard)


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
library(Seurat)
library(TSCAN)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)



Sys.setenv(R_MAX_VSIZE = 16e9)

#cdScFiltAnnot <- loadHDF5SummarizedExperiment(dir="cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnotK <- loadHDF5SummarizedExperiment(dir="cdScFiltAnnotHDF5", prefix="")


cdScFiltAnnot <-  as(cdScFiltAnnotK, "SingleCellExperiment")


cellLabels <- cdScFiltAnnot$cellType
deng <- counts(cdScFiltAnnot)
colnames(deng) <- cellLabels




c30 <- c("dodgerblue2",#1
         "#E31A1C", #2 red
         "green4", #3
         "#FF7F00", #4 orange
         "green1",#5
         "purple",#6
         "blue1",#7
         "deeppink1",#8
         "darkorange4",#9
         "black",#10
         "gold1",#11
         "darkturquoise",#12
         "#6A3D9A", #13 purple
         "orchid1",#14
         "gray70",#15
         "maroon",#16
         "palegreen2",#17
         "#333333",#18
         "#CAB2D6", #19 lt purple
         "#FDBF6F", #20 lt orange
         "khaki2",#21
         "skyblue2",#22
         "steelblue4",#23
         "green1",#24
         "yellow4",#25
         "yellow3",#26
         "#FB9A99", #27 lt pink
         "brown",#28
         "#000099",#29
         "#CC3300"#30
)



c_sample_col <- c30[c(1,3,23,19,30)]
c_clust_col <- c30[c(1,2,3,4,5,6,7,8,9,11,12,14,19,22,24,25)]

my.clusters <- cdScFiltAnnot$Clusters
#options(bitmapType='cairo')


df_shiny <<- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
df_shiny[,'Clusters'] <- as.factor(cdScFiltAnnot$Clusters)
df_shiny$Group <- "gray"
df_shiny$key <- row.names(df_shiny)


df_shiny_ForTwoClust <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
rownames(df_shiny_ForTwoClust) <- cdScFiltAnnot$Barcode
df_shiny_ForTwoClust$Group <- "gray88"
df_shiny_ForTwoClust$key <- row.names(df_shiny_ForTwoClust)


df_shiny_ForTwoSample <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
rownames(df_shiny_ForTwoSample) <- cdScFiltAnnot$Barcode
df_shiny_ForTwoSample$Group <- "gray88"
df_shiny_ForTwoSample$key <- row.names(df_shiny_ForTwoSample)

#cbPalette <- paletteer_c(package="pals", palette="glasbey", n=12)

colData(cdScFiltAnnot)[,'Clust_Sample'] <- paste0(colData(cdScFiltAnnot)$Clusters,'_',colData(cdScFiltAnnot)$Sample)

colData(cdScFiltAnnot)[,'Clust_SampleWithClust'] <- paste0(colData(cdScFiltAnnot)$Sample,'_',colData(cdScFiltAnnot)$Clusters)


GeneNameSorted <- sort(rownames(logcounts(cdScFiltAnnot)))
df = data.frame(Condition=colData(cdScFiltAnnot)$Sample)
#df <- df[order(df$Class),]
ha = HeatmapAnnotation(df = df)

flag <<- 0

## centre parametrer in trajectory


scatter_plot_dot_size <- list(
  min = 1,
  max = 20,
  step = 1,
  default = 5
)

scatter_plot_dot_opacity <- list (
  min = 0.1,
  max = 1.0,
  step = 0.1,
  default = 1.0
)
scatter_plot_point_opacity <- list(

  min = 0.1,
  max = 1.0,
  step = 0.1,
  default = 1.0
)

scatter_plot_percentage_cells_to_show <- list(
  min = 10,
  max = 100,
  step = 10,
  default = 100
)



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
               menuSubItem('Tscan', tabName = 'trajectory_Tscan'),
               menuSubItem('Slingshot', tabName = 'trajectory_slingshot'),
               menuSubItem('Monocle', tabName = 'trajectory_monocle'),
               menuSubItem('Monocle3', tabName = 'trajectory_monocle3'),
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
      
      tabItem(tabName = 'trajectory_Tscan',
              box(
                title = "TSCAN", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_TscanOT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      tabItem(tabName = 'trajectory_slingshot',
              box(
                title = "Slingshot", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_slingshotOT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
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
                title = "Monocle3", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_monocle3OT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )

      ),
      
      tabItem(tabName = 'trajectory_slicer',
              box(
                title = "Slicer", status = "primary", solidHeader = TRUE,
                collapsible = TRUE, width = 12,
                plotOutput("trajectory_slicerOT", width = "100%")%>% withSpinner(type = getOption("spinner.type", default = 8))
              )
              
      ),
      
      tabItem(tabName = 'analysisInfo',
              h2('Analysis details')
              #includeMarkdown("../snRNA_seq_dataset_Hiseq127_JenniferScott/snRNA_seq_dataset_Hiseq127_JenniferScott.md")
      )
      
    )
  )
)

server <- function(input, output, session) { 
  
  #################################
  #
  # For overview panel
  #
  #################################
  
 
  
  output$violinPlot <- renderPlot({
    
    geneListToTest <- input$geneName
    #cdScFiltAnnotTmp <- cdScFiltAnnot
    dfViolin <- data.frame(Sample=colData(cdScFiltAnnot)$Sample, logCounts=logcounts(cdScFiltAnnot)[geneListToTest,], title=geneListToTest)
    p <- ggplot(dfViolin, aes(factor(Sample), logCounts)) +
      geom_violin(aes(fill=factor(Sample)), scale="width", trim=TRUE) +
      xlab('Sample') +
      ylab('Expression (logcounts)') +
      scale_fill_manual(values = c_sample_col) +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) + 
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Samples")) + 
      facet_grid(. ~ title)
    p
    
    
  })
  
  
  output$violinPlotClusterOrig <- renderPlot({
    
    geneListToTest <- input$geneName
    
    
    dfViolin <- data.frame(Cluster=colData(cdScFiltAnnot)$Clusters, logCounts=logcounts(cdScFiltAnnot)[geneListToTest,], title=geneListToTest)
    p <- ggplot(dfViolin, aes(factor(Cluster), logCounts)) +
      geom_violin(aes(fill=factor(Cluster)), scale="width", trim=TRUE) +
      xlab('Cluster') +
      ylab('Expression (logcounts)') +
      scale_fill_manual(values = c_clust_col) +
      #      scale_fill_manual(values = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF"))+
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Clusters")) + 
      facet_grid(. ~ title)
    p
    
  })
  
  
  plot_tsnePlotCluster <- function(){ 	
    
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_tsnePlotCellType <- function(){ 	
    
    
  
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~cellType, colors = c30, type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_tsnePlotSample <- function(){ 	
    
  
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Cell']=as.factor(colData(cdScFiltAnnot)$Barcode)
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'cellType'] <- as.factor(colData(cdScFiltAnnot)$cellType)
    df %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Sample, colors = c_sample_col[c(1:4)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'tSNE with Samples',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
  }
  plot_tsnePlotnUMI <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nUMI = log10(cdScFiltAnnot$total)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nUMI, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with total UMI in Log10',showlegend = FALSE)
    
  }
  plot_tsnePlotnGene <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nGenes = log10(cdScFiltAnnot$detected)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nGenes, type="scatter", mode="markers", hoverinfo = 'text' , 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with total gene in log10',showlegend = FALSE)
    
  }
  
  plot_tsnePlotPct_MT <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(pct_mt = cdScFiltAnnot$subsets_Mt_percent) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~pct_mt, type="scatter", mode="markers", hoverinfo = 'text' , 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'tSNE with pct_mt',showlegend = FALSE)
    
  }
  
  
  
  plot_UMAPPlotCluster <- function(){    
    
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~cellType, colors = c_clust_col[c(1:8)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_UMAPPlotCellType <- function(){    
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Clusters, colors = c_clust_col[c(1:8)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  }
  
  plot_UMAPPlotSample <- function(){    
    
    

    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      filter(Clusters %in% input$checkboxProjectionCluster) %>%
      filter(Sample %in% input$checkboxProjectionSample) %>%
      filter(cellType %in% input$checkboxProjectionCellType) %>%
      plot_ly(color = ~Sample, colors = c_sample_col[c(1:4)], type="scatter", mode="markers", hoverinfo = 'text',
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Sample,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with Samples',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
  }
  plot_UMAPPlotnUMI <- function(){    
    
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nUMI = log10(cdScFiltAnnot$total)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nUMI, type="scatter", mode="markers", hoverinfo = 'text' ,
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>%
      layout(title = 'UMAP with total UMI in log10',showlegend = FALSE)
  }
  
  plot_UMAPPlotnGene <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(log10_nGenes = log10(cdScFiltAnnot$detected)) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~log10_nGenes, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'UMAP with total gene in log10',showlegend = FALSE)
    
  }
  
  plot_UMAPPlotPercentMt <- function(){ 	
    
    as.tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Samples = colData(cdScFiltAnnot)$Sample) %>%
      mutate(Cell = colData(cdScFiltAnnot)$Barcode) %>%
      mutate(Clusters = as.factor(my.clusters)) %>%
      mutate(pct_mt = cdScFiltAnnot$subsets_Mt_percent) %>%
      mutate(cellType = colData(cdScFiltAnnot)$cellType) %>%
      plot_ly(color = ~pct_mt, type="scatter", mode="markers", hoverinfo = 'text', 
              marker = list(size = input$plotOverviewDotSize, opacity = input$plotOverviewDotOpacity),
              text = ~paste('</br> Cell: ', Cell,
                            '</br> Clusters: ', Clusters,
                            '</br> Samples: ', Samples,
                            '</br> CellType: ', cellType)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2) %>% 
      layout(title = 'UMAP with pct_mt',showlegend = FALSE)
    
  }
  
  
  output$PlotGeneExpr <- renderPlot({
    if(input$geneExprProjection == 'tSNE')
      plot_tsnePlotGene()
    else if(input$geneExprProjection == 'UMAP')
      plot_UMAPPlotGene()
  })
  
  
  output$tsnePlotCluster <- renderPlotly({  
    if(input$projection == 'tSNE' & input$colorCellsBy == 'Cluster')
      plot_tsnePlotCluster()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'Sample')
      plot_tsnePlotSample()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'CellType')
      plot_tsnePlotCellType()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'nUMI')
      plot_tsnePlotnUMI()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'nGene')
      plot_tsnePlotnGene()
    else if(input$projection == 'tSNE' & input$colorCellsBy == 'percent_mt')
      plot_tsnePlotPct_MT()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'Cluster')
      plot_UMAPPlotCluster()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'Sample')
      plot_UMAPPlotSample()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'CellType')
      plot_UMAPPlotCellType()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'nUMI')
      plot_UMAPPlotnUMI()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'nGene')
      plot_UMAPPlotnGene()
    else if(input$projection == 'UMAP' & input$colorCellsBy == 'percent_mt')
      plot_UMAPPlotPercentMt()
  })
  
  output$tsnePlotSample <- renderPlotly({   
    if(input$projection == 'tSNE')
      plot_tsnePlotSample()
    else
      plot_UMAPPlotSample()
  })
  
  
  output$checkboxProjectionSample <- renderUI({
    choice <-  unique(levels(as.factor(cdScFiltAnnot$Sample)))
    pickerInput(
      inputId = "checkboxProjectionSample", 
      label = "Select/deselect samples for display", 
      choices = choice, 
      options = list(
        `actions-box` = TRUE, 
        size = 10,
        `selected-text-format` = "count > 3"
      ), 
      multiple = TRUE
    )
    
    
  })
  
  output$checkboxProjectionCluster <- renderUI({
    choice <-  unique(levels(cdScFiltAnnot$Clusters))
    pickerInput(
      inputId = "checkboxProjectionCluster", 
      label = "Select/deselect clusters for display", 
      choices = choice, 
      options = list(
        `actions-box` = TRUE, 
        size = 10,
        `selected-text-format` = "count > 3"
      ), 
      multiple = TRUE
    )
    
  })
  
  ############################################
  #
  # For Gene expression panel
  #
  ############################################
  
  plot_tsnePlotGene <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(GeneExp = as.vector(logcounts(cdScFiltAnnot[input$geneName,]))) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      filter(Clusters %in% input$checkboxGeneExpressionCluster) %>%
      filter(Sample %in% input$checkboxGeneExpressionSample) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellType) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSize,alpha=input$geneExpressionplotOverviewDotOpacity, aes(colour = GeneExp)) +
      #scale_colour_viridis_c()+
      scale_colour_gradientn(colours=c(input$colmingeneExp, input$colmaxgeneExp),
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('tSNE Gene Exp:',input$geneName))+
      theme_classic(base_size=14) 
    # theme(strip.background = element_blank(),
    #       strip.text.x     = element_blank(),
    #       axis.text.x      = element_blank(),
    #       axis.text.y      = element_blank(),
    #       axis.ticks       = element_blank(),
    #       axis.line        = element_blank(),
    #       panel.border     = element_blank())
    
  }
  
  
  plot_tsnePlotGeneFacetGene <- function(){
    
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    GeneExp=logcounts(cdScFiltAnnot)[input$geneName,]
    Sample = as.factor(colData(cdScFiltAnnot)$Sample)
    
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_gradientn(colours=c("#FDBB84","#FC8D59","#E34A33","#B30000"))
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Sample, ncol=4,nrow=3)
    
  }
  
  
  plot_tsnePlotGeneFacetSample <- function(){
    
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    GeneExp=logcounts(cdScFiltAnnot)[input$geneName,]
    Sample = as.factor(colData(cdScFiltAnnot)$Sample)
    
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_gradientn(colours=c("#FDBB84","#FC8D59","#E34A33","#B30000"))
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Sample, ncol=4,nrow=3)
    
  }
  
  plot_tsnePlotGeneFacetCluster <- function(){
    
    #GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    #GeneName = 'SPN'
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    #df[,'GeneExp']=logcounts(cdScFiltAnnot)[input$geneName,]
    #df[,'Clusters'] <- as.factor(my.clusters)
    GeneExp <- logcounts(cdScFiltAnnot)[input$geneName,]
    as.tibble(df) %>%
      dplyr::mutate("GeneExp" = GeneExp) %>% 
      dplyr::mutate("Clusters" = as.factor(my.clusters)) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=0.6,alpha=1, aes(colour = GeneExp)) +
      scale_colour_gradientn(colours=c("gray88", "red"))+
      #scale_colour_viridis_c()+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('Gene Exp:',input$geneName)) +
      theme_classic(base_size=14) +  
      theme(#strip.background = element_blank(),
        #strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
        panel.spacing.x=unit(0, "lines"), 
        panel.spacing.y=unit(0,"lines"),
        #legend.position = ("none")) +
        panel.background =  element_blank(), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank()) + 
      facet_wrap(~Clusters, ncol=3,nrow=4)
    
  }
  
  
  plot_UMAPPlotGene <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(GeneExp = as.vector(logcounts(cdScFiltAnnot[input$geneName,]))) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      filter(Clusters %in% input$checkboxGeneExpressionCluster) %>%
      filter(Sample %in% input$checkboxGeneExpressionSample) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellType) %>%
      ggplot(aes(x=V1, y=V2, GeneExp = GeneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSize,alpha=input$geneExpressionplotOverviewDotOpacity, aes(colour = GeneExp)) +
      #scale_colour_viridis_c()+
      scale_colour_gradientn(colours=c(input$colmingeneExp, input$colmaxgeneExp))+
      #scale_colour_gradientn(colours=c("grey90", "orangered", "orangered4"))+
      #guides(colour = guide_legend(override.aes = list(size=4))) +
      xlab("") + ylab("") +
      ggtitle(paste0('UMAP Gene Exp:',input$geneName))+
      theme_classic(base_size=14) 
    # theme(strip.background = element_blank(),
    #       strip.text.x     = element_blank(),
    #       axis.text.x      = element_blank(),
    #       axis.text.y      = element_blank(),
    #       axis.ticks       = element_blank(),
    #       axis.line        = element_blank(),
    #       panel.border     = element_blank())
    
  }
  
  
  #######################################
  #
  # Gene expression for multiple genes
  #
  ######################################
  
  plot_tsnePlotGeneMultiple <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'tSNE')) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      bind_cols(as_tibble(t(logcounts(cdScFiltAnnot[input$geneNameMultiple,])))) %>%
      filter(Clusters %in% input$checkboxGeneExpressionClusterMultiple) %>%
      filter(Sample %in% input$checkboxGeneExpressionSampleMultiple) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellTypeMultiple) %>%
      gather(.,key="geneName", value="geneExp", -c(V1,V2,Sample,Clusters,cellType) ) %>%
      ggplot(aes(x=V1, y=V2, geneExp = geneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSizeMultiple,
                 alpha=input$geneExpressionplotOverviewDotOpacityMultiple, 
                 aes(colour = geneExp)) + 
      scale_colour_gradientn(colours=c(input$colmingeneExpMultiple, input$colmaxgeneExpMultiple),
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      xlab("") + ylab("") +
      theme_classic(base_size=14) + 
      facet_wrap(~geneName)
    
  }
  
  
  plot_UMAPPlotGeneMultiple <- function(){
    
    as_tibble(reducedDim(cdScFiltAnnot,'UMAP')) %>%
      mutate(Sample = as.factor(cdScFiltAnnot$Sample)) %>%
      mutate(Clusters = as.factor(cdScFiltAnnot$Clusters)) %>%
      mutate(cellType = as.factor(cdScFiltAnnot$cellType)) %>%
      bind_cols(as_tibble(t(logcounts(cdScFiltAnnot[input$geneNameMultiple,])))) %>%
      filter(Clusters %in% input$checkboxGeneExpressionClusterMultiple) %>%
      filter(Sample %in% input$checkboxGeneExpressionSampleMultiple) %>%
      filter(cellType %in% input$checkboxGeneExpressionCellTypeMultiple) %>%
      gather(.,key="geneName", value="geneExp", -c(V1,V2,Sample,Clusters,cellType) ) %>%
      ggplot(aes(x=V1, y=V2, geneExp = geneExp)) +
      geom_point(size=input$geneExpressionplotOverviewDotSizeMultiple,
                 alpha=input$geneExpressionplotOverviewDotOpacityMultiple, 
                 aes(colour = geneExp)) + 
      scale_colour_gradientn(colours=c(input$colmingeneExpMultiple, input$colmaxgeneExpMultiple),
                             guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      xlab("") + ylab("") +
      theme_classic(base_size=14) + 
      facet_wrap(~geneName)
    
  }
  PlotGeneExprMultipleFunction <- eventReactive(input$buttonForMultipleGeneExpresion, {
    if(input$geneExprProjectionMultiple == 'tSNE')
      plot_tsnePlotGeneMultiple()
    else
      plot_UMAPPlotGeneMultiple()
  })
  
  output$PlotGeneExprMultiple <- renderPlot({   
    PlotGeneExprMultipleFunction()
  })
  
  
  #####################################
  
  
  ####################################
  #
  # Summary panel plots
  #
  ####################################
  
  
  output$sampleClusterProp <- renderPlot({   
    table_samples_by_clusters <- as_tibble(colData(cdScFiltAnnot)) %>%
      group_by(Sample, Clusters) %>%
      summarize(count = n()) %>%
      spread(Clusters, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c('Sample', 'total_cell_count', everything())) 
    
    p1 <- table_samples_by_clusters %>%
      dplyr::select(-c('total_cell_count')) %>%
      reshape2::melt(id.vars = 'Sample') %>%
      mutate(Sample = factor(Sample)) %>%
      ggplot(aes(Sample, value, fill = variable)) +
      geom_bar(position = 'fill', stat = 'identity', show.legend = TRUE) +
      scale_fill_manual(name = 'Cluster', values = c_clust_col) +
      scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
      theme_bw() +
      theme(
        legend.position = 'left',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
    p1
  })
  
  
  
  
  output$samplenUMI <- renderPlot({   
    dfViolin <- data.frame(Sample = cdScFiltAnnot$Sample, nUMI = log10(cdScFiltAnnot$total))
    p <- ggplot(dfViolin, aes(factor(Sample), nUMI)) +
      geom_violin(aes(fill=factor(Sample)), scale="width", trim=TRUE) +
      xlab('Sample') +
      ylab('log10(number of UMIs)') +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Samples")) + scale_fill_manual(values = c_sample_col)
    p
  })
  
  
  output$samplenGene <- renderPlot({   
    dfViolin <- data.frame(Sample = cdScFiltAnnot$Sample, nGene = (cdScFiltAnnot$detected))
    p <- ggplot(dfViolin, aes(factor(Sample), nGene)) +
      geom_violin(aes(fill=factor(Sample)), scale="width", trim=TRUE) +
      xlab('Sample') +
      ylab('Number of geness expressed') +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Samples")) + scale_fill_manual(values = c_sample_col)
    p
  })
  
  
  
  output$clusternUMI <- renderPlot({   
    dfViolin <- data.frame(Cluster = cdScFiltAnnot$Clusters, nUMI = log10(cdScFiltAnnot$total))
    p <- ggplot(dfViolin, aes(factor(Cluster), nUMI)) +
      geom_violin(aes(fill=factor(Cluster)), scale="width", trim=TRUE) +
      xlab('Cluster') +
      ylab('log10(number of UMIs)') +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Cluster")) + scale_fill_manual(values = c_clust_col)
    p
  })
  
  
  output$clusternGene <- renderPlot({   
    dfViolin <- data.frame(Cluster = cdScFiltAnnot$Clusters, nGene = (cdScFiltAnnot$detected))
    p <- ggplot(dfViolin, aes(factor(Cluster), nGene)) +
      geom_violin(aes(fill=factor(Cluster)), scale="width", trim=TRUE) +
      xlab('Cluster') +
      ylab('Number of geness expressed') +
      #      scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
      theme_classic(base_size=14) +
      guides(fill=guide_legend(title="Clusters")) + scale_fill_manual(values = c_clust_col)
    p
  })
  
  
  tableSampleClustRendering<- function(){ 	
    
    table_samples_by_clusters <- as_tibble(colData(cdScFiltAnnot)) %>%
      group_by(Sample, Clusters) %>%
      summarize(count = n()) %>%
      spread(Clusters, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c('Sample', 'total_cell_count', everything())) %>%
      dplyr::rename('Total cell' = 'total_cell_count') %>% datatable(rownames = FALSE)
    
    table_samples_by_clusters
    
  }
  
  
  
  output$tableSampleClust = DT::renderDataTable(tableSampleClustRendering(), server = FALSE, extensions = 'Buttons', 
                                                options = list(dom = 'lBfrtip',
                                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                               pageLength = 5, autoWidth = TRUE))
  
  ####################
  
  #####################################
  #
  # Highly expressed genes
  #
  ####################################
  
  
  
  tableRenderinghegChooseSample <- function(){
    
    
    sampleID <- as.character(input$hegChooseSample)
    dfColName <- paste0('PercentSample',sampleID)
    
    df <- data.frame("Sample"=sampleID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "PercentInSample" = rowData(cdScFiltAnnot)[,dfColName],
                     "PercentTotal"= rowData(cdScFiltAnnot)$PercentTotal)
    df <- df[order(df$PercentInSample, decreasing = TRUE),]
    df <- as.data.frame(df[1:200,])
    row.names(df) <- NULL
    
    sampleIndex <- which(levels(as.factor(cdScFiltAnnot$Sample)) == sampleID)
    as.datatable(formattable(df, list(
      Sample = color_tile('white',c30[sampleIndex]),
      area(col = c(PercentInSample)) ~ normalize_bar("pink", 0.2),
      area(col = c(PercentTotal)) ~ normalize_bar("yellow", 0.2)
    )), colnames = c('% of cells in this sample' = 'PercentInSample',
                     '% of cells in total dataset' = 'PercentTotal'), rownames=FALSE)
    
  }
  
  
  
  output$hegTableSample = DT::renderDataTable(tableRenderinghegChooseSample(), server = FALSE, extensions = 'Buttons', 
                                              options = list(dom = 'lBfrtip',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                             pageLength = 5, autoWidth = TRUE))
  
  
  
  
  tableRenderinghegChooseSCluster <- function(){
    
    clustID <- as.numeric(input$hegChooseCluster)
    dfColName <- paste0('PercentClust',clustID)
    df <- data.frame("Cluster"=clustID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "PercentInClust" = rowData(cdScFiltAnnot)[,dfColName],
                     "PercentTotal"= rowData(cdScFiltAnnot)$PercentTotal)
    df <- df[order(df$PercentInClust, decreasing = TRUE),]
    df <- df[1:200,]
    
    
    as.datatable(formattable(df, list(
      Cluster = color_tile('white',c30[clustID+4]),
      area(col = c(PercentInClust)) ~ normalize_bar("pink", 0.2),
      area(col = c(PercentTotal)) ~ normalize_bar("yellow", 0.2)
    )), colnames = c('% of cells in this cluster' = 'PercentInClust',
                     '% of cells in total dataset' = 'PercentTotal'), rownames=FALSE)
    
    
  }
  
  
  
  output$hegTableCluster = DT::renderDataTable(tableRenderinghegChooseSCluster (), server = FALSE, extensions = 'Buttons', 
                                               options = list(dom = 'lBfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                              pageLength = 5, autoWidth = TRUE))
  
  ####################################
  
  
  
  
  ####################################
  #
  # Marker genes
  #
  ####################################
  
  
  
  
  tableRenderingmarkerTableSample <- function(){
    
    
    sampleID <- as.character(input$markerChooseSample)
    dfColName <- paste0('PercentSample',sampleID)
    df <- data.frame("Sample"=sampleID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "FDR" = metadata(cdScFiltAnnot)[['Sample']][[1]][[sampleID]][rownames(cdScFiltAnnot),'FDR'],
                     "PercentInClust" = rowData(cdScFiltAnnot)[,dfColName]/100)
    df <- as.data.frame(cbind(df, metadata(cdScFiltAnnot)[['Sample']][[1]][[sampleID]][rownames(cdScFiltAnnot),3:dim(metadata(cdScFiltAnnot)[['Sample']][[1]][[sampleID]])[2]]))
    df <- df[order(df$FDR, decreasing = FALSE),]
    df$FDR <- formatC(df$FDR, format = "E", digits = 2)
    
    #df <- df[1:100,]
    #rownames(df) <- NULL
    df %>%
      datatable(colnames = c('% cells in this sample' = 'PercentInClust'), 
                rownames = FALSE,
                caption = 'Table : Marker genes for cluster.') %>%
      #formatRound(columns = 'FDR', digits = 3) %>%
      formatRound(columns=c(colnames(df)[grep('logFC',colnames(df))]), digits=3) %>%
      formatPercentage(columns = '% cells in this sample') %>% 
      formatStyle(names(df[,5:dim(df)[2]]),
                  background = styleColorBar(range(df[,5:dim(df)[2]]), 'lightblue'),
                  backgroundSize = '98% 58%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center')
    
  }
  
  
  
  output$markerTableSample = DT::renderDataTable(tableRenderingmarkerTableSample(), server = FALSE, extensions = 'Buttons', 
                                                 options = list(dom = 'lBfrtip',
                                                                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                                pageLength = 5, autoWidth = TRUE))
  
  
  
  
  tableRenderingmarkerTableCluster <- function(){
    
    
    clustID <- input$markerChooseCluster
    dfColName <- paste0('PercentClust',clustID)
    df <- data.frame("Cluster"=clustID,
                     "Gene"= rownames(cdScFiltAnnot),
                     "FDR" = metadata(cdScFiltAnnot)[['Cluster']][[1]][[clustID]][rownames(cdScFiltAnnot),'FDR'],
                     "PercentInClust" = rowData(cdScFiltAnnot)[,dfColName]/100)
    df <- as.data.frame(cbind(df, metadata(cdScFiltAnnot)[['Cluster']][[1]][[clustID]][rownames(cdScFiltAnnot),3:(max(as.numeric(cdScFiltAnnot$Clusters))+1)]))
    df <- df[order(df$FDR, decreasing = FALSE),]
    df$FDR <- formatC(df$FDR, format = "E", digits = 2)
    
    #df <- df[1:100,]
    #rownames(df) <- NULL
    df %>%
      datatable(colnames = c('% cells in this cluster' = 'PercentInClust'), 
                rownames = FALSE,
                caption = 'Table : Marker genes for cluster.') %>%
      #formatRound(columns = 'FDR', digits = 3) %>%
      formatRound(columns=c(colnames(df)[grep('logFC',colnames(df))]), digits=3) %>%
      formatPercentage(columns = '% cells in this cluster') %>% 
      formatStyle(names(df[,5:dim(df)[2]]),
                  background = styleColorBar(range(df[,5:dim(df)[2]]), 'lightblue'),
                  backgroundSize = '98% 58%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center')
    
    
    
  }
  
  
  
  output$markerTableCluster = DT::renderDataTable(tableRenderingmarkerTableCluster (), server = FALSE, extensions = 'Buttons', 
                                                  options = list(dom = 'lBfrtip',
                                                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                                 pageLength = 5, autoWidth = TRUE))
  
  
  ################################################
  
  
  #####################################
  #
  # Pathway enrichment
  #
  #####################################
  
  observe({
    sampleID <- input$ENChooseSample
    
    # Can also set the label and select items
    updateSelectInput(session, "ENChooseSampleGO",
                      label = "Choose sample enrichr",
                      choices = names(which(sapply(metadata(cdScFiltAnnot)[['enrichSample']][[1]][[sampleID]], NROW) != 0)),
                      selected = names(which(sapply(metadata(cdScFiltAnnot)[['enrichSample']][[1]][[sampleID]], NROW) != 0))[1])
    
    
    clusterID <- as.numeric(input$ENChooseCluster)
    
    # Can also set the label and select items
    updateSelectInput(session, "ENChooseClusterGO",
                      label = "Choose cluster enrichr",
                      choices = names(which(sapply(metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]], NROW) != 0)),
                      selected = names(which(sapply(metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]], NROW) != 0))[1])
  })
  
  
  
  tableRenderingENTableSample <- function(){
    
    
    sampleID <- input$ENChooseSample
    sampleEN <- input$ENChooseSampleGO
    dfColName <- paste0('PercentClust',sampleID)
    df <- metadata(cdScFiltAnnot)[['enrichSample']][[1]][[sampleID]][[sampleEN]]  
    if(df$P.value <0.001)
      df$P.value <- formatC(df$P.value, format = "E", digits = 2)
    if(df$Adjusted.P.value < 0.01)
      df$Adjusted.P.value <- formatC(df$Adjusted.P.value, format = "E", digits = 2)
    df
  }
  
  
  
  output$ENTableSample = DT::renderDataTable(tableRenderingENTableSample(), server = FALSE, extensions = 'Buttons', 
                                             options = list(dom = 'lBfrtip',
                                                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                            pageLength = 5, autoWidth = TRUE))
  
  
  tableRenderingENTableCluster <- function(){
    
    
    clusterID <- as.numeric(input$ENChooseCluster)
    clusterEN <- input$ENChooseClusterGO
    df <- metadata(cdScFiltAnnot)[['enrichCluster']][[1]][[clusterID]][[clusterEN]]  
    if(df$P.value <0.001)
      df$P.value <- formatC(df$P.value, format = "E", digits = 2)
    if(df$Adjusted.P.value < 0.01)
      df$Adjusted.P.value <- formatC(df$Adjusted.P.value, format = "E", digits = 2)
    df
  }
  
  
  
  output$ENTableCluster = DT::renderDataTable(tableRenderingENTableCluster(), server = FALSE, extensions = 'Buttons', 
                                              options = list(dom = 'lBfrtip',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                             pageLength = 5, autoWidth = TRUE))
  
  
  
  #####################################
  
  ######################################
  #
  # Multiple cluster heatmap
  #
  ####################################
  
  HeatmapRendering <- eventReactive(input$buttonForClusterHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    id <- c(input$ChooseClustersHeatmap)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$ClusterGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Clusters %in% id])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Clusters %in% id]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(clustTemp)]
    
    # heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
    #   q10 = quantile(x, 0.1)
    #   q97 = quantile(x, 0.97)
    #   x[x < q10] = q10
    #   x[x > q97] = q97
    #   scale(x)
    # }))
    
    heat.vals <- as.matrix(chosenGene.exprs)
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Clusters = as.factor(clustTemp), Sample = as.factor(cdScFiltAnnot[,names(clustTemp)]$Sample))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(id)], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    ht_list <- Heatmap(heat.vals,
                       col = colorRamp2(c(-1.5,0,1.5), c(input$colminClustHeatmap, input$colmidClustHeatmap, input$colmaxClustHeatmap)),
                       heatmap_legend_param = list(
                         color_bar = "continuous",
                         title = "Scaled expr",
                         direction = "horizontal"
                       ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
                       cluster_columns  = FALSE, show_column_dend = FALSE,
                       row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    
    
  })
  
  
  
  
  output$plotClusterHeatmap <- renderPlot({
    
    HeatmapRendering()
    
  })
  
  
  
  ######################################
  #
  # Multiple sample heatmap
  #
  ####################################
  
  HeatmapRenderingSample <- eventReactive(input$buttonForSampleHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    id <- c(input$ChooseSampleHeatmap)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$SampleGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Sample %in% id])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Sample %in% id]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    sampleTemp <- cdScFiltAnnot[,names(clustTemp)]$Sample
    names(sampleTemp) <- cdScFiltAnnot[,names(clustTemp)]$Barcode
    
    
    #sampleTemp <- sampleTemp[order(sampleTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(sampleTemp)]
    
    # heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
    #   q10 = quantile(x, 0.1)
    #   q97 = quantile(x, 0.97)
    #   x[x < q10] = q10
    #   x[x > q97] = q97
    #   scale(x)
    # }))
    
    heat.vals <- as.matrix(chosenGene.exprs)
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Sample = as.factor(cdScFiltAnnot[,names(sampleTemp)]$Sample), Clusters = as.factor(clustTemp))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(1:nlevels(cdScFiltAnnot$Clusters))], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    ht_list <- Heatmap(heat.vals,
                       col = colorRamp2(c(-1.5,0,1.5), c(input$colminSampleHeatmap, input$colmidSampleHeatmap, input$colmaxSampleHeatmap)),
                       heatmap_legend_param = list(
                         color_bar = "continuous",
                         title = "Scaled expr",
                         direction = "horizontal"
                       ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
                       cluster_columns  = FALSE, show_column_dend = FALSE,
                       row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    
    
  })
  
  
  
  
  output$plotSampleHeatmap <- renderPlot({
    
    HeatmapRenderingSample()
    
  })
  
  
  
  
  ######################################
  #
  # Multiple sample & cluster heatmap
  #
  ####################################
  
  HeatmapRenderingSampleCluster <- eventReactive(input$buttonForSampleClusterHeatmap, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    
    
    # cdScFiltAnnotTmp <- cdScFiltAnnot
    # cdScFiltAnnotTmp$Cluster <- as.numeric(cdScFiltAnnotTmp$Clusters)
    
    #id <- c("1","2","3","4","5","6","7" ,"8","9","10","11")
    #id <- c("1","2","4","10","3","7","11" ,"6","8","9","5")
    idSample <- c(input$ChooseSampleClusterHeatmapSample)
    idCluster <- c(input$ChooseSampleClusterHeatmapCluster)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$SampleClusterGeneList, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    # 
    # #chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,colData(cdScFiltAnnotTmp)$Cluster %in% c(input$ChooseClustersHeatmap)]
    # chosenGene.exprs <- logcounts(cdScFiltAnnotTmp)[geneListTmp,]
    clustTemp <- as.numeric(cdScFiltAnnot$Clusters[cdScFiltAnnot$Sample %in% idSample &
                                                     cdScFiltAnnot$Clusters %in% idCluster])
    names(clustTemp) <- cdScFiltAnnot$Barcode[cdScFiltAnnot$Sample %in% idSample &
                                                cdScFiltAnnot$Clusters %in% idCluster]
    clustTemp <- clustTemp[order(clustTemp, decreasing = FALSE)]
    
    sampleTemp <- cdScFiltAnnot[,names(clustTemp)]$Sample
    names(sampleTemp) <- cdScFiltAnnot[,names(clustTemp)]$Barcode
    
    
    #sampleTemp <- sampleTemp[order(sampleTemp, decreasing = FALSE)]
    
    chosenGene.exprs <- logcounts(cdScFiltAnnot)[geneListTmp,names(sampleTemp)]
    
    heat.vals = t(apply(as.matrix(chosenGene.exprs), 1, function(x) {
      q10 = quantile(x, 0.1)
      q97 = quantile(x, 0.97)
      x[x < q10] = q10
      x[x > q97] = q97
      scale(x)
    }))
    
    #df = data.frame(Cluster = colData(cdScFiltAnnot)[colData(cdScFiltAnnot)$Clusters %in% c(input$ChooseClustersHeatmap),'Clusters'])
    df = data.frame(Sample = as.factor(cdScFiltAnnot[,names(sampleTemp)]$Sample), Clusters = as.factor(clustTemp))
    clusterClass <- c_clust_col
    names(clusterClass) <- c(1:length(c_clust_col))
    
    sampleClass <- c_sample_col[1:nlevels(as.factor(cdScFiltAnnot$Sample))]
    names(sampleClass) <- levels(as.factor(cdScFiltAnnot$Sample))
    
    #ha = HeatmapAnnotation(df = df)
    ha = HeatmapAnnotation(df = df, col = list(Clusters = clusterClass[as.numeric(1:nlevels(cdScFiltAnnot$Clusters))], Sample = sampleClass), 
                           annotation_legend_param = list(
                             Clusters = list(nrow = 1),
                             Sample=list(nrow = 1)))
    # rm(cdScFiltAnnotTmp)
    # gc()
    
    
    ht_list <- Heatmap(heat.vals,
                       col = colorRamp2(c(-1.5,0,1.5), c(input$colminSampleClusterHeatmap, input$colmidSampleClusterHeatmap, input$colmaxSampleClusterHeatmap)),
                       heatmap_legend_param = list(
                         color_bar = "continuous",
                         title = "Scaled expr",
                         direction = "horizontal"
                       ),        top_annotation = ha, show_column_names=FALSE, cluster_rows = FALSE, show_row_dend = FALSE,
                       cluster_columns  = FALSE, show_column_dend = FALSE,
                       row_names_gp = gpar(fontsize = 10))
    draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    
    
  })
  
  
  
  
  output$plotSampleClusterHeatmap<- renderPlot({
    
    HeatmapRenderingSampleCluster()
    
  })
  
  ##################################
  
  
  
  
  ######################################
  #
  # Multiple cluster Bubbleplot
  #
  ####################################
  
  BubblePlotClusterRendering <- eventReactive(input$buttonForClusterBubbleplot, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    id <- c(input$ChooseClustersBubblePlot)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$ClusterGeneListClusterBubblePlot, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    
    group_ids <- id
    
    genes_to_show <- geneListTmp
    
    expression_levels_per_group <- vapply(
      group_ids, FUN.VALUE = numeric(length(genes_to_show)), function(x) {
        cells_in_current_group <- which(cdScFiltAnnot$Clusters %in% x)
        rowMeans(logcounts(cdScFiltAnnot)[genes_to_show,cells_in_current_group])
      }
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(cluster = rownames(.)) %>%
      select(cluster, everything()) %>%
      pivot_longer(
        cols = c(2:ncol(.)),
        names_to = 'gene'
      ) %>%
      dplyr::rename(expression = value) %>%
      mutate(id_to_merge = paste0(cluster, '_', gene))
    
    percentage_of_cells_expressing_gene <- vapply(
      group_ids, FUN.VALUE = numeric(length(genes_to_show)), function(x) {
        cells_in_current_group <- which(cdScFiltAnnot$Clusters %in% x)
        rowSums(logcounts(cdScFiltAnnot)[genes_to_show,cells_in_current_group] != 0)
      }
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(cluster = rownames(.)) %>%
      select(cluster, everything()) %>%
      pivot_longer(
        cols = c(2:ncol(.)),
        names_to = 'gene'
      ) %>%
      dplyr::rename(cell_count = value) %>%
      left_join(
        .,
        as_tibble(colData(cdScFiltAnnot)) %>%
          group_by(Clusters) %>%
          tally() %>%
          dplyr::rename(cluster = Clusters),
        by = 'cluster') %>%
      mutate(
        id_to_merge = paste0(cluster, '_', gene),
        percent_cells = cell_count / n
      )
    
    p <- left_join(
      expression_levels_per_group,
      percentage_of_cells_expressing_gene %>% select(id_to_merge, percent_cells),
      by = 'id_to_merge'
    ) %>%
      mutate(
        cluster = factor(cluster, levels = rev(group_ids)),
        gene = factor(gene, levels = genes_to_show)
      ) %>%
      arrange(gene) %>%
      ggplot(aes(gene, cluster)) +
      geom_point(aes(color = expression, size = percent_cells)) +
      scale_color_distiller(
        palette = 'Reds',
        direction = 1,
        name = 'Log-normalised\nexpression',
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      scale_size(name = 'Percent\nof cells', labels = scales::percent) +
      labs(y = 'Cluster', color = 'Expression') +
      coord_fixed() +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    p
    
  })
  
  output$plotClusterBubbleplot <- renderPlot({
    
    BubblePlotClusterRendering()
    
  })
  
  
  #############################################
  
  
  
  
  ######################################
  #
  # Multiple sample Bubbleplot
  #
  ####################################
  
  BubblePlotSampleRendering <- eventReactive(input$buttonForSampleBubbleplot, {
    
    # validate(
    #   need(input$ChooseClusters != "All Clusters", "Please select a single valid cluster")
    # )
    
    id <- c(input$ChooseSampleBubblePlot)
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[, colData(cdScFiltAnnotTmp)$Cluster %in% id]
    # cdScFiltAnnotTmp <-  cdScFiltAnnotTmp[,order(match(colData(cdScFiltAnnotTmp)$Cluster, id))]
    # 
    # 
    # #cdScFiltAnnotTmp <- cdScFiltAnnotTmp[,order(colData(cdScFiltAnnotTmp)$Cluster)]
    geneListTmp <- strsplit(input$sampleGeneListSampleBubblePlot, split = '\n')
    geneListTmp <- geneListTmp[[1]]
    
    group_ids <- id
    
    genes_to_show <- geneListTmp
    
    expression_levels_per_group <- vapply(
      group_ids, FUN.VALUE = numeric(length(genes_to_show)), function(x) {
        cells_in_current_group <- which(cdScFiltAnnot$Sample %in% x)
        #print(length(rowMeans(logcounts(cdScFiltAnnot)[genes_to_show,cells_in_current_group])))
        rowMeans(logcounts(cdScFiltAnnot)[genes_to_show,cells_in_current_group])
      }
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(sample = rownames(.)) %>%
      select(sample, everything()) %>%
      pivot_longer(
        cols = c(2:ncol(.)),
        names_to = 'gene'
      ) %>%
      dplyr::rename(expression = value) %>%
      mutate(id_to_merge = paste0(sample, '_', gene))
    
    percentage_of_cells_expressing_gene <- vapply(
      group_ids, FUN.VALUE = numeric(length(genes_to_show)), function(x) {
        cells_in_current_group <- which(cdScFiltAnnot$Sample %in% x)
        rowSums(logcounts(cdScFiltAnnot)[genes_to_show,cells_in_current_group] != 0)
      }
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(sample = rownames(.)) %>%
      select(sample, everything()) %>%
      pivot_longer(
        cols = c(2:ncol(.)),
        names_to = 'gene'
      ) %>%
      dplyr::rename(cell_count = value) %>%
      left_join(
        .,
        as_tibble(colData(cdScFiltAnnot)) %>%
          group_by(Sample) %>%
          tally() %>%
          dplyr::rename(sample = Sample),
        by = 'sample') %>%
      mutate(
        id_to_merge = paste0(sample, '_', gene),
        percent_cells = cell_count / n
      )
    
    p <- left_join(
      expression_levels_per_group,
      percentage_of_cells_expressing_gene %>% select(id_to_merge, percent_cells),
      by = 'id_to_merge'
    ) %>%
      mutate(
        sample = factor(sample, levels = rev(group_ids)),
        gene = factor(gene, levels = genes_to_show)
      ) %>%
      arrange(gene) %>%
      ggplot(aes(gene, sample)) +
      geom_point(aes(color = expression, size = percent_cells)) +
      scale_color_distiller(
        palette = 'Reds',
        direction = 1,
        name = 'Log-normalised\nexpression',
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      scale_size(name = 'Percent\nof cells', labels = scales::percent) +
      labs(y = 'Sample', color = 'Expression') +
      coord_fixed() +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    p
    
  })
  
  
  output$plotSampleBubbleplot <- renderPlot({
    
    BubblePlotSampleRendering()
    
  })
  
  
  #############################################
  
  
  
  #####################################
  #
  # DE between two clusters
  #
  ################################
  
  
  output$AllClustTwoClustComp <- renderPlotly({
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df %>% plot_ly(color = ~Clusters, colors = c_clust_col[c(1:nlevels(cdScFiltAnnot$Clusters))], type="scatter", mode="markers", hoverinfo = 'text',
                   text = ~paste('</br> Clusters: ', Clusters,
                                 '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2,
                marker = list(
                  size = 3)) %>% 
      layout(title = 'tSNE with Clusters',showlegend = TRUE, legend = list(title='Cluster',font = list(size = 10), itemsizing='constant'))
    
  })
  
  
  output$plotSelectClust <- renderPlot({
    #df$Group <- "All Group"
    
    df_shiny_ForTwoClust[colData(cdScFiltAnnot)$Clusters %in% input$clust1, "Group"] <- "blue"
    df_shiny_ForTwoClust[colData(cdScFiltAnnot)$Clusters %in% input$clust2, "Group"] <- "red"
    
    
    
    p <- ggplot(df_shiny_ForTwoClust, aes(V1, V2, col = I(Group))) +
      geom_point(size=0.5)  +
      theme_classic(base_size=10) +
      theme(strip.background = element_blank(),
            strip.text.x     = element_blank(),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            panel.border     = element_blank())
    
    p
  })
  
  
  tableRenderingTwoClust <- eventReactive(input$buttonForTwoClustDE, {
    
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust1])
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust2])
    
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    print(table(conds.Sel))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    
    normCounts1 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust1])
    normCounts2 <-  as.matrix(logcounts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Clusters %in% input$clust2])
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
  }
  )
  
  
  output$mytableTwoClust = DT::renderDataTable(tableRenderingTwoClust(), server = FALSE, extensions = 'Buttons',
                                               
                                               options = list(dom = 'lBfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                              pageLength = 5, autoWidth = TRUE)
  )
  
  
  ####################################
  
  
  #####################################
  #
  # DE between two Samples
  #
  ################################
  
  
  output$AllClustTwoSampleComp <- renderPlotly({
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df %>% plot_ly(color = ~Sample, colors = c_sample_col[c(1:nlevels(as.factor(cdScFiltAnnot$Sample)))], type="scatter", mode="markers", hoverinfo = 'text',
                   text = ~paste('</br> Sample: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      add_trace(x=~V1,y=~V2,
                marker = list(
                  size = 3)) %>% 
      layout(title = 'tSNE with Sample',showlegend = TRUE, legend = list(title='Sample',font = list(size = 10), itemsizing='constant'))
    
    
    # Sample <- as.factor(colData(cdScFiltAnnot)[,'Sample'])
    # 
    # ggplotly(ggplot(as.data.frame(reducedDim(cdScFiltAnnot,'tSNE')), aes(x=V1, y=V2, color=Sample)) +
    #            geom_point(size=1) +
    #            guides(colour = guide_legend(override.aes = list(size=1))) +
    #            xlab("") + ylab("") +
    #            ggtitle("t-SNE 2D coloured by Sample") +
    #            scale_colour_manual(values = c_sample_col) + 
    #            theme_bw(base_size=10))
    #            # theme(strip.background = element_blank(),
    #            #       strip.text.x     = element_blank(),
    #            #       axis.text.x      = element_blank(),
    #            #       axis.text.y      = element_blank(),
    #            #       axis.ticks       = element_blank(),
    #            #       axis.line        = element_blank(),
    #            #       panel.border     = element_blank()))
    
  })
  
  
  output$plotSelectSample <- renderPlot({
    #df$Group <- "All Group"
    
    df_shiny_ForTwoSample[colData(cdScFiltAnnot)$Sample %in% c(input$sample1), "Group"] <- "blue"
    df_shiny_ForTwoSample[colData(cdScFiltAnnot)$Sample %in% c(input$sample2), "Group"] <- "red"
    
    
    
    p <- ggplot(df_shiny_ForTwoSample, aes(V1, V2, col = I(Group))) +
      geom_point(size=0.8, alpha = 0.7)  +
      theme_classic(base_size=10)
    # theme(strip.background = element_blank(),
    #       strip.text.x     = element_blank(),
    #       axis.text.x      = element_blank(),
    #       axis.text.y      = element_blank(),
    #       axis.ticks       = element_blank(),
    #       axis.line        = element_blank(),
    #       panel.border     = element_blank())
    
    p
  })
  
  
  tableRenderingTwoSample <- eventReactive(input$buttonForTwoSampleDE, {
    
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample1])
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample2])
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    
    normCounts1 <- exprs(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample1]
    normCounts2 <- exprs(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample2]
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
    
  }
  )
  
  
  output$mytableTwoSample = DT::renderDataTable(tableRenderingTwoSample(), server = FALSE, extensions = 'Buttons',
                                                
                                                options = list(dom = 'lBfrtip',
                                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                               pageLength = 5, autoWidth = TRUE)
  )
  
  
  ####################################
  
  
  
  
  #####################################
  #
  # DE between sample & cluster
  #
  ################################
  
  
  output$AllClustMixedSelection <- renderPlot({
    
    if(input$selectProjectionMixtureSelection == 'tSNE'){
      df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
      projectionType <- 'tSNE'
    }
    else{
      df <- as.data.frame(reducedDim(cdScFiltAnnot,'UMAP'))
      projectionType <- 'UMAP'
    }
    
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    # df <- df %>%
    #   filter(Clusters %in% input$checkboxMultiSelectionCluster) %>%
    #   filter(Sample %in% input$checkboxMultiSelectionSample)
    
    if(input$colorCellsByMixtureSelection == 'Sample'){
      # p <- df %>% plot_ly(color = ~Sample, colors = c_sample_col[c(1:nlevels(as.factor(cdScFiltAnnot$Sample)))], type="scatter", mode="markers", hoverinfo = 'text',
      #                         text = ~paste('</br> Sample: ', Sample))
      p <- ggplot(df, aes(V1, V2, color=Sample)) +
        geom_point(size=1)  +
        theme_classic(base_size=10) + 
        scale_colour_manual(values = c_sample_col)
      projectionCatagory <- 'Sample'
    }
    else{
      # p <- df %>% plot_ly(color = ~Cluster, colors = c_sample_col[c(1:nlevels(as.factor(cdScFiltAnnot$Cluster)))], type="scatter", mode="markers", hoverinfo = 'text',
      #                         text = ~paste('</br> Cluster: ', Cluster))
      p <- ggplot(df, aes(V1, V2, color=Clusters)) +
        geom_point(size=1)  +
        theme_classic(base_size=10) + 
        scale_colour_manual(values = c_clust_col)
      projectionCatagory <- 'Cluster'
    }
    # p %>% layout(legend= list(font=list(size=8))) %>%
    # add_trace(x=~V1,y=~V2,
    #           marker = list(
    #             size = 3)) %>% 
    # layout(title = paste0(projectionType,' with Sample'),showlegend = TRUE, legend = list(title='Sample',font = list(size = 10), itemsizing='constant'))
    # 
    
    p <- p + ggtitle(paste(projectionType,'with',projectionCatagory))
    p
    
    # Sample <- as.factor(colData(cdScFiltAnnot)[,'Sample'])
    # 
    # ggplotly(ggplot(as.data.frame(reducedDim(cdScFiltAnnot,'tSNE')), aes(x=V1, y=V2, color=Sample)) +
    #            geom_point(size=1) +
    #            guides(colour = guide_legend(override.aes = list(size=1))) +
    #            xlab("") + ylab("") +
    #            ggtitle("t-SNE 2D coloured by Sample") +
    #            scale_colour_manual(values = c_sample_col) + 
    #            theme_bw(base_size=10))
    #            # theme(strip.background = element_blank(),
    #            #       strip.text.x     = element_blank(),
    #            #       axis.text.x      = element_blank(),
    #            #       axis.text.y      = element_blank(),
    #            #       axis.ticks       = element_blank(),
    #            #       axis.line        = element_blank(),
    #            #       panel.border     = element_blank()))
    
  })
  
  
  output$selectedCellsClustMixedSelection <- renderPlot({
    
    if(input$selectProjectionMixtureSelection == 'tSNE'){
      df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
      projectionType <- 'tSNE'
    }
    else{
      df <- as.data.frame(reducedDim(cdScFiltAnnot,'UMAP'))
      projectionType <- 'UMAP'
    }
    
    df[,'Sample']=as.factor(colData(cdScFiltAnnot)$Sample)
    df[,'Clusters'] <- as.factor(my.clusters)
    df[,'Group'] <-'gray88'
    #print(head(df))
    # df <- df %>%
    #    filter(Clusters %in% input$checkboxMultiSelectionClusterSel1) %>%
    #    filter(Sample %in% input$checkboxMultiSelectionSampleSel1) %>%
    #    filter(Clusters %in% input$checkboxMultiSelectionClusterSel2) %>%
    #    filter(Sample %in% input$checkboxMultiSelectionSampleSel2)
    
    cells_in_current_group1 <- which(cdScFiltAnnot$Clusters %in% input$checkboxMultiSelectionClusterSel1 &
                                       cdScFiltAnnot$Sample %in% input$checkboxMultiSelectionSampleSel1)
    
    cells_in_current_group2 <- which(cdScFiltAnnot$Clusters %in% input$checkboxMultiSelectionClusterSel2 &
                                       cdScFiltAnnot$Sample %in% input$checkboxMultiSelectionSampleSel2)
    if(!is.null(cells_in_current_group1))
      df[cells_in_current_group1,'Group'] <- 'blue'
    if(!is.null(cells_in_current_group2))
      df[cells_in_current_group2,'Group'] <- 'red'
    
    # p <- df %>% plot_ly(color = ~Sample, colors = c_sample_col[c(1:nlevels(as.factor(cdScFiltAnnot$Sample)))], type="scatter", mode="markers", hoverinfo = 'text',
    #                         text = ~paste('</br> Sample: ', Sample))
    p <- ggplot(df, aes(V1, V2, col=I(Group))) +
      geom_point(size=1, alpha=0.5)  +
      theme_classic(base_size=10) 
    
    
    
    p <- p + ggtitle(paste(projectionType))
    p
    
    
    
  })
  
  
  tableRenderingMixedSampleCluster <- eventReactive(input$buttonForDEMixedSelection, {
    
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,cdScFiltAnnot$Clusters %in% input$checkboxMultiSelectionClusterSel1 &
                                                 cdScFiltAnnot$Sample %in% input$checkboxMultiSelectionSampleSel1])
    
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,cdScFiltAnnot$Clusters %in% input$checkboxMultiSelectionClusterSel2 &
                                                 cdScFiltAnnot$Sample %in% input$checkboxMultiSelectionSampleSel2])
    
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    print(table(conds.Sel))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    
    normCounts1 <- exprs(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample1]
    normCounts2 <- exprs(cdScFiltAnnot)[,colData(cdScFiltAnnot)$Sample %in% input$sample2]
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
    
  }
  )
  
  
  output$mytableMixedSelection = DT::renderDataTable(tableRenderingMixedSampleCluster(), server = FALSE, extensions = 'Buttons',
                                                     
                                                     options = list(dom = 'lBfrtip',
                                                                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                                    pageLength = 5, autoWidth = TRUE)
  )
  
  
  ####################################
  
  
  
  ########################################
  #
  # Manual selection for DE
  #
  ########################################
  
  
  
  plot_tsnePlotSampleManualSelection <- function(){ 	
    
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df$Clusters <- as.factor(cdScFiltAnnot$Clusters)
    df$Sample <- as.factor(cdScFiltAnnot$Sample)
    df$Barcode <- cdScFiltAnnot$Barcode
    df$key <- cdScFiltAnnot$Barcode
    df$Group <- rev(c_sample_col[as.integer(as.factor(cdScFiltAnnot$Sample))]) 
    rownames(df) <- cdScFiltAnnot$Barcode
    df <- df %>%
      filter(Clusters %in% input$checkboxMultiSelectionCluster) %>%
      filter(Sample %in% input$checkboxMultiSelectionSample)
    
    click_data <- event_data("plotly_click")
    select_data <- event_data("plotly_selected")
    
    # df %>%
    #   filter(Clusters %in% input$checkboxMultiSelectionCluster) %>%
    #   filter(Sample %in% input$checkboxMultiSelectionSample) %>%
    #   plot_ly(color = ~Sample, colors = c_sample_col[c(1:4)], type="scatter", mode="markers", hoverinfo = 'text',
    #           marker = list(size = 1.5, opacity = 1),
    #           text = ~paste('</br> Cell: ', Cell,
    #                         '</br> Clusters: ', Clusters,
    #                         '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
    #   add_trace(x=~V1,y=~V2) %>%
    #   layout(title = 'tSNE with Samples',showlegend = TRUE, legend = list(font = list(size = 10), itemsizing='constant'))
    
    click_data <- event_data("plotly_click", priority = 'event')
    select_data <- event_data("plotly_selected")
    
    #print(head(click_data))
    #print(head(select_data))
    
    
    if(is.null(select_data)){
      # df %>%
      # plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      
      p <- ggplot(df, aes(V1, V2, col = Sample, key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        scale_colour_manual(values = c_sample_col) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    else if (!is.null(select_data) && flag == 0) {
      saved_key <<- select_data$key
      df[df$key %in% select_data$key, "Group"] <- "gray5"
      #print(select_data)
      saveClust1 <<- select_data$key  
      flag <<- 1
      # df %>%
      # plot_ly(color = ~I(Group), colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      
      p <- ggplot(df, aes(V1, V2, col = I(Group), key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    else if (!is.null(select_data) && flag == 1) {
      df$Group <- 'gray88'
      df[df$key %in% saved_key, "Group"] <- "blue"
      #df_shiny[df_shiny$key %in% select_data$key, "Group"] <- "red"
      df[df$key %in% select_data$key, "Group"] <- "red"
      saveClust2 <<- select_data$key  
      flag <<- 0
      # df %>%
      # plot_ly(color = ~I(Group), colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      select_data <- NULL
      click_data <- NULL
      p <- ggplot(df, aes(V1, V2, col = I(Group), key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    
  }
  
  
  
  plot_tsnePlotClusterManualSelection <- function(){ 	
    
    df <- as.data.frame(reducedDim(cdScFiltAnnot,'tSNE'))
    df$Clusters <- as.factor(cdScFiltAnnot$Clusters)
    df$Sample <- as.factor(cdScFiltAnnot$Sample)
    df$Barcode <- cdScFiltAnnot$Barcode
    df$key <- cdScFiltAnnot$Barcode
    df$Group <- c_clust_col[cdScFiltAnnot$Clusters] 
    rownames(df) <- cdScFiltAnnot$Barcode
    df <- df %>%
      filter(Clusters %in% input$checkboxMultiSelectionCluster) %>%
      filter(Sample %in% input$checkboxMultiSelectionSample)
    
    click_data <- event_data("plotly_click")
    select_data <- event_data("plotly_selected")
    
    if(is.null(select_data)){
      # df %>%
      # plot_ly(color = ~Clusters, colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      #print(input$checkboxMultiSelectionCluster)
      p <- ggplot(df, aes(V1, V2, col = Clusters, key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        scale_colour_manual(values = c_clust_col[as.numeric(input$checkboxMultiSelectionCluster)]) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    else if (!is.null(select_data) && flag == 0) {
      saved_key <<- select_data$key
      df[df$key %in% select_data$key, "Group"] <- "gray5"
      #print(select_data)
      saveClust1 <<- select_data$key  
      flag <<- 1
      # df %>%
      # plot_ly(color = ~I(Group), colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      
      p <- ggplot(df, aes(V1, V2, col = I(Group), key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    
    else if (!is.null(select_data) && flag == 1) {
      df$Group <- 'gray88'
      df[df$key %in% saved_key, "Group"] <- "blue"
      #df_shiny[df_shiny$key %in% select_data$key, "Group"] <- "red"
      df[df$key %in% select_data$key, "Group"] <- "red"
      saveClust2 <<- select_data$key  
      flag <<- 0
      # df %>%
      # plot_ly(color = ~I(Group), colors = c_clust_col[c(1:9)], type="scatter", mode="markers", hoverinfo = 'text',
      #         marker = list(size = 1, opacity = 1),
      #         text = ~paste('</br> Clusters: ', Clusters,
      #                       '</br> Samples: ', Sample)) %>% layout(legend= list(font=list(size=8))) %>%
      #   add_trace(x=~V1,y=~V2) %>%
      #   layout(title = 'tSNE with Clusters',showlegend = TRUE, 
      #          legend = list(font = list(size = 10), itemsizing='constant'),
      #          dragmode = "lasso")
      p <- ggplot(df, aes(V1, V2, col = I(Group), key = Barcode)) + 
        geom_point(size=0.1)  +
        theme_classic(base_size=10) +
        theme(strip.background = element_blank(),
              strip.text.x     = element_blank(),
              axis.text.x      = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks       = element_blank(),
              axis.line        = element_blank(),
              panel.border     = element_blank(),
              legend.position = "none")
      
      ggplotly(p) %>% layout(dragmode = "lasso")
    }
    
    
  }
  
  
  output$AllClustManualSelection <- renderPlotly({
    if(input$colorCellsByManualSelection == 'Sample'){
      plot_tsnePlotSampleManualSelection()
    }
    else{
      plot_tsnePlotClusterManualSelection()
    }
    
  })
  
  
  tableRendering <- eventReactive(input$buttonForDEManualSelection, {
    
    
    # clusterMod <- as.factor(my.clusters)
    # chosen.clust <- which(levels(clusterMod)==input$ChooseClusters)
    # for (clust in 1:length(clusterMod)) {
    #   if (clusterMod[clust]==chosen.clust[1]) { 
    #     clusterMod[clust] <- "1"
    #   }
    #   else {
    #     clusterMod[clust] <- "2"
    #   }
    #   
    # }
    # 
    counts1 <- as.matrix(counts(cdScFiltAnnot)[,saveClust1])
    counts2 <- as.matrix(counts(cdScFiltAnnot)[,saveClust2])
    countsTable <- cbind(counts1,counts2)
    rownames(countsTable) <- make.names(rownames(countsTable), unique = TRUE)
    
    conds.Sel <- c(rep("1", dim(counts1)[2]), rep("2", dim(counts2)[2]))
    print(table(conds.Sel))
    
    res1 <- nbTestSH(countsTable, conds.Sel, "1", "2")
    res1[,'FDR'] <- p.adjust(res1[,7], method = "bonferroni")
    
    res1$Mean <- round(res1$Mean,3)
    res1$rawMeanA <- round(res1$rawMeanA,3)
    res1$rawMeanB <- round(res1$rawMeanB,3)
    res1$rawLog2FoldChange <- round(res1$rawLog2FoldChange,3)
    
    res1$pval <- formatC(res1$pval, format = "E", digits = 2)
    res1$FDR <- formatC(res1$FDR, format = "E", digits = 2)
    
    #res1[,'Z-score_1'] <- colMeans(scale(t(normCounts1)))
    #res1[,'Z-score_2'] <- colMeans(scale(t(normCounts2)))
    res1[,c('Mean','rawMeanA','rawMeanB','rawLog2FoldChange','pval','FDR')]
  }
  )
  
  
  output$mytableManualSelection = DT::renderDataTable(tableRendering(), server = FALSE, extensions = 'Buttons', 
                                                      options = list(dom = 'lBfrtip',
                                                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                                                                     pageLength = 5, autoWidth = TRUE)
  )
  
  observeEvent(input[["reset"]], {
    runjs("Shiny.setInputValue('plotly_selected-A', null);")
  })
  ##########################################
  
  
  
  
  ####################################
  #
  # All info messages
  #
  ####################################
  
  observeEvent(input$projectionInfo, {
    
    showModal(modalDialog(
      title = "Projection",
      tags$div(
        tags$ul(
          tags$li("Projection of cells into 2D space. The projection can be toggled between tSNE and UMAP from the Projection dropdown menu"),
          tags$li("Cells can be coloured based on samples, clusters, number of UMIs each cell have and number of genes each cell express"),
          tags$li("Initial size of the dots are 1.5 which can be changed by moving the slide bar"),
          tags$li("Dot opacity controls the visual transperancy of the cells. It can be lowered by moving the slide bar so that cells 
                  underneath one eather can be can be viewed"),
          tags$li("Samples or clusters can be removed by clicking on each of the lables in the plot. One click will remove the cells associated with 
                  that feature and double click will remove all but the one that is clicked"),
          tags$li("Hovering mouse over each cell would show the cell barcode, cluster, sample and CellType"),
          tags$li("If you want to see a specific cell type, just double click on the legend of the image, so other cell
                  types would be deactivated from visualization."),
          tags$li("You can select/deselect samples or clusters from the dropdown list bottom of the figure. This would show only the selected cells.")
        )
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$sampleByclustersInfo, {
    
    showModal(modalDialog(
      title = "Samples into clusters",
      tags$div(
        "This table shows the number of samples distributed into different clusters. It also shows the total number of cells in a sample."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  
  observeEvent(input$percentClusterInSampleInfo, {
    
    showModal(modalDialog(
      title = "Percent of cells from cluster going to sample",
      tags$div(
        "This figures shows the % of cells going into a cluster. All the cells in a specific sample are divided into different clusters and the percentage
        is calcualted based on cells only this sample."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$umiinSampleInfo, {
    
    showModal(modalDialog(
      title = "UMI in sample",
      tags$div(
        "This figures shows the total number of Unique Molecular Identifier (UMI)s for each individul sample in log10 scale. Unique molecular identifiers (UMI) are molecular tags that are used to detect and 
         quantify unique mRNA transcripts. In this method, mRNA libraries are generated by fragmentation and reverse-transcribed to cDNA. 
         Oligo(dT) primers with specific sequencing linkers are added to the cDNA."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$numberOfGenesExpressedInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This violin plot shows the total number of genes expressed in each sample, i.e. genes having a UMI count > 0"
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$umiInClustersInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This figures shows the total number of Unique Molecular Identifier (UMI)s for each individul cluster in log10 scale. Unique molecular identifiers (UMI) are molecular tags that are used to detect and 
         quantify unique mRNA transcripts. In this method, mRNA libraries are generated by fragmentation and reverse-transcribed to cDNA. 
        Oligo(dT) primers with specific sequencing linkers are added to the cDNA. "
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$numberOfGenesExpressedInClustersInfo, {
    
    showModal(modalDialog(
      title = "Total genes expressed in sample",
      tags$div(
        "This violin plot shows the total number of genes expressed in each cluster, i.e. genes having a UMI count > 0"
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$DE_between_sample_and_clustersSelectionInfo, {
    
    showModal(modalDialog(
      title = "Differential Expression between multiple sample and cluster",
      tags$div(
        "This input panel allows you to choose samples and clusters within a sample. DE will be calculated on cells that
        belong to the selected samples and clusters. This panel helps to select the cells that belong to a specific sample
        and cluster. For eg. cluster 1 is shared between Sample 1 and Sample 2. If you select cluster 1 and sample 1 then
        cells only in cluster 1 and sample 1 would be selected for the DE.",
        tags$ul(
          tags$li("The first selection is for the first group of cells."),
          tags$li("Second selection is for the second group of cells.")
        )
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  
  observeEvent(input$DE_between_sample_and_clustersProjectionInfo, {
    showModal(modalDialog(
      title = "Projection of cells",
      tags$div(
        "This panel projects the cells into a 2D dimension. Users can choose the projection of tSNE or UMAP. 
        They can also choose to colour the cells based on sample or cluster. This panel helps the users to identify
        the cells they want to choose for differential expression."
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  observeEvent(input$DE_in_manual_selection_input, {
    showModal(modalDialog(
      title = "Manual selection of cells for DE",
      tags$div(
        tags$ul(
          tags$li("This input panel configures the projection panel based on user selection of the inputs."),
          tags$li("In order to select cells manually on which the users want to run DE, they first need to select the first group of cells
                  by dragging with mouse pointer over the cells. The second group of cells are then selected by dragging again on a 
                  second group of cells. The first group of cells would then be coloured blue and the second group would be coloured blue."),
          tags$li("The user would then click DE between selected cells to run the Differential Expression analysis between the 
                  two groups of cells"),
          tags$li("In order to clear all the selection the uer needs to click Reset all selection")
        )
        
      ),
      
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  ############################################
  #
  # For Trajectory panel
  #
  ############################################
  
  
  ## Slingshot
  output$trajectory_slingshotOT <- renderPlot({
    
    # Run PCA on Deng data. Use the runPCA function from the SingleCellExperiment package.
    cdScFiltAnnot <- runPCA(cdScFiltAnnot, ncomponents = 50)
    
    # Use the reducedDim function to access the PCA and store the results. 
    pca <- reducedDim(cdScFiltAnnot, "PCA")
    
    # Describe how the PCA is stored in a matrix. Why does it have this structure?
    head(pca)
    dim(pca)
    
    
    # Add PCA data to the deng_SCE object.
    cdScFiltAnnot$PC1 <- pca[, 1]
    cdScFiltAnnot$PC2 <- pca[, 2]
    
    
    # Plot PC biplot with cells colored by cell_type2. 
    # colData(deng_SCE) accesses the cell metadata DataFrame object for deng_SCE.
    # Look at Figure 1A of the paper as a comparison to your PC biplot.
    ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = PC1, y = PC2, color = cellType)) + geom_quasirandom(groupOnX = FALSE) +
      scale_color_tableau() + theme_classic() +
      xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")
    
    
    
    
    
    
    
    
    
    
    
    
    # Read the Slingshot documentation (?slingshot) and then run Slingshot below. 
    # Given your understanding of the algorithm and the documentation, what is one 
    # major set of parameters we omitted here when running Slingshot?
    sce <- slingshot(cdScFiltAnnot, reducedDim = 'PCA')  # no clusters
    
    # Plot PC1 vs PC2 colored by Slingshot pseudotime.
    colors <- rainbow(50, alpha = 1)
    plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2)
    
    # Plot Slingshot pseudotime vs cell stage. 
    ggplot(as.data.frame(colData(cdScFiltAnnot)), aes(x = sce$slingPseudotime_1, y = cellType, 
                                                      colour = cellType)) +
      geom_quasirandom(groupOnX = FALSE) +
      scale_color_tableau() + theme_classic() +
      xlab("Slingshot pseudotime") + ylab("Timepoint") +
      ggtitle("Cells ordered by Slingshot pseudotime")
    plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2)
    
    
   })                                                                                                                                                                                           
  ##----------------------------------------------------------------------------##
  ## Projection.
  ##----------------------------------------------------------------------------##
  output$trajectory_projection <- plotly::renderPlotly({
    # don't do anything before these inputs are selected
    req(
      input[["trajectory_to_display"]],
      input[["trajectory_samples_to_display"]],
      input[["trajectory_clusters_to_display"]],
      input[["trajectory_percentage_cells_to_show"]],
      input[["trajectory_dot_color"]],
      input[["trajectory_dot_size"]],
      input[["trajectory_dot_opacity"]]
    )
    
    trajectory_to_display <- input[["trajectory_to_display"]]
    samples_to_display <- input[["trajectory_samples_to_display"]]
    clusters_to_display <- input[["trajectory_clusters_to_display"]]
    cells_to_display <- which(
      (cdScFiltAnnot$cellType$sample %in% samples_to_display) &
        (cdScFiltAnnot$cellType$cluster %in% clusters_to_display)
    )
    
    # randomly remove cells
    if ( input[["trajectory_percentage_cells_to_show"]] < 100 ) {
      number_of_cells_to_plot <- ceiling(
        input[["trajectory_percentage_cells_to_show"]] / 100 * length(cells_to_display)
      )
      cells_to_display <- cells_to_display[ sample(1:length(cells_to_display), number_of_cells_to_plot) ]
    }
    
    # extract cells to plot
    to_plot <- cbind(
      cdScFiltAnnot$trajectory$monocle2[[ trajectory_to_display ]][["meta"]][ cells_to_display , ],
      cdScFiltAnnot$cellType[ cells_to_display , ]
    ) %>%
      dplyr::filter(!is.na(pseudotime))
    to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]
    
    color_variable <- input[["trajectory_dot_color"]]
  
    # convert edges of trajectory into list format to plot with plotly
    trajectory_edges <- cdScFiltAnnot$trajectory$monocle2[[trajectory_to_display]][["edges"]]
    trajectory_lines <- list()
    for (i in 1:nrow(trajectory_edges) ) {
      line = list(
        type = "line",
        line = list(color = "black"),
        xref = "x",
        yref = "y",
        x0 = trajectory_edges$source_dim_1[i],
        y0 = trajectory_edges$source_dim_2[i],
        x1 = trajectory_edges$target_dim_1[i],
        y1 = trajectory_edges$target_dim_2[i]
      )
      trajectory_lines <- c(trajectory_lines, list(line))
    }
    
    if ( is.factor(to_plot[[ color_variable ]]) || is.character(to_plot[[ color_variable ]]) ) {
      if ( color_variable == "sample" ) {
        colors_this_plot <- reactive_colors()$sampless
      } else if ( color_variable == "cluster" ) {
        colors_this_plot <- reactive_colors()$clusters
      } else if ( color_variable %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
        colors_this_plot <- cell_cycle_colorset
      } else if ( is.factor(to_plot[[ color_variable ]]) ) {
        colors_this_plot <- setNames(
          default_colorset[1:length(levels(to_plot[[ color_variable ]]))],
          levels(to_plot[[ color_variable ]])
        )
      } else {
        colors_this_plot <- default_colorset
      }
      plot <- plotly::plot_ly(
        to_plot,
        x = ~DR_1,
        y = ~DR_2,
        color = ~to_plot[[ color_variable ]],
        colors = colors_this_plot,
        type = "scatter",
        mode = "markers",
        marker = list(
          opacity = input[["trajectory_dot_opacity"]],
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["trajectory_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>State</b>: ", to_plot[ , "state" ], "<br>",
          "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
        )
      ) %>%
        plotly::layout(
          shapes = trajectory_lines,
          xaxis = list(
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE,
            range = range(to_plot$DR_1) * 1.1
          ),
          yaxis = list(
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE,
            range = range(to_plot$DR_2) * 1.1
          ),
          hoverlabel = list(font = list(size = 11))
        )
      if ( preferences[["use_webgl"]] == TRUE ) {
        plot %>% plotly::toWebGL()
      } else {
        plot
      }
    } else {
      plot <- plotly::plot_ly(
        data = to_plot,
        x = ~DR_1,
        y = ~DR_2,
        type = "scatter",
        mode = "markers",
        marker = list(
          colorbar = list(
            title = colnames(to_plot)[which(colnames(to_plot) == color_variable)]
          ),
          color = ~to_plot[[ color_variable ]],
          opacity = input[["trajectory_dot_opacity"]],
          colorscale = "YlGnBu",
          reversescale = TRUE,
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["trajectory_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>State</b>: ", to_plot[ , "state" ], "<br>",
          "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
        )
      ) %>%
        plotly::layout(
          shapes = trajectory_lines,
          xaxis = list(
            title = colnames(to_plot)[1],
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE,
            range = range(to_plot$DR_1) * 1.1
          ),
          yaxis = list(
            title = colnames(to_plot)[2],
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE,
            range = range(to_plot$DR_2) * 1.1
          ),
          hoverlabel = list(font = list(size = 11))
        )
      if ( preferences$use_webgl == TRUE ) {
        plotly::toWebGL(plot)
      } else {
        plot
      }
    }
  })
  ##----------------------------------------------------------------------------##
  ## Distribution along pseudotime.
  ##----------------------------------------------------------------------------##
  
  output[["trajectory_density_plot"]] <- plotly::renderPlotly({
    # don't do anything before these inputs are selected
    req(
      input[["trajectory_to_display"]],
      input[["trajectory_samples_to_display"]],
      input[["trajectory_clusters_to_display"]],
      input[["trajectory_dot_color"]]
    )
    
    trajectory_to_display <- input[["trajectory_to_display"]]
    samples_to_display <- input[["trajectory_samples_to_display"]]
    clusters_to_display <- input[["trajectory_clusters_to_display"]]
    cells_to_display <- which(
      (cdScFiltAnnot$cells$sample %in% samples_to_display) &
        (cdScFiltAnnot$cells$cluster %in% clusters_to_display)
    )
    
    # extract cells to plot
    to_plot <- cbind(
      input$trajectory_projection_info[[ trajectory_to_display ]][["meta"]][ cells_to_display , ],
      cdScFiltAnnot$cells[ cells_to_display , ]
    ) %>%
      dplyr::filter(!is.na(pseudotime))
    to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]
    
    color_variable <- input[["trajectory_dot_color"]]
    
    if ( is.factor(to_plot[[ color_variable ]]) || is.character(to_plot[[ color_variable ]]) ) {
      if ( color_variable == "sample" ) {
        colors_this_plot <- reactive_colors()$samples
      } else if ( color_variable == "cluster" ) {
        colors_this_plot <- reactive_colors()$clusters
      } else if ( color_variable %in% c("cell_cycle_seurat","cell_cycle_cyclone") ) {
        colors_this_plot <- cell_cycle_colorset
      } else if ( is.factor(to_plot[[ color_variable ]]) ) {
        colors_this_plot <- setNames(
          default_colorset[1:length(levels(to_plot[[ color_variable ]]))],
          levels(to_plot[[ color_variable ]])
        )
      } else {
        colors_this_plot <- default_colorset
      }
      p <- ggplot(to_plot, aes_string(x = "pseudotime", fill = color_variable)) +
        geom_density(alpha = 0.4, color = "black") +
        theme_bw() +
        labs(x = "Pseudotime", y = "Density") +
        scale_fill_manual(values = colors_this_plot) +
        guides(fill = guide_legend(override.aes = list(alpha = 1)))
      plotly::ggplotly(p, tooltip = "text") %>%
        plotly::style(
          hoveron = "fill"
        )
    } else {
      colors_this_plot <- setNames(
        default_colorset[1:length(levels(cdScFiltAnnot$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state))],
        levels(cdScFiltAnnot$trajectory$monocle2[[ input[["trajectory_to_display"]] ]][["meta"]]$state)
      )
      plot <- plotly::plot_ly(
        data = to_plot,
        x = ~pseudotime,
        y = ~to_plot[[ color_variable ]],
        type = "scatter",
        mode = "markers",
        color = ~state,
        colors = colors_this_plot,
        marker = list(
          opacity = input[["trajectory_dot_opacity"]],
          line = list(
            color = "rgb(196,196,196)",
            width = 1
          ),
          size = input[["trajectory_dot_size"]]
        ),
        hoverinfo = "text",
        text = ~paste(
          "<b>Cell</b>: ", to_plot[ , "cell_barcode" ], "<br>",
          "<b>Sample</b>: ", to_plot[ , "sample" ], "<br>",
          "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br>",
          "<b>Transcripts</b>: ", formatC(to_plot[ , "nUMI" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>Expressed genes</b>: ", formatC(to_plot[ , "nGene" ], format = "f", big.mark = ",", digits = 0), "<br>",
          "<b>State</b>: ", to_plot[ , "state" ], "<br>",
          "<b>Pseudotime</b>: ", round(to_plot[ , "pseudotime" ], 3)
        )
      ) %>%
        plotly::layout(
          xaxis = list(
            title = "Pseudotime",
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE
          ),
          yaxis = list(
            title = color_variable,
            mirror = TRUE,
            showline = TRUE,
            zeroline = FALSE
          ),
          hoverlabel = list(font = list(size = 11))
        )
      if ( preferences$use_webgl == TRUE ) {
        plotly::toWebGL(plot)
      } else {
        plot
      }
    }
  })
  
  
  ###################################
  
  
  
  
  
  
}

shinyApp(ui, server)