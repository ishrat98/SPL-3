library(shiny)
library(shinydashboard)


library(scater)
library(plotly)
library(reshape2)
library(circlize)
library(sSeq)
library(shinycssloaders)
library(ComplexHeatmap)
library(colourpicker)
library(shinyWidgets)
library(shinyjs)

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

ui <- dashboardPage(
  #skin = "purple",
  dashboardHeader(
    title = "UoM Single cell"
  )
  shinyApp(ui, server)
  