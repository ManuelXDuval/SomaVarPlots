#GeneCreek
#Nov. 24th, 2017
#application scope: visualize somatic variants onto an Errb2 gene diagram.
#SomaVarPlot ui
library(shiny)
#library(shinyjs)

# Define UI ----
ui <- fluidPage(
  titlePanel("Oncogene Gene1 somatic variants visualization"),
  sidebarLayout(
    sidebarPanel(
      h2("tumor origin"),
      selectInput("tumorSource", "Tumor Source", c("biliary", "bladder", "cervical")),
      h2("plot"),
      selectInput("plotType", "plot type", c("", "Lolliplot", "Oncoplot", "Waterfall and swim plots"))
    ),
    mainPanel(
      #h1("Repartition of the somatic variants w.r.t. the Best Overall Response"),
      h1("Repartition of the somatic variants w.r.t. the clinical end point"),
      h2(textOutput("result")),
      h3(textOutput("plotselec")),
      #plotOutput("lolliplot", height = "500px"),
      plotOutput("plot", height = "500px")
      
    )
  )
)