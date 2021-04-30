library(shinythemes) # adding the shinythemese package
library(shiny)
library(DESeq2)
library(dplyr)
library(cluster)
library(dendextend)
library(colorRamps)
library(RColorBrewer)
library(gplots)
library(vegan)
library(apcluster)
library(dynamicTreeCut)
library(e1071)
library(tidyr)
library("fpc")
library(reshape2)
library(ggplot2)
library(DT)

ui <- fluidPage(theme = shinytheme("darkly"), 
         titlePanel(strong("Clustering of RNAseq data")), # using strong as a direct tag
                  mainPanel(
                    tabsetPanel(
                      tabPanel("introduction",  
                               fluidRow( 
                               
                               column(12,
                                      
                               htmlOutput("intro1")
                               ),
                               column(12,
                                      wellPanel(
                                        htmlOutput("intro2"),
                                        tableOutput("counts1") 
                               )),
                               column(12,
                                      htmlOutput("tab1")
                               ),
                               column(9,
                                      wellPanel(
                                        tableOutput("counts2")
                               )),
                               column(12,
                                      
                                      htmlOutput("tab2")
                               ),
                               column(6,
                                      wellPanel(
                                        tableOutput("sample") 
                               )),
                               column(12,
                                      htmlOutput("intro3")
                                      
                              )
                               
                           )
                        ), 
                      tabPanel("Input",
                        fluidRow(
                          column(4,
                            wellPanel(
                                 br(),
                                 fileInput("file1", "Counts File",
                                           multiple = TRUE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 radioButtons("filetype", label ="file type", 
                                             choices = list("featureCounts" = 0 ,"preprocessed" = 1), 
                                             selected = 0),
                                 radioButtons("disp", "Display",
                                              choices = c(Head = "head",
                                                          All = "all"),
                                              selected = "head"),
                                 actionButton("goButton", "submit", 
                                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),

                                 br(),
                                 br(),
                                 tags$hr(),
                                 #br(),
                                 fileInput("file2", "Choose sample setup file",
                                           multiple = TRUE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 #tags$br(),
                                 radioButtons("disp2", "Display",
                                              choices = c(Head = "head",
                                                          All = "all"),
                                              selected = "head"),
                                 actionButton("goButton2", " submit",
                                              style="color: #fff; background-color: #1ddc5b; border-color: #19cf54"),
                                 br(),
                                 br()
                                 )
                          ),
                          column(8, 
                                 br(),
                                 h4(htmlOutput("text1")), 
                                 #h4(textOutput(p("text1", style = "color:#19cf54"))),
                                 br(),
                                 br(),
                                 tableOutput("contents"),
                                 tags$hr(),
                                 br(),
                                 br(),
                                 h4(htmlOutput("text2")),
                                 tableOutput("contents2")
                                 )
                        )),
                      tabPanel("Normalization",
                         fluidRow(
                            column(4,
                              wellPanel(
                                  sliderInput("f", "Foldchange", min = 2, max = 20, value = 2, step = 1),
                                  br()
                                )),
                            
                            column(4,
                                   wellPanel(
                                     sliderInput("p", "P-value", min = 0.001, max = 0.5, value = 0.01, step = 0.001),
                                     br()
                                   )),
                            
                            column(2,
                                   br(),
                                   br(),
                                   radioButtons("disp3", "Display",
                                                choices = c(Head = "head",
                                                            All = "all"),
                                                selected = "head"),
                                   br()
                            ),
                            column(2,
                                   br(),
                                   actionButton("goNorm", "submit",
                                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4; height:100px"),
                                   br()
                            ),

                            column(12,
                                  br(),
                                  htmlOutput("text3"), 
                                  plotOutput(outputId = "heatmap1"),
                                  br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                  tableOutput("contents3"),
                                  br()
                            )
                          )), 
                      tabPanel("Cluster number",
                               #h4(htmlOutput("text4")), 
                               selectInput("opt_K", label =" ", 
                                           choices = list("Choose_method" = 0 ,"elbow" = 1, "average_silhouette_width" = 2, "Calinsky_criterion" = 3, "gap_statistic" = 4, "affinity_propogation" = 5 ), 
                                           selected = 0),
                               
                               br(),
                               actionButton("goButton3", "submit", 
                                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                               br(),
                               br(),
                               h4(htmlOutput("text9")),
                               br(),
                               br(),
                               plotOutput(outputId = "opt_K_plot")
                               
 
                               ),
                      tabPanel("Clustering",
                               fluidRow(
                                 column(9,         
                               wellPanel(
                        tabsetPanel(id = "ClusMet",
                          tabPanel("hierarchical",
                                   br(),
                                   selectInput("selectH", label =h4("Number of clusters"), 
                                               choices = list("2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"= 8,"9"= 9,"10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,"21"=21,"22"=22,"23"=23,"24"=24,"25"=25), 
                                               selected = 2)
                            
                                  ),
                          tabPanel("Dynamic Hr",
                                   br(),
                                   h4(htmlOutput("text10"))
                                   

                                   
                          ),
                          tabPanel("K-mean",
                                   br(),
                                   selectInput("selectK", label =h4("Number of clusters"), 
                                               choices = list("2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"= 8,"9"= 9,"10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,"21"=21,"22"=22,"23"=23,"24"=24,"25"=25), 
                                               selected = 2)
                                   
                                  ),
                          tabPanel("Fuzzy K-mean",
                                   br(),
                                   selectInput("selectF", label =h4("Number of clusters"), 
                                               choices = list("2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"= 8,"9"= 9,"10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,"21"=21,"22"=22,"23"=23,"24"=24,"25"=25), 
                                               selected = 2)
                            
                                  ),
                          tabPanel("PAM",
                                   br(),
                                   selectInput("selectP", label =h4("Number of clusters"), 
                                               choices = list("2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=7,"8"= 8,"9"= 9,"10"=10,"11"=11,"12"=12,"13"=13,"14"=14,"15"=15,"16"=16,"17"=17,"18"=18,"19"=19,"20"=20,"21"=21,"22"=22,"23"=23,"24"=24,"25"=25), 
                                               selected = 2)
                                  ), 
                          tabPanel("dpscan",
                                   br(),
                                   sliderInput("eps", "eps", min = 0.001, max = 1.0, value = 0.1, step = 0.001),
                                   sliderInput("minPts", "minPts", min = 1, max = 25, value = 5, step = 1)
                          )
                        )
                      )
                                  ),
                        column(2,
                               br(),
                               br(),
                               hr(),
                               radioButtons("Avarage", "Average to conditions",
                                            choices = c(Yes = "yes",
                                                        No = "no"),
                                            selected = "no"),
                               br()
                        ),
                        column(1,
                               br(),
                               br(),
                               hr(),
                               actionButton("goClus", "submit",
                                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4; height:100px"),
                               br()
                        ),
                        column(12,  
                               br(),
                               h4(htmlOutput("text5")),
                               plotOutput(outputId = "heatmap2"),
                               br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                               plotOutput(outputId = "cluCen")
                           )
                         )     
                       ),
                      tabPanel("individual clusters",                          
                             fluidRow(
                              column(4,
                               br(),
                               selectInput('cluster', 'cluster', ""),
                               actionButton("sinClus", "submit",
                                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                              ),

                              column(12,
                               h4(htmlOutput("text11")),
                               plotOutput(outputId = "heatmap3"),
                               br(),
                               plotOutput(outputId = "SinClus"),
                               br(),
                               h4(htmlOutput("text12")),
                               dataTableOutput("clusTab"),
                               br(),
                               br()
                              ))
                      ),
                      tabPanel("individual genes",
                          fluidRow(
                              column(4, 
                                     br(),
                                     textAreaInput('gene', 'Input gene name(s) separated by comma:', value = "", placeholder = 'e.g. ENSMUSG00000000028,ENSMUSG00000000001')
                              ),
                              column(1,
                                     br(),
                                     br(),
                                     actionButton("do", "submit",style="color: #fff; background-color: #337ab7; border-color: #2e6da4; height:70px")
                              ),
                              column(5,
                                     br(),
                                     htmlOutput("text16")
                                     
                              ),       
                              column(2,      
                                     br(),
                                     downloadButton('download',HTML("Download"),
                                                    style="color: #0a0800; background-color: #f7c00a; border-color: #fadb73; padding-top: 30px; height:90px; width:130px")  
                              ),
                              column(12,
                                      h4(htmlOutput("text14")),
                                      plotOutput("plot1"),
                                      br(),    
                                      br(),
                                      h4(htmlOutput("text15")),
                                      plotOutput("plot2")
                             
                                 
                          )
                        )
                      )
                    )
                  )
                )
