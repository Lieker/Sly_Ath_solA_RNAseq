server <- function(input, output, session) {

#########################################
### set output texts to starting position
#########################################
  
  source("introduction.R")
  source("Plots.R")
  output$intro1 <- renderText(introductie1())
  output$intro2 <- renderText(introductie2())
  output$counts1 <- renderTable(
    exampleFeatureCounts()
  )
  output$tab1 <- renderText(header1())
  output$counts2 <- renderTable(
    exampleProcessedCounts()
  )
  output$tab2 <- renderText(header2())
  output$sample <- renderTable(
    exampleSample()
  )
  output$intro3 <- renderText(introductie3())
output$text1 <- renderText(
  {"select a featureCounts raw counts file and press the <font color=\"#2da5fa\">blue</font> submit."})
output$text3 <- renderText(
  {"<h4><font color=\"#fa461e\"><b>First select counts and sample files in the \"Input\" tab.</b></font></h4>"})
output$text9 <- renderText(
  {"<font color=\"#fa461e\"><b>First select counts and sample files in the \"Input\" tab.</b></font>"})
output$text5 <- renderText(
  {"<font color=\"#fa461e\"><b>First select counts and sample files in the \" Input\" tab.</b></font>"})
output$text6 <- renderText(
  {"<font color=\"#fa461e\"><b>First select counts and sample files in the \"Raw Counts Input\" tab.</b></font>"})
output$text10 <- renderText(
  "When using Dynamic Hierarchical Clustering, the number of clusters is determined by the algorithm, and can not be manually fixed.<br><br>")
output$text11 <- renderText(
  {"<font color=\"#fa461e\"><b>First select counts and sample files in the \"Input\" tab.</b></font>"})
output$text14 <- renderText(
  {"<font color=\"#fa461e\"><b><br>First select counts and sample files in the \"Input\" tab.</b></font>"})
output$text15 <- renderText(
  {" "})

rv <- reactiveValues(df = NULL, ft = NULL, df2 = NULL, select = NULL, samples = NULL, scaledat = NULL, mydat=NULL,
                     avarages = NULL, scaledata = NULL, method = NULL, clusts = NULL, k = 1,
                     hr = NULL, my_palette = NULL, Kmolten = NULL, selection = NULL, centroids_long= NULL,
                     step = 0
                     )


########################################  
### upload counts file
########################################

observeEvent(input$goButton, {   #submit input counts.txt
  output$text1 <- renderText(
    "data table showing the raw counts values")
  output$text2 <- renderText(
    "Choose sample file containing sample names and conditions and press the <font color=\"#1ddc5b\">green</font> submit button")
  output$text3 <- renderText(
    {"<h4><font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font></h4>"})
  output$text9 <- renderText(
    {"<font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font>"})
  output$text5 <- renderText(
    {"<font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font>"})
  output$text6 <- renderText(
    {"<font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font>"})
  output$text11 <- renderText(
    {"<font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font>"})
  output$text14 <- renderText(
    {"<font color=\"#fa461e\"><b>First select a sample file in the \"Input\" tab.</b></font>"})
  print(rv$step)
  req(input$file1)
  rv$ft <- input$filetype
  if(rv$ft == 0) {
    df <- read.table(input$file1$datapath, 
                     header=T,
                     skip=1, 
                     row.names=1,
                     stringsAsFactors = F, 
                     check.names = F)
    rv$df <- df[ ,6:ncol(df)]
    rv$step=1
  } else {
    rv$df <- read.table(input$file1$datapath, 
                     header=T,
                     row.names=1,
                     stringsAsFactors = F, 
                     check.names = F)
  
    rv$step=1
  }
  })
  
  output$contents <- renderTable({   
    if(input$disp == "head") {
      return(head(rv$df))
    } else {
      return(rv$df)
    }
  }, rownames = T)

  
######################################
### upload samles file
######################################
  
env2 <- observeEvent(input$goButton2,{     #submit input sample.csv


    output$text3 <- renderText(
      "<h4>Please select minimum   fold change and maximum P-value and next press 
      <font color=\"#2da5fa\">submit</font>.
      <br><br><font color=\"#fa782d\">
      Please take note that the normalization takes up to a couple of minutes. So please have some patience.
      </font></h4>")
    output$text9 <- renderText(
      {"<font color=\"#fa461e\"><b>Data still needs to be normalized, go back to the \"Normalization\" tab.</b></font>"})
    output$text5 <- renderText(
      {"<font color=\"#fa461e\"><b>Data still needs to be normalized, go back to the \"Normalization\" tab.</b></font>"})
    output$text6 <- renderText(
      {"<font color=\"#fa461e\"><b>Data still needs to be normalized, go back to the \"Normalization\" tab.</b></font>"})
    output$text11 <- renderText(
      {"<font color=\"#fa461e\"><b>Data still needs to be normalized, go back to the \"Normalization\" tab.</b></font>"})
    output$text14 <- renderText(
      {"<font color=\"#fa461e\"><b>Data still needs to be normalized, go back to the \"Normalization\" tab.</b></font>"})

    if(rv$step == 0){
      step = 0
      output$text2 <- renderText(
        {"<font color=\"#fa461e\">First select a counts and sample file in the \"Raw Counts Input\" tab.</font>"})
    } else {
      output$text2 <- renderText(
        "data table showing names of samples and their condition")
      step = 2
    req(input$file2)
    rv$df2 <- read.delim(input$file2$datapath, 
                      header=T,
                      stringsAsFactors = F, 
                      check.names = F)
    }
    rv$step=step
    output$contents2 <- renderTable({
      if(input$disp2 == "head") {
        return(head(rv$df2))
      }
      else {
        return(rv$df2)
      } 
    }, rownames = T)
})
  

############################################################
### normalize counts data and create heatmap
############################################################
  
env3 <- observeEvent(input$goNorm,{      # submit to start normalisation
  if(rv$step == 0){
    step = 0
    output$text3 <- renderText(
      {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a counts and sample file in the \"Input\" tab.</center></font></h2>"})
  } else if(rv$step == 1){
    step = 1
    output$text3 <- renderText(
      {"<h2><font color=\"#fa461e\"><center>No really, You first have to select a sample file in the \"Input\" tab.</center></font></h2>"})
  } else {
    output$text3 <- renderText(
      "<h4>Scroll down for table containing normalized values (head or full table)</h4>")
    output$text9 <- renderText(
      "This step is optional. It might help you to choose a number of clusters to be used in the clustering.But as you will see the different included algorithms can predict a different number.
      In the end it is best to try different methods of clustering and clusternumbers to see wich combination fits your data best.")
    output$text5 <- renderText(
      "The clustering can be done with different methods. As this step runs pretty fast it is possible to try any of the given methods with different numbers of clusters, either based on the optimal number algorithms or based on your own choose.")
    output$text11 <- renderText(
      "<font color=\"#fa461e\"><b>First the clustering needs to be done in the tab \"clustering\" tab.</b></font>")
    output$text14 <- renderText(
      "Choose a single gene or multiple genes (comma separated) to get expression plot")
    step = 3
  foldchange <- log2(input$f)
  pvalue <- input$p 
  source("DESeqNormalize.R")
  normData <- normalize(rv$df, rv$df2, foldchange, pvalue)
  rv$select <- normData$selection
  rv$samples <- normData$samples
  rv$scaledat <- normData$scaledata
  rv$avarages <- normData$avarages
  rv$scaledata <- rv$scaledat
  output$text3 <- renderText(
    paste0("<h4><font color=\"#1ddc5b\">", nrow(rv$select), "</font> out of a total of <font color=\"#fc0000\">",nrow(rv$df), "</font> genes were found to be significantly differientially expressed.<br>
           Scroll down for table containing normalized values (head or full table)</h4>"))
  
  
  rv$mydat <- rv$select
  rv$mydat$ID <- rownames(rv$mydat)
  
  rv$hr <- hclust(as.dist(1-cor(t(rv$scaledat), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
  TreeR = as.dendrogram(rv$hr, method="average")
  # colors of the heatmap
  rv$my_palette <- colorRampPalette(c("white","green","green4","violet","purple"))(100)
  #(c("magenta", "black", "green"))(n = 299)
  
  
  output$heatmap1 <- renderPlot({heatmap.2(rv$scaledat,
                                           na.rm = TRUE,
                                           Rowv=as.dendrogram(rv$hr), 
                                           Colv=NA,
                                           col=rv$my_palette,
                                           scale="row",
                                           margins = c(7, 7),
                                           cexCol = 0.7,
                                           labRow = F,
                                           dendrogram = "row",
                                           main = "Heatmap.2",
                                           trace = "none")}, height=800)
  
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(head(rv$select))
    }
    else {
      return(rv$select)
    }
  }, rownames = T)
  }
  rv$step=step
}) 
#########################################################
### determin optimal number of clusters (optinal feature)
#########################################################
 
  env4 <- observeEvent(input$goButton3,{      # submit to start normalisation
    
    if(rv$step == 0){
      step = 0
      output$text9 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a counts and sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 1){
      step = 1
      output$text9 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 2){
      step = 2
      output$text9 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> Your data still needs to be normalized, go back to the \"Normalization\" tab.</center></font></h2>"})

      
     } else {
      
    if(input$opt_K == 1){
      output$text4 <- renderText(
        "<br>")
      wss <- (nrow(rv$scaledata)-1)*sum(apply(rv$scaledata,2,var))
      for (i in 2:20) wss[i] <- sum(kmeans(rv$scaledata,
                                           centers=i)$withinss)
      output$opt_K_plot <- renderPlot({plot(1:20, wss, type="b", xlab="Number of Clusters",
           ylab="Within groups sum of squares")}, height=600)
      output$text9 <- renderText(
        "Optimal number of clusters based on the <font color=\"#fab834\">elbow</font> method is where line has the greatest bend")
    } else if(input$opt_K == 2){
      output$text4 <- renderText(
        "<br>")
      sil <- rep(0, 20)
      #repeat k-means for 1:20 and extract silhouette:
      for(i in 2:20){
        k1to20 <- kmeans(rv$scaledata, centers = i, nstart = 25, iter.max = 20)
        ss <- silhouette(k1to20$cluster, dist(rv$scaledata))
        sil[i] <- mean(ss[, 3])
      }
      # Plot the  average silhouette width
      output$opt_K_plot <- renderPlot({plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
      abline(v = which.max(sil), lty = 2)
      }, height=600) 
      output$text9 <- renderText({
        paste("Optimal number of clusters based on the <font color=\"#fcf800\">Average_silhouette_width</font> method is <font color=\"#fcf800\">", which.max(sil), "</font>")})

      
    } else if(input$opt_K == 3){
      output$text4 <- renderText(
        "<br>") 
      fit <- cascadeKM(rv$scaledata, 1, 20, iter = 100)
      output$opt_K_plot <- renderPlot({plot(fit, sortg = TRUE, grpmts.plot = TRUE)
      }, height=600)
      calinski.best <- as.numeric(which.max(fit$results[2,]))
      output$text9 <- renderText({
        paste("Optimal number of clusters based on the <font color=\"#fa5ff5\">Calinsky_criterion</font> method is <font color=\"#fa5ff5\">", calinski.best, "</font>")})

    } else if(input$opt_K == 4){
      output$text4 <- renderText(
        "<br>") 
      set.seed(13)
      gap <- clusGap(rv$scaledata, kmeans, 20, B = 100, verbose = interactive())
      output$opt_K_plot <- renderPlot({plot(gap, main = "Gap statistic")
      abline(v=which.max(gap$Tab[,3]), lty = 2)}, height=600)
      output$text9 <- renderText({
        paste("Optimal number of clusters based on the <font color=\"#17fc03\">gap_statistic</font> method is <font color=\"#17fc03\">", which.max(gap$Tab[,3]), "</font>")})
      
    } else if(input$opt_K == 5){
      output$text4 <- renderText(
        "<br>")   
      d.apclus <- apcluster(negDistMat(r=2), rv$scaledata)
      cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
      #uncomment the next line for the heatmap, it takes a long time to run
      #heatmap(d.apclus,cexRow=0, cexCol=0)
      output$opt_K_plot <- renderPlot(NULL)
      output$text9 <- renderText({
        paste("Optimal number of clusters based on the <font color=\"#5fd8fa\">affinity_propogation</font> method is <font color=\"#5fd8fa\">", length(d.apclus@clusters), "</font>")})
 
      } 
    }
  })

  
#############################################
### The Clustering
#############################################
  
  env4 <- observeEvent(input$goClus,{
    
    if(input$Avarage == "no"){
      rv$scaledata = rv$scaledat
    } else if(input$Avarage == "yes"){
      rv$scaledata = rv$avarages
    }
    
    if(rv$step == 0){
      step = 0
      output$text5 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a counts and sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 1){
      step = 1
      output$text5 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 2){
      step = 2
      output$text5 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> Your data still needs to be normalized, go back to the \"Normalization\" tab.</center></font></h2>"})
    } else {
      step = 4
    rv$method <- input$ClusMet
    output$text11 <- renderText("Choose cluster number you are interested in."
       )
    output$text16 <- renderText("<center>The list containing the differentially expressed genes and the cluster number they have been assigned to can be downloaded by pressing the <font color=\"#f7c00a\"><b>download</b></font> button.</center>"
       )
    
       if(rv$method == "hierarchical"){
          output$text5 <- renderText({paste("Clustering has been performed using", rv$method, " with ", input$selectH, "clusters" )})
          rv$clusts = cutree(rv$hr, k=input$selectH)
          rv$k = input$selectH
          
          mycolhc1 <- rainbow(length(unique(rv$clusts)), start=0.1, end=0.9)
          mycolhc <- mycolhc1[as.vector(rv$clusts)]
          
          output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                    Rowv=as.dendrogram(rv$hr), 
                    Colv=NA,
                    col=rv$my_palette,
                    scale="row",
                    margins = c(7, 7),
                    cexCol = 0.7,
                    labRow = F,
                    dendrogram = "row",
                    main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                    trace = "none",
                    RowSideColors=mycolhc,
                    key = FALSE)
          }, height=800)
          clust.centroid = function(i, dat, clusters) {
            ind = (clusters == i)
            colMeans(dat[ind,])
          }
          clusUniq <- unique(rv$clusts)
          kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
          
          rv$Kmolten <- melt(kClustcentroids)
          colnames(rv$Kmolten) <- c('sample','cluster','value')
          
          output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) +
            scale_colour_manual(values = mycolhc1) +
            geom_point() + 
            geom_line() +
            xlab("Time") +
            ylab("Expression") +
            labs(title= "Cluster Expression of the samples",color = "Cluster")
          })

          
          
       }else if(rv$method == "Dynamic Hr"){  
         output$text5 <- renderText({paste("Cluster will run with", rv$method)})
         clusts <- cutreeDynamic(rv$hr, distM = as.matrix(as.dist(1-cor(t(rv$scaledata)))), method = "hybrid")
         clusts <- clusts +1

         names(clusts) <- rownames(rv$scaledata)
         #clusts <- as.data.frame(clusts)
         #colnames(clusts) <- "cluster"
         rv$k = length(unique(clusts))
         rv$clusts <- clusts

         mycolhc1 <- rainbow(length(unique(clusts)), start=0.1, end=0.9)
         mycolhc <- mycolhc1[as.vector(clusts)]

         output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                                      Rowv=as.dendrogram(rv$hr), 
                                      Colv=NA,
                                      col=rv$my_palette,
                                      scale="row",
                                      margins = c(7, 7),
                                      cexCol = 0.7,
                                      labRow = F,
                                      dendrogram = "row",
                                      main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                                      trace = "none",
                                      RowSideColors=mycolhc,
                                      key = FALSE)
         }, height=800)
         clust.centroid = function(i, dat, clusters) {
           ind = (clusters == i)
           colMeans(dat[ind,])
         }
         clusUniq <- unique(rv$clusts)
         kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
         
         rv$Kmolten <- melt(kClustcentroids)
         colnames(rv$Kmolten) <- c('sample','cluster','value')
         
         output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
             scale_colour_manual(values = mycolhc1) +
             geom_point() + 
             geom_line() +
             xlab("Time") +
             ylab("Expression") +
             labs(title= "Cluster Expression of the samples",color = "Cluster")
         })
          
       }else if(rv$method == "K-mean"){
         rv$k = as.integer(input$selectK)
         output$text5 <- renderText({paste("Clustering has been performed using", rv$method, " with ", input$selectK, "clusters")})
         clusts <- kmeans(rv$scaledata, centers = rv$k, nstart = 1000, iter.max = 20)
         rv$clusts <- clusts$cluster
         
         mycolhc1 <- rainbow(length(unique(rv$clusts)), start=0.1, end=0.9)
         mycolhc <- mycolhc1[as.vector(rv$clusts)]

         output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                                                  Rowv=as.dendrogram(rv$hr), 
                                                  Colv=NA,
                                                  col=rv$my_palette,
                                                  scale="row",
                                                  margins = c(7, 7),
                                                  cexCol = 0.7,
                                                  labRow = F,
                                                  dendrogram = "row",
                                                  main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                                                  trace = "none",
                                                  RowSideColors=mycolhc,
                                                  key = FALSE)
         }, height=800)
         clust.centroid = function(i, dat, clusters) {
           ind = (clusters == i)
           colMeans(dat[ind,])
         }
         clusUniq <- unique(rv$clusts)
         kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
         
         rv$Kmolten <- melt(kClustcentroids)
         colnames(rv$Kmolten) <- c('sample','cluster','value')
         
         output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
             scale_colour_manual(values = mycolhc1) +
             geom_point() + 
             geom_line() +
             xlab("Time") +
             ylab("Expression") +
             labs(title= "Cluster Expression of the samples",color = "Cluster")
         })

               
         
       }else if(rv$method == "Fuzzy K-mean"){  
         output$text5 <- renderText({paste("Cluster will run with", rv$method)})
         mestimate<- function(df){
           N <-  dim(df)[[1]]
           D <- dim(df)[[2]]
           m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
           return(m.sj)
         }
         rv$k <- as.integer(input$selectF)
         mes <- mestimate(rv$scaledata)
         fcm_results <- cmeans(rv$scaledata,centers=rv$k,m=mes)
         fcm_plotting_df <- data.frame(rv$scaledata)
         fcm_plotting_df$gene <- row.names(fcm_plotting_df)
         fcm_plotting_df$cluster <- fcm_results$cluster
         fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster),function(row){
           clust <- fcm_plotting_df$cluster[row]
           fcm_results$membership[row,clust]
         })
         clusts <- fcm_plotting_df$cluster
         names(clusts) <- fcm_plotting_df$gene
         # get the core data
         fcm_centroids <- fcm_results$centers
         fcm_centroids_df <- data.frame(fcm_centroids)
         fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
         rv$centroids_long <- tidyr::gather(fcm_centroids_df,"sample",'value', 1:ncol(rv$scaledata))
         #start with the input data
         fcm_plotting_df <- data.frame(rv$scaledata)
         #add genes
         fcm_plotting_df$gene <- row.names(fcm_plotting_df)
         #bind cluster assinment
         fcm_plotting_df$cluster <- fcm_results$cluster
         #fetch the membership for each gene/top scoring cluster
         fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster),function(row){
           clust <- fcm_plotting_df$cluster[row]
           fcm_results$membership[row,clust]
         })
         # filter out genes don't really belong to any of the clusters
         rv$selection <- fcm_plotting_df %>% filter(membership > 0.2)
         rv$clusts <- clusts

         mycolhc1 <- rainbow(length(unique(rv$clusts)), start=0.1, end=0.9)
         mycolhc <- mycolhc1[as.vector(rv$clusts)]

         output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                                                  Rowv=as.dendrogram(rv$hr), 
                                                  Colv=NA,
                                                  col=rv$my_palette,
                                                  scale="row",
                                                  margins = c(7, 7),
                                                  cexCol = 0.7,
                                                  labRow = F,
                                                  dendrogram = "row",
                                                  main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                                                  trace = "none",
                                                  RowSideColors=mycolhc,
                                                  key = FALSE)
         }, height=800)
         clust.centroid = function(i, dat, clusters) {
           ind = (clusters == i)
           colMeans(dat[ind,])
         }
         clusUniq <- unique(rv$clusts)
         kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
         
         rv$Kmolten <- melt(kClustcentroids)
         colnames(rv$Kmolten) <- c('sample','cluster','value')
         
         output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) +
             scale_colour_manual(values = mycolhc1) +
             geom_point() + 
             geom_line() +
             xlab("Time") +
             ylab("Expression") +
             labs(title= "Cluster Expression of the samples",color = "Cluster")
         })
         

       }else if(rv$method == "PAM"){
         output$text5 <- renderText({paste("Clustering has been performed using", rv$method, " with ", input$selectP, "clusters")})
         pam.Pclusts <- pam(rv$scaledata, k=input$selectP)
         rv$clusts <- pam.Pclusts$'clustering'
         rv$k = input$selectP
         
         mycolhc1 <- rainbow(length(unique(rv$clusts)), start=0.1, end=0.9)
         mycolhc <- mycolhc1[as.vector(rv$clusts)]
         
         output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                                                  Rowv=as.dendrogram(rv$hr), 
                                                  Colv=NA,
                                                  col=rv$my_palette,
                                                  scale="row",
                                                  margins = c(7, 7),
                                                  cexCol = 0.7,
                                                  labRow = F,
                                                  dendrogram = "row",
                                                  main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                                                  trace = "none",
                                                  RowSideColors=mycolhc,
                                                  key = FALSE)
         }, height=800)
         clust.centroid = function(i, dat, clusters) {
           ind = (clusters == i)
           colMeans(dat[ind,])
         }
         clusUniq <- unique(rv$clusts)
         kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
         
         rv$Kmolten <- melt(kClustcentroids)
         colnames(rv$Kmolten) <- c('sample','cluster','value')
         
         output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
             scale_colour_manual(values = mycolhc1) +
             geom_point() + 
             geom_line() +
             xlab("Time") +
             ylab("Expression") +
             labs(title= "Cluster Expression of the samples",color = "Cluster")
         }, height=600)
         
       }else if(rv$method == "dpscan"){
         output$text5 <- renderText({paste("Clustering has been performed using", rv$method, " with ", input$selectK, "clusters")})
         db <- fpc::dbscan(rv$scaledata, eps = input$eps, MinPts = input$minPts)
         rv$k <- length(unique(db$cluster))
         clusts <- db$cluster
         rv$clusts = clusts + 1
         names(rv$clusts) <- rownames(rv$scaledata)
         
         mycolhc1 <- rainbow(length(unique(rv$clusts)), start=0.1, end=0.9)
         mycolhc <- mycolhc1[as.vector(rv$clusts)]
         
         output$heatmap2 <- renderPlot({heatmap.2(rv$scaledata,
                                                  Rowv=as.dendrogram(rv$hr), 
                                                  Colv=NA,
                                                  col=rv$my_palette,
                                                  scale="row",
                                                  margins = c(7, 7),
                                                  cexCol = 0.7,
                                                  labRow = F,
                                                  dendrogram = "row",
                                                  main = paste0("Heatmap.2 with colour bar indicating the clusters (",rv$k,")"),
                                                  trace = "none",
                                                  RowSideColors=mycolhc,
                                                  key = FALSE)
           
         }, height=800)
         clust.centroid = function(i, dat, clusters) {
           ind = (clusters == i)
           colMeans(dat[ind,])
         }
         clusUniq <- unique(rv$clusts)
         kClustcentroids <- sapply(levels(factor(rv$clusts)), clust.centroid, rv$scaledata, rv$clusts)
         
         rv$Kmolten <- melt(kClustcentroids)
         colnames(rv$Kmolten) <- c('sample','cluster','value')
         
         output$cluCen <- renderPlot({ggplot(rv$Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
             scale_colour_manual(values = mycolhc1) +
             geom_point() + 
             geom_line() +
             xlab("Time") +
             ylab("Expression") +
             labs(title= "Cluster Expression of the samples",color = "Cluster")
         }, height=600)
       }
    clus1 <- seq(1,rv$k,1)
    updateSelectInput(session, "cluster",
                      choices = clus1)
    }
    rv$step <- step
  })
  
  
  
  env4 <- observeEvent(input$sinClus,{
    
    if(rv$step == 0){
      step = 0
      output$text11 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a counts and sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 1){
      step = 1
      output$text11 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 2){
      step = 2
      output$text11 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> Your data still needs to be normalized, go back to the \"Normalization\" tab.</center></font></h2>"})
    } else if(rv$step == 3){
      step = 3
      output$text11 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>Please,<br> Your data still needs to be clustered, go back to the \"Clustering\" tab.</center></font></h2>"})
    } else {
      output$text11 <- renderText({paste0(input$cluster,  " were found while clustering with ", rv$method)})
  
      rv$i <- input$cluster
      output$text11 <- renderText({paste0("Cluster ", rv$i,  " were found while clustering with ", rv$method)})
      if(rv$method == "Fuzzy K-mean"){


          cluster_plot_df <- dplyr::filter(rv$selection, cluster == rv$i) %>%
            dplyr::select(.,1:ncol(rv$scaledata),membership,gene) %>%
            tidyr::gather(.,"sample",'value',1:ncol(rv$scaledata)) 
          
          K2 <- cluster_plot_df
          K1 <- (rv$scaledata[rv$clusts==rv$i,])
          clusterTable <- (rv$select[row.names(K1),])
          #order the dataframe by score
          cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership),]
          #set the order by setting the factors using forcats
          cluster_plot_df$gene = forcats::fct_inorder(cluster_plot_df$gene)
          
          #subset the cores by cluster
          core <- dplyr::filter(rv$centroids_long, cluster == rv$i)
          
          output$SinClus  <- renderPlot({ggplot(cluster_plot_df, aes(x=sample,y=value)) + 
            geom_line(aes(colour=membership, group=gene)) +
            scale_colour_gradientn(colours=c('blue1','red2')) +
            #this adds the core 
            geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
            geom_point(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
            xlab("Sample") +
            ylab("Expression") +
            labs(title= paste0("Cluster ",rv$i," Expression of the Samples"),color = "Score")
          })
          
          group <- rv$selection %>% filter(rv$selection$cluster==rv$i)
          rownames(group) <- group$gene
          group <- as.matrix(group[,1:ncol(rv$scaledata)])
          if(nrow(group) > 1){
            hrc <- hclust(as.dist(1-cor(t(group), method="pearson")), method="complete")
            #pheatmap(group, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = F)  #cellheight = 8
            output$heatmap3 <- renderPlot({heatmap.2(group,
                      Rowv=as.dendrogram(hrc), 
                      Colv=NA,
                      col= rv$my_palette,
                      labRow=rownames(group),
                      scale="row",
                      margins = c(7, 7),
                      cexCol = 0.7,
                      dendrogram = "row",
                      main = paste0("Cluster ",rv$i, " consisting ", nrow(group), " genes"),
                      trace = "none",key = F)
            })
          }

      }else if(rv$method == "K-mean" || rv$method == "hierarchical" || rv$method == "PAM" || rv$method == "Dynamic Hr" || rv$method == "dpscan"){

          K2 <- (rv$scaledata[rv$clusts==rv$i,])
          clusterTable <- (rv$select[row.names(K2),])
          if(nrow(K2) > 1){
            Kmolten <- rv$Kmolten
            #calculate the correlation with the core
            core <- Kmolten[Kmolten$cluster==rv$i,]
            corscore <- function(x){cor(x,core$value)}
            score <- apply(K2, 1, corscore)
            # add the correlations to the list
            #cors <- c(cors, score)
            #get the data frame into long format for plotting
            K2molten <- melt(K2)
            colnames(K2molten) <- c('gene','sample','value')
            #add the score
            K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
            colnames(K2molten) <- c('gene','sample','value','score')
            #order the dataframe by score
            #to do this first create an ordering factor
            K2molten$order_factor <- 1:length(K2molten$gene)
            #order the dataframe by score
            K2molten <- K2molten[order(K2molten$score),]
            #set the order by setting the factors
            K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
            
            # Everything on the same plot
            if(rv$method == "K-mean" || rv$method =="PAM"){
              output$SinClus <- renderPlot({ggplot(K2molten, aes(x=sample,y=value)) + 
                geom_line(aes(colour=score, group=gene)) +
                scale_colour_gradientn(colours=c('blue1','red2')) +
                #this adds the core 
                geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
                geom_point(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
                xlab("Samples") +
                ylab("Expression") +
                labs(title= paste0("Cluster ",rv$i, " consisting ", nrow(K2), " genes"),color = "Score")
              })
            }
            else if(rv$method == "hierarchical" || rv$method == "Dynamic Hr" || rv$method == "dpscan"){
              output$SinClus <- renderPlot({ggplot(K2molten, aes(x=sample,y=value)) + 
                geom_line(color="grey", aes(color="grey", group=gene)) +
                #this adds the core 
                geom_line(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
                geom_point(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
                xlab("Samples") +
                ylab("Expression") +
                labs(title= paste0("Cluster ",rv$i, " consisting ", nrow(K2), " genes"),color = "Score")
              })
            }

            hrc <- hclust(as.dist(1-cor(t(K2), method="pearson")), method="complete")
            #pheatmap(group, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = F)  #cellheight = 8
            output$heatmap3 <- renderPlot({heatmap.2(K2,
                      Rowv=as.dendrogram(hrc), 
                      Colv=NA,
                      col= rv$my_palette,
                      labRow=rownames(K2),
                      scale="row",
                      margins = c(7, 7),
                      cexCol = 0.7,
                      dendrogram = "row",
                      main = paste0("Cluster ",rv$i, " of ", rv$k, " containing ", nrow(K2), " genes"),
                      trace = "none",key = F)
            })
          }

      }  
      clusterTable <- tibble::rownames_to_column(clusterTable, "gene")
      output$text12 <- renderText(
        "Table showing normalized read counts")
      output$clusTab <- renderDataTable(datatable(clusterTable, 
                                                  options = list(pageLength = 10,lengthChange=FALSE), rownames= FALSE) %>%
                                          formatRound(c(2:ncol(clusterTable)), 1) %>% 
                                          formatStyle(c(1:ncol(clusterTable)), 'text-align' = 'center',  color = 'white', backgroundColor = 'black'))
        
        
        
        
        
        #renderDataTable({        clusterTable      })
    }
 })
  
  gene <- eventReactive(input$do, {
    if(rv$step == 0){
      step = 0
      output$text14 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a counts and sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 1){
      step = 1
      output$text14 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> You first have to select a sample file in the \"Input\" tab.</center></font></h2>"})
    } else if(rv$step == 2){
      step = 2
      output$text14 <- renderText(
        {"<h2><font color=\"#fa461e\"><center>No really,<br> Your data still needs to be normalized, go back to the \"Normalization\" tab.</center></font></h2>"})
    } else {
    output$text14 <- renderText(
      "plot showing normalized counts values in the different samples.")
    output$text15 <- renderText(
      "plot showing normalized counts values averaged to the conditions")
    unlist(strsplit(as.character(input$gene), ',', fixed=TRUE))  }
  }, ignoreNULL= T)
  geneCounts <- reactive({
    rv$mydat[rownames(rv$mydat) %in% c(gene()),]
    })
  AvGene <- reactive({
    Avar(row.names(geneCounts()),rv$select,rv$samples)
    })
  geneCountsLong <- reactive({
    melt(geneCounts())
    })
  p1 <- reactive({
    ggplot(geneCountsLong(), aes(x=geneCountsLong()$variable, y=geneCountsLong()$value, group=ID, fill=ID)) +
      geom_bar(position="dodge", stat="identity") +
      xlab("Sample") +
      ylab("counts")
  })
  output$plot1 <- renderPlot({
    print(p1())
  })
  p2 <- reactive({
    ggplot(AvGene()$AV, aes(x=AvGene()$AV$Var2, beside=T, y=AvGene()$AV$value, group=Var1 , fill=Var1 )) + 
      geom_bar(position="dodge", stat="identity") +
      geom_errorbar(aes(ymin=value-AvGene()$ST$value, ymax=value+AvGene()$ST$value),position="dodge", stat="identity",width=0.9,alpha=0.9, colour="black") +
      xlab("Sample") +
      ylab("counts")
  })
  output$plot2 <- renderPlot({
    print(p2())
  })
  thedata <- reactive(merge(as.data.frame(rv$select), as.data.frame(rv$clusts), by="row.names", sort=FALSE))
  output$download <- downloadHandler(
    filename = function(){"DE_genes.csv"}, 
    content = function(fname){
      write.csv(thedata(), fname)
    }
  )
}