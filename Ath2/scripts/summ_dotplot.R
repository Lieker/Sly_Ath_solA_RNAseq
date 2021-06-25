sumplotdot <- function(x = "Ath2/output/panther/",
                       revigopath = "Ath2/output/panther/revigo/",
                       nutrcond = c("all","Pstarv","Nstarv","replete"),
                       updown = c("all", "up","down"),
                       sign = 0.05,
                       double = "no",
                       show = "all")
{
  files = list.files(path = "Ath2/output/panther/", pattern=".csv")
  sites <- files
  setwd(x)
  
  ans <- map2(files, sites, ~read_csv(.x) %>% mutate(id = .y))
  colunames <- c("GOterm","ref","ngenes","expected","over/under","foldenrich","rawPvalue","fdr","group")
  ans <- lapply(ans, setNames, colunames)
  setwd("~/github/Ath_RNAseq")
  
  eureka <- do.call(rbind, ans)
  eureka$group <- gsub("GOann-pantherfdr_","", eureka$group)
  eureka$group <- gsub(".csv","",eureka$group)
  eureka$foldenrich <- as.numeric(as.character(eureka$foldenrich))
  
  out <- strsplit(as.character(eureka$group),'_')
  
  out <- sapply(out, '[', seq(max(sapply(out, length)))) %>% t() %>% as.data.frame()
  names(out) <- c("nutritioncondition", "whichgenes")
  out[is.na(out)] <- "all"
  eureka <- cbind(eureka, out)
  
  eureka <- eureka %>% dplyr::filter(nutritioncondition %in% nutrcond)
  eureka <- eureka %>% dplyr::filter(whichgenes %in% updown)
  eureka <- eureka %>% dplyr::filter(eureka$fdr < sign)
  
  files = list.files(path = "Ath2/output/panther/revigo/", pattern=".csv")
  sites <- files
  setwd(revigopath)
  ans <- map2(files, sites, ~read_csv(.x) %>% mutate(id = .y))
  setwd("~/github/Ath_RNAseq")
  revigo <- do.call(rbind, ans)
  revigo <- revigo[revigo$Eliminated == FALSE, ]
  revigo$id <- gsub("Selfcr_", "", revigo$id)
  revigo$id <- gsub(".csv", "", revigo$id)
  revigo_out <- strsplit(as.character(revigo$id), "_")
  revigo_out <- sapply(revigo_out, '[', seq(max(sapply(revigo_out, length)))) %>% t() %>% as.data.frame()
  names(revigo_out) <- c("nutritioncondition", "whichgenes")
  revigo_out[is.na(revigo_out)] <- "all"
  revigo <- cbind(revigo, revigo_out)
  

  if(double == "no"){
    eureka <- eureka
  } else if(double == "yes") {
    e <- unique(eureka[duplicated(eureka$GOterm), ])
    eureka <- eureka[eureka$GOterm %in% e$GOterm, ]
  } else if(double == "revigo") {
    eureka_go <- strsplit(as.character(eureka$GOterm), "GO" )
    eureka_go <- sapply(eureka_go, '[', seq(max(sapply(eureka_go, length)))) %>% t() %>% as.data.frame()
    eureka_go$V2 <- gsub(":", "GO:", eureka_go$V2)
    eureka_go$V2 <- substr(eureka_go$V2, start = 1, stop = 10)
    eureka <- cbind(eureka, eureka_go$V2)
    names(eureka)[12] <- "TermID"
    eureka <- inner_join(revigo, eureka, by = c("TermID", "whichgenes", "nutritioncondition"))
  } else if(double == "triple") {
    tt <- table(eureka$GOterm)
    eureka <- subset(eureka, GOterm %in% names(tt[tt >= 3]))
  }
  
  eureka <- eureka[eureka$ref < 1000, ]

  
  if(show == "all"){
    eureka <- eureka
  } else if(class(show) == "numeric" | class(show) == "integer") {
    l <- unique(eureka$group)
    if(length(l) == 1) {
    eureka <- head(eureka, show)
    } else if(length(l) == 2) {
      eureka1 <- eureka[eureka$group == l[1], ]
      eureka2 <- eureka[eureka$group == l[2], ]
      eureka <- rbind(head(eureka1, show), head(eureka2, show))
    } else if(length(l) == 3) {
      eureka1 <- eureka[eureka$group == l[1], ]
      eureka2 <- eureka[eureka$group == l[2], ]
      eureka3 <- eureka[eureka$group == l[3], ]
      eureka <- rbind(head(eureka1, show), head(eureka2, show), head(eureka3, show))
    }
  }
  
  eureka$group <- factor(eureka$group, levels = c("replete_down", "Pstarv_down", "NPstarv_down"), ordered = TRUE)
  eureka <- eureka[order(eureka$group), ]

  p <- ggplot(eureka, 
              aes(x=foldenrich, 
                  y=reorder(GOterm, 
                            foldenrich), 
                  size=ngenes, 
                  col=fdr)) + 
    geom_point(alpha = 0.8) + 
    theme_minimal() +
    scale_colour_gradient(low = "firebrick1",
                          high = "blue4",
                          space = "Lab", 
                          limit = c(0, sign),
                          breaks = c(0,(sign/5), ((sign/5)*2), ((sign/5)*3), ((sign/5)*4), sign)) +
    theme(axis.title.y=element_blank(),
          text = element_text(size=14))+ 
    facet_wrap(~group) +
    xlab("Enrichment factor")
  
  return(p)
  
}
