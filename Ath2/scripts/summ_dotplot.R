sumplotdot <- function(x = "Ath2/output/panther/",
                       nutritioncondition = c("all","Pstarv","Nstarv","replete"),
                       updown = c("all", "up","down"),
                       sign = 0.05)
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
  
  eureka <- eureka %>% dplyr::filter(nutritioncondition %in% nutritioncondition)
  eureka <- eureka %>% dplyr::filter(whichgenes %in% updown)
  eureka <- eureka %>% dplyr::filter(eureka$fdr < sign)
  
  p <- ggplot(eureka, 
              aes(x=group, 
                  y=GOterm, 
                  size=foldenrich, 
                  color=fdr)) + 
    geom_point(alpha = 0.8) + 
    theme_minimal() +
    scale_color_gradient(low = "firebrick1",
                         high = "blue4",
                         space = "Lab", 
                         limit = c(0, sign)) +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
          text = element_text(size=8)) +
    coord_fixed(ratio = 0.4)
  return(p)
  
}
