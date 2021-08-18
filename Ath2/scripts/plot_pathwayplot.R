plot_pathwayplot <- function(t){
  tl <- length(unique(t$variable))
  t$solA <- gsub("yes", "SolA", t$solA)
  t$solA <- gsub("no", "cntrl", t$solA)
  t <- data.frame(lapply(t, function(x) {gsub("yes", "+", x)}))
  t <- data.frame(lapply(t, function(x) {gsub("no", "-", x)}))
  t$group <- paste0("N",t$N,"P",t$P)
  t$value <- as.numeric(as.character(t$value))
  t <- t %>% group_by(variable,group,solA) %>% summarise(sd = sd(value, na.rm = TRUE),
                                                         meancounts = mean(value))
  t$group2 <- paste0(t$variable,"_",t$solA)
  t <- t[!(t$group == "N-P+"),]
  
  t1 <- t[t$group == "N-P-",]
  tsum1 <- t1 %>% group_by(solA) %>% summarise(sd = sd(meancounts, na.rm = TRUE),
                                               meancounts = mean(meancounts))
  delta1 <- tsum1[2,3] - tsum1[1,3]
  t2 <- t[t$group == "N+P-",]
  tsum2 <- t2 %>% group_by(solA) %>% summarise(sd = sd(meancounts, na.rm = TRUE),
                                               meancounts = mean(meancounts))
  delta2 <- tsum2[2,3] - tsum2[1,3]
  t3 <- t[t$group == "N+P+",]
  tsum3 <- t3 %>% group_by(solA) %>% summarise(sd = sd(meancounts, na.rm = TRUE),
                                               meancounts = mean(meancounts))
  delta3 <- tsum3[2,3] - tsum3[1,3]
  windowsFonts(Times=windowsFont("TT Times New Roman"))
  g1 <- ggplot(data = t1, aes(x = solA, y = meancounts)) +
    geom_line(aes(group = variable),
              size = 0.5, 
              alpha = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank()) +
    scale_y_continuous(limits = c(1,13000),
                       trans = "log10") +
    geom_line(data = tsum1,
              aes(x = solA, y = meancounts, group = 1),
              size = 1.5,
              colour = "red") +
    ggtitle(paste0(t1$group)) +
    geom_point(data = tsum1,
               aes(x = solA, y = meancounts),
               colour = "red",
               size = 3) +
    geom_text(x = 1.5, y = 0.4, size = 4,
              mapping = aes(label = paste0(expression("\u0394"),"= ",round(delta1))),
              family = "Arial") +
    geom_text(x=1.5, y = 0.1, size = 4,
              mapping = aes(label = paste0("n=",tl)),
              family = "Arial")
  
  g2 <- ggplot(data = t2, aes(x = solA, y = meancounts)) +
    geom_line(aes(group = variable),
              size = 0.5, 
              alpha = 0.2) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(limits = c(1,13000),
                       trans = "log10") +
    geom_line(data = tsum2,
              aes(x = solA, y = meancounts, group = 1),
              size = 1.5,
              colour = "red") +
    ggtitle(paste0(t2$group)) +
    geom_point(data = tsum2,
               aes(x = solA, y = meancounts),
               colour = "red",
               size = 3) +
    geom_text(x = 1.5, y = 0.4, size = 4,
              mapping = aes(label = paste0(expression("\u0394"),"= ",round(delta2))),
              family = "Arial")  +
    geom_text(x=1.5, y = 0.1, size = 4,
              mapping = aes(label = paste0("n=",tl)),
              family = "Arial")
  
  g3 <- ggplot(data = t3, aes(x = solA, y = meancounts)) +
    geom_line(aes(group = variable),
              size = 0.5, 
              alpha = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(limits = c(1,13000),
                       trans = "log10") +
    geom_line(data = tsum3,
              aes(x = solA, y = meancounts, group = 1),
              size = 1.5,
              colour = "red") +
    ggtitle(paste0(t3$group)) +
    geom_point(data = tsum3,
               aes(x = solA, y = meancounts),
               colour = "red",
               size = 3)+
    geom_text(x = 1.5, y = 0.4, size = 4,
              mapping = aes(label = paste0(expression("\u0394"),"= ",round(delta3))),
              family = "Arial") +
    geom_text(x=1.5, y = 0.1, size = 4,
              mapping = aes(label = paste0("n=",tl)),
              family = "Arial")
  g <- g1 + g2 + g3
  return(g)
}
