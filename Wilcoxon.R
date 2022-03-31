library(ggplot2)
library(ggthemes)
library(latex2exp)
library(gridExtra)
wilcox.simple <- function(MSE, miss.no){
  RMSE <- sqrt(MSE)
  method.to.compare <- colnames(RMSE)
  prop.win <- matrix(0, nrow = 1, ncol = ncol(RMSE))
  colnames(prop.win) <- method.to.compare
  prop.temp <- table(apply(RMSE,1,which.min))/nrow(RMSE)
  for(i in names(prop.temp)){
    prop.win[1,as.numeric(i)] <- prop.temp[i]
  }
  df.min <- data.frame(method = method.to.compare,
                       Proportion = round(c(prop.win), 2))
  
  
  m.RMSE <- RMSE
  median.methods   <- apply(m.RMSE, 2, median)
  m.RMSE <- m.RMSE[,names(sort(median.methods))]
  
  df.min$y <- max(median.methods)
  p.value <- c()
  compare.name <- c()
  compare.jia <- c()
  compare.yi  <- c()
  order <- c()
  for(i in 1:(ncol(RMSE)-1)){
    for(j in (i+1):ncol(RMSE)){
      p.value <- append(p.value, wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                             paired = TRUE)$p.value)
      if(p.value[length(p.value)] < 0.05){
        direction <- ifelse(wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                        paired = TRUE, conf.int = TRUE)$estimate>0, 1, 0)
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        if(direction == 1){
          compare.jia <- append(compare.jia, name.1)
          compare.yi  <- append(compare.yi, name.2)
          order       <- append(order, 0)
        } else{
          compare.jia <- append(compare.jia, name.2)
          compare.yi  <- append(compare.yi, name.1)
          order       <- append(order, 1)
        }
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      } else{
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        compare.jia <- append(compare.jia, name.1)
        compare.yi  <- append(compare.yi, name.2)
        order       <- append(order, 0)
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      }
    }
  }
  df.test <- data.frame(p.value = p.value,
                        comparison = compare.name,
                        compare.1 = compare.jia,
                        compare.2 = compare.yi,
                        order = order)
  df.pp <- data.frame(x = factor(colnames(m.RMSE),
                                 levels = colnames(m.RMSE)),
                      y = sort(median.methods))
  df.test$p.value.ind <- ifelse(df.test$p.value<0.01, "very_sig",
                                ifelse(df.test$p.value>0.05, "insig_dif", "moderate_sig"))
  df.test$y1 <- median.methods[compare.jia]
  df.test$y2 <- median.methods[compare.yi]
  
  p <- ggplot(df.pp, aes(x=x, y=y))+
    geom_point()
  
  can.color <- c("#F8766D", "#7CAE00","#00BFC4")
  curvature <- c(0.5, 0.5, -0.5)
  names(can.color) <- c("very_sig", "moderate_sig", "insig_dif")
  names(curvature) <- c("very_sig", "moderate_sig", "insig_dif")
  for(i in 1:nrow(df.test)){
    for(sig.level in c("very_sig", "moderate_sig", "insig_dif")){
      if(df.test[i,]$p.value.ind==sig.level){
        if(sig.level != "insig_dif"){
          if(df.test[i,]$order==1){
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = -curvature[sig.level], arrow = arrow(length = unit(0.15, "cm"),
                                                             type = "closed"), show.legend = TRUE)
          } else{
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = curvature[sig.level], arrow = arrow(length = unit(0.15, "cm"),
                                                            type = "closed"), show.legend = TRUE)
          }
        } else{
          p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
          ), col = can.color[sig.level],  data = df.test[i,], 
          curvature = curvature[sig.level], show.legend = TRUE)
        }
      }
    }
  }
  
  p <- p + ylim(min(median.methods)-1.5, max(median.methods)+0.8)+
    theme_bw()+
    ylab("Median RMSE")+
    theme(axis.title.x=element_blank(),
          legend.position = "none")+
    ggtitle(paste0("# Missing variables: ", miss.no))+
    geom_text(aes(x = method, y = y+0.6, label = Proportion),
              size = 3.5, data = df.min)
  return(p)
}
wilcox.sparse <- function(MSE, miss.x, miss.n){
  RMSE <- sqrt(MSE)
  method.to.compare <- colnames(RMSE)
  prop.win <- matrix(0, nrow = 1, ncol = ncol(RMSE))
  colnames(prop.win) <- method.to.compare
  prop.temp <- table(apply(RMSE,1,which.min))/nrow(RMSE)
  for(i in names(prop.temp)){
    prop.win[1,as.numeric(i)] <- prop.temp[i]
  }
  df.min <- data.frame(method = method.to.compare,
                       Proportion = round(c(prop.win), 2))
  
  
  m.RMSE <- RMSE
  median.methods   <- apply(m.RMSE, 2, median)
  m.RMSE <- m.RMSE[,names(sort(median.methods))]
  
  df.min$y <- max(median.methods)
  p.value <- c()
  compare.name <- c()
  compare.jia <- c()
  compare.yi  <- c()
  order <- c()
  for(i in 1:(ncol(RMSE)-1)){
    for(j in (i+1):ncol(RMSE)){
      p.value <- append(p.value, wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                             paired = TRUE)$p.value)
      if(p.value[length(p.value)] < 0.05){
        direction <- ifelse(wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                        paired = TRUE, conf.int = TRUE)$estimate>0, 1, 0)
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        if(direction == 1){
          compare.jia <- append(compare.jia, name.1)
          compare.yi  <- append(compare.yi, name.2)
          order       <- append(order, 0)
        } else{
          compare.jia <- append(compare.jia, name.2)
          compare.yi  <- append(compare.yi, name.1)
          order       <- append(order, 1)
        }
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      } else{
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        compare.jia <- append(compare.jia, name.1)
        compare.yi  <- append(compare.yi, name.2)
        order       <- append(order, 0)
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      }
    }
  }
  df.test <- data.frame(p.value = p.value,
                        comparison = compare.name,
                        compare.1 = compare.jia,
                        compare.2 = compare.yi,
                        order = order)
  df.pp <- data.frame(x = factor(colnames(m.RMSE),
                                 levels = colnames(m.RMSE)),
                      y = sort(median.methods))
  df.test$p.value.ind <- ifelse(df.test$p.value<0.01, "very_sig",
                                ifelse(df.test$p.value>0.05, "insig_dif", "moderate_sig"))
  df.test$y1 <- median.methods[compare.jia]
  df.test$y2 <- median.methods[compare.yi]
  
  p <- ggplot(df.pp, aes(x=x, y=y))+
    geom_point()
  
  can.color <- c("#F8766D", "#7CAE00","#00BFC4")
  curvature <- c(0.5, 0.5, -0.5)
  names(can.color) <- c("very_sig", "moderate_sig", "insig_dif")
  names(curvature) <- c("very_sig", "moderate_sig", "insig_dif")
  for(i in 1:nrow(df.test)){
    for(sig.level in c("very_sig", "moderate_sig", "insig_dif")){
      if(df.test[i,]$p.value.ind==sig.level){
        if(sig.level != "insig_dif"){
          if(df.test[i,]$order==1){
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = -curvature[sig.level], arrow = arrow(length = unit(0.15, "cm"),
                                                             type = "closed"))
          } else{
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = curvature[sig.level], arrow = arrow(length = unit(0.15, "cm"),
                                                            type = "closed"))
          }
        } else{
          p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
          ), col = can.color[sig.level],  data = df.test[i,], 
          curvature = curvature[sig.level])
        }
      }
    }
  }
  
  p <- p + ylim(min(median.methods)-1.5, max(median.methods)+0.8)+
    theme_bw()+
    ylab("Median RMSE")+
    theme(axis.title.x=element_blank(),
          legend.position = "none")+
    ggtitle(paste0("# Miss re. & irre. Var.= \n", miss.x, " and ", miss.n))+
    geom_text(aes(x = method, y = y+0.6, label = Proportion),
              size = 3.5, data = df.min)
  return(p)
}
wilcox.realistic <- function(MSE, topn){
  RMSE <- sqrt(MSE)
  method.to.compare <- colnames(RMSE)
  prop.win <- matrix(0, nrow = 1, ncol = ncol(RMSE))
  colnames(prop.win) <- method.to.compare
  prop.temp <- table(apply(RMSE,1,which.min))/nrow(RMSE)
  for(i in names(prop.temp)){
    prop.win[1,as.numeric(i)] <- prop.temp[i]
  }
  df.min <- data.frame(method = method.to.compare,
                       Proportion = round(c(prop.win), 2))
  
  
  m.RMSE <- RMSE
  median.methods   <- apply(m.RMSE, 2, median)
  m.RMSE <- m.RMSE[,names(sort(median.methods))]
  
  df.min$y <- max(median.methods)
  p.value <- c()
  compare.name <- c()
  compare.jia <- c()
  compare.yi  <- c()
  order <- c()
  for(i in 1:(ncol(RMSE)-1)){
    for(j in (i+1):ncol(RMSE)){
      p.value <- append(p.value, wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                             paired = TRUE)$p.value)
      if(p.value[length(p.value)] < 0.05){
        direction <- ifelse(wilcox.test(m.RMSE[,i], m.RMSE[,j],
                                        paired = TRUE, conf.int = TRUE)$estimate>0, 1, 0)
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        if(direction == 1){
          compare.jia <- append(compare.jia, name.1)
          compare.yi  <- append(compare.yi, name.2)
          order       <- append(order, 0)
        } else{
          compare.jia <- append(compare.jia, name.2)
          compare.yi  <- append(compare.yi, name.1)
          order       <- append(order, 1)
        }
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      } else{
        name.1 <- colnames(m.RMSE)[i]
        name.2 <- colnames(m.RMSE)[j]
        compare.jia <- append(compare.jia, name.1)
        compare.yi  <- append(compare.yi, name.2)
        order       <- append(order, 0)
        compare.name <- append(compare.name, paste0(name.1,"_vs_",name.2) )
      }
    }
  }
  df.test <- data.frame(p.value = p.value,
                        comparison = compare.name,
                        compare.1 = compare.jia,
                        compare.2 = compare.yi,
                        order = order)
  df.pp <- data.frame(x = factor(colnames(m.RMSE),
                                 levels = colnames(m.RMSE)),
                      y = sort(median.methods))
  df.test$p.value.ind <- ifelse(df.test$p.value<0.01, "very_sig",
                                ifelse(df.test$p.value>0.05, "insig_dif", "moderate_sig"))
  df.test$y1 <- median.methods[compare.jia]
  df.test$y2 <- median.methods[compare.yi]
  
  p <- ggplot(df.pp, aes(x=x, y=y))+
    geom_point()
  
  can.color <- c("#F8766D", "#7CAE00","#00BFC4")
  curvature <- c(0.5, 0.5, -0.5)
  names(can.color) <- c("very_sig", "moderate_sig", "insig_dif")
  names(curvature) <- c("very_sig", "moderate_sig", "insig_dif")
  for(i in 1:nrow(df.test)){
    for(sig.level in c("very_sig", "moderate_sig", "insig_dif")){
      if(df.test[i,]$p.value.ind==sig.level){
        if(sig.level != "insig_dif"){
          if(df.test[i,]$order==1){
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = -curvature[sig.level], arrow = arrow(length = unit(0.2, "cm"),
                                                             type = "closed"))
          } else{
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = curvature[sig.level],arrow = arrow(length = unit(0.2, "cm"),
                                                           type = "closed"))
          }
        } else{
            p <- p + geom_curve(aes(x = compare.1, y = y1, xend = compare.2, yend = y2,
            ), col = can.color[sig.level],  data = df.test[i,], 
            curvature = curvature[sig.level])
        }
      }
    }
  }
  
  p <- p + ylim(min(median.methods)-1.5, max(median.methods)+3)+
    theme_bw()+
    ylab("Median RMSE")+
    theme(axis.title.x=element_blank(),
          legend.position = "none")+
    ggtitle(paste0("# Top variables: ", topn))+
    geom_text(aes(x = method, y = y+2, label = Proportion),
              size = 3.5, data = df.min)
  return(p)
}

