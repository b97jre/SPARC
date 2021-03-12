


plotExpressionViolinPlot  <- function(expressionInfo2,protein, moleculePalette){ 
  expressionInfo3 = expressionInfo2 %>% select(geneID, sample, RPKM, Cq, time, type)
  expressionInfo4 = expressionInfo3 %>% gather(RPKM, Cq , key = molecule, value = expression)
  expressionInfoGene = expressionInfo4[expressionInfo4$geneID == protein,]
  plot = ggplot(expressionInfoGene[expressionInfoGene$type == "sc",], aes (x = time, y = expression, fill = molecule))+ 
    geom_violin(size = 0.5)+
    scale_fill_manual(values=moleculePalette)+
    geom_point(data =expressionInfoGene[expressionInfoGene$type == "FACS",],aes (x = time, y = expression),shape = 21,colour = "white", fill = "black", size = 1, stroke = 1) +
    facet_grid(.~molecule) +  
    theme(legend.position="none")+ 
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  return(plot)
  
  
}

plotExpressionViolinPlotTitle  <- function(expressionInfo2,protein, moleculePalette){ 
  expressionInfo3 = expressionInfo2 %>% select(geneID, sample, RPKM, Cq, time, type)
  expressionInfo4 = expressionInfo3 %>% gather(RPKM, Cq , key = molecule, value = expression)
  expressionInfoGene = expressionInfo4[expressionInfo4$geneID == protein,]
  plot = ggplot(expressionInfoGene[expressionInfoGene$type == "sc",], aes (x = time, y = expression, fill = molecule))+ 
    geom_violin(size = 0.5)+
    scale_fill_manual(values=moleculePalette)+
    geom_point(data =expressionInfoGene[expressionInfoGene$type == "FACS",],aes (x = time, y = expression),shape = 21,colour = "white", fill = "black", size = 1, stroke = 1) +
    facet_grid(.~molecule) +  
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    ggtitle(protein)+
    theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
  
  
}





scatterDensityPlot = function(expressionInfoSpec, geneName,cbPalette){
  expressionInfoSpec = expressionInfoSpec %>% 
    filter(type == "sc") %>%
    filter(geneID == geneName )
  scatterPlot = ggplot(expressionInfoSpec, 
                       aes(RPKM,Cq, color = time)) + 
    geom_point(size = 0.5) + 
    theme(legend.position = "none")+
    ggtitle(geneName)+ 
    scale_colour_manual(values=cbPalette)+ 
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(size = 7, angle = 0),  
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size=10))
  
  
  scatterPlot2 =  axis_canvas(scatterPlot, axis = "x")+ 
    geom_density(data = expressionInfoSpec, aes(RPKM, color = time, fill = time),alpha = 0.7, size = 0.5,adjust = 3/4) + scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
  
  scatterPlot3 =  axis_canvas(scatterPlot, axis = "y", coord_flip = TRUE)+ 
    geom_density(data = expressionInfoSpec, aes(Cq, color = time,size = 4, fill = time),alpha = 0.7, size = 0.5,adjust = 3/4) +coord_flip()+ scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
  
  p1 <- insert_xaxis_grob(scatterPlot, scatterPlot2, grid::unit(.2, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, scatterPlot3, grid::unit(.2, "null"), position = "right")
  
  return(p2)  
  
}

scatterDensityPlotFigure = function(expressionInfoSpec, geneName,cbPalette){
  expressionInfoSpec = expressionInfoSpec %>% 
    filter(type == "sc") %>%
    filter(geneID == geneName )
  scatterPlot = ggplot(expressionInfoSpec, 
                       aes(RPKM,Cq, color = time)) + 
    geom_point(size = 0.5) + 
    scale_colour_manual(values=cbPalette)+ 
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = .5),
          axis.text.y = element_text(size = 7, angle = 0),  
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size=10))+
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
    
  
  
  scatterPlot2 =  axis_canvas(scatterPlot, axis = "x")+ 
    geom_density(data = expressionInfoSpec, aes(RPKM, color = time, fill = time),alpha = 0.7, size = 0.5,adjust = 3/4) + scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
  
  scatterPlot3 =  axis_canvas(scatterPlot, axis = "y", coord_flip = TRUE)+ 
    geom_density(data = expressionInfoSpec, aes(Cq, color = time,size = 4, fill = time),alpha = 0.7, size = 0.5,adjust = 3/4) +coord_flip()+ scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
  
  p1 <- insert_xaxis_grob(scatterPlot, scatterPlot2, grid::unit(.2, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, scatterPlot3, grid::unit(.2, "null"), position = "right")
  
  return(p2)  
  
}





  

plotPseudotime  <- function(expressionInfo2,protein, timePalette){ 
  
  expressionInfo3 = expressionInfo2 %>% select(geneID, sample, RPKM, Cq, pseudoTime, type, time)
  expressionInfo4 = expressionInfo3 %>% gather(RPKM, Cq , key = molecule, value = expression)
  expressionInfoSC = expressionInfo4[expressionInfo4$type == "sc",]
  expressionInfoGene = expressionInfoSC[expressionInfoSC$geneID == protein,]
  expressionInfoGene = expressionInfoGene[!is.na(expressionInfoGene$molecule),]
  
  plot = ggplot(expressionInfoGene, aes (x = pseudoTime, y = expression, color = time))+ 
    geom_point(size = 0.5)+
    scale_color_manual(values=timePalette)+ 
    geom_smooth(aes(colour = "Linear model"), method = "lm") +
    ggtitle(protein) + facet_grid(rows = vars(molecule), scales = "free_y")+
    theme(legend.position="none")
  
  return(plot)
  
  
}


plotPseudotimeFigure  <- function(expressionInfo2,protein, timePalette){ 
  
  expressionInfo3 = expressionInfo2 %>% select(geneID, sample, RPKM, Cq, pseudoTime, type, time)
  expressionInfo4 = expressionInfo3 %>% gather(RPKM, Cq , key = molecule, value = expression)
  expressionInfoSC = expressionInfo4[expressionInfo4$type == "sc",]
  expressionInfoGene = expressionInfoSC[expressionInfoSC$geneID == protein,]
  expressionInfoGene = expressionInfoGene[!is.na(expressionInfoGene$molecule),]
  
  plot = ggplot(expressionInfoGene, aes (x = pseudoTime, y = expression, color = time))+ 
    geom_point(size = 0.5)+
    scale_color_manual(values=timePalette)+ 
    geom_smooth(aes(colour = "Linear model"), method = "lm") +
    facet_grid(rows = vars(molecule), scales = "free_y")+
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
  theme(legend.position="none")
    
  
  return(plot)
  
  
}

getROCinfo <- function(genie3DF ,TF = "POU5F1", molecule = "RNA", time = "All", FPR = 0.1, regulation = "cisRegulated_POU5F1"){
  ROCit_obj <- rocit( score =genie3DF[["score"]],
                      class=genie3DF[[regulation]])
  
  ROCplotInfo = data.frame(time = time, molecule = molecule, TPR = ROCit_obj$TPR, FPR = ROCit_obj$FPR)
  
  InfoPlot = plot(ROCit_obj)
  ksInfo = ksplot(ROCit_obj)
  
  
  Youden  = data.frame(time = time, molecule = molecule,  
                       TPR = InfoPlot$`optimal Youden Index point`[3] ,
                       FPR = InfoPlot$`optimal Youden Index point`[2],
                       YoudenCutoff = InfoPlot$`optimal Youden Index point`[4],
                       ksCutoff = ksInfo$`KS Cutoff`,
                       AUC = ROCit_obj$AUC,
                       Neg_count = ROCit_obj$neg_count,
                       hardCutOff = ROCit_obj$Cutoff[sum(ROCit_obj$FPR <FPR)],
                       hardTPR = ROCit_obj$TPR[sum(ROCit_obj$FPR < FPR)], 
                       FalsePositive = sum(ROCit_obj$FPR < FPR)
  )
  ROCInfo = list(ROCplotInfo=ROCplotInfo, Youden =  Youden) 
  return (ROCInfo)
  
}


getShuffleAUC <- function(TF_target_scores , TF1 =  "POU5F1" ,
                          molecule1 = "RNA", time1 = "All",
                          n =1000, regulation = "POU5F1_targets" ){
  AUCs = list()
  genie3DF =  TF_target_scores %>% 
    filter(method == method1, 
           molecule == molecule1,
           time == time1,
           TF == TF1)
  genie3DF =  TF_target_scores %>% 
    filter(method == method1, 
           molecule == molecule1,
           time == time1,
           TF == TF1)
  
  
  for(i in 1:n){
    rows <- sample(nrow(genie3DF), replace = FALSE)
    shuffleCis=genie3DF[[regulation]][rows]
    scores = genie3DF[["score"]]

    ROCit_obj <- rocit( score =scores,
                        class=shuffleCis)
    AUCs[i] = ROCit_obj$AUC
  }    
  
  DF = t(data.frame(AUCs))
  permuted = data.frame(time = time1 ,
                        molecule = molecule1,
                        TF = TF1,
                        AUC = DF[,1],
                        type = "permuted"
  )
  
  
  
  return (permuted)
}


plotScatterPlot <- function(ExpressionInfo, TF1, target1 , timePalette = "#009E73"){
  
  TF_target = ExpressionInfo %>% filter(target == target1 & TF == TF1) 
  TDGF1_protein = ggscatter(TF_target, x ="TF_expression", y = "Target_expression", 
                            color =timePalette  ,
                            add = "reg.line",
                            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                            conf.int = TRUE, 
                            cor.coef = TRUE, 
                            cor.coeff.args = list(method = "pearson",label.sep = "\n",
                                                  label.x.npc = "left", label.y.npc = "top"),
                            cor.method = "pearson",
                            xlab = TF1, ylab = paste(target1,"RNA", sep=" "))
  
}



own_residuals <- function(x, y, log=F, p=1, own_weights=F, return_labels=F, labels=c())
{
  filter = is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y) | (x==0); x=x[!filter]; y=y[!filter]; if (own_weights!=F){own_weights=own_weights[!filter]}
  if (log==T){x = log2(x); y=log2(y)}
  if (p==1){if (own_weights==F){fit = lm(y ~ x)} else {fit = lm(y ~ x, weights=own_weights)}}
  if (p==2){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2), weights=own_weights)}}
  if (p==3){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3), weights=own_weights)}}
  if (p==4){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4), weights=own_weights)}}
  if (p==5){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5), weights=own_weights)}}
  if (p==6){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6), weights=own_weights)}}
  if (p==7){if (own_weights==F){fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7))} else {fit = lm(y ~ I(x^0) + I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7), weights=own_weights)}}
  if (!return_labels)
  {
    if (log==T){return(list(2**x, residuals(fit)))}
    else{return(list(x, residuals(fit)))}
  }
  else
  {
    names(x) = labels[!filter]
    z = residuals(fit)
    names(z) = labels[!filter]
    if (log==T){return(list(2**x, z))}
    else{return(list(x, z))}    
  }
}

least_total_squares <- function(x, y)
{
  filter = (!is.na(x)) & (!is.na(y)); x = x[filter]; y = y[filter]
  vector <- prcomp(cbind(x, y))$rotation
  beta <- (-vector[(-ncol(vector)),ncol(vector)] / vector[ncol(vector),ncol(vector)])
  inter <- mean(y) - beta * mean(x)
  return(c(inter, beta))
}



