


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
    theme(legend.position="none")+ 
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


