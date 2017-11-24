#SomaVarPlot server
library(shiny)
library(GenomicRanges)
library(rtracklayer)
library(trackViewer)
library(plyr)
library(dplyr)
library(DT)
library(stringr)
library(tidyr)
library(magrittr)
library(ComplexHeatmap)
library(GetoptLong)
library(gridExtra)
library(cowplot)
library(grid)

#Dependencies.  
#sapply(c("rtracklayer", "trackViewer", "plyr", "dplyr", "DT", "stringr", "tidyr", "magrittr","ComplexHeatmap", "GetoptLong", "gridExtra", "cowplot", "grid"), library, character.only = TRUE)

#the Erbb2Features' GRanges object used to draw lolliplots with the Erbb2 domains.  
#instantiating the Erbb2Features GenomicRanges object which holds ErBb2 protein domains. 
Erbb2Features <- GRanges("ERBB2", IRanges(c(1, 652, 675, 720, 1003),  width = c(651, 23, 45, 285, 405), names = c("Extracellular", "TM", "JM", "Kinase", "Tail")))
Erbb2Features$fill <- c("bisque", "coral2", "darkseagreen1","skyblue", "chartreuse2")
Erbb2Features$height <- c(0.02, 0.05, 0.035, 0.04, 0.03)

#setting the Best Overall Response legend for the lolliplots
response <- c("ND", "PD", "SD", "PR", "CR")
response.color.set <- as.list(as.data.frame(rbind(c("gray", "black", "yellow", "green", "red"), "#FFFFFFFF"), stringsAsFactors = FALSE))
names(response.color.set) <- response
###~~~~

#function for the oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "green", col = NA))
  },
  DeepDel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Promoter = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "purple", col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "black", col = NA))
  },    
  MissenseAmp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "cyan", col = NA))
  },
  fusionAmp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "cadetblue1", col = NA))
  },
  splice = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "coral3", col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "darksalmon", col = NA))
  }, 
  rearrangement = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "khaki", col = NA))
  }, 
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.1, "mm"), gp = gpar(fill = "darkgreen", col = NA))
  }
)
col = c("Missense" = "green", "AMP" = "red", "Promoter" = "purple", "DeepDel" = "blue", "Nonsense" = "black", "MissenseAmp" = "cyan", "Indel" = "darkgreen", "splice" = "coral3", "Gain" = "darksalmon", "rearrangement" = "khaki", "fusionAmp" = "cadetblue1")


# Define server logic ----
shinyServer(function(input, output) {
  observe({
  output$result <- renderText({input$tumorSource})
  ts <- input$tumorSource
  output$plotselec <- renderText({input$plotType})
  ps <- input$plotType
  
 if(ps == "Lolliplot") {
   
  tumorSource4lolliplot <- readRDS(paste(ts,".RDS", sep = ""))
#####
  tsDf1 <- data.frame(tumorSource4lolliplot$ErBb2SomaticVarCoord, paste0(tumorSource4lolliplot$ErBb2WildTypeAA, tumorSource4lolliplot$ErBb2SomaticVarCoord), tumorSource4lolliplot$ErBb2SubsAA, tumorSource4lolliplot$Best.Overall.Response)
  #the tsDf1 df Erbb2Features 4 columns, the coordinates where the variations map, the WT aa with the coordinate, the variant AA and the clinical response.
  colnames(tsDf1) <- c("coord", "WTaaCoord", "Var", "Response")
  #creating a 2nd df in order to concatenate the substituions mapping to the same coordinate
  tsDf2 <- as.data.frame(unique(tsDf1[,-4]) %>% group_by(WTaaCoord, coord) %>% do(Varalleles = paste(.$Var, collapse = "/")) %>% ungroup() %>% mutate(Varalleles = unlist(Varalleles)))
  #adding one variable referred to as SNV which reads the WT aa, the coordinate and the substituted aa
  tsDf2$SNV <- paste0(tsDf2$WTaaCoord, tsDf2$Varalleles)
  #adding it to the initial tsDf1 df
  tsDf1 <- merge(tsDf1, tsDf2[,-c(2,3)], by = "WTaaCoord")
  #adding one attribute: concatenation of the clinical response and the variant
  tsDf1$VarRes <- paste0(tsDf1$Response, tsDf1$SNV)
  
  # set order for clinical response value to bottom CR, PR, SD, PD, ND
  tsDf1$Response <- str_replace_all(str_c(tsDf1$Response), c(CR = "E-CR", PR = "D-PR", SD = "C-SD", PD = "B-PD", ND = "A-ND"))
  # order the df
  tsDf1 <- tsDf1[order(tsDf1$coord, tsDf1$Response), ]
  
  ## 
  tsDf3 <- tsDf1
  ## 
  
  tsDf1 <- unique(ddply(tsDf1[,-c(3)], "VarRes", mutate, score = length(VarRes)))
  tsDf1 <- tsDf1[order(tsDf1$coord, tsDf1$VarRes),]
  
  #instantiating the GRanges object with the coordinates and names of somatic variations
  ts.gr <- GRanges("ERBB2", IRanges(unlist(sapply(1:dim(tsDf1)[1], function(i) rep(tsDf1$coord[i], tsDf1$score[i]))), width = 1, names = unlist(sapply(1:dim(tsDf1)[1], function(i) rep(tsDf1$SNV[i], tsDf1$score[i])))))
  
  #adding the stack.factor attribute
  ts.gr$stack.factor <- unlist(sapply(1:dim(tsDf1)[1], function(i) paste0(tsDf1$Response[i], "_", seq(1:tsDf1$score[i]))))
  ts.gr$value1 <- 100
  ts.gr$value2 <- 100 - ts.gr$value1
  
  ts.gr$color <- response.color.set[gsub("_\\d*|\\w*-", "", ts.gr$stack.factor)]
  legend <- list(labels = response, col = "gray80", fill = sapply(response.color.set, `[`, 1))
  
  #editing the gr object's stack.factor variable
  tsDf3 <- ddply(tsDf3, .(WTaaCoord), mutate, idx = seq_along(coord))
  tsDf3 <- tsDf3[order(tsDf3$coord, tsDf3$Response),]
  
  tsDf3$idx <- str_replace_all(str_c(tsDf3$idx), c("10" = "J", "11" = "K", "1" = "A", "2" = "B", "3" = "C", "4" = "D", "5" = "E", "6" = "F", "7" = "G", "8" = "H", "9" = "I"))
  ts.gr$stack.factor <- as.character(tsDf3$idx)
  
  #plotting
  output$plot <- renderPlot({
    lolliplot(ts.gr, Erbb2Features, type = "pie.stack", legend = legend, dashline.col = "gray", cex = .7, ylab = c(ts))
  })
 }
  else if (ps == "Oncoplot"){
    tumorSource4Oncoplot <- readRDS(paste(ts,"Somatics.RDS", sep = ""))
    #Nov3: replacing the cfIMPACT string of the PUMA_ID column by NA
    tumorSource4Oncoplot[tumorSource4Oncoplot == "cfIMPACT"] <- NA
    #Nov3: filling out blank cells.
    tumorSource4Oncoplot  <- tumorSource4Oncoplot  %>% fill(PUMA_ID, Gene)
    tumorSource4Oncoplot$Subject.Identifier.for.the.Study <- gsub("\\d+-\\d+\\-", "", tumorSource4Oncoplot$PUMA_ID)
    #IMPACT data
    tumorSource4OncoplotImpact <- tumorSource4Oncoplot[(tumorSource4Oncoplot$Gene != "not done"),]
    #subsetting
    tumorSource4OncoplotImpact <- tumorSource4OncoplotImpact[,c(7,2,3)]
    #filtering out duplication rows
    tumorSource4OncoplotImpact <- unique(tumorSource4OncoplotImpact)
    
    #filtering it out rows where missing values
    tumorSource4OncoplotImpact <- na.omit(tumorSource4OncoplotImpact)
    
    #assigning variants to their respective type
    tumorSource4OncoplotImpact$Alteration <- sub("[[:alnum:]]+\\*[[:alnum:]]*", "Nonsense;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("[[:alnum:]]+_[[:alnum:]]+dup", "Indel;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("AMP", "AMP;", tumorSource4OncoplotImpact$Alteration, ignore.case = T)
    tumorSource4OncoplotImpact$Alteration <- sub("DeepDel", "DeepDel;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("promoter", "Promoter;", tumorSource4OncoplotImpact$Alteration, ignore.case = T)
    tumorSource4OncoplotImpact$Alteration <- sub("splice", "splice;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("deletion", "Indel;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("gain", "Gain;", tumorSource4OncoplotImpact$Alteration, ignore.case = T)
    tumorSource4OncoplotImpact$Alteration <- sub("[[:alnum:]]+_[[:alnum:]]+ins[[:alnum:]]*", "Indel;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("[[:alnum:]]+_[[:alnum:]]del", "Indel;", tumorSource4OncoplotImpact$Alteration)
    tumorSource4OncoplotImpact$Alteration <- sub("\\w+\\d+\\w+", "Missense;", trimws(tumorSource4OncoplotImpact$Alteration, which = c("both")))
    
    
    #concatenate variations' types whenever occuring in a given gene for a given subject
    tumorSource4OncoplotImpact <- tumorSource4OncoplotImpact %>% group_by(Subject.Identifier.for.the.Study, Gene) %>% do(Alteration = paste(.$Alteration, collapse = '')) %>% ungroup() %>% mutate(Alteration = unlist(Alteration))
    mat <- as.data.frame(spread(tumorSource4OncoplotImpact, Subject.Identifier.for.the.Study, Alteration))
    row.names(mat) <- mat$Gene
    mat <- mat[,-1]
    
    #Non IMPACT data
    tumorSource4OncoplotNonImpact <- tumorSource4Oncoplot[(tumorSource4Oncoplot$Gene == "not done"),]
    tumorSource4OncoplotNonImpact <- tumorSource4OncoplotNonImpact[,c(7,5,6)]
    tumorSource4OncoplotNonImpact <- unique(tumorSource4OncoplotNonImpact)
    tumorSource4OncoplotNonImpact$Subject.Identifier.for.the.Study <- as.character(paste(tumorSource4OncoplotNonImpact$Subject.Identifier.for.the.Study, "*", sep = ""))
    tumorSource4OncoplotNonImpact <- na.omit(tumorSource4OncoplotNonImpact)
    
    #assigning variants to their respective type
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("[[:alnum:]]+\\*[[:alnum:]]*", "Nonsense;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("[[:alnum:]]+_[[:alnum:]]+dup", "Indel;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("AMP", "AMP;", tumorSource4OncoplotNonImpact$Alteration.1, ignore.case = T)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("DeepDel", "DeepDel;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("promoter", "Promoter;", tumorSource4OncoplotNonImpact$Alteration.1, ignore.case = T)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("splice", "splice;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("deletion", "Indel;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("rearrangement", "rearrangement;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("[[:alnum:]]+_[[:alnum:]]+ins[[:alnum:]]*", "Indel;", tumorSource4OncoplotNonImpact$Alteration.1)
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("[[:alnum:]]+_[[:alnum:]]del", "Indel;", tumorSource4OncoplotNonImpact$Alteration.1)
    
    tumorSource4OncoplotNonImpact$Alteration.1 <- sub("\\w+\\d+\\w+", "Missense;", trimws(tumorSource4OncoplotNonImpact$Alteration.1, which = c("both")))
   
     tumorSource4OncoplotNonImpact <- tumorSource4OncoplotNonImpact %>% group_by(Subject.Identifier.for.the.Study, Gene.1) %>% do(Alteration.1 = paste(.$Alteration.1, collapse = '')) %>% ungroup() %>% mutate(Alteration.1 = unlist(Alteration.1))
    
    #creating the mutations matrix for oncoprint
    mat2 <- as.data.frame(spread(tumorSource4OncoplotNonImpact, Subject.Identifier.for.the.Study, Alteration.1))
    mat2 <- subset(mat2, select = colnames(mat2) != "NA*")
    row.names(mat2) <- mat2$Gene.1
    mat2 <- mat2[,-1]
    
    mat <- merge(mat, mat2, by = "row.names", all = T)
    row.names(mat) <- mat$Row.names
    mat <- mat[,-1]
    
    output$plot <- renderPlot({ oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 6), pct_gp = gpar(fontsize = 6), column_title = paste("Oncoplot for ",ts," tumor", sep = ""), heatmap_legend_param = list(title = "Alterations", at = c("Missense", "AMP", "Nonsense", "Indel", "splice", "Gain", "Promoter", "DeepDel", "rearrangement", "MissenseAmp"), labels = c("Missense", "Amplification", "Nonsense", "Indel", "splice", "Gain", "Promoter", "Deep Del", "rearrangement","Missense & Amp")), show_column_names = TRUE)})
    
  }
  else{
    output$plot <- renderPlot({})
  }
#####
  })
})