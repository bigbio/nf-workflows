library(MsTats)
comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)

#---------------

# qualtity control
qcplot <- dataProcessPlots(processed.quant, type="QCplot", 
                 ylimDown=0, 
                 which.Protein = 'allonly',
                 width=7, height=7, address=F)


#---------------

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)

#---------------

test.MSstats.cr <- test.MSstats$ComparisonResult

# Rename spiked ins to A,B,C....
pnames <- c("A", "B", "C", "D", "E", "F")
names(pnames) <- c(
  "sp|P44015|VAC2_YEAST",
  "sp|P55752|ISCB_YEAST",
  "sp|P44374|SFG2_YEAST",
  "sp|P44983|UTR6_YEAST",
  "sp|P44683|PGA4_YEAST",
  "sp|P55249|ZRT4_YEAST"
)

test.MSstats.cr.spikedins <- bind_rows(
  test.MSstats.cr[grep("P44015", test.MSstats.cr$Protein),],
  test.MSstats.cr[grep("P55752", test.MSstats.cr$Protein),],
  test.MSstats.cr[grep("P44374", test.MSstats.cr$Protein),],
  test.MSstats.cr[grep("P44683", test.MSstats.cr$Protein),],
  test.MSstats.cr[grep("P44983", test.MSstats.cr$Protein),],
  test.MSstats.cr[grep("P55249", test.MSstats.cr$Protein),]
)
# Rename Proteins
test.MSstats.cr.spikedins$Protein <- sapply(test.MSstats.cr.spikedins$Protein, function(x) {pnames[as.character(x)]})
test.MSstats.cr$Protein <- sapply(test.MSstats.cr$Protein, function(x) {
  
  x <- as.character(x)
  
  if (x %in% names(pnames)) {
    
    return(pnames[as.character(x)]) 
  } else {
      return("")
  }
})

#---------------

groupComparisonPlots(data=test.MSstats.cr, type="VolcanoPlot",
                     width=12, height=12,dot.size = 2,ylimUp = 7,
                     which.Comparison = "C2-C1",
                     address=F)
