



viaEdgeR <- function(longTitle, exptsetup, countFile, countsToAdd, conditionVector, prefix, contrasts, detest_script, categoryFile){
  if ( ! file.exists( countFile)) {
    cat( "\nNot found: ", countFile)
    next
  }
  #rm(list=ls())

  #for yaml
  yamlFileName <- paste(prefix, ".yml", sep="")
  yamlFileName
  metaList <- c(longtitle = longTitle)
  exptfilesList <- c(detestscript_file = detest_script)
  countfilesList <- c(rawcounts_file = countFile)
  
  library(rtracklayer)
  library(edgeR)
  getwd()
  
  cat( "\nReading in raw count file :", countFile )
  data <- read.delim(countFile, header=TRUE,  row.names=1, stringsAsFactors=FALSE)
  print(head(data))
  libNames <- colnames(data)
  
  cat("\n\n Ploting Readcount Distribution plot:\n")
  categoryData <- read.table(categoryFile, header=TRUE)
  categoryNames <- colnames(categoryData)
  categoryNames
  
  tmpData <- mergeTabedFiles(data, categoryData)
  for(libName in libNames){
    print(libName)
 
    tmpData <- tmpData[with(tmpData, order(-get(libName))),]
    head(tmpData)
    for(categoryName in categoryNames){
      print(categoryName)
      
      cat <- tapply(1:nrow(tmpData), factor(tmpData[[categoryName]]), FUN=NULL)
      cat
      mycolors <- c("CDS"="black", "ribosomalp"="green", "rRNA"="red","snoRNA"="cyan", "snRNA"="cyan", "tRNA"="orange", "tubulin"="blue")
      mycolors.cat = mycolors[cat]
      mycolors.cat
      
      mypch <- c("CDS"=".", "ribosomalp"="*", "rRNA"="*","snoRNA"="*", "snRNA"="*", "tRNA"="+", "tubulin"="*")
      mypch.cat = mypch[cat]
      mypch.cat
      
      plotName <- paste("ReadCountDistribution.Raw", libName, "png", sep=".")
      yamlKey = paste("ReadCountDistribution.Raw", libName, sep=".")
      exptfilesList[[yamlKey]] <- plotName #for yaml
      plotTitle <- paste("Read count distribution(Raw Reads)\n", longTitle, "\nLibrary:", libName, "  Category: ", categoryName, sep="")
      png(plotName)
      plot(tmpData[[libName]], log="y", pch=mypch.cat, col=mycolors.cat, las=2, xlab="Genes", ylab="Readcounts in log2",  yaxt='n')
      legend("topright", legend=names(mypch), pch=mypch, col=mycolors)
      title(main=plotTitle, sub="Gene categories are highlighted in color" )
      axis(2, at=(c(1,10,100,1000,10000,100000)), las=2)
      dev.off()
    }
  }

  

  
  cat("\n\n Ploting box plot:\n")
  plotName <- paste(prefix, ".boxplot.RawReadcounts.png", sep="")
  exptfilesList[['boxplot_raw']] <- plotName #for yaml
  plotTitle <- paste("Read count distribution(Raw Reads)\n", longTitle, "\nData Source:", dataSource, sep="")
  png(plotName)
  boxplot(log2(data), use.cols = TRUE, main=plotTitle, ylab="log2 of readcounts", las=2)
  dev.off()
  
  cat( "\nAdding 5 to all reads:")
  data <- data + as.integer(countsToAdd)
  print(head(data))
  
  cat("\n\n Calculating all. vs all correlation:\n")
  cor(log2(data))
  all.vs.all_correlation <- cor(log2(data))
  all.vs.all_correlation.filename <- paste(prefix, ".All.vs.All.correlation.tab", sep="")
  exptfilesList[['allvsall_correlation_file']] <- plotName #for yaml
  write.table(all.vs.all_correlation, all.vs.all_correlation.filename, sep="\t")
  
  cat("\n\n Ploting correlation among samples:\n")
  plotName <- paste(prefix, ".correlation.RawReadcounts.png", sep="")
  exptfilesList[['allvsall_correlation_plot']] <- plotName #for yaml
  plotTitle <- paste("Correlation of genelevel counts(Raw Reads) between samples\n", longTitle, "\nData Source:", dataSource, sep="")
  png(plotName)
  boxplot(log2(data), use.cols = TRUE, main=plotTitle, ylab="log2 of readcounts", las=2)
  pairs(log2(data), pch='.', lower.panel = panel.cor, cex = 2, main=plotTitle)
  dev.off()
  
  #setup the design
  group <- conditionVector
  setDGE <- DGEList(count=data, group=conditionVector)
  cat("\n\nPrinting Samples details before Normalization:\n")
  print(setDGE$samples)
  
  #calculate Normalization factors
  setDGE <- calcNormFactors(setDGE)
  cat("\n\nPrinting Samples details after Normalization:\n")
  print(setDGE$samples)
  
  # Filtering reads
  # at least three of the libraris should have one CPM
  cat("\n\n Dimention of data before filtering\n")
  print(dim(setDGE))
  keep <- rowSums(cpm(setDGE)>0) >= 3
  keptDGE <- setDGE[keep,]
  cat("\n\n Dimention of data after filtering\n")
  print(dim(keptDGE))
  
  # adjust the library size
  cat("\n\nAdjusting the library size after filtering\n")
  keptDGE$samples$lib.size <- colSums(keptDGE$counts)
  print(keptDGE$samples)
  
  # recalculate normalization factors
  cat("\n\nRecalculating normalization factors for the new library size\n")
  keptDGE <- calcNormFactors(keptDGE)
  print(keptDGE$samples)
  cat("\n\nCounts after filtering\n")
  print(head(keptDGE$counts))
  cat("\n\nCPM values after filtering\n")
  print(head(cpm(keptDGE)))
  
  # Data exploration
  plotName <- paste(prefix, ".MDSplot.BCVdistance.png", sep="")
  exptfilesList[['MDSbcv_plot']] <- plotName #for yaml
  plotTitle <- paste("MDS plot (BCVdistance)\n", longTitle, "\nData Source:", dataSource, sep="")
  png(plotName)
  plotMDS(keptDGE, main=plotTitle)
  dev.off()
  
  logCPM <- predFC(keptDGE, prior.count=2*ncol(keptDGE))
  plotName <- paste(prefix, ".MDSplot.logFCdistance.png", sep="")
  exptfilesList[['MDSfc_plot']] <- plotName #for yaml
  plotTitle <- paste("MDS plot (logFCdistance)\n", longTitle, "\nData Source:", dataSource, sep="")
  png(plotName)
  plotMDS(logCPM, main=plotTitle)
  dev.off()
  
  
  # Estimate dispersions
  cat("\n\n Estimating common, tagwise dispersions\n")
  keptDGE <- estimateCommonDisp(keptDGE, verbose=TRUE)
  keptDGE <- estimateTagwiseDisp(keptDGE)
  
  plotName <- paste(prefix, ".BCVplot.png", sep="")
  plotTitle <- paste("BCV plot\n", longTitle, "\nData Source:", dataSource, sep="")
  exptfilesList[['BCV_plot']] <- plotName #for yaml
  png(plotName)
  plotBCV(keptDGE, main=plotTitle)
  dev.off()
  
  # Print adjusted RAW and normalized read counts
  cat("\nPrinting Raw counts from DGE:\n")
  print(head(keptDGE$counts))
  adjraw.count.filename <- paste(prefix, ".adjustedRawReadCounts.tab", sep="")
  countfilesList[['adjustedrawcounts_file']] <- adjraw.count.filename #for yaml
  cat("\nPrinting adjusted Raw counts to a file:",adjraw.count.filename, "\n")
  write.table(keptDGE$counts, adjraw.count.filename, sep="\t")
  cat("\nPrinting Normalized counts from DGE:\n")
  print(head(keptDGE$pseudo.counts))
  norm.count.filename <- paste(prefix, ".normalizedReadCounts.tab", sep="")
  countfilesList[['normalizedcounts_file']] <- norm.count.filename #for yaml
  cat("\nPrinting Normalized counts to a file:",norm.count.filename, "\n")
  write.table(keptDGE$pseudo.counts, norm.count.filename, sep="\t")
  cat("\nPrinting log CPMvalues from DGE:\n")
  print(head(keptDGE$AveLogCPM))
  norm.count.filename <- paste(prefix, ".LogCPMvalues.tab", sep="")
  cat("\nPrinting log CPMvalues to a file:",norm.count.filename, "\n")
  write.table(keptDGE$AveLogCPM, norm.count.filename, sep="\t")
  
  cat("\nPrinting raw CPMvalues from DGE:\n")
  print(head(cpm(keptDGE)))
  norm.count.filename <- paste(prefix, ".CPMvalues.tab", sep="")
  cat("\nPrinting raw CPMvalues to a file:",norm.count.filename, "\n")
  write.table(cpm(keptDGE), norm.count.filename, sep="\t")
  
  cat("\n\n Ploting box plot after normalization:\n")
  plotName <- paste(prefix, ".boxplot.NormReadcounts.png", sep="")
  exptfilesList[['boxplot_normalized']] <- plotName #for yaml
  plotTitle <- paste("Read count distribution (Normalized)\n", longTitle, "\nData Source:", dataSource, sep="")
  png(plotName)
  boxplot(log2(keptDGE$pseudo.counts), use.cols = TRUE, main=plotTitle, ylab="log2 of readcounts", las=2)
  dev.off()
  
  

  
  cat("\n\n#----------------------------------------------#")
  cat("\n# Doing a GLM ananlysis:")
  cat("\n#----------------------------------------------#\n\n")
  

  keptDGE <- estimateGLMCommonDisp(keptDGE, design, verbose=TRUE)
  keptDGE <- estimateGLMTrendedDisp(keptDGE, design)
  keptDGE <- estimateGLMTagwiseDisp(keptDGE, design)
  keptDGE.fit <- glmFit(keptDGE, design)
  
  cat(contrasts)
  contrasts.filename <- paste(prefix, ".all.contrasts.tab", sep="")
  exptfilesList[['contrastset_file']] <- contrasts.filename #for yaml
  cat("\nPrinting all the contrast to a file for reference:", contrasts.filename, "\n")
  write.table(contrasts, contrasts.filename, sep="\t")
  colnames(contrasts)
  
  contrastList <- list() #for yaml
  
  #Do GLMtest for each contrast
  for (contrastName in colnames(contrasts)){
    contrastStr<-contrasts[,contrastName]
    cat("\nCurrent contrast is:", contrastName, "\n")
    print(contrastStr)
    
    #setup lists for yaml
    compfilesList <- list()
    resultfilesList <- list()
    groupName = strsplit(contrastName, "_")
    contrastList[[contrastName]]$contrast_string = paste(contrastStr, sep=" ")
    contrastList[[contrastName]]$basegroup = toupper(groupName[[1]][3])
    contrastList[[contrastName]]$querygroup = toupper(groupName[[1]][1])

    prefix <- paste(contrastName, basePrefix, sep=".")
    longTitle <- paste(baseLongTitle, "\n(", contrastName ,")")
    print(longTitle)
    
    #perform the DE test
    keptDGE.lrt <- glmLRT(keptDGE.fit, contrast=contrastStr)
    topTags(keptDGE.lrt, n=100)
    keptDGE$samples
  
  
    cat("\n\nPrinting summary of GLM analysis:\n")
    de.summary <- summary(de <- decideTestsDGE(keptDGE.lrt))
    print(de.summary)
    detags <- rownames(keptDGE)[as.logical(de)]
    
    
    plotName <- paste(contrastName, ".SMEARplot.raw.png", sep="")
    compfilesList$smearplot_raw = plotName # for yaml
    plotTitle <- paste("Smear plot - GLM Analysis\n", longTitle, "\nData Source:", dataSource, sep="")
    png(plotName)
    plotSmear(keptDGE.lrt, de.tags=detags, main=plotTitle)
    abline(h=c(-1, 1), col="blue")
    dev.off()
    
    plotName <- paste(contrastName, ".SMEARplot.scaled.png", sep="")
    compfilesList$smearplot_scaled = plotName # for yaml
    compfilesList
    plotTitle <- paste("Smear plot - GLM Analysis\n", longTitle, "\nData Source:", dataSource, sep="")
    png(plotName)
    plotSmear(keptDGE.lrt, de.tags=detags, main=plotTitle, ylim=c(-8,8))
    abline(h=c(-1, 1), col="blue")
    dev.off()
    
    cat("\nPrinting Trended Dispersion (GLM) from DGE:\n")
    print(head(keptDGE$trended.dispersion))
    norm.count.filename <- paste(prefix, ".trendedDispersion.tab", sep="")
    cat("\nPrinting Trended Dispersion (GLM) to a file:",norm.count.filename, "\n")
    write.table(keptDGE$trended.dispersion, norm.count.filename, sep="\t")
    
    cat("\nPrinting Tagwise Dispersion (GLM) from DGE:\n")
    print(head(keptDGE$tagwise.dispersion))
    norm.count.filename <- paste(prefix, ".tagwiseDispersion.tab", sep="")
    cat("\nPrinting Tagwise Dispersion (GLM) to a file:",norm.count.filename, "\n")
    write.table(keptDGE$tagwise.dispersion, norm.count.filename, sep="\t")
    
    cat("\nPrinting Master data table from DGE:\n")
    master.data.table <- cbind(keptDGE$counts, keptDGE$pseudo.counts, cpm(keptDGE), keptDGE$tagwise.dispersion, keptDGE$trended.dispersion, keptDGE$AveLogCPM)
    print(head(master.data.table))
    z<- colnames(keptDGE$counts)
    colnames(master.data.table) <- c(paste(z, "rawcount", sep="_"), paste(z, "pseudocount", sep="_"), paste(z, "cpm", sep="_"), "TagwiseDispersion", "TrendedDispersion", "AveLogCPM")
    norm.count.filename <- paste(prefix, ".masterdata.fromDGE.counts.psuedocounts.cpm.tagwisedisp.trenddisp.avelogcpm.tab", sep="")
    cat("\nPrinting Master data table (Counts, PseudoCounts, CPM, TagwiseDisp, TrendedDisp, AveLogCPM) from DGE to a file:",norm.count.filename, "\n")
    write.table(master.data.table, norm.count.filename, sep="\t")
    
    master.data.table.product <- addProductName(master.data.table, productTable)
    norm.count.filename <- paste(prefix, ".masterdata.fromDGE.counts.psuedocounts.cpm.tagwisedisp.trenddisp.avelogcpm.PRODUCT.tab", sep="")
    cat("\nPrinting Master data table (Counts, PseudoCounts, CPM, TagwiseDisp, TrendedDisp, AveLogCPM) from DGE to a file:",norm.count.filename, "\n")
    print(head(master.data.table.product))
    write.table(master.data.table.product, norm.count.filename, sep="\t")
    
    #write all tags
    maxRows <- nrow(keptDGE$counts)
    outfileName <- writeResultTable(keptDGE.lrt, prefix=prefix, suffix="allGenes", maxRows=maxRows)
    resultfilesList$DEresult = outfileName
    
    #write only significant tags
    sig.down.genes <- de.summary[1]
    sig.up.genes <- de.summary[3]
    nonsig.genes <- de.summary[2]
    sig.genes = sig.down.genes + sig.up.genes
    maxRows <- sig.genes
    writeResultTable(keptDGE.lrt, prefix=prefix, suffix="significantGenes", maxRows=maxRows)
    
    
    #manually set maxRows
    if (sigTopCount){
      maxRows <- sigTopCount
      writeResultTable(keptDGE.lrt, prefix=prefix, suffix="manualMaxRows", maxRows=maxRows)
    }
    
    contrastList[[contrastName]]$compfiles = compfilesList
    contrastList[[contrastName]]$resultfiles = resultfilesList
    
  }

# Print YAML file
#yaml file name


write("---", file = yamlFileName)
write("Global:", file = yamlFileName, append = TRUE)
write(paste("    ", 'meta',':', sep=""), file = yamlFileName, append = TRUE)
for(name in names(metaList)){
  print(name)
  print(metaList[[name]])
  write(paste("      ", name,': ', metaList[[name]], sep=""), file = yamlFileName, append = TRUE)
}

write(paste("      ", 'exptsetup',':', sep=""), file = yamlFileName, append = TRUE)
for(name in names(exptsetup)){
  print(name)
  print(exptsetup[[name]]$notes)
  write(paste("        ", name, ':', sep=""), file = yamlFileName, append = TRUE)
  write(paste("          ", 'libs', ': ', exptsetup[[name]]$libs, sep=""), file = yamlFileName, append = TRUE)
  write(paste("          ", 'notes', ': \"', exptsetup[[name]]$notes, "\"",  sep=""), file = yamlFileName, append = TRUE)
}

write(paste("      ", 'countfiles',':', sep=""), file = yamlFileName, append = TRUE)
for(name in names(countfilesList)){
  print(name)
  print(countfilesList[[name]])
  write(paste("          ", name, ': ', countfilesList[[name]], sep=""), file = yamlFileName, append = TRUE)
}
write(paste("      ", 'exptfiles',':', sep=""), file = yamlFileName, append = TRUE)
for(name in names(exptfilesList)){
  print(name)
  print(exptfilesList[[name]])
  write(paste("          ", name, ': ', exptfilesList[[name]], sep=""), file = yamlFileName, append = TRUE)
}

print(contrastList)
#compfilesList$smearplot_scaled = plotName
write(paste('Contrasts', ':', sep=""), file = yamlFileName, append = TRUE)
for(contrastName in names(contrastList)){
  write(paste("    ", contrastName,':', sep=""), file = yamlFileName, append = TRUE)
  for(name in names(contrastList[[contrastName]])){
    if(name == 'compfiles'){
      write(paste("      ", 'compfiles',':', sep=""), file = yamlFileName, append = TRUE)
      for(key in names(contrastList[[contrastName]][[name]])){
         write(paste("        ", key,': ', contrastList[[contrastName]][[name]][[key]], sep=""), file = yamlFileName, append = TRUE)
       }
    } else if (name == 'resultfiles'){
      write(paste("      ", 'resultfiles',':', sep=""), file = yamlFileName, append = TRUE)
      for(key in names(contrastList[[contrastName]][[name]])){
        write(paste("        ", key,': ', contrastList[[contrastName]][[name]][[key]], sep=""), file = yamlFileName, append = TRUE)
      }
    } else {
      write(paste("      ", name,': ', contrastList[[contrastName]][[name]], sep=""), file = yamlFileName, append = TRUE)
    }
  }
}
}

writeResultTable <- function(keptDGE.lrt, prefix="", suffix="", maxRows=""){
  
  # If maxRows not specified print full table
  
  if ( is.null( maxRows) ) {
    maxRows <- 100
  }
  resultTable <- topTags(keptDGE.lrt, n=maxRows, sort.by="logFC")
  resultTable <- addProductName(resultTable, productTable)
  #Up regulated at the top
  fileName <- paste(prefix, ".DEtestResult.", "top",maxRows,"rows.", suffix, ".tab",  sep="")
  write.table(resultTable, fileName)
  
  allResults <- read.table(fileName, header=TRUE)
  cnames <- colnames(allResults)
  cnames <- sub("logFC", "LOG_2_FC", cnames)
  cnames <- sub("PRODUCT", "PRODUCT", cnames)
  cnames <- sub("Row.names", "GENEID", cnames)
  cnames <- sub("PValue", "PVALUE", cnames)
  cnames <- sub("FDR", "FDR", cnames)
  
  colnames(allResults) <- cnames
  #Down regulated at the top
  fileName <- paste(prefix, ".DEtestResult.", "top",maxRows,"rows.", suffix, ".DOWNonTOP.tab", sep="")
  write.table(allResults[order(allResults$LOG_2_FC),], fileName, sep="\t")
  
  return(fileName)
}



addProductName <- function(resultTable, productTable){
  cat("\nReading product name file in\n")
  if ( ! file.exists( productTable)) {
    cat( "\nNot found: ", productTable)
    next
  }
  productNames=read.delim(productTable,  sep="\t", header=TRUE)
  print(head(productNames))
  print(dim(productNames))
  
  #print(head(resultTable))
  print(dim(resultTable))
  
  resultTable <- merge(resultTable, productNames, by=0)
  print(head(resultTable))
  print(dim(resultTable))
  
  
  return(resultTable)
  
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


mergeTabedFiles <- function(table1, table2){
  
  cat("\nPasting/merging two tabed files side by side using first column as key. It should not have header so set as rowname\n")
  
  colsToKeep <- colnames(table1)
  colsToKeep <- c(colsToKeep, colnames(table2))
  tableMerged <- merge(table1, table2, by=0)
  rownames(tableMerged) <- tableMerged$Row.names
  tableMerged <- tableMerged[colsToKeep]
  
  return(tableMerged)
}


storeAll <- function(){
  sessionInfo()
  historyFileName <- paste(prefix, ".Rhistory", sep="")
  rimageFileName <- paste(prefix, ".Rdata", sep="")
  save.history(file=historyFileName)
  save.image(file=rimageFileName)
}
#rm(list=ls())
