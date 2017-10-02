############################################
# Original Author: Jason Cochran           #
############################################

stats = FALSE
if ( stats == TRUE) {
  idCols <- 5 + sampleStartCol + numValues
  outputFile <- output
  
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
  library(ggplot2)
  if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel")}
  library("ggrepel")
  if("multcompView" %in% rownames(installed.packages()) == FALSE) {install.packages("multcompView")}
  library("multcompView")
  if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}
  library("plyr")
  
  quantifiedAmounts <- subset.data.frame(quantifiedAmounts, quantifiedAmounts[, sampleStartCol + numValues + 3] == TRUE )
  
  # Preparing the data
  quantifiedAmounts_t <- quantifiedAmounts
  colnames(grouping) <- colnames(quantifiedAmounts_t)
  quantifiedAmounts_t <- rbind( grouping, quantifiedAmounts_t)
  quantifiedAmounts_t <- t(quantifiedAmounts_t)
  
  #Removing the garbage... need to automate this
  graphData <- quantifiedAmounts_t[-c(1:(sampleStartCol - 1), (sampleStartCol + numValues):ncol(quantifiedAmounts_t)) , ]
  
  ####################  Graphing Function  ####################
  
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  graphData_df <- as.data.frame(graphData, stringsAsFactors = F, as.is = T)
  
  # Box plots with mean difference
  for(i in 2:ncol(graphData)) {
    y <- graphData[,i]
    lev <- graphData[,1]
    a <- aov( y~lev , data = graphData_df)
    tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)
    
    generate_label_df <- function(HSD, flev){
      # Extract labels and factor levels from Tukey post-hoc 
      Tukey.levels <- HSD[[flev]][,4]
      Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
      plot.labels <- names(Tukey.labels[['Letters']])
      # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
      # upper quantile and label placement
      boxplot.df <- ddply(graphData_df, flev, function (x) max(fivenum(as.numeric(x[,i]))) + 0.2)
      # Create a data frame out of the factor levels and Tukey's homogenous group letters
      plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']], stringsAsFactors = FALSE)
      # Merge it with the labels
      labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
      return(labels.df)
    }
    p_base <- ggplot(graphData_df, aes(x=lev, y=as.numeric(y) )) + geom_boxplot(aes(fill=lev))
    if( length(unique(lev)) > 2 ) {
      dataTitles <- generate_label_df(tHSD, 'lev')
      dataTitles$V1 <- ggplot_build(p_base)$panel$ranges[[1]]$y.range[[1]]
      p_base <- p_base + geom_text(data = dataTitles, aes(x = plot.labels, y = V1, label = labels))
    }
    p_base <- p_base + theme( plot.title = element_text(size=15, face="plain", margin = margin(10,0,10,0)), panel.background = element_rect(fill = "white"), legend.position = "right", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.ticks = element_blank(),axis.text.y = element_text(size = rel(1.2)), axis.text.x = element_text(size = rel(1.2) ) ) +
      stat_summary(fun.y=mean, geom="point", size=.1) + scale_fill_manual( guide = guide_legend(title = NULL) ,values=cbPalette)
    p_base <- p_base + labs( x = "Sample Group", y = "Concentration")
    ggsave(file = paste0(output, "/graphs", "/test", i, ".png"), plot = p_base, width = 3.5, height = 3.5)
  }
  
  # Master Pie Plot by Class
  #Automote the selection here
  # count <- table(quantifiedAmounts[,25])
  count <- table(quantifiedAmounts[,sampleStartCol+numValues+1])
  # count <- table(quantifiedAmounts[,27])
  count <- as.data.frame(count, header = FALSE, stringsAsFactors = FALSE, as.is = TRUE )
  count <- count[-1,]
  y.breaks <- cumsum(count[,2]) - count[,2]/2
  p <- ggplot(data=count,  aes(x=factor(''), y = count[,2], fill = count[,1], label= count[,1] )) + geom_bar(stat = "identity", width = 1)
  p <- p + coord_polar(theta= "y" ) + theme_minimal() + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), panel.grid=element_blank())
  p <- p + guides(fill = guide_legend(title = "", title.position = "top"))
  ggsave(file = paste0(output, "/graphs", "/","totalIdentified" , ".png"), plot = p, width = 3.5, height = 3.5)
  
  # Master Pie Plot by Amount (average across all samples (C1-C18), and than sum across all rows (all ceramides))
  for(i in 1:numValues) {
    quantifiedAmounts[,sampleStartCol+i-1] <- as.numeric(quantifiedAmounts[,sampleStartCol+i-1])
  }
  averaged <- data.frame( rowMeans( quantifiedAmounts[,c(sampleStartCol:numValues-1)], na.rm = T), stringsAsFactors = F )
  averaged <- cbind( quantifiedAmounts$Class, averaged)
  averaged <- aggregate.data.frame(averaged[,2], list(averaged[,1]), FUN = sum )
  y.breaks <- cumsum(averaged[,2]) - averaged[,2]/2
  p <- ggplot(data=averaged,  aes(x=factor(''), y = averaged[,2], fill = averaged[,1] )) + geom_bar(stat = "identity", width = 1)
  p <- p + coord_polar(theta= "y" ) + theme_minimal() + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), axis.ticks = element_blank() ) #,legend.position="none")
  p <- p + guides(fill = guide_legend(title = "", title.position = "top"))
  # p <- p + scale_y_continuous( breaks=y.breaks, labels= paste0( averaged[,1], " - ", ceiling(averaged[,2]), " ug" ) )
  ggsave(file = paste0(output, "/graphs", "/","totalIdentified_Amount" , ".png"), plot = p, width = 3.5, height = 3.5)
  
  # Master Pie Plot by Amount (average across all samples (C1-C18), and than sum across all rows (all ceramides))
  for(i in 1:length(unique(lev))) {
    averaged <- data.frame( rowMeans(quantifiedAmounts[,c( ((6*i)+1):(6*(1+i)) ) ], na.rm = T), stringsAsFactors = F )
    averaged <- cbind( quantifiedAmounts$Class, averaged)
    averaged <- aggregate.data.frame(averaged[,2], list(averaged[,1]), FUN = sum )
    y.breaks <- cumsum(averaged[,2]) - averaged[,2]/2
    p <- ggplot(data=averaged,  aes(x=factor(''), y = averaged[,2], fill = averaged[,1] )) + geom_bar(stat = "identity", width = 1)
    p <- p + coord_polar(theta= "y" ) +
      theme_minimal() + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), panel.grid=element_blank()) # , legend.position="none")
    # p <- p + scale_y_continuous( breaks=y.breaks, labels= paste0(averaged[,1], " - ", ceiling(averaged[,2])," ug" ) )
    p <- p + guides(fill = guide_legend(title = "", title.position = "top"))
    ggsave(file = paste0(output, "/graphs", "/","totalIdentified_Amount","_", i , ".png"), plot = p, width = 3.5, height = 3.5)
  }
  
  # Color by Class Dot plot
  p <- ggplot(quantifiedAmounts, aes(x = quantifiedAmounts[,2], y = quantifiedAmounts[,3]))
  p <- p + theme_minimal() + theme(legend.position= "none") + geom_point( aes(color = factor(quantifiedAmounts[,sampleStartCol+numValues+1])) )
  ggsave(file = paste0(output, "/graphs/" , "TotalLipidomeClassDistribution.png"), plot = p, width = 3.5, height = 3.5)
  
  #################################
  # Setting up the excel notebook #
  #################################
  if("rJava" %in% rownames(installed.packages()) == FALSE) {install.packages("rJava")}
  if("xlsxjars" %in% rownames(installed.packages()) == FALSE) {install.packages("xlsxjars")}
  if("xlsx" %in% rownames(installed.packages()) == FALSE) {install.packages("xlsx")}
  library(rJava)
  library(xlsxjars)
  library(xlsx)
  
  PostQuantAnalysis.wb <- createWorkbook()
  sheet.0 <- createSheet(PostQuantAnalysis.wb, sheetName = "Lipid Names" )
  sheet.1 <- createSheet(PostQuantAnalysis.wb, sheetName = "By Class (count)" )
  sheet.2 <- createSheet(PostQuantAnalysis.wb, sheetName = "By Sample Group and Class (ug)" )
  sheet.3 <- createSheet(PostQuantAnalysis.wb, sheetName = "By Class (ug)" )
  sheet.4 <- createSheet(PostQuantAnalysis.wb, sheetName = "By Feature Box Plots" )
  sheet.5 <- createSheet(PostQuantAnalysis.wb, sheetName = "All features mz by rt" )
  # sheet.5 <- createSheet(PostQuantAnalysis.wb, sheetName = "ANOVA" )
  
  # Styles
  TITLE_STYLE <- CellStyle(PostQuantAnalysis.wb)+ Font(PostQuantAnalysis.wb,  heightInPoints=16, isBold=TRUE, underline=1)
  SUB_TITLE_STYLE <- CellStyle(PostQuantAnalysis.wb) + Font(PostQuantAnalysis.wb,  heightInPoints=12, isItalic=FALSE, isBold=FALSE)
  
  # Title making function
  xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle, colIndex)	{
    rows <- createRow(sheet,rowIndex=rowIndex)
    sheetTitle <-createCell(rows, colIndex=colIndex)
    setCellValue(sheetTitle[[1,1]], title)
    setCellStyle(sheetTitle[[1,1]], titleStyle) # Could be: TITLE_STYLE or SUB_TITLE_STYLE
  }
  
  # Making multiple titles per row... otherwise they overwrite each other
  xlsx.addTitleMultiple<-function(sheet, rowIndex, titles, titleStyle, colIndexs)	{
    rows <- createRow(sheet,rowIndex=rowIndex)
    for(i in 1:titles) {
      sheetTitle <-createCell(rows, colIndex=colIndex[i] )
      setCellValue(sheetTitle[[1,1]], title)
      setCellStyle(sheetTitle[[1,1]], titleStyle) # Could be: TITLE_STYLE or SUB_TITLE_STYLE
    }
  }
  
  # Put lipid names in first table
  xlsx.addTitle(sheet.0, 1, "Names of All Lipids Included in Report", TITLE_STYLE, 1)
  # quantifiedAmounts_final <- quantifiedAmounts[,-27]
  quantifiedAmounts_final <- quantifiedAmounts[,-(sampleStartCol+numValues+2)]
  addDataFrame(x = quantifiedAmounts_final, sheet = sheet.0, row.names = FALSE, startColumn = 1, startRow = 2) # write data to sheet starting on line 1, column 4
  
  # Master table - shows how many of each lipid class we identified in our sample
  xlsx.addTitle(sheet.1, rowIndex=1, title="Total Distribution of Lipids Identified in Sample", titleStyle = TITLE_STYLE, colIndex=2)
  addPicture(file =paste0(output, "/graphs", "/","totalIdentified" , ".png"), startRow = 2, startColumn = 2, scale = .6, sheet = sheet.1)
  
  # By group and class (ug)
  xlsx.addTitle(sheet.2, rowIndex=1, title="Distribution of Lipids Identified by Class and Sample Grouping", titleStyle = TITLE_STYLE, colIndex=2)
  groups <- LETTERS[seq( from = 1, to = 3 )]
  x <- 1 # Horizontal location
  y <- 2 # Vertical location
  groups <- unique(graphData[,1])
  # for(i in 1:length(unique(lev)) ) {
  for(i in 1:length(LETTERS[seq( from = 1, to = 3 )]) ) {
    if( x == 1 ) { # Create our row for titling
      rows <- createRow(sheet.2, rowIndex = y)
    }
    sheetTitle <-createCell(rows, colIndex = x )
    setCellValue(sheetTitle[[1,1]], groups[i] )
    setCellStyle(sheetTitle[[1,1]], SUB_TITLE_STYLE)
    # xlsx.addTitle(sheet.2, rowIndex = 2, title = groups[i], titleStyle = TITLE_STYLE, colIndex=i)
    # xlsx.addTitleMultiple(sheet.2, rowIndex = y, titles = , SUB_TITLE_STYLE, colIndexs = x)
    addPicture(file = paste0(output, "/graphs", "/","totalIdentified_Amount","_", i , ".png"), startRow = (y + 1), startColumn = x, scale = .4, sheet = sheet.2)
    x <- x + 7
    if( x == 22 ) {
      x <- 1
      y <- y + 23
    }
  }
  
  # Feature (row)
  xlsx.addTitle(sheet.4, rowIndex=1, title="Distribution by Feature (also shows ANOVA by sample group)", titleStyle = TITLE_STYLE, colIndex=2)
  x <- 1 # Horizontal location
  y <- 2 # Vertical location
  for(i in 2:(ncol(graphData)) ) {
    if( x == 1 ) { # Create our row for titling
      rows <- createRow(sheet.4, rowIndex = y)
    }
    sheetTitle <-createCell(rows, colIndex = x )
    setCellValue(sheetTitle[[1,1]], quantifiedAmounts[i,sampleStartCol + numValues + 1] )
    setCellStyle(sheetTitle[[1,1]], SUB_TITLE_STYLE)
    # xlsx.addTitle(sheet.2, rowIndex = 2, title = groups[i], titleStyle = TITLE_STYLE, colIndex=i)
    # xlsx.addTitleMultiple(sheet.2, rowIndex = y, titles = , SUB_TITLE_STYLE, colIndexs = x)
    addPicture(file = paste0(output, "/graphs", "/test", i, ".png"), startRow = (y + 1), startColumn = x, scale = .4, sheet = sheet.4)
    x <- x + 7
    if( x == 22 ) {
      x <- 1
      y <- y + 23
    }
  }
  
  # Master table - amount by group
  xlsx.addTitle(sheet.3, rowIndex=1, title="Total Distribution of Lipids Identified in Sample (ug)", titleStyle = TITLE_STYLE, colIndex=2)
  addPicture(file =paste0(output, "/graphs", "/","totalIdentified_Amount" , ".png"), startRow = 2, startColumn = 2, scale = .6, sheet = sheet.3)
  
  xlsx.addTitle(sheet.5, rowIndex=1, title="All Features Visualized Using RT and m/z", titleStyle = TITLE_STYLE, colIndex=2)
  addPicture(file =paste0(output, "/graphs/TotalLipidomeClassDistribution.png"), startRow = 2, startColumn = 2, scale = .6, sheet = sheet.5)
  
}
####################   Closing up   ####################

saveWorkbook( PostQuantAnalysis.wb, paste( output, "/QuantStatistics", ".xlsx",sep = "") )

# (NMFData, k, ExperimentalDesign, lookup, output)
# NMF(quantifiedAmounts, k, ExpDesign, lookup, output)