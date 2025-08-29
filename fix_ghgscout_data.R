#######################################################
# Extract linear relationship from Picarro scout files
#
# Gustaf Granath 
#######################################################

# Batch processing all of .dat files in a folder.
# Function is interactive and for each file you define the measurement window (by clicking)
# where you think the 'gas ~ time' relationship is linear. A loess smooth is plotted
# to help you.
# A csv file with data of the defined interval is saved for each gas (co2, ch4)

# First load this function
plot.locator <- function(num.points, col=1){
  xs.list<-c(NULL)
  ys.list<-c(NULL)
  for( I in 1:num.points){
    location<-locator(1,type="n") #, par(pch=16,col="blue"))
    xs.list<-c(xs.list,location$'x')
    ys.list<-c(ys.list,location$'y')
    points(location$'x',location$'y', cex=2, col=col, pch=19)
  }
  locations<-list(x=xs.list)
  return(unlist(locations))
}

library(lubridate) # to easy get minutes

# Now choose folder and FIRST remove junk files (eg data recorded between measurements) 
# It may be easier to use "Set working directory" in RStudio (under Session) then  
# to choose the folder you want to work with.
# If you write the directory path, then note that on windows machine you should use "\" instead...I think..
# The csv files are saved in the same folder.
#fold <- "./" # if getwd() shows the same folder as the files you want to work with
# ALL DATES
fold <- "C:/Picarro/ghg/" # path to the folder where you have the files

filenames = list.files(fold, pattern="*.dat")
filenames

# Run loop over all files
counter=0
save.plots = list()
for (i in 1:length(filenames)) {
  counter = counter+1
  x<-filenames[i]
  hdr <- read.table(paste(fold,x,sep=""), header = TRUE)
  #select columns
  hdr <- hdr[,c("DATE", "TIME", "EPOCH_TIME", "CH4_dry", "CO2_dry", "H2O")]
  # make a column with seconds since start
  hdr$Time_s <- hdr$EPOCH_TIME-hdr$EPOCH_TIME[1]
  hdr$TIME <- sapply(strsplit(as.character(hdr$TIME), "[.]"), "[",1) # remove milliseconds and keep a nicer format
  hdr <- hdr[,-which(colnames(hdr) == "EPOCH_TIME")] # remove column
  
  # Not needed I think but often easier with a separate window
  x11() 
  #windows() # maybe this instead on windowes machines?
  
  #CO2 file
  plot(CO2_dry~Time_s, data = hdr, main=paste("CO2 file: ", unlist(strsplit(filenames[i], "[-]"))[2],
                                              unlist(strsplit(filenames[i], "[-]"))[3], sep=" - ") , cex.main=1, col.main=2)
  
  # add axis with minutes
  mins <- minute(as.POSIXct(paste(hdr$DATE,hdr$TIME,sep=" ")))
  mins.un <- unique(mins)
  axis(1,labels = mins.un, at=c(0, cumsum(rle(mins)$length)[1:(length(mins.un)-1)]), 
       line=3,col="red",col.ticks="red",col.axis="red", xpd=T)
  mtext(side=1, text="min", line = 2, col="red", at=c(0,0))
  
  # add loess smooth
  ls <- loess(CO2_dry~Time_s, data = hdr, span = 0.85)
  lines(predict(ls), col='red', lwd=2)
  #ss = identifyPch(hdr$Time_s, hdr$CO2_dry)
  ss = plot.locator(2, col=2)
  save.plots[[counter]] <- recordPlot()
  legend("top", "CLICK THE MOUSE TO CONTINUE", text.col="deeppink")
  locator(1)
  #keyPressed = readkeygraph("[press any key to continue]")
  hdr.co2 <- hdr[hdr$Time_s>ss[1] & hdr$Time_s<ss[2], -which(colnames(hdr) == "CH4_dry")] # subset range of interest and remove ch4 column
  
  #CH4 file
  counter = counter+1
  #x11()
  #windows() # maybe needed on windowes machines?
  plot(CH4_dry~Time_s, data = hdr, main=paste("CH4 file:",unlist(strsplit(filenames[i], "[-]"))[2],
                                              unlist(strsplit(filenames[i], "[-]"))[3], sep=" - "), cex.main=1, col.main=4)
  
  # add axis with minutes
  mins <- minute(as.POSIXct(paste(hdr$DATE,hdr$TIME,sep=" ")))
  mins.un <- unique(mins)
  axis(1,labels = mins.un, at=c(0, cumsum(rle(mins)$length)[1:(length(mins.un)-1)]), 
       line=3,col="blue",col.ticks="blue",col.axis="blue", xpd=T)
  mtext(side=1, text="min", line = 2, col="blue", at=c(0,0))
  
  # add loess smooth
  ls <- loess(CH4_dry~Time_s, data = hdr, span = 0.85)
  lines(predict(ls), col='blue', lwd=2)
  #ss = identifyPch(hdr$Time_s, hdr$CH4_dry, col=4)
  ss = plot.locator(2, col=4) # blue points
  save.plots[[counter]] <- recordPlot()
  legend("top", "CLICK THE MOUSE TO CONTINUE", text.col="deeppink")
  locator(1)
  dev.off()
  hdr.ch4 <- hdr[hdr$Time_s>ss[1] & hdr$Time_s<ss[2], -which(colnames(hdr) == "CO2_dry")] # subset range of interest and remove co2 column

  # fix filename to save as csv. Puts file in the same folder as the raw data
  nn <- unlist(strsplit(filenames[i], "[.]"))[1]
  nn.co2 <- paste(unlist(strsplit(nn, "[-]"))[2], unlist(strsplit(nn, "[-]"))[3], "CO2", sep="_")
  nn.ch4 <- paste(unlist(strsplit(nn, "[-]"))[2], unlist(strsplit(nn, "[-]"))[3], "CH4", sep="_")
  nn.co2 <- paste(fold,nn.co2, sep="")
  nn.ch4 <- paste(fold,nn.ch4, sep="")
  write.csv(hdr.co2, file=paste(nn.co2,"csv", sep="."), row.names = FALSE)
  write.csv(hdr.ch4, file=paste(nn.ch4,"csv", sep="."), row.names = FALSE)
}

# save plots to pdf
pdf("raw_picarro_plots.pdf")
save.plots
dev.off()
