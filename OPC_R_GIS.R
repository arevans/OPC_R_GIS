### PBDB Workshop
### Phenotypic evolution
### Practical 1

### I. GIS in R

### 1. 

InputG <- read.table("Felis_silvestris-1.780-uw_lrm1bemnewor150b.dat", header=FALSE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
RasterG <- rasterFromXYZ(InputG)
projection(RasterG) <- "+proj=eqc +lat_ts=0"

### 2. 

plot(x = RasterG, col = grey(0:100/100))
plot(x = RasterG, col = terrain.colors(255))

### 3.

RasterS <- terrain(RasterG, opt='slope')
RasterA <- terrain(RasterG, opt='aspect')
RasterH <- hillShade(RasterS, RasterA, 35, 335)
plot(RasterH, col=grey(0:100/100), legend=FALSE)
plot(RasterG, col=rev(terrain.colors(255, alpha=0.35)), add=TRUE)

### II. Dental Complexity

### 1. Define aspect reclassify matrix

Define_rcl <- function (NoCat) {
  rcl <- matrix(nrow = NoCat+1, ncol = 3)
  CatSize <- 360 / NoCat
  rcl[,1] <- c(0,seq(from = CatSize/2, to = 360-(CatSize/2), by = CatSize))
  rcl[,2] <- c(seq(from = CatSize/2, to = 360-(CatSize/2), by = CatSize), 360)
  rcl[,3] <- c(1:NoCat,1)
  result <- rcl
}

rcl <- Define_rcl(8)
rcl

### 2. 

RasterAD <- terrain(RasterG, opt='aspect', unit='degrees')

### 3. 

RasterADC <- reclassify(RasterAD,rcl)

aspect.colors <- colorRampPalette(c("yellow","lightblue","blue","lightgreen","darkgreen","purple", "red", "brown"), space = "rgb")
plot(RasterAspectDegClass, col=aspect.colors(255))

### 4.

RasterC1 <- clump(RasterADC==1,directions=4,gaps=FALSE)
freq(RasterC1)

### 5. Clump all aspect categories
ClumpAspect <- function(RasterADC, NoCat) {
  result <- 0
  for (i in 1:NoCat) {
    RasterClump <- clump(RasterADC==i,directions=4,gaps=FALSE)
    
    if(i==1) { RasterClumpStack <- RasterClump }
    else { RasterClumpStack <- stack(RasterClumpStack,RasterClump) }
  }
  names(RasterClumpStack) <- c(paste(rep("clumps",NoCat),1:NoCat, sep = ""))
  result <- RasterClumpStack
}

NoCat <- 8
RasterC <- ClumpAspect(RasterADC, NoCat)
plot(RasterC)

### 6. Compile clump data into dataframe
CompileClumpData <- function(ClumpsOutput) {
  result <- 0
  ClumpsFreq <- freq(ClumpsOutput)
  NoCat <- length(ClumpsFreq)
  
  temp1 <- as.data.frame(ClumpsFreq$clumps1)
  temp1$Cat <- rep(1, length(temp1$value))
  temp2 <- as.data.frame(ClumpsFreq$clumps2)
  temp2$Cat <- rep(2, length(temp2$value))
  temp3 <- as.data.frame(ClumpsFreq$clumps3)
  temp3$Cat <- rep(3, length(temp3$value))
  temp4 <- as.data.frame(ClumpsFreq$clumps4)
  temp4$Cat <- rep(4, length(temp4$value))
  temp5 <- as.data.frame(ClumpsFreq$clumps5)
  temp5$Cat <- rep(5, length(temp5$value))
  temp6 <- as.data.frame(ClumpsFreq$clumps6)
  temp6$Cat <- rep(6, length(temp6$value))
  temp7 <- as.data.frame(ClumpsFreq$clumps7)
  temp7$Cat <- rep(7, length(temp7$value))
  temp8 <- as.data.frame(ClumpsFreq$clumps8)
  temp8$Cat <- rep(8, length(temp8$value))
  
  tempa <- merge(temp1,temp2, all = TRUE)
  tempa <- merge(tempa,temp3, all = TRUE)
  tempa <- merge(tempa,temp4, all = TRUE)
  tempa <- merge(tempa,temp5, all = TRUE)
  tempa <- merge(tempa,temp6, all = TRUE)
  tempa <- merge(tempa,temp7, all = TRUE)
  tempa <- merge(tempa,temp8, all = TRUE)
  
  result <- tempa
}

ClumpResults <- CompileClumpData(RasterC)

### 7. 

MinClumpSize <- 2 #Only count clumps larger than MinClumpSize pixels in size
OPCClumps <- subset(ClumpResults, ClumpResults[,2] > MinClumpSize)
OPCClumps <- na.omit(OPCClumps)
OPC <- nrow(OPCClumps)

OPCStats <- t(as.matrix(summary(OPCClumps$count)))
OPCStats <- as.data.frame(OPCStats)
OPCStats$SD <- sd(OPCClumps[,2])                          
OPCStats <- rbind(c(OPC,OPCStats))
colnames(OPCStats)[1] <- "OPC"

### 8.

OPCCalc <- function(RasterG, NoCat=8, MinClumpSize=2) {
  RasterAD <- terrain(RasterG, opt='aspect',unit='degrees')
  rcl <- Define_rcl(NoCat)
  RasterADC <- reclassify(RasterAD,rcl)
  RasterC <- ClumpAspect(RasterADC,NoCat)
  ClumpResults <- CompileClumpData(RasterC)

  OPCClumps <- subset(ClumpResults, ClumpResults[,2] > MinClumpSize)
  OPCClumps <- na.omit(OPCClumps)
  OPC <- nrow(OPCClumps)
  OPCData <- as.data.frame(OPC)
  
  OPCStats <- t(as.matrix(summary(OPCClumps$count)))
  OPCStats <- as.data.frame(OPCStats)
  OPCStats$SD <- sd(OPCClumps[,2])                          
  OPCStats <- merge(OPCData,OPCStats)
  return(OPCStats)
}

OPCCalc(RasterG)
OPCCalc(RasterG,MinClumpSize=4)

RunOPCDir <- function(NoCat=8,MinClumpSize=2) {
  Filenames <- list.files(pattern = "\\.dat$")
  NoFiles <- length(Filenames)
  OPCResultsFiles <- data.frame(Filenames, "OPC"=NA, "Min."=NA ,"1st Qu."=NA,"Median"=NA,"Mean"=NA,"3rd Qu."=NA,"Max."=NA,"SD"=NA)
  result <- 0
  for (i in 1:4) { #NoFiles
    InputGrid <- read.table(Filenames[i], header=FALSE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
    RasterGrid <- rasterFromXYZ(InputGrid)
    projection(RasterGrid) <- "+proj=eqc +lat_ts=0"
    
    OPCStats <- OPCCalc(RasterGrid)
    
    OPCResultsFiles[i,2:9] <- OPCStats
  }
  
  result <- OPCResultsFiles  
}

OPCResults <- RunOPCDir()
OPCResults


### III. OPCR

Define_rclR <- function (NoCat, NoRot, CurrRot) {
  rcl <- matrix(nrow = NoCat+1, ncol = 3)
  CatSize <- 360 / NoCat
  OffsetSize <- CatSize / NoRot
  Offset <- (CurrRot - 1) * OffsetSize
  
  if (CurrRot <= NoRot/2 + 1) {
    rcl[,1] <- c(0,seq(from = (CatSize/2) + Offset, to = 360-(CatSize/2) + Offset, by = CatSize))
    rcl[,2] <- c(seq(from = (CatSize/2) + Offset, to = 360-(CatSize/2) + Offset, by = CatSize), 360)
    rcl[,3] <- c(1:NoCat,1) 
  }
  else {
    rcl[,1] <- c(seq(from = Offset - (CatSize/2), to = 360, by = CatSize),0)
    rcl[,2] <- c(seq(from = (CatSize/2) + Offset, to = 360, by = CatSize), 360, Offset - CatSize/2)
    rcl[,3] <- c(1:NoCat,8)
  }
  result <- rcl
}

OPCRCalc <- function(RasterGrid, NoCat=8, MinClumpSize=2, NoRot=8, Plot=FALSE) {
  RasterAspectDeg <- terrain(RasterGrid, opt='aspect',unit='degrees')

  Rots <- c(1:NoRot)
  OPCResultsRot <- data.frame(Rots, "OPC"=NA, "Min."=NA ,"1st Qu."=NA,"Median"=NA,"Mean"=NA,"3rd Qu."=NA,"Max."=NA,"SD"=NA)
  
  CatSize <- 360 / NoCat
  
  for (i in 1:NoRot) {
    rclR <- Define_rclR(NoCat,NoRot,i)
    RasterAspectDegClass <- reclassify(RasterAspectDeg,rclR)
    if (Plot==TRUE) { plot(RasterAspectDegClass) }
    RasterClumped <- ClumpAspect(RasterAspectDegClass,NoCat)
  
    ClumpResults <- CompileClumpData(RasterClumped)
  
    OPCClumps <- subset(ClumpResults, ClumpResults[,2] > MinClumpSize)
    OPCClumps <- na.omit(OPCClumps)
    OPC <- nrow(OPCClumps)
    OPCData <- as.data.frame(OPC)
    
    OPCStats <- t(as.matrix(summary(OPCClumps$count)))
    OPCStats <- as.data.frame(OPCStats)
    OPCStats$SD <- sd(OPCClumps[,2])                          
    OPCStats <- merge(OPCData,OPCStats)

    OPCResultsRot[i,1] <- i
    OPCResultsRot[i,2:9] <- OPCStats
    
  }
  OPCRStats <- apply(OPCResultsRot, 2, mean)
  OPCRStats <- OPCRStats[2:9]
  return(OPCRStats)
}

OPCRCalc(RasterGrid)


RunOPCRDir <- function(NoCat=8, MinClumpSize=2, NoRot=8, MaxFiles=0, Plot=FALSE) {
  Filenames <- list.files(pattern = "\\.dat$")
  if (MaxFiles == 0) {
    NoFiles <- length(Filenames)
  }
  else {
    NoFiles <- MaxFiles
  }
  OPCRResultsFiles <- data.frame(Filenames, "OPCR"=NA, "Min."=NA ,"1st Qu."=NA,"Median"=NA,"Mean"=NA,"3rd Qu."=NA,"Max."=NA,"SD"=NA)
  result <- 0
  for (i in 1:NoFiles) {
    InputGrid <- read.table(Filenames[i], header=FALSE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
    RasterGrid <- rasterFromXYZ(InputGrid)
    projection(RasterGrid) <- "+proj=eqc +lat_ts=0"
    
    OPCRStats <- OPCRCalc(RasterGrid, Plot=Plot)
    
    OPCRResultsFiles[i,2:9] <- OPCRStats
  }
  
  result <- OPCRResultsFiles  
}

OPCRResults <- RunOPCRDir()
OPCRResults

