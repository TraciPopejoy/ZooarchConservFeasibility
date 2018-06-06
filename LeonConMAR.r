###PRIORITIZING MUSSEL BEDS TO CONSERVE BASED ON SIMILARITY TO THE PAST AND FEASIBILITY OF PROTECTION
#code by Traci Popejoy

#libraries
library(maps)
library(maptools)
library(sp)
library(rgdal)
library(rgeos)  
library(xlsx)
library(raster)
library(ggplot2)

######Map the Leon River to Begin######
#Shapefiles
HUC1207River <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LowerBrazos", layer="NHDFlowline")
HUC1207Reser<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/LowerBrazos", layer="NHDWaterbody")
TX<-map("state", "texas")
#bringing in the raster files earlier to use their projection
TXLandCover <- raster("C:/Users/Owner/Documents/GISfile/NLCD2011_LC_Texas/NLCD2011_LC_Texas.tif")
projection<-proj4string(TXLandCover)
#land use is in a different projection; reprojecting the raster to use later; works but loses data
TXLandUseutm<-raster("C:/Users/Owner/Documents/GISfile/TXLandUse/NASS_TX/cdl_30m_r_tx_2015_utm14.tif")
TXproj<-map2SpatialLines(TX, proj4string = CRS(as.character(projection)))

#isolate the Leon River (and associate tributaries) out of the HUC1207 file
leonmap <- HUC1207River[grep('Leon River', HUC1207River$GNIS_NAME),]
BeltonLake <- grep('Belton', HUC1207Reser$GNIS_NAME)
ProctorLake <-grep('Proctor', HUC1207Reser$GNIS_NAME)
LeonLake<-grep('Leon', HUC1207Reser$GNIS_NAME)
reservoirmap<- c(BeltonLake, ProctorLake, LeonLake)
reservoirs<-HUC1207Reser[reservoirmap,]

library(hydroMap)
library(maps)
#def.par<-par()
par(bg="white", mar=c(.5,.5,.5,.5))
mapRange1<-as.vector(bbox(leonmap))
mapRange<-c(mapRange1[1]-.1,mapRange1[3]+.1,mapRange1[2]-.1,mapRange1[4]+.1)
#flowLines <- getFlowLines(mapRange, 3)

#####Load the Mussel Data#####
#Late Holocene; from Popejoy et al (2016) and Randklev (2010)
LeonHData<- read.xlsx("LeonData.xlsx", sheetIndex=1)
#Modern; 2011 survey from Randklev et al (2013)
RandklevData<-read.xlsx("LeonData.xlsx", sheetIndex=2)
modernSpTotal<-apply(RandklevData[,6:22], 2, function(x)sum(x,na.rm=T))
RandklevData<-RandklevData[-53,]

#### Make Spatial Data Frames so can do Spatial Analysis ####
coordinates(RandklevData)<-~Easting+Northing
coordinates(LeonHData)<-~Easting+Northing
proj4string(RandklevData)<-CRS("+proj=utm +zone=14R +datum=WGS84")
proj4string(LeonHData)<-CRS("+proj=utm +zone=14R +datum=WGS84")

Rmussel<-SpatialPoints(RandklevData[1:9,], proj4string=CRS("+proj=utm +zone=14S +datum=WGS84"))
Smussel<-SpatialPoints(RandklevData[10:52,], proj4string=CRS("+proj=utm +zone=14R +datum=WGS84"))
RmusselT<-spTransform(Rmussel, proj4string(TXproj))
SmusselT<-spTransform(Smussel, proj4string(TXproj))
LModernMussel<-spRbind(RmusselT, SmusselT)

proj4string(LeonHData)<-CRS("+proj=utm +zone=14R +datum=WGS84")
LHoloceneMussel<-spTransform(LeonHData, proj4string(TXproj))

#### Build rarefaction curve to test sample adequacy ####
head(LeonHData)
HistComM<-as.data.frame(LeonHData)
rownames(HistComM)<-HistComM[,1]
HistComM<-HistComM[,-c(1:5,23)]
rarecurve(HistComM)
#minimum richness: 5?
LeonHData$richness<-as.vector(specnumber(HistComM))
LeonHDataMOD<-LeonHData[LeonHData$richness>=5,]
HistComMOD<-as.data.frame(LeonHDataMOD)

head(RandklevData)
ModComM<-as.data.frame(RandklevData)
rownames(ModComM)<-ModComM[,1]
ModComM<-ModComM[,-c(1:5,23:27)]
rarecurve(ModComM)
ModComM$LRKM<-RandklevData$LRKM

#####Fuzzy ordination to assess similarity of mussel community#####
library(fso)
library(vegan)
#need data in one dataframe
rownames(HistComMOD)<-HistComMOD[,1]
LeonTData<-rbind(HistComMOD[,c(2,6:22)],ModComM)
#eliminate sites with no mussels
LeonTDataFIN<-LeonTData[rowSums(LeonTData[,-1])!=0,]
#taxa in relative abundance
LeonRAData<-LeonTDataFIN[,-1]/rowSums(LeonTDataFIN[,-1])
LeonRAData$LRKM<-LeonTDataFIN$LRKM

leondis<-vegdist(LeonRAData[,-18], "bray") #distance matrix
leon.fso<-fso(LeonRAData$LRKM, leondis)
plot(leon.fso)
summary(leon.fso)

LRKMtotalstand<-leon.fso$mu*max(leon.fso$data)#un-standardization of ordinated RKM (mu)
LeonRAData$standLRKM<-as.factor(LRKMtotalstand)#making labels for graphs
LeonRAData$LRKMlabel<-as.factor(LeonRAData$LRKM)
LeonRAData$Site.Letter<-as.factor(rownames(LeonRAData))

#make a pretty graph to get the point that ordination shows which beds are most similar in terms of sp. comp.
library(reshape2)
Leongraph<-melt(LeonRAData[,-18])
ordinatedlabels<-LeonRAData[order(LeonRAData$standLRKM),]
unordinatedlabels<-LeonRAData[order(LeonRAData$LRKMlabel),]
greyC<-gray.colors(length(unique(Leongraph$variable)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)

Leongraph$var2<-factor(Leongraph$variable, levels=c("QH","QV","AP","LF","LT","MN","QA","PG","PP","UI","AC","UT","CT","FM","LH","TX","TM"))

Leonplot<-ggplot(data = Leongraph, aes(x = standLRKM, y = value, fill=var2)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyC, name="Species")+
  scale_x_discrete(labels=rownames(ordinatedlabels))+
  coord_flip()+
  labs(x="Ordinated Sites", y="Relative Abundance") +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
Leonplot 

##### Identify Conservation Stream Segments #####
#Find which stream units were surveyed for mussels
LeonMainstem<-spTransform(leonmap, proj4string(TXproj))

snap<-snapPointsToLines(LModernMussel, LeonMainstem)
LeonMainstem$LineID<-row.names(LeonMainstem)
ConLReaches<-LeonMainstem[match(snap$nearest_line_id, LeonMainstem$LineID),]

bufbeds<-gBuffer(LModernMussel, width=300, byid=T) #this makes circles 300m radius around the bed
reaches<-ConLReaches[0,]
for (i in 1:length(bufbeds)){ 
   ex<-bbox(bufbeds[i,])
  test<-crop(LeonMainstem, ex)
  row.names(test)<-as.vector(paste(i,"N", rep(1:length(test))))
  reaches<-rbind(reaches,test)
}
plot(reaches)
reaches$lengthact<-gLength(reaches, byid=T)

summary(reaches$lengthact)

snapD<-snapPointsToLines(LModernMussel, reaches)
reaches$rowid<-row.names(reaches)
ConDReaches<-reaches[match(snapD$nearest_line_id, reaches$rowid),] ####these are the reaches I want to use

summary(ConDReaches$lengthact)

#### START THE REGRESSION TABLE
#rank beds and identify line segments those beds fall on 
#based on similarity to the past; using the ordination graph as my ranking system
sites<-rownames(LeonRAData)
zooranks<-cbind(sites, c(1-leon.fso$mu))
conranks<-data.frame(LineID=row.names(ConLReaches), 
					 siteID=as.character(RandklevData$Site.Letter), 
           reachlength=ConDReaches[match(ConLReaches$LineID, ConDReaches$LineID), 14], 
					 QHtotal=RandklevData$QH, 
					 density=RandklevData$CPUE,
                     zoorank=rep(NA, length(ConLReaches)), 
					 Houses=rep(NA))
conranks$zoorank<- zooranks[match(conranks$siteID, zooranks[,1]),2]
conranks[is.na(conranks$zoorank),6] <- 0
conranks$zoorank<-as.numeric(conranks$zoorank)
conranks$QHtotal<-as.numeric(conranks$QHtotal)
conranks$density[conranks$density=="*"] <- 0
conranks$density<-as.numeric(as.character(conranks$density))

cor.test(conranks$zoorank, conranks$QHtotal, method="spearman")
cor.test(conranks$zoorank,conranks$density, method="spearman")

##### Evaluate Housing Next to Stream Reaches #####
#Using census data to get housing density in areas close to the stream
#census data: http://osd.texas.gov/Data/Decennial/2010/
BellData<- read.csv("LeonCensusData/Bell/Bell/tx00042_sf1.csv")
CorData<- read.csv("LeonCensusData/Coryell/Coryell/tx00042_sf1.csv")
HamData<- read.csv("LeonCensusData/Hamilton/Hamilton/tx00042_sf1.csv")
ComData<- read.csv("LeonCensusData/Comanche/Comanche/tx00042_sf1.csv")
EastData<- read.csv("LeonCensusData/Eastland/Eastland/tx00042_sf1.csv")

#combine it all into one dataframe for a quick for loop below
HousingData<-rbind(BellData,CorData,HamData,ComData,EastData)

#census block shapefiles:http://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2010&layergroup=Blocks
Bellblock <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LeonBlockShapefiles", layer="tl_2010_48027_tabblock10")
Corblock <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LeonBlockShapefiles", layer="tl_2010_48099_tabblock10")
Hamblock <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LeonBlockShapefiles", layer="tl_2010_48193_tabblock10")
Comblock <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LeonBlockShapefiles", layer="tl_2010_48093_tabblock10")
Eastblock <- readOGR(dsn="C:/Users/Owner/Documents/GISfile/LeonBlockShapefiles", layer="tl_2010_48133_tabblock10")

#need to find the blocks that intersect with the river
#have to add planar projection to the data to use the buffer function to ensure you get all the blocks, 
#not those just bisected by the river
Bellblock1<-spTransform(Bellblock, proj4string(LeonMainstem))
leonmapbuff<-gBuffer(LeonMainstem, width=1000, byid=T)

#Bell County
BellInd<-gIntersects(Bellblock1, leonmapbuff, byid=T)
Bellclipped<-apply(BellInd==F, MARGIN=2, FUN=all)
intbellblock<-Bellblock1[which(!Bellclipped),]
#Coryell County
Corblock1<-spTransform(Corblock, proj4string(LeonMainstem))
CorInd<-gIntersects(leonmapbuff, Corblock1, byid=T)
Corclipped<-apply(CorInd==F, MARGIN=1, FUN=all)
intcorblock<-Corblock1[which(!Corclipped),]
#Comanche County
Comblock1<-spTransform(Comblock, proj4string(LeonMainstem))
ComInd<-gIntersects(leonmapbuff, Comblock1, byid=T)
Comclipped<-apply(ComInd==F, MARGIN=1, FUN=all)
intcomblock<-Comblock1[which(!Comclipped),]
#Eastland County
Eastblock1<-spTransform(Eastblock, proj4string(LeonMainstem))
EastInd<-gIntersects(leonmapbuff, Eastblock1, byid=T)
Eastclipped<-apply(EastInd==F, MARGIN=1, FUN=all)
inteastblock<-Eastblock1[which(!Eastclipped),]
#Hamilton County
Hamblock1<-spTransform(Hamblock, proj4string(LeonMainstem))
HamInd<-gIntersects(leonmapbuff, Hamblock1, byid=T)
Hamclipped<-apply(HamInd==F, MARGIN=1, FUN=all)
inthamblock<-Hamblock1[which(!Hamclipped),]

#plot to show progress
plot(LeonMainstem)
plot(inthamblock, add=T, col="blue")
plot(inteastblock, add=T, col="green")
plot(intcorblock, add=T, col="red")
plot(intcomblock, add=T, col="purple")
plot(intbellblock, col="yellow")
plot(reservoirs, add=T)

#get all the block datasets in one file
#first have to give rownames to different Spatial Data Frames
#can't get it all into one dataset because "non-unique polygon IDs" found this answer online
makeUniform<-function(SPDF){
  pref<-substitute(SPDF)  #just putting the file name in front.
  newSPDF<-spChFIDs(SPDF,as.character(paste(pref,rownames(as(SPDF,"data.frame")),sep="_")))
  return(newSPDF)
}

inthamblock<-makeUniform(inthamblock)
intbellblock<-makeUniform(intbellblock)
intcomblock<-makeUniform(intcomblock)
intcorblock<-makeUniform(intcorblock)
inteastblock<-makeUniform(inteastblock)

SuperBlock<-spRbind(inthamblock, intbellblock)
SuperBlock<-spRbind(SuperBlock, intcomblock)
SuperBlock<-spRbind(SuperBlock, intcorblock)
SuperBlock<-spRbind(SuperBlock, inteastblock)
SuperBlock$PolyID<-rep.int(1:length(SuperBlock), 1)

#could eliminate blocks that intersect with reservoirs; focused on conserving river habitat
#NOT doing this because some conservation reaches need these blocks (had problems later)
#ResInd<-gIntersects(reservoirs, SuperBlock, byid=T)
#Resclipped<-apply(ResInd==F, MARGIN=1, FUN=all)
#SuperBlockwoR<-SuperBlock[which(Resclipped),]

#so now I need to use the geoID's in the block polygon dataset to pull the housing numbers out
InterestingBlocks<-data.frame(County=rep(NA, length(SuperBlock)), GeoID=rep(NA, length(SuperBlock)), 
                              Block=rep(NA, length(SuperBlock)), LandArea=rep(NA, length(SuperBlock)), 
                              WaterArea=rep(NA, length(SuperBlock)), TotalArea=rep(NA, length(SuperBlock)), 
                              Housing=rep(NA, length(SuperBlock)),  HpA=rep(NA, length(SuperBlock)), 
                              PolyID=rep(NA, length(SuperBlock)))

InterestingBlocks$County<-SuperBlock$COUNTYFP10
InterestingBlocks$GeoID<-SuperBlock$GEOID10
InterestingBlocks$Block<-SuperBlock$NAME10
InterestingBlocks$LandArea<- as.numeric(paste(SuperBlock$ALAND10))/1000000
InterestingBlocks$WaterArea<-as.numeric(paste(SuperBlock$AWATER10))/1000000
InterestingBlocks$TotalArea<- InterestingBlocks$LandArea + InterestingBlocks$WaterArea
InterestingBlocks$Housing<- HousingData[match(SuperBlock$GEOID10, HousingData$GEOID10),7]
InterestingBlocks$HpA<-InterestingBlocks$Housing/InterestingBlocks$TotalArea
InterestingBlocks$PolyID<-SuperBlock$PolyID

#Only need the blocks next to the reaches
conblock<-makeUniform(ConLReaches)
conblock2<-gBuffer(conblock,width=100, byid=T)
ConInd<-gIntersects(conblock2, SuperBlock, byid=T)
Conclipped<-apply(ConInd==F, MARGIN=1, FUN=all)
SuperBlockCon<-SuperBlock[which(!Conclipped),]
plot(SuperBlockCon)
plot(conblock, add=T, col="blue", lwd=4)
conLineIDs<-apply(ConInd==F, MARGIN = 2, FUN=all)

###so at this point - SuperBlockCon has the blocks next to mussel reaches, 
###InterestingBlocks has the housing information, and conranks has the LineID and mussel data

#need to merge these tables so I can do a regression or a gss regression
#want to merge based on blocks?
#conblock and conraknk has line segments, blocks are polygons

#merge conblock and conrank
namesL<-t(t(apply(ConInd, 2, function(u) paste(names(which(u)), collapse=", " ))))
namesLIST<-strsplit(namesL, ", ")
overlap<-NULL
for(u in 1:length(namesLIST)){
  blT<-data.frame(ConIndLineID=colnames(ConInd)[u], 
             Block=namesLIST[[u]],
             LineID=unlist(strsplit(as.character(colnames(ConInd)[u]), "_"))[2])
  overlap<-rbind(overlap, blT)
}

#adding housing data to overlap
houseoverlap<-data.frame(LineID=overlap$LineID, Blocks=overlap$Block, GeoID=rep(NA), HpA=rep(NA), Houses=rep(NA),LandArea=rep(NA))
houseoverlap$GeoID<-as.matrix(as.data.frame(SuperBlockCon[match(overlap$Block,row.names(SuperBlockCon)),5]))
houseoverlap$HpA<-InterestingBlocks[match(houseoverlap$GeoID, InterestingBlocks$GeoID), 8]
houseoverlap$Houses<-InterestingBlocks[match(houseoverlap$GeoID, InterestingBlocks$GeoID), 7]
houseoverlap$LandArea<-InterestingBlocks[match(houseoverlap$GeoID, InterestingBlocks$GeoID),4]
#get the mean housing for each reach
reachhouseland<-aggregate(houseoverlap[,4:6], list(houseoverlap$LineID),FUN=mean, na.rm=T)
#add this data to conranks
for(i in 1:length(conranks$LineID)){
  conranks$HpA[i]<-reachhouseland[match(conranks$LineID[i], reachhouseland$Group.1),2]
  conranks$Houses[i]<-reachhouseland[match(conranks$LineID[i], reachhouseland$Group.1), 3]
  conranks$LandArea[i]<-reachhouseland[match(conranks$LineID[i], reachhouseland$Group.1), 4]
}

SuperBlockCon$Housing<-InterestingBlocks[match(SuperBlockCon$PolyID, InterestingBlocks$PolyID),7]
library(GISTools)
shades=auto.shading(SuperBlockCon$Housing, cutter = rangeCuts, cols=brewer.pal(5, 'Blues'))
choropleth(SuperBlockCon,SuperBlockCon$Housing, shades)
choro.legend(-256800, 950000,shades)

#####Evaluating Land Use Close to the Conservation Reaches ##### 
blocktemp<-spTransform(SuperBlock, crs(TXLandUseutm))
e<-extent(bbox(blocktemp))
CenTXLU<-crop(TXLandUseutm, e)
CenTXLUsp<-as(CenTXLU, 'SpatialGridDataFrame')

LandUseID<-as.data.frame(CenTXLU@data@attributes)
utmmussel<-spTransform(LModernMussel, crs(CenTXLUsp))

LUlist<-extract(CenTXLU, utmmussel, buffer=1000, factors=T)#L
LUlist
#making a shorter list of landuse types
for(n in 1:length(LUlist)) {
  print(table(LUlist[[n]]))
}
LUIDs<-as.factor(c(1,2,4,21,23,24,27,28,29,36,37,50,57,61,67,74,111,121,122,123,124,
                   131,141,142,143,152,176,190,195,205,211,236,238))
LUdf<-data.frame(landusetypes=LUIDs)
for(i in 1:length(LUlist)) {
  x<-as.vector(LUlist[[i]])
  y<-as.data.frame(table(x))
  w<-match(y[,1], LUdf$landusetypes)
  LUdf[w,i+1]<-y[,2]
}
rownames(LUdf)<-LUdf[,1]
LUdf<-LUdf[,-1]
LUrelab<-LUdf
for(j in 1:52){
  LUrelab[,j]<-100*(LUdf[,j]/colSums(LUdf, na.rm=T)[j])
}
rownames(LUrelab)<-LandUseID[match(rownames(LUrelab),LandUseID$ID),3]
LUrelab<-t(LUrelab)

#should combine into major types
LUgrouped<-LUrelab
LUgrouped[is.na(LUgrouped)] <- 0
grain<-LUgrouped[,3]+LUgrouped[,4]+LUgrouped[,5]+LUgrouped[,6]+LUgrouped[,7]+
  LUgrouped[,8]+LUgrouped[,9]+LUgrouped[,10]+LUgrouped[,11]+LUgrouped[,30]+
  LUgrouped[,32]+LUgrouped[,33]
developed<-LUgrouped[,18]+LUgrouped[,19]+LUgrouped[,20]+LUgrouped[,21]
crops<-LUgrouped[,1]+LUgrouped[,2]+LUgrouped[,12]+LUgrouped[,13]+LUgrouped[,15]+
       LUgrouped[,16] + LUgrouped[,31]
LUfinalG<-data.frame(LUgrouped)
LUfinalG<-LUfinalG[, -c(1:13,15,16,18:21,30:33)]
LUfinalG$grain<-grain
LUfinalG$developed<-developed
LUfinalG$crops<-crops

#adding LUdata to the regression table
for(k in 1:length(colnames(LUfinalG))) {
  conranks[,9+k]<-LUfinalG[,k]
}
colnames(conranks)[10:22]<-colnames(LUfinalG)
rowSums(conranks[,10:22])

LandUseID<-as.data.frame(TXLandUseutm@data@attributes)
LandUseID$LUcolors<-TXLandUseutm@legend@colortable
LandUseID$Class.Names<-as.character(LandUseID$Class.Names)
LandUseID[257,"Class.Names"]<-"grain"
LandUseID[257,"LUcolors"]<-"#D8b56B"
#grain will be #D8B56B

LANDmeltedranks<-melt(conranks[,c(2,10:22)])
LANDmeltedranks$siteF<-factor(unique(LANDmeltedranks$siteID), levels=conranks[order(conranks$zoorank),2])
LANDmeltedranks$color<-LandUseID[match(as.character(LANDmeltedranks$variable), LandUseID$Class.Names), 5]

(Leonplot<-ggplot(data = LANDmeltedranks, aes(x = siteID, y = value, fill = variable)) +
      geom_bar(stat="identity")+coord_flip()) #+ 
      scale_fill_manual(values=as.character(LANDmeltedranks$color), 
                                               name="Land Cover Type"))

#putting coordinates on conranks data.frame
LeonMTaxa$siteID<-row.names(LeonMTaxa)
conrankPoints<-merge(LeonMTaxa, conranks, by="siteID")
conrankPoints<-conrankPoints[, -c(2:18)] #eliminates superflous species counts; not necessary
conrankPoints<-spTransform(conrankPoints, proj4string(CenTXLUsp))

################# graphs ##################
#property lines;
library(ggplot2)
middlehouse<- mean(conranks$HpA)
HousesZ<-ggplot(conranks, aes(y=zoorank, x=log10(HpA)))+ 
  geom_text(label=conranks$siteID, size=4) +scale_y_continuous(limits=c(0,1))+
  scale_x_reverse(c(1,-.5), 
                  breaks=c(1.25, 1,.75, .5, .25, 0, -.25), 
                  labels = c("17.8","10.0","5.6","3.2","1.8","1.0", "0.6"))+
  geom_hline(yintercept = 0.5)+geom_vline(xintercept = log10(middlehouse))+
  labs(x="Houses per Square Kilometer", y="Similarity to the Past") + 
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
HousesZ


#land use; index = forest:other index
conranks$landuseindex<-(conranks$Deciduous.Forest + 
                          conranks$Mixed.Forest + 
                          conranks$Evergreen.Forest + 
                          conranks$Woody.Wetlands) /100 

LandUseZ<-ggplot(conranks, aes(y=zoorank, x=landuseindex)) + 
  geom_text(label=conranks$siteID, size=4)+
  labs(x="Land Use Index", y="Similarity to the Past") +
  scale_y_continuous(limits=c(0,1))+
  geom_hline(yintercept = 0.5)+geom_vline(xintercept = .3)+
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
LandUseZ



#### Making a River Map ####
par(mar=c(0,0,0,0))
LModernMap<-spTransform(LModernMussel, crs(leonmap))
LHoloMap<-spTransform(LeonHDataMOD, crs(leonmap))
plot(flowLines, col="lightgrey")
plot(leonmap, add=T,lwd=3)
plot(reservoirs, bg="black", col="black",add=T)
plot(LModernMap[-c(27:30,40),], add=T, pch=21, bg="lightgrey", cex=1)
plot(LHoloMap, add=T, pch=22, bg="lightgrey", cex=1)
plot(LModernMap[c(27:30,40),], add=T, pch=21, bg="white", cex=1.7)
maps::map.scale(-99,31, relwidth = .17, ratio=F)
points(-98.120278,31.703333,pch=15)
points(-97.49,31.058889, pch=15)
legend(-99.2,31.4,
       legend=c("Mussel Beds","Zooarch. sites","Conservation Reaches", "Leon River","3rd order streams"),
       col=c("black","black","black","black","lightgrey"), lty=c(0,0,0,1,1), 
       pch=c(21,22,21,NA,NA),lwd=1.2)
compassRose(-99,31.65)
TX<-c(-106,-93.5,25.5,37)
map('state',"Texas", fill=T, col="white")
TXLines<- getFlowLines(TX, 6)
map('state','.',col="lightgrey")

plot(TXLines, add=T, col="darkgrey", lwd=2)
map('state',"Oklahoma",col="white",fill=T,add=T, border="lightgrey")
map('state',"Louisiana", col="white", add=T,fill=T, border="lightgrey")
map('state',"Arkansas", col="white", add=T, fill=T, border="lightgrey")
map('state',"New Mexico",col="white", add=T, fill=T,border="lightgrey")
map('state',"Kansas", col="white",add=T,fill=T, border="lightgrey")
map('state',"Missouri",col="white",add=T, fill=T,border="lightgrey")
map('state','Texas')
plot(leonmap, add=T, col="black", lwd=2)



plot(LeonMainstem)
legend(-270000, 925000, legend=c("2011 mussel beds","late Holocene samples","Best Conserv. Opportunities"), col=col, pch=c(19,17,19))
map.scale(grconvertX(0.3,"npc"), grconvertY(0.095, "npc"), metric = T, ratio=F, relwidth = 0.25)
plot(LModernMussel, pch=19, col=col[1], add=TRUE, cex=.7)
plot(LHoloceneMussel, pch=17, col=col[2], add=T, cex=1.3)
plot(LModernMussel[29:30,], pch=19, col=col[3], add=T, cex=1.3)
plot(res, add=T)
compassRose(-240000, 942000, cex=.7)

#### identify reaches by type ####
for(k in 1:nrow(conranks)){
  if(conranks$landuseindex[k]<0.247 & conranks$zoorank[k]<.5){
      CF[k]<-"Type4"} else{
        if(conranks$landuseindex[k]>0.247 & conranks$zoorank[k]>.5){
          CF[k]<-"Type1"} else{
            if(conranks$landuseindex[k]>0.247 & conranks$zoorank[k]<.5){
              CF[k]<-"Type2"}else{CF[k]<-"Type3"}}
          }
        
      }
  
}
conranks$LUtypes<-CF

for(k in 1:nrow(conranks)){
  if(conranks$HpA[k]<3.158 & conranks$zoorank[k]<.5){
    CF[k]<-"Type2"} else{
      if(conranks$HpA[k]>3.158 & conranks$zoorank[k]>.5){
        CF[k]<-"Type4"} else{
          if(conranks$HpA[k]>3.158 & conranks$zoorank[k]<.5){
            CF[k]<-"Type3"}else{CF[k]<-"Type1"}}
    }
  
}

conranks$Htypes<-CF


LModernMap$Htype<-conranks[match(LeonMTaxa$siteID, conranks$siteID),25]
LModernMap$LUtype<-conranks[match(LeonMTaxa$siteID, conranks$siteID),24]

shadesT<-c("#c2e699","#78c679","#31a354","#006837")
par(mfrow=c(1,2))
plot(leonmap, col="grey", cex=250, main="Conservation Ranking based on Housing Abundance")
plot(LModernMap[LModernMap$Htype=="Type1",], col=shadesT[4],add=T, pch=19)
plot(LModernMap[LModernMap$Htype=="Type2",], col=shadesT[3],add=T, pch=19)
plot(LModernMap[LModernMap$Htype=="Type3",], col=shadesT[2],add=T, pch=19)
plot(LModernMap[LModernMap$Htype=="Type4",], col=shadesT[1],add=T, pch=19)
plot(leonmap, col="grey", cex=250, main="Conservation Ranking based on Forest Abundance")
plot(LModernMap[LModernMap$LUtype=="Type1",], col=shadesT[4],add=T, pch=19)
plot(LModernMap[LModernMap$LUtype=="Type2",], col=shadesT[3],add=T, pch=19)
plot(LModernMap[LModernMap$LUtype=="Type3",], col=shadesT[2],add=T, pch=19)
plot(LModernMap[LModernMap$LUtype=="Type4",], col=shadesT[1],add=T, pch=19)
maps::map.scale(-99,31, relwidth = .5, ratio=F)
legend(-99.3,31.4, 
       legend=c("Type 1","Type 2","Type 3", "Type 4"),
       col=shadesT[4:1], 
       pch=19)

##### write important tables out #####
write.csv(conranks, "Mar11Conblocks.csv")

summary(conranks$Houses)
summary(conranks$HpA)
summary(conranks$)