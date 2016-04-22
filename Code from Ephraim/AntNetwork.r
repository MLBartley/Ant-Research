##############################
##
## Read in data
##
##   antlist.rg2 = A "list" where each object is an individual ant
##     ID = ant ID
##     type = forager, nest, or Queen
##     cells = vector of cell IDs the ant was in, in temporal order
##     t.in.cell = time spent in the cell
##     time = actual time the ant moved into the cell
##            format: MonthDay.Fraction of day
##     day = day of observation (8 days total)
##
##############################

load("antlist.rg2.Rdata")
str(antlist.rg2)

##############################
##
## Create raster of nest with cell labels
##
##############################

library(raster)
nx=11
ny=7
extent.mat=matrix(c(1-.5,nx+.5,1-.5,ny+.5),nrow=2,byrow=T)
ee=extent(extent.mat)
ee
rast=raster(ee,nrows=ny,ncol=nx,crs="+proj=longlat +datum=WGS84")
values(rast) <- 0
wall.pts=cbind(6,1:6)
values(rast)[cellFromXY(rast,wall.pts)] <- NA
wall.pts=cbind(2:10,4)
values(rast)[cellFromXY(rast,wall.pts)] <- NA
image(rast,main="Nest Raster with Cell Labels")
xy=xyFromCell(rast,1:77)
text(xy[,1],xy[,2])

##############################
##
## Get "contacts" when ants are in the same cell
##
## Columns:
##   i = ID of first ant
##   j = ID of second ant
##   start = time stamp when contact started
##   end = time stamp when contact ended
##   length = end-start
##
##############################

source("get.contacts.r")
C=get.contacts(antlist.rg2)
## remove NA's
na.idx=which(is.na(C[,3]))
c=C[-na.idx,]
c.rg2=data.frame(c)

str(c.rg2)
summary(c.rg2)

## histogram of length of contact
hist(c.rg2[,5],breaks=200)

## table of the absolute number of contacts between each pair
table(c.rg2[,1:2])


###############################
##
## Plotting observed contact networks aggregated over each day
##
###############################

library(network)

days=sort(unique(floor(c.rg2$start)))
days

network.day=list()
for(i in 1:length(days)){
    day=days[i]
    idx.day=which(c.rg2$start>day &c.rg2$start<day+1)
    network.day[[i]]=network(c.rg2[idx.day,],vertex.attr=5,vertex.attrnames="length",directed=FALSE)
}

par(mfrow=c(2,4))
for(i in 1:8){
    plot(network.day[[i]],main=days[i])
}

