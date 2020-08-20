## funcoes necessarias

## exploracao dos dados

## fazer um histograma pra saber o numero de eventIDS por ano
barplot_function <- function (df1,df2) { 
  
  barplot(
    table(as.numeric(substr (unique(df1$eventID), 
                             nchar (as.character(unique(df1$eventID)))-3,
                             nchar(as.character(unique(df1$eventID)))))),
    xlab="Year",col="green",main="Number of eventIDs per year",
    space=1.2,ylim=c(0,60
    ))
  
  barplot(
    table(as.numeric(substr (as.character(unique(df2$eventID)), 
                             nchar (as.character(unique(df2$eventID)))-3,
                             nchar(as.character(unique(df2$eventID)))))),
    col="green4",main="Number of eventIDs per year",add=T,
    space=c(2,1.2,1.2,1.2,1.2),xaxt="n")
  
  legend ("topleft",legend= c("Benthos", "Fishes"),
          pch=15,bty="n",
          col=c("green","green4"),
          cex=1.5)
  
}

## mapa exploratorio dos sitios

initial_map_function <- function (df1, df2) { 
  
  # mapa mundi
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  require(ggplot2)
  library(ggrepel)
  
  # cortar o mapa para ver a america do Sul e parte da central
  wm <- ggplot() + 
    geom_sf (data=world, size = 0.1, 
             fill= "gray90",colour="gray90") +
    coord_sf (xlim = c(-50, -30),  ylim = c(-25, -1), expand = FALSE) +
    theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),axis.ticks.y=element_blank(),
          title = element_text(size=8)) 
  
  jitter <- position_jitter(width = 0.2, height = 0.5)
  
  bentos_coord <- wm + geom_point(data=df1,aes (x= Lon, y=Lat),
                                  stroke=1,shape=1, size=1, 
                                  position = jitter,col="red") 
  
  peixes_coord <- bentos_coord + geom_point(data=df2,aes (x= Lon, y=Lat),
                                            shape=19, size=0.1, 
                                            position = jitter)
  
  peixes_coord ## mostre o mapa
}


### 

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}


