# imputeTSRFsp

imputeTSRFsp is a technique that uses imputeTS and missForest algorithms, and interpolates the data interpolated by 
considering the temporal characteristics of imputeTS by estimating multivariate through RandomForest like the SpatialRandomForest algorithm. 
Based on this interpolated data, it is provided in grid data format like the SpatialRandomForest algorithm.




    devtools::install_github('qkdrk7777775/imputeTSRFsp')
    library(GSIF)
    library(imputeTSRFsp)
    load("C:/Users/cj/Desktop/imputeTSRFsp/df.rda")
    GADM=raster::getData('GADM', country='KOR', level=1)
    srtm1 <- raster::getData ( 'SRTM', lon = 128, lat = 36)
    srtm2 <- raster::getData ( 'SRTM', lon = 128, lat = 34)
    srtm3=raster::merge(srtm1,srtm2)

    #df2=df

    for(year in 2010:2019){
      for(mon in c(paste0('0',1:9),10:12)){
        message(year,mon)
        df=df2
        df=na.omit(df[substr(df$date,1,7)==paste0(year,'-',mon),c('station','date','temp','ws','rh','lat','lon','alt')])
        sp::coordinates(df)=~lon+lat
        sp::proj4string(df) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")
        full_df=df
        output=imputeTSRFsp(full_df, grid.dist,station_unique_name=c('alt','station'))

        df.grid0=expand.grid(seq(min(station@bbox[1,1]),max(station@bbox[1,2]),len=100),
                             seq(min(station@bbox[2,1]),max(station@bbox[2,2]),len=100))
        colnames(df.grid0)=c('lon','lat')

        df_pixel0<- sp::SpatialPixelsDataFrame(points = df.grid0[,c('lon','lat')], data = df.grid0)
        pred.grid=GSIF::buffer.dist(station,df_pixel0,classes=as.factor(1:nrow(station)))

        for(i in colnames(pred.grid@data)){
          pred.grid@data[,i]<-NULL
        }

        pred.grid@data=output$Grid$ws@data

        # library(animation)
        # dir.create("examples")
        # saveGIF({
        #   for(i in names(output$Grid$ws)[1:100]){
        #     message(i)
        #     raster::plot(srtm3,legend=F,col=0,xlim=c(125,131),ylim=c(33.2,38),main=i,cex=2)
        #     raster::plot(raster::raster(pred.grid[i]),main=i,cex.main=2,zlim=c(0,10),
        #                  alpha=.9,add=T,col= rev(cm.colors(255)))
        #   }
        # }, movie.name = paste0(getwd(),"/examples/ws2.gif"))
        for(t in c('temp','rh','ws')){
          pred.grid@data=output$Grid$rh@data
          write.csv(data.frame(pred.grid@coords,pred.grid@data),paste0('C:/Users/cj/Desktop/imputeTSRF/',t,'_',year,'_',mon,'.csv'))
        }
      }
    }
