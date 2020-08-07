#' Spatio-temporal interpolation method
#'
#' This function is a package created based on the function of 'imputeTS', 'ranger' packages algorithm
#' After creating all the time data for each point by using the imputeTS package, create the randomForest
#' model except when the target variable is missing. Then, based on the data interpolated by b KalmanFilter method,
#' it is used as test data of the RandomForest model to finally produce the data.
#' Based on the produced data, grid data is generated using the Spatial RandomForest model.
#'
#' @param 'full_df' is a missing material with location and time information.
#' @export
#' @examples
#'
#' output=imputeTSRFsp(full_df, grid.dist,station_unique_name=c('alt','station'))
#'
#' df.grid0=expand.grid(seq(min(station@bbox[1,1]),max(station@bbox[1,2]),len=100),
#' seq(min(station@bbox[2,1]),max(station@bbox[2,2]),len=100))
#' colnames(df.grid0)=c('lon','lat')
#'
#' df_pixel0<- sp::SpatialPixelsDataFrame(points = df.grid0[,c('lon','lat')], data = df.grid0)
#' pred.grid=GSIF::buffer.dist(station,df_pixel0,classes=as.factor(1:nrow(station)))
#'
#' for(i in colnames(pred.grid@data)){
#'   pred.grid@data[,i]<-NULL
#'   }
#'
#'   pred.grid@data=output$Grid$rh@data
#'
#'   GADM=raster::getData('GADM', country='KOR', level=1)
#'   srtm1 <- raster::getData ( 'SRTM', lon = 128, lat = 36)
#'   srtm2 <- raster::getData ( 'SRTM', lon = 128, lat = 34)
#'   srtm3=raster::merge(srtm1,srtm2)
#'
#'   dir.create("examples")
#'   animation::saveGIF({
#'       for(i in names(output$Grid$ws)[1:100]){
#'           message(i)
#'           raster::plot(srtm3,legend=F,col=0,xlim=c(125,131),ylim=c(33.2,38),main=i,cex=2)
#'           raster::plot(raster::raster(pred.grid[i]),main=i,cex.main=2,zlim=c(4,100),
#'           alpha=.9,add=T,col= rev(cm.colors(255)))
#'         }
#'    }, movie.name = paste0(getwd(),"/examples/rh.gif"), interval = 0.1)
#'
#'
#'
imputeTSRFsp=function(full_df, grid.dist,station_unique_name=c('alt','station'),date_col='date'){
  full_df_temp=list();full_df_temp2=list()
  for(stat in unique(full_df$station)){

    stat=as.character(stat)
    message(stat)
    unique_name=colnames(full_df@coords)
    unique_name=c(unique_name,station_unique_name)
    full_df_temp[[stat]]=data.frame(full_df[full_df$station==stat,])
    na_df=data.frame(date=seq(min(full_df_temp[[stat]]$date),max(full_df_temp[[stat]]$date),3600))

    for(i in setdiff(1:ncol(full_df),which(colnames(full_df_temp[[stat]])%in%date_col))){
      na_df[,colnames(full_df_temp[[stat]])[i]]=NA
    }

    for(i in unique_name){
      na_df[,i]=full_df_temp[[stat]][1,i]
    }
    full_df_temp[[stat]]=dplyr::bind_rows(full_df_temp[[stat]],na_df)

    #full_df_temp[['full_data']][[stat]]
    numeric_col=names(which(sapply(full_df_temp[[stat]],is.numeric)))
    full_df_temp2[[stat]]=full_df_temp[[stat]]
    for(i in numeric_col){
      #message(paste0(i,' impute using Kalman'))
      idx=is.na(full_df_temp[[stat]][,i])
      if(sum(idx)!=0){
        tryCatch({
          full_df_temp2[[stat]][,i]<-imputeTS::na_kalman(full_df_temp[[stat]][,i])
        }
        ,error=function(x){
          full_df_temp2[[stat]][idx,i]<-mean(full_df_temp[[stat]][-idx,i],na.rm=T)
        })
      }
    }
  }
  raw_df2 =dplyr::bind_rows(full_df_temp)
  full_df_temp2=dplyr::bind_rows(full_df_temp2)
  full_df2=full_df_temp2
  pred.grid=grid.dist
  for(i in colnames(pred.grid@data)){
    pred.grid@data[,i]<-NULL
  }

  #Spatial RandomForest model fit
  output=list()
  target_names=names(which(sort(apply(is.na(raw_df2),2,sum),decreasing = T)!=0))
  target_names=target_names[target_names%in%numeric_col]
  for(target in target_names){
    train=full_df2[!is.na(raw_df2[,target]),]
    test=full_df2[is.na(raw_df2[,target]),]
    #temp_df2$station=NULL
    sp::coordinates(train)=~lon+lat
    sp::coordinates(test )=~lon+lat
    from_crs = sp::CRS("+proj=tmerc +lat_0=38 +lon_0=127 +k=1 +x_0=200000 +y_0=500000 +ellps=bessel +units=m")
    sp::proj4string(train) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")
    sp::proj4string(test ) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")
    train <- sp::spTransform(train, from_crs)
    test  <- sp::spTransform(test, from_crs)
    dn0 <- paste(names(grid.dist), collapse="+")

    fm0 <- as.formula(paste(target," ~ ", dn0,'+',
                            paste0(setdiff(names(which(sapply(train@data,is.numeric))),target),collapse = '+')))
    ov.zinc <- sp::over(train, grid.dist)
    rm.zinc=cbind(train@data,ov.zinc)

    date_list=as.character(unique(test@data$date))
    for(date_temp in date_list){
      message(paste0(target,' ', date_temp))
      train2=rm.zinc[rm.zinc$date==date_temp,]
      test2=test[test@data$date==date_temp,]
      m.zinc <- ranger::ranger(fm0, train2, quantreg=TRUE, num.trees=150, seed=1)
      full_df2[intersect(which(is.na(raw_df2[,target])),which(full_df2$date==date_temp)),target]=
        predict(m.zinc,cbind(test2@data,sp::over(test2,grid.dist)))$prediction
    }

    full_df3=full_df2
    sp::coordinates(full_df3)=~lon+lat
    from_crs = sp::CRS("+proj=tmerc +lat_0=38 +lon_0=127 +k=1 +x_0=200000 +y_0=500000 +ellps=bessel +units=m")
    sp::proj4string(full_df3) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")
    full_df3 <- sp::spTransform(full_df3, from_crs)
    ov.zinc3 <- sp::over(full_df3, grid.dist)
    rm.zinc3=cbind(full_df3@data,ov.zinc3)
    fm1 <- as.formula(paste(target," ~ ", dn0))

    pred_list=list()
    date_list2=as.character(unique(full_df2$date))
    for(date_temp in date_list2){
      message(paste0(target,' ', date_temp))
      train3=rm.zinc3[rm.zinc3$date==date_temp,]
      m.zinc2 <- ranger::ranger(fm1, train3, quantreg=TRUE, num.trees=150, seed=1)
      pred_list[[date_temp]]=predict(m.zinc2,grid.dist@data)$prediction
    }
    pred.grid@data=dplyr::bind_cols(pred_list)
    output[[target]]=pred.grid
  }

  return(list(Points=full_df2,Grid=output,kalmanPoints=full_df_temp2))
}

