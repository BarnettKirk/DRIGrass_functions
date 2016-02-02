library(doBy)
library(car)
PCbiplot <- function(PC, x="PC1", y="PC2", txtcol=NULL, arrowcol=NULL, ital=NULL, colbyfact=NULL, main=NULL, maxlength=NULL, symbol=NULL, colors=NULL) {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- if(is.null(colbyfact)){
    ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  }else {levels <- colbyfact
         data <- cbind(data, levels)
         ggplot(data=data, aes(x=PC1, y=PC2)) + geom_point(aes(colour=levels), size=4, shape=symbol)+theme(legend.title=element_blank())+
           scale_colour_manual(values=colors)}
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  datapc <- datapc[sqrt((datapc$v2^2+datapc$v1^2)) > (maxlength),]
  datapc$varnames <- gsub("[.]", " ",datapc$varnames)
  plot <- if(is.null(ital)){
    plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5.1, vjust=1, color="white") + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color=txtcol)
  }else {plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5.1, vjust=1, color="white", fontface=3) + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color=txtcol, fontface=3)}
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color=arrowcol) + ggtitle(main)+theme(plot.title = element_text(lineheight=.8, face="bold"))
  plot
}


se <- function(x)
{
  sqrt(var(x)/length(x))
}

############################Biomass completing dataset script
complete.biomass <- function(x){
  
  columns <- c("plot","spp","mass")
  
  plotnum <- 60
  
  rows <- plotnum * length(unique(x$spp))
  
  cell <- matrix(ncol=length(columns), nrow=0)
  
  for(g in 1:plotnum){
    for(i in 1:length(unique(x$spp)))
    {
      cell<-rbind(cell, c(g,levels(factor(x$spp))[i],0))
    }
  }
  
  cell <- data.frame(cell)
  names(cell) <- columns
  mergeddata <- cell
  
  mergeddata <- merge(x[,-4], cell, all=TRUE, na.omit=TRUE)
  mergeddata$mass <- as.numeric(mergeddata$mass)
  mergeddata <- with(mergeddata, aggregate(mass, by=list(plot,spp), FUN=sum, na.rm=FALSE))
  mergeddata$type <- "S"
  names(mergeddata) <- names(x)
  match <- x[with(x, which(type == 'W')),]
  for(i in 1:nrow(match)){
    mergeddata$type[mergeddata$plot == match$plot[i] & mergeddata$spp == levels(factor(match$spp[i]))] <- "W"
  }
  mergeddata <- mergeddata[order(mergeddata$plot,mergeddata$spp),]
  
  biomass <<- mergeddata
}

#######################################################################

biomass.harvest.scale.all <- function(x){
  
  ####Biomass cleaning script
  
  x <- biomass
  
  x$mass <- as.numeric(as.character(x$mass)) ##turn mass numeric
  
  x$spp <- recode(x$spp, "c('Paspalum dilitatum','Paspalum notatum','Paspalum spp','Paspalum.dilitatum','Paspalum.spp')= 'Paspalum spp'") #Rename all paspalums species
  
  x <- aggregate(mass ~ plot+type+spp, data=x, FUN=sum) # add multiple paspalum together
  
  x$spp <- factor(x$spp) # change factors
  
  splist<- levels(factor((x$spp[x$spp != 'Dead' & x$spp != 'Unsorted'& x$spp != 'Dead.raked'& x$spp != 'vortis']))) ##create list of species
  
  
  x <- x[with(x, order(plot,spp)),] ## reorder data frame by plot
  
  ## loop to calculate live mass
  if(sum(grepl("vortis",levels(x$spp)))>0){
    for (j in 1:length(unique(x$plot))){
      x$live[x$plot == j] <- sum(x$mass[x$plot==j])-
        x$mass[x$plot ==j & x$spp == 'Unsorted']-
        x$mass[x$plot ==j & x$spp == 'Dead.raked']-
        x$mass[x$plot ==j & x$spp == 'vortis']
    }}else{
      if('Dead.raked' %in% levels(x$spp)){
        for (j in 1:length(unique(x$plot))){
          x$live[x$plot == j] <- sum(x$mass[x$plot==j])-
            x$mass[x$plot ==j & x$spp == 'Unsorted']-
            x$mass[x$plot ==j & x$spp == 'Dead.raked']
        }
      }else{
        for (j in 1:length(unique(x$plot))){
          x$live[x$plot == j] <- sum(x$mass[x$plot==j])-
            x$mass[x$plot ==j & x$spp == 'Unsorted']
        }
      }
    }
  
  
  ## loop to calculate proportion of species in live mass
  for (j in 1:length(unique(x$plot))){
    x$proportion[x$plot == j] <- x$mass[x$plot ==j]/x$live[x$plot ==j]
  }
  
  x$proportion[x$type == "W"|x$spp=='Unsorted'|x$spp=='Dead.raked'|x$spp == 'vortis'] <- 1 #change whole plot forbs to 1
  
  ## loop to determine presence counts of species in plots
  for (j in 1:length(unique(x$plot))){
    for (i in 1:length(unique(x$spp))){
      x$pres[x$plot == j][i] <- ifelse(x$proportion[x$plot == j][i] > 0.05, 1, 0)
    }}
  x$pres <- as.numeric(as.character(x$pres))
  
  ## loop to find abundant species and make a list
  #   splistabund <- vector()
  #   for (i in 1:length(splist)){
  #     if (sum(x$pres[x$spp == splist[i]]) >= 15){
  #       splistabund[i] <-splist[i] 
  #     }}
  #   splistabund <- splistabund[!is.na(splistabund)]  
  #   levels(x$spp) <- c(levels(x$spp), "Other species")
  #   ###################################################################
  #   
  #   x$spp[!x$spp %in% splistabund & x$spp != 'Dead' & x$spp != 'Unsorted'] <- "Other species"
  
  x_agg <- aggregate(mass ~ plot+spp+type, data=x, FUN=sum)
  
  x_agg <- x_agg[with(x_agg, order(plot,spp)),] ## reorder data frame by plot
  
  
  ###############################################################
  ## loop to calculate live mass
  if(sum(grepl("vortis",levels(x$spp)))>0){
    for (j in 1:length(unique(x_agg$plot))){
      x_agg$live[x_agg$plot == j] <- sum(x_agg$mass[x_agg$plot==j])-
        x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] -
        x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Dead.raked'] -
        x_agg$mass[x_agg$plot ==j & x_agg$spp == 'vortis']
    }}else{
      if('Dead.raked' %in% levels(x$spp)){
        for (j in 1:length(unique(x_agg$plot))){
          x_agg$live[x_agg$plot == j] <- sum(x_agg$mass[x_agg$plot==j])-
            x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] -
            x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Dead.raked']
        }
      }else{
        for (j in 1:length(unique(x_agg$plot))){
          x_agg$live[x_agg$plot == j] <- sum(x_agg$mass[x_agg$plot==j])-
            x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted']
        }
      }
    }
  ## loop to calculate proportion of species in live mass
  for (j in 1:length(unique(x_agg$plot))){
    x_agg$proportion[x_agg$plot == j] <- 
      x_agg$mass[x_agg$plot ==j]/x_agg$live[x_agg$plot ==j]
  }
  x_agg$proportion[x_agg$type == "W"|x_agg$spp=='Unsorted'|x_agg$spp=='Dead.raked'|x_agg$spp == 'vortis'] <- 1 #change whole plot forbs to 1
  ###################################################################
  
  ## create scaled values for g/m^2
  if(sum(grepl("vortis",levels(x$spp)))>0){
    for (j in 1:length(unique(x_agg$plot))){
      for (i in 1:length(unique(x_agg$spp[x_agg$plot ==j]))){
        x_agg$scaled[x_agg$plot == j][i] <- 
          (x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] + 
             x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Dead.raked'] +
             x_agg$mass[x_agg$plot ==j & x_agg$spp == 'vortis'] +
             x_agg$live[x_agg$plot == j][i])*
          x_agg$proportion[x_agg$plot == j][i]
      }}}else{
        if('Dead.raked' %in% levels(x$spp)){
          for (j in 1:length(unique(x_agg$plot))){
            for (i in 1:length(unique(x_agg$spp[x_agg$plot ==j]))){
              x_agg$scaled[x_agg$plot == j][i] <- 
                (x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] + 
                   x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Dead.raked'] +
                   x_agg$live[x_agg$plot == j][i])*
                x_agg$proportion[x_agg$plot == j][i]
            }}}else{
              for (j in 1:length(unique(x_agg$plot))){
                for (i in 1:length(unique(x_agg$spp[x_agg$plot ==j]))){
                  x_agg$scaled[x_agg$plot == j][i] <- 
                    (x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] + 
                       x_agg$live[x_agg$plot == j][i])*
                    x_agg$proportion[x_agg$plot == j][i]
                }}}
      }
  x_agg$scaled <- round(x_agg$scaled, digits=4)
  x_agg <- x_agg[x_agg$spp != "Unsorted",]
  x_agg$scaled[x_agg$spp == 'Dead.raked'] <- x_agg$mass[x_agg$spp == 'Dead.raked']
  x_agg <- x_agg[x_agg$spp != 'vortis',]
  if('Dead.raked' %in% levels(x$spp)){
    x_agg$scaled[x_agg$spp == 'Dead'] <- x_agg$scaled[x_agg$spp == 'Dead'] + x_agg$scaled[x_agg$spp == 'Dead.raked'] # make a whole plot expected dead
    x_agg <- x_agg[x_agg$spp != "Dead.raked",]
  }
  x_agg$scaled[x_agg$type == "W"] <- x_agg$mass[x_agg$type == "W"]
  biomass <<- x_agg
}
#######################################################################

biomass.harvest.scale.other <- function(x){
  
  ####Biomass cleaning script
  
  x <- biomass
  
  x$mass <- as.numeric(as.character(x$mass)) ##turn mass numeric
  
  x$spp <- recode(x$spp, "c('Paspalum dilitatum','Paspalum notatum','Paspalum spp')= 'Paspalum spp'") #Rename all paspalums species
  
  x <- aggregate(mass ~ plot+type+spp, data=x, FUN=sum) # add multiple paspalum together
  
  x$spp <- factor(x$spp) # change factors
  
  splist<- levels(factor((x$spp[x$spp != 'Dead' & x$spp != 'Unsorted']))) ##create list of species
  
  
  x <- x[with(x, order(plot,spp)),] ## reorder data frame by plot
  
  ## loop to calculate live mass
  for (j in 1:length(unique(x$plot))){
    x$live[x$plot == j] <- sum(x$mass[x$plot==j])-
      x$mass[x$plot==j & x$spp == 'Dead'] - 
      x$mass[x$plot ==j & x$spp == 'Unsorted']
  }
  
  
  ## loop to calculate proportion of species in live mass
  for (j in 1:length(unique(x$plot))){
    x$proportion[x$plot == j] <- x$mass[x$plot ==j]/x$live[x$plot ==j]
  }
  
  x$proportion[x$type == "W"] <- 1 #change whole plot forbs to 1
  
  ## loop to determine presence counts of species in plots
  for (j in 1:length(unique(x$plot))){
    for (i in 1:length(unique(x$spp))){
      x$pres[x$plot == j][i] <- ifelse(x$proportion[x$plot == j][i] > 0.05, 1, 0)
    }}
  x$pres <- as.numeric(as.character(x$pres))
  
  ## loop to find abundant species and make a list
  splistabund <- vector()
  for (i in 1:length(splist)){
    if (sum(x$pres[x$spp == splist[i]]) >= 15){
      splistabund[i] <-splist[i] 
    }}
  splistabund <- splistabund[!is.na(splistabund)]  
  levels(x$spp) <- c(levels(x$spp), "Other species")
  ###################################################################
  
  x$spp[!x$spp %in% splistabund & x$spp != 'Dead' & x$spp != 'Unsorted'] <- "Other species"
  
  x_agg <- aggregate(mass ~ plot+spp, data=x, FUN=sum)
  
  x_agg <- x_agg[with(x_agg, order(plot,spp)),] ## reorder data frame by plot
  
  
  ###############################################################
  ## loop to calculate live mass
  for (j in 1:length(unique(x_agg$plot))){
    x_agg$live[x_agg$plot == j] <- sum(x_agg$mass[x_agg$plot==j])-
      x_agg$mass[x_agg$plot==j & x_agg$spp == 'Dead'] - 
      x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted']
  }
  
  
  ## loop to calculate proportion of species in live mass
  for (j in 1:length(unique(x_agg$plot))){
    x_agg$proportion[x_agg$plot == j] <- 
      x_agg$mass[x_agg$plot ==j]/x_agg$live[x_agg$plot ==j]
  }
  x_agg$proportion[x_agg$type == "W"] <- 1 #change whole plot forbs to 1
  ###################################################################
  
  ## create scaled values for g/m^2
  
  for (j in 1:length(unique(x_agg$plot))){
    for (i in 1:length(unique(x_agg$spp[x_agg$plot ==j]))){
      x_agg$scaled[x_agg$plot == j][i] <- 
        (x_agg$mass[x_agg$plot ==j & x_agg$spp == 'Unsorted'] + 
           x_agg$live[x_agg$plot == j][i])*
        x_agg$proportion[x_agg$plot == j][i]
    }}
  x_agg$scaled <- round(x_agg$scaled, digits=4)
  x_agg <- x_agg[x_agg$spp != "Unsorted",]
  
  biomass <<- x_agg
}

##################################Add factors to DRIGRass

add.factors <- function(x){
  treatment<-x$plot
  treatment<-recode(treatment, "c(2,10,12,18,24,27,36,39,44,46,53,55)='Ambient'")
  treatment<-recode(treatment, "c(1,8,13,16,26,29,31,33,47,50,52,54)='Pulsed drought'")
  treatment<-recode(treatment, "c(6,7,15,20,21,23,34,38,42,49,51,58)='Ambient(no shelter)'")
  treatment<-recode(treatment, "c(3,17,22,40,45,59)='Seasonal'")
  treatment<-recode(treatment, "c(4,5,11,19,28,30,35,37,41,43,56,60)='Drought'")
  treatment<-recode(treatment, "c(9,14,25,32,48,57)='Increased'")
  x$treatment<-treatment
  x$plot <- as.factor(x$plot)
  
  herb <- x$plot
  herb<-recode(herb, "c(4,6,8,10,12,13,19,20,23,26,27,28,31,35,36,38,41,42,44,50,51,54,55,60)='Added'")
  herb<- recode(herb, "c(1,11,14,15,16,17,18,2,21,22,24,25,29,3,30,32,33,34,37,39,40,43,45,46,47,48,49,5,52,53,56,57,58,59,7,9)='Ambient'")
  x$herb<-herb
  
  side <- x$plot
  side<-recode(side, "c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)='Vineyard'")
  side<-recode(side, "c(31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60)='Fence'")
  
  x$side <- side
  x$side <- factor(x$side)
  x$herb <- factor(x$herb, levels=c("Ambient","Added"))
  x$treatment <- factor(x$treatment, levels = c("Ambient(no shelter)", "Ambient", "Increased","Drought","Pulsed drought", "Seasonal"))
  
  x
}
#Standard Error
se <- function(x)
{
  sqrt(var(x)/length(x))
}

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summaryasdataframe <- function(x) {
  if(length(x) == 1) {
    as.data.frame(x[[1]])
  } else {
    lapply(unlist(x, FALSE), as.data.frame)
  }
}

########Function to recode based on a dataframe########################
recoderFunc <- function(data, oldvalue, newvalue) {
  #from: https://susanejohnston.wordpress.com/2012/10/01/find-and-replace-in-r-part-2-how-to-recode-many-values-simultaneously/
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}
###############Splitting for Cover community
coversplit_herb <- function(x){
  x <- x[x$treatment!='Ambient(no shelter)'&x$treatment!='Increased'&x$treatment!='Seasonal',]
  x$treatment <- factor(x$treatment, levels=c("Ambient", "Drought","Pulsed drought"))
  x$herb <- factor(x$herb, levels=c("Ambient","Added"))
  filtercolumns <- c("plot","treatment","herb","side")
  colskeep <- colSums(x[,!names(x)%in%filtercolumns]) > 0
  colskeep <- colskeep[colskeep == TRUE]
  colskeep <- c(filtercolumns,names(colskeep))
  x <- x[,names(x)%in%colskeep]
  x
}
coversplit_noherb <- function(x){
  x <- x[x$treatment!='Ambient(no shelter)'&x$herb!='Added',]
  x$treatment <- factor(x$treatment, levels=c("Ambient","Increased", "Drought","Pulsed drought","Seasonal"))
  filtercolumns <- c("plot","treatment","herb","side")
  colskeep <- colSums(x[,!names(x)%in%filtercolumns]) > 0
  colskeep <- colskeep[colskeep == TRUE]
  colskeep <- c(filtercolumns,names(colskeep))
  x <- x[,names(x)%in%colskeep]
  x
}