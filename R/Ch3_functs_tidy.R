#====================
# Chapter 3 functions
#====================
# A place to keep long functions and objects for all chapter 3 R scripts


#--------------------------------------
# plot theme
plot_theme <- theme(legend.position="none",
                    panel.background =  element_rect(fill = NA),
                    panel.grid.major = element_line(),
                    panel.grid.minor = element_line(),
                    panel.border = element_blank(),
                    axis.line = element_line(colour = "black",size = 1,lineend = "butt"),
                    axis.title = element_text(colour = "black", size = 11),
                    axis.text = element_text(colour = "black"))

#-------------------------------
# Plot a histogram with quantiles 

disthistwquants<- function(distance, data, xlabname, title, binsize) {
  # calculate a distance matrix
  dist<- distance(data)
  # Calculate some quantiles & max and min vals
  mean <- quantile(dist, probs=c(0.5), na.rm = T)
  quants <- quantile(dist, probs=c(0.05,0.95),na.rm = T)
  max <- quantile(dist, probs=c(1),na.rm = T)
  #convert dist matrix into a df for plotting 
  distdf<- as.data.frame(dist[1:length(dist)]) %>%
    rename(., distance=1)
  # Plot in a histogram
  ggplot(distdf,aes(x=distance))+ geom_histogram(bins = binsize) +
    geom_vline(xintercept=mean, color="blue", linetype="solid", size=1) +
    geom_vline(xintercept=quants, color="blue", linetype="dashed", size=0.7) +
    geom_vline(xintercept=max, color="red", linetype="dotted", size=0.7) + 
    ylab("Count") + xlab(xlabname) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(title)+
    plot_theme
}

# seperate one for bruvos cause I can't work out how to optionally remove the replen arguenent
disthistwquantsbruvo<- function(distance, data, xlabname, title, binsize) {
  # calculate a distance matrix
  dist<- distance(data, replen=sareplen)
  # Calculate some quantiles & max and min vals
  mean <- quantile(dist, probs=c(0.5))
  quants <- quantile(dist, probs=c(0.05,0.95))
  max <- quantile(dist, probs=c(1))
  #convert dist matrix into a df for plotting 
  distdf<- as.data.frame(dist[1:length(dist)]) %>%
    rename(., distance=1)
  # Plot in a histogram
  ggplot(distdf,aes(x=distance))+ geom_histogram(bins = binsize) +
    geom_vline(xintercept=mean, color="blue", linetype="solid", size=1) +
    geom_vline(xintercept=quants, color="blue", linetype="dashed", size=0.7) +
    geom_vline(xintercept=max, color="red", linetype="dotted", size=0.7) + 
    ylab("Count") + xlab(xlabname) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(title)+
    plot_theme
}


#--------------------------------------------------------------------------------------
# Plot nj tree from a distance metric

plotnjtree<- function(distmatname){
            x<- nj(distmatname) #make the object
            treename <-paste(distmatname, "_njtree")
             x$tip.label<-IDS$FullNames[match(x$tip.label,IDS$IDs)]
             assign(treename, x, envir=.GlobalEnv)
}

#---------------------------------------------------------------
# Function to to plot a tree with resistance on it 
Restreeplotter <- function(tree, layout){
  traittree<- ggtree(tree,layout =  layout) %<+% trait_data
  traittree + geom_tippoint(aes(colour=Res2,shape=Res2), size=2)+ 
    theme(legend.position = "NULL",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
}


#-------------------------------------------------------------------
# extract a subsets of 50 samples from our genind object
extract50subset<- function(genind) {
  set.seed(100)
  fiftysamps <- sample(nInd(genind), 50)
  genind[fiftysamps] %>% 
    popsub(., exclude = c("Unique", "nolabel"))
}

extract50subset2<- function(genind) {
  set.seed(34)
  fiftysamps <- sample(nInd(genind), 50)
  genind[fiftysamps] %>% 
    popsub(., exclude = c("Unique", "nolabel"))
}

extract50subset3<- function(genind) {
  set.seed(617)
  fiftysamps <- sample(nInd(genind), 50)
  genind[fiftysamps] %>% 
    popsub(., exclude = c("Unique", "nolabel"))
}

#------------------------------------------------------------------------------------------------
# FUNCTION TO MAKE AN ID2MLG DF FROM A GENCLONE OBJ

ID2MLG<- function(genclone.obj) {
  mlg.list<-mlg.id(genclone.obj)
  genclone.obj<-as.data.frame(do.call(rbind, mlg.list))
  genclone.obj$MLG<-rownames(genclone.obj)
  genclone.obj$MLG<- paste("MLG", genclone.obj$MLG, sep=".")
  genclone.obj<- pivot_longer(genclone.obj, 1:(length(genclone.obj)-1),values_to = "sampleIDs")
  genclone.obj$name<-NULL
  genclone.obj<-unique(genclone.obj)
}

# Values to "IDs" not "SampleIDs"
ID2MLG2<- function(genclone.obj) {
  mlg.list<-poppr::mlg.id(genclone.obj)
  genclone.obj<-as.data.frame(do.call(rbind, mlg.list))
  genclone.obj<-as.data.frame(do.call(rbind, mlg.list))
  genclone.obj$MLG<-rownames(genclone.obj)
  genclone.obj$MLG<- paste("MLG", genclone.obj$MLG, sep=".")
  genclone.obj<-tidyr::pivot_longer(genclone.obj, 1:(length(genclone.obj)-1),values_to = "IDs")
  genclone.obj$name<-NULL
  genclone.obj<-unique(genclone.obj)
}

#-------------------------------------------------------------------------------------------------------
# Calculate the Fullnames dataframe from a genind object

CalcFullNames <- function(genind) {
  # create a full names list to cross reference against
  mll.custom(genind) <- paste(genind@strata$Country, 
                              genind@strata$Year,
                              genind@strata$Location, 
                              genind@strata$Res,
                              genind@strata$CloneRaw, sep="_") 
  mll(ukwogcc) <- "custom" # set the mll slot to be this custom value
  # Calculate a list and merge to a dataframe
  mlglist <- mlg.id(genind)
  mlgdf<- data.frame(lapply(mlglist, `length<-`, max(lengths(mlglist)))) # pads different lengths with NAs
  FullNames<- mlgdf %>%
    pivot_longer(., 1:length(.),names_to = "FullNames", values_to = "IDs") %>% drop_na()
  # Add the information from strata
  FullNames$Country<-str_split_fixed(FullNames$FullNames, "_", n=6)[,1]
  FullNames$Year<- str_split_fixed(FullNames$FullNames, "_", n=6)[,2]
  FullNames$Location <- str_split_fixed(FullNames$FullNames, "_", n=6)[,3]
  FullNames$Res <-str_split_fixed(FullNames$FullNames, "_", n=6)[,4]
  FullNames$CloneRaw <-str_split_fixed(FullNames$FullNames, "_", n=6)[,5]
  FullNames$TraitMatName<-paste(FullNames$Country,
                                FullNames$Year,
                                FullNames$Location,
                                FullNames$CloneRaw,
                                paste("ID", FullNames$IDs, sep=""),
                                sep="_")
  return(FullNames)
}

#------------------------------------------------------------------------------
# FUNCTION TO CALCULATE A PROPORTIONAL COMMUNITY MATRIX FROM MLG TABLE

CalcComMat <-function(df)  {
  setPop(df) <- ~Location # set the location
  commtab <- mlg.table(df) # get a community table
  temp<-commtab/rowSums(commtab) #make into props
  temp
  }

#==============================================================================
# SESMPD FUNCTIONS 

#--------------------------------------------------------------------------------
# main function creates and output objects needed for sesmpd from a genind object
#-------------------------------------------------------------------------------

SesmpdObjectsLat<- function(pop){
  #1) set the strata to location
  setPop(pop)<-~Location
  #2) create ID table for cross referencing
  IDs <- ID2MLG(pop)
  # 3) extract a community matrix 
  ComMat<-as.matrix(mlg.table(pop))
  ComMat[ComMat > 0] = 1
  # 4) Calculate a dsitance matrix
  MolDist <- bruvo.dist(pop, replen = sareplen)
  bruvomat<- as.matrix(MolDist)
  #5) Rename the distance matrix to match the comm matrix
  rownames(bruvomat)<- IDs$MLG[match(rownames(bruvomat),IDs$sampleIDs)]
  colnames(bruvomat)<- IDs$MLG[match(colnames(bruvomat),IDs$sampleIDs)]
  # 6) Ensure both data sets are ordered the same
  ordering <-colnames(ComMat)
  bruvomat<- bruvomat[ordering,ordering]
  #bruvomat<- bruvomat[,match(colnames(ComMat), colnames(bruvomat))]
  bruvomat<-as.matrix(bruvomat)
  # 7) Calculate ses.mpd
  sesmpdstats<- ses.mpd(ComMat, bruvomat, runs=9999, null.model = "taxa.labels")
  # Add columns of identifiers to aid subsetting
  setPop(pop)<-~Year
  sesmpdstats$Year <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  setPop(pop)<-~Lat
  sesmpdstats$Lat <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  setPop(pop)<-~Long
  sesmpdstats$Long <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  #return the results
  temp<-list(sesmpdstats=sesmpdstats, IDs=IDs, ComMat=ComMat,
             MolDist=MolDist, bruvomat=bruvomat, ordering=ordering)
  return(temp)
}

SesmpdObjectsYear<- function(pop){
  #1) set the strata to location
  setPop(pop)<-~Year
  #2) create ID table for cross referencing
  IDs <- ID2MLG(pop)
  # 3) extract a community matrix 
  ComMat<-as.matrix(mlg.table(pop))
  ComMat[ComMat > 0] = 1
  # 4) Calculate a dsitance matrix
  MolDist <- bruvo.dist(pop, replen = sareplen)
  bruvomat<- as.matrix(MolDist)
  #5) Rename the distance matrix to match the comm matrix
  rownames(bruvomat)<- IDs$MLG[match(rownames(bruvomat),IDs$sampleIDs)]
  colnames(bruvomat)<- IDs$MLG[match(colnames(bruvomat),IDs$sampleIDs)]
  # 6) Ensure both data sets are ordered the same
  ordering <-colnames(ComMat)
  bruvomat<- bruvomat[ordering,ordering]
  #bruvomat<- bruvomat[,match(colnames(ComMat), colnames(bruvomat))]
  bruvomat<-as.matrix(bruvomat)
  # 7) Calculate ses.mpd
  sesmpdstats<- ses.mpd(ComMat, bruvomat, runs=9999, null.model = "taxa.labels")
  # Add columns of identifiers to aid subsetting
  setPop(pop)<-~Year
  sesmpdstats$Year <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  setPop(pop)<-~Lat
  sesmpdstats$Lat <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  setPop(pop)<-~Long
  sesmpdstats$Long <- rep(unique(pop$pop), length=nrow(sesmpdstats))
  #return the results
  temp<-list(sesmpdstats=sesmpdstats, IDs=IDs, ComMat=ComMat,
             MolDist=MolDist, bruvomat=bruvomat, ordering=ordering)
  return(temp)
}


#-----------------------
# Just get the stats
#-----------------------

CalcSesmpdLat<- function(pop){
  temp<- SesmpdObjectsLat(pop) # runs sesmpdobjects funct
  stats<-temp[["sesmpdstats"]] #extracts just the stats
  #return the stats
  return(stats)
}


CalcSesmpdYear<- function(pop){
  temp<- SesmpdObjectsYear(pop) # runs sesmpdobjects funct
  stats<-temp[["sesmpdstats"]] #extracts just the stats
  #return the stats
  return(stats)
}

#-----------------------
# Get some plots
#-----------------------

PlotSesmpd<- function(pop) {
  # Run sesmpd function to get relevant objects
  temp<- SesmpdObjects(pop)
  
  # Calculate a tree
  NJtree<-njs(temp[["MolDist"]])
  NJtree$tip.label <- temp[["IDs"]]$MLG[match(NJtree$tip.label,temp[["IDs"]]$sampleIDs)]
  
  # From the ComMat make a trait dataset
  traitdat<-as.data.frame(temp[["ComMat"]])
  traitdat$Location <- rownames(traitdat) # make the rownames a col
  traitdat<-pivot_longer(traitdat,
                         1:(ncol(traitdat)-1),
                         names_to ="tip.labs", 
                         values_to = "Presence") # make the data long
  traitdat<-traitdat[,c("tip.labs",names(traitdat)[names(traitdat) != "tip.labs"])] #moves the tip labels to the front col
  
  # NOW LETS MAKE OUR PLOTS!!!
  par(mfrow = c(1, 1))
  plotlist<-list(Basetree=NJtree) # create a plotlist to store the objects in
  for (i in row.names(temp[["ComMat"]])) {
    # filter the traits to only include loc and prescence
    traits<-mutate(traitdat, Presence=(ifelse(Presence>0, 1, 0))) %>% filter(., Location==i, Presence >0)
    traittree<-ggtree(NJtree) %<+% traits
    traittree<-traittree + geom_tippoint(aes(colour=Presence), size=3)+ 
      theme(legend.position = "NULL",
            plot.margin = unit(c(0, 0, 0, 0), "cm"))
    plotlist[[i]]<-traittree # adds the tree to the plotlist
  }
  return(plotlist) # return the plotlist
}

#=========================================================================
# Map theme

# Define the map theme
maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))


# Network graph function => 

netgraph <-function(gencloneobj,longcoords,latcoords,gendistrange,graphtitle){
  
  # STEP 1 # Derive the nodes from the genclone df
  mll(gencloneobj)<- "original"
  setPop(gencloneobj) <- ~Location
  IDs<- ID2MLG2(gencloneobj) # get IDs and MLGs
  IDs<-left_join(IDs,FullNames)
  IDs$TrapNameMS<- str_extract(IDs$FullNames, "[^_]+") # Get just the locations
  nodes<- merge(IDs,trapdat) %>%
    dplyr::select(MLG, IDs, Lat, Long)# merge with geodat
  
  # STEP 2 #  Derive the edges from the distance matrix
  bdist<-bruvo.dist(gencloneobj, replen=sareplen) #Calculate a pairwise distance matrix using a distance function
  matdat<- as.matrix(bdist) # convert it to a matrix ...
  dist.df<- as.data.frame(matdat)  #then to a dataframe...
  dist.df$from<- rownames(dist.df)  # Take the rownames and make them a column...
  dist.df<- pivot_longer(dist.df, cols=1:(ncol(dist.df)-1), names_to = "to", values_to = "weights") # and pivot to get mapping of IDs to distances 
  edges<- filter(dist.df, weights>gendistrange[1] & weights<=gendistrange[2])
  
  # STEP 3 # Add jitter to coords so lines shown for each PWD
  newcoords<- as.data.frame(cbind(Lat=nodes$Lat,Long= nodes$Long, IDs=nodes$ID))
  newcoords$Lat <- as.numeric(newcoords$Lat)
  newcoords$Long <- as.numeric(newcoords$Long)
  newcoords$newLat<- jitter(newcoords$Lat, amount= 0.1)
  newcoords$newLong<-jitter(newcoords$Long, amount= 0.1)
  nodes$origLat <-nodes$Lat
  nodes$origLong <-nodes$Long
  nodes$Lat <- newcoords$newLat[match(nodes$IDs, newcoords$IDs)]  # do a matchy thing 
  nodes$Long <- newcoords$newLong[match(nodes$IDs, newcoords$IDs)]
  
  # STEP 4 # Create a dataframe for plotting
  edges_for_plot <- edges %>%
    inner_join(nodes %>% dplyr::select(IDs, Long, Lat, origLat, origLong), by = c('from' = 'IDs')) %>%
    rename(x = Long, y = Lat) %>%
    inner_join(nodes %>% dplyr::select(IDs, Long, Lat,origLat, origLong), by = c('to' = 'IDs')) %>%
    rename(xend = Long, yend = Lat)
  
  # STEP 5 # Define the area to plot over 
  worldmap = map_data('world')
  country_shapes<-geom_polygon(data = worldmap, 
                               aes(x = long, 
                                   y = lat, 
                                   group = group,
                               ),
                               fill = 'gray90', 
                               color = 'black')
  map_coords<- coord_fixed(xlim = c(longcoords[1], longcoords[2]), 
                           ylim = c(latcoords[1], latcoords[2])) # create a country shape for plotting over
  
  # STEP 6 # And make the plot
  networkplot<- ggplot(nodes) + country_shapes +
    geom_curve(aes(x = x, y = y, xend=xend, yend=yend,    # draw edges as arcs
                   colour=weights),
               data = edges_for_plot,
               alpha = 0.5) +  scale_color_viridis() +
    map_coords +
    ggtitle(graphtitle) 
 
   # Create a list of objects to return
  temp<- list(networkplot=networkplot, nodes=nodes, edges=edges, plottingdf=edges_for_plot, distmat=matdat)
  return(temp)
}

#===========================================================================
# PLOT A SPATIAL NETWORK PER MLG FUNCTION
  
plotsepClones<- function(id) {
    edges_for_plot<-dplyr::filter(pdf, from==id & to==id)
    nodesvec <- as.data.frame(unique(c(edges_for_plot$from, edges_for_plot$to))) %>% dplyr::rename(., nodeID=1)
     merged<- merge(nodesvec,trapdat)
     nodesdf<-dplyr::select(merged, 1, Lat, Long)
    worldmap = map_data('world')
    country_shapes<-geom_polygon(data = worldmap, 
                                 aes(x = long, 
                                     y = lat, 
                                     group = group,
                                 ),
                                 fill = 'gray90', 
                                 color = 'black')
    map_coords<- coord_fixed(xlim = c(-10,3),
                             ylim = c(50.3, 59))
    
    plot<- ggplot(nodesdf) + country_shapes +
      geom_curve(aes(x = x, y = y, xend=xend, yend=yend,    # draw edges as arcs
                     colour=weights),
                 data = edges_for_plot,
                 alpha = 0.5) +  scale_color_viridis() +
      map_coords 
    return(plot)
  }
  

#=========================================================================
# Calculate Mantel test

CalcMantelTest <- function(genind){
  gendistmat<-bruvo.dist(genind, replen=sareplen)
  IDs<-unique(attributes(gendistmat)$Labels)
  geodf<- FullNames[FullNames$IDs %in% IDs,]
  geodf<- merge(geodf, trapdat)
  geodistmat<- dist(cbind(geodf$Lat, geodf$Long))
  
  results<- mantel.rtest(gendistmat, geodistmat, nrepet = 9999)
  return(results)
}


#=========================================================================
CalcNoClones <- function(gid, threshold) {
  gidcc<- clonecorrect(gid)
  mlg.filter(gid, 
             distance = bruvo.dist, 
             replen = sareplen) <- threshold # Recalculates mlgs
  mlgidvec <- ID2MLG(gid) 
  IDs<-left_join(mlgidvec,FullNames)
  return(IDs)
}

