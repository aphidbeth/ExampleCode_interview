# Plotting an ML tree from the Aphrad19 field samples with ggtree # 

# Authors: BM
# Updated on: 03 March 2023

# DATA: 
# Phylip file generated from field samples in stacks => 
# ML tree generated in raxml using GTGAMMA + lewis correction 
# ascertainment bias 

#===================================================================

# housekeeping

library(phytools)
library(treeio)
library(treemap)
library(ape)
library(pegas)
library(phangorn)
library(ggplot2)
library(ggtree)
library(dplyr)
library(ggtreeExtra)

#setwd("C:/Users/beth-/OneDrive - University of Aberdeen/Data/AphRad19")
setwd("C:/Users/r03bm17/OneDrive - University of Aberdeen/Data/AphRad19")

rm(list=ls())

#-----------------#
# Import tree
#-----------------#

# Load read tree in that was created with RAxML: 
tree<- read.raxml("./Data/FilteringOutputs/FinalTreeOutputs/RAxML_bipartitionsBranchLabels.FINAL1000.tre")

tree@phylo$tip.label
tree@data$bootstrap

#Root by the midpoint
midrootML<-midpoint.root(tree@phylo)
midrootML

tree@phylo <- midrootML
tree

ggtree(midrootML, layout="circular") + geom_text2(aes(label=label))


# Add the bs data to the phylo object so easier for plotting 
tree@data
midrootML$edge # bootstraps are in the same order as the edges 
midrootML$node.label <- as.matrix(tree@data)
midrootML["node.label"] <- (as.matrix(tree@data))
midrootML$node.label<- (as.matrix(tree@data))


# -------------------------------#
 # Import and curate trait data   #
#--------------------------------#

sampleinfo <- read.csv("./Data/Sample_Metadata.csv")
sampleinfo <- dplyr::rename(sampleinfo, Plate=3, DNAConc=5)
str(sampleinfo)

fielddat <-read.csv("./Data/Site_info_2019.csv")
fielddat$Pop<- sub("*19_", "", fielddat$Site_name) # add the populations
fielddat$Pop[fielddat$Pop=="X_A"]<-"X"
fielddat$Pop[fielddat$Pop=="X_B"]<-"XX" 
fielddat <- dplyr::rename(fielddat, FieldName=3, SiteName=8)
str(fielddat)

# Read missingness info: 
IDmiss <- read.csv("./Data/PerSampleMissingness.csv")
IDmiss

# select just the variables we want: 
fieldsel <- dplyr::select(fielddat, Date,FieldName, SiteName, Lat, Long, Crop_type, Sampling_strategy, Insecticide_sprayed, Pop)
sampsel <- dplyr::select(sampleinfo, Sample.ID, FieldName, Plate, Well, DNAConc)

traitdat<- left_join(sampsel, fieldsel)
traitdat <- filter(traitdat, Sample.ID %in% tree@phylo$tip.label)
traitdat$fmiss <- IDmiss$f_miss[match(IDmiss$Sample, traitdat$Sample.ID)]
traitdat<- cbind(label = traitdat$Sample.ID,         # Append new column with ID labels to first col of df
                   traitdat)
str(traitdat)

# make latitude a rounded character for plotting 
traitdat$LatC <- as.character(round(traitdat$Lat, digits=3))
IDs<- traitdat$Sample.ID[is.na(traitdat$LatC)==TRUE]
traitdat$LatC[is.na(traitdat$LatC)==TRUE] <- IDs
traitdat$LatC  
traitdat$LatF <- as.factor(traitdat$LatC)

#-------------------------------#
  # Plot tree with trait data  #
#-------------------------------#

treeplot<- ggtree(midrootML, layout="circular") %<+% traitdat

# Plot tree + host crop
treeplot + geom_tippoint(aes(color=Crop_type), size=2,alpha=.5)+ 
  scale_color_brewer("Host crop", palette="Set1")

# Plot tree + seq plate
treeplot+ geom_tippoint(aes(color=Plate), size=2,alpha=.5)+ 
  scale_color_brewer("Plate", palette="Set1") 

# Plot tree + latitude
treeplot + geom_tippoint(aes(color=Lat), size=3,alpha=.5)+ 
  scale_color_continuous(type="viridis") 

# Make bootstraps a category
bs<-tree@data
bs$support <- ifelse(bs$bootstrap>70, "yes", "no")

# CHANGE THE SA1 AND SA6 NAMES to SA3_B and SA3_A respectively:
treeplot$data$label<- sub("SA1-", "SA3_B-", treeplot$data$label)
treeplot$data$label<- sub("SA6-", "SA3_A-", treeplot$data$label)

# Rotate so SA3 clade on outer part of tree for clarity: 
treeplot %<+% bs + geom_nodepoint(aes(shape=support), size=4) + 
  geom_nodelab(aes(label=node, vjust=2))

treeplot<- rotate(treeplot, 181)


# Plot tree + latitude + bootstrap values
treeplot1 <-treeplot %<+% bs + geom_nodepoint(aes(shape=support), size=4) + 
  geom_tiplab(aes(label=label, color=Lat), size=5,align=T) +
  scale_color_gradient("Latitude", high= "darkblue", low = "darkred") + xlim(-0.1, NA) + 
  scale_shape_manual("Boostrap support", #legend name
                     values=c(20,8),
                     na.translate = F,#remove na,
                     labels = c(levels(bs$support), "0-70","70-100"))  + 
  theme(legend.text = element_text(size=12),  legend.title = element_text(size=14),
        legend.position = "bottom")
treeplot1


#--------------------------------------------#
# Add an outer ring of colours showing MLGs
#--------------------------------------------#

# import MLG groupings created in defineMLGs.R
MLGgrps <- read.csv("./Data/MLGlist_tree.csv")
str(MLGgrps)
fieldclones <-c("Unique", "MLG.142", "MLG.138", "MLG.162", "MLG.119", "MLG.6", "MLG.110",
                "MLG.169", "MLG.114", "MLG.103", "MLG.93", "MLG.144", "MLG.140", "MLG.2", "MLG.165")
MLG.grps <- filter(MLGgrps, CloneCat %in% fieldclones)

treeplot_MLG <- treeplot1 + geom_fruit(data=MLG.grps, geom=geom_tile,
                    mapping=aes( x=CloneCat, y=sampleIDs, fill=CloneCat), width=0.1,
                    color="white",
                    pwidth=0.05,
                    offset=0.3
)+
  scale_fill_manual(
    name="MLG",
    values=c("#114B5F", "#1A936F", "#88D498", "#C6DABF", "#F3E9D2",
             "#0A369D", "#4472CA", "#5E7CE2", "#92B4F4", "#F45B69",
             "#FFCDB2", "#FFB4A2", "#66101F", "#B5838D", "#6D6875"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=2,
                       keyheight=2,
                       order=3
    )
  ) 
 


treeplot_MLG
png("./Trees/MLGtree.png", width=1400, height=2500)
print(treeplot_MLG)
dev.off()

#===================================================================#
