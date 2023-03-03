# Network Metrics and analysis for the chapter 3 microsatellite data # 

# Created by BM
# Updated on 03/03/2023

# Data: 
# Original data collected and scored by Gaynor Malloch at JHI
# Curated and tidied by BM. Details of curation reported in chapter three
# of thesis 

#=====================================================================#

setwd("C:/Users/beth-/OneDrive - University of Aberdeen/Data/Chapter_3")
#setwd("C:/Users/r03bm17/OneDrive - University of Aberdeen/Data/Chapter_3")


#remotes::install_github("dyerlab/popgraph")
#install.packages('devtools')
library(devtools)
devtools::install_github('thomasp85/tidygraph')
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(adegenet)
library(poppr)
library(tidyr)
library(ape)
library(stringr)
library(ggraph)
library(tidygraph)
library(treeio)
library(ape)
library(ggtree)
library(ggpubr)

rm(list=ls())

# source functions
source("./Rscripts/Ch3_functs.R")

#============================
# DATA IMPORT AND CURATION
#============================

# Read in data
#sa.complex<- read.genalex("./Genotype_data/genalex_micro.csv" , sep=",", ploidy=2)
sa.complex <- read.genalex("./Genotype_data/ALL_rawmicrosat_NArm_genalex_Treesamples.csv" , 
                           sep=",", ploidy=2)

#====================
# Data manipulation
#====================
# Split the strata
splitStrata(sa.complex) <- ~Country/Year/Location/CloneRaw/ID
sa.complex

# Add new strata (resistance, latitude, long, alt)
sa.complex@strata$Res <- rep("SS", times=nrow(sa.complex@strata)) #adding resistance strata
sa.complex@strata$Res[sa.complex@strata$CloneRaw=="SA3"] <- "SR"

trapdat<- read.csv("RIS_trap_info.csv")
sa.complex@strata$Lat <- rep(NA, times=nrow(sa.complex@strata)) #adding trap info
sa.complex@strata$Long <- rep(NA, times=nrow(sa.complex@strata)) 
sa.complex@strata$Alt <- rep(NA, times=nrow(sa.complex@strata)) 
sa.complex@strata$Lat<-trapdat$Lat[match(sa.complex@strata$Location,trapdat$TrapNameMS)]
sa.complex@strata$Long<-trapdat$Long[match(sa.complex@strata$Location,trapdat$TrapNameMS)]
sa.complex@strata$Alt<-trapdat$Alt[match(sa.complex@strata$Location,trapdat$TrapNameMS)]

# Define the rep length of our microsats: 
sareplen = c(2,2,2,2,2)

# All tree samples samples:
treesamps<-c("Eng", "Scot", "China")
setPop(sa.complex) <- ~Country
ukwog<-popsub(sa.complex, sublist=treesamps) 
setPop(ukwog) <- ~Location
ukwog<- popsub(ukwog, exclude="SRcontrol") 

# CLONE CORRECT 
# add a strata that covers all samples
ukwog$strata$ALL <- rep("ALL", times = nrow(ukwog$strata))
ukwogcc<- clonecorrect(ukwog, strata = ~ALL)

# create a full names list to cross reference against
FullNames <- CalcFullNames(ukwogcc)

# Create own colour shape scale to deal with the number of colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_pch_cross = expand.grid(
  pch = 1:8,
  colors = gg_color_hue(6)
  ,stringsAsFactors = F)


#================================================================
# Whole graph:
setPop(ukwogcc) <- ~Res
set.seed(109)
# collapse by MLL
msn<- bruvo.msn(ukwogcc, replen=sareplen, threshold = 0.18, mlg.compute="original", include.ties =TRUE)
# original MLGs
msn2<- bruvo.msn(ukwogcc, replen=sareplen, threshold = 0, mlg.compute="original", include.ties =TRUE)

g <- msn$graph 
g2<- msn2$graph


# plot with ggraph
ggraph.df<-tidygraph::as_tbl_graph(g)
ggraph.df<-mutate(ggraph.df, 
                       Location= FullNames$Location[match(name, FullNames$IDs)],
                       Res= FullNames$Res[match(name, FullNames$IDs)],
                       Year= FullNames$Year[match(name, FullNames$IDs)])

ggraph(ggraph.df, layout = "kk") + 
  geom_edge_link(aes(alpha=weight),edge_width=1) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res, size=size)) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))

# plot with ggraph
ggraph.df<-tidygraph::as_tbl_graph(g2)
ggraph.df<-mutate(ggraph.df, 
                  Location= FullNames$Location[match(name, FullNames$IDs)],
                  Res= FullNames$Res[match(name, FullNames$IDs)],
                  Year= FullNames$Year[match(name, FullNames$IDs)])

completenet<- ggraph(ggraph.df, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=2) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=3) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))


 # png(file = "./Graphs/CompleteNetwork.png", width = 1200, height = 1000, units = 'px')
 # print(completenet)
 # dev.off()


#=============================================================
# Subset the data to the SR CLADE: 
#==============================================================

graph_tidy <- as_tbl_graph(g2) 

# plot with ggraph
ggraph.df<-mutate(graph_tidy, 
                  Location= FullNames$Location[match(name, FullNames$IDs)],
                  Res= FullNames$Res[match(name, FullNames$IDs)],
                  Year= FullNames$Year[match(name, FullNames$IDs)])

resggraph<- ggraph.df %>% 
  activate(nodes) %>% 
  filter(Res=="SR")

resggraph

ggraph(resggraph, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=1) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=1) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))


#================================================================
# ...but we want to see those adjacent aswell
#=================================================================
resvertices<- V(resggraph)$name

adjvs<-adjacent_vertices(g2, resvertices)
N = unlist(adjvs)
resadj<- names(V(g2)[N])

resadjggraph<- ggraph.df %>% 
  activate(nodes) %>% 
  filter(name %in% resadj)

resadjggraph

ggraph(resadjggraph, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=1) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=1) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))


# and second order adjacency: 

resvertices2<- V(resadjggraph)$name

adjvs2<-adjacent_vertices(g2, resvertices2)
N2 = unlist(adjvs2)
resadj2<- names(V(g2)[N2])

resadjggraph2<- ggraph.df %>% 
  activate(nodes) %>% 
  filter(name %in% resadj2)

resadjggraph2

ggraph(resadjggraph2, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=1) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=1) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))


# and third order: 

resvertices3<- V(resadjggraph2)$name

adjvs3<-adjacent_vertices(g2, resvertices3)
N3 = unlist(adjvs3)
resadj3<- names(V(g2)[N3])

resadjggraph3<- ggraph.df %>% 
  activate(nodes) %>% 
  filter(name %in% resadj3)

resadjggraph3

ggraph(resadjggraph3, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=1) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res, shape=Res),size=1) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))

resubgraph<-ggraph(resadjggraph3, layout = "fr") + 
  geom_edge_link(aes(alpha=weight),edge_width=2) + 
  scale_edge_width(range = c(0.05, 2))  + 
  geom_node_point(aes(colour=Res, shape=Res),size=3, alpha=0.7) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))

# png(file = "./Graphs/ResNet_subgraph.png", width = 800, height = 600, units = 'px')
# print(resubgraph)
# dev.off()


#===============================================================================
# Looking at network connectivity of in the RES SUBCLADE (AS EXTRACTED FROM THE PHYLOGENY)


# Import the list of IDs from the SR clade on the tree:

tree<- read.tree("./Trees/njtree_oneoutgroup.tre")

# get the node labels
tree$tip.label<-FullNames$TraitMatName[match(tree$tip.label, FullNames$FullNames)]
# Increasing levels back shows 6 about is the right number!
subtree<-tree_subset(tree, node= "Eng_2015_Writtle_SA3_ID856",
                     levels_back=6)
#write.tree(subtree, "./Trees/SRSubtree.tre")
trait_data <- read.csv("Ch3_TraitData.csv")
str(trait_data)

dend<-ggtree(subtree,layout = "dendrogram") %<+% trait_data
dend + geom_tippoint(aes(color=Resistance), size=2,alpha=.5)+ 
  scale_color_brewer("Resistance", palette="Set1") 

# Extract list of names # should be 346
SRCladetips<- subtree$tip.label
setPop(uksa) <- ~Country/Year/Location/CloneRaw/ID
SRClade<-popsub(uksa, sublist= SRCladetips) 
SSClade<- popsub(uksa, exclude=SRCladetips)

setPop(SRClade) <- ~Year
setPop(uksa) <-~Year

SRClade<-popsub(SRClade, sublist=c("2013", "2014", "2015","2016", "2017"))
SR2013<-popsub(SRClade, sublist="2013")
SR2014<-popsub(SRClade, sublist="2014") 
SR2015<-popsub(SRClade, sublist="2015")
SR2016<-popsub(SRClade, sublist="2016") 
SR2017<-popsub(SRClade, sublist="2017")
ALL2013<-popsub(uksa, sublist="2013")
ALL2014<-popsub(uksa, sublist="2014") 
ALL2015<-popsub(uksa, sublist="2015")
ALL2016<-popsub(uksa, sublist="2016") 
ALL2017<-popsub(uksa, sublist="2017")
SRClade<-popsub(SRClade, sublist=c("2013", "2014", "2015"))
SR2015cc<-clonecorrect(SR2015, strata = ~Location)
SR2014cc<-clonecorrect(SR2014, strata = ~Location)
SR2013cc<-clonecorrect(SR2013, strata = ~Location)
SR2016cc<-clonecorrect(SR2016, strata = ~Location)
SR2017cc<-clonecorrect(SR2017, strata = ~Location)
ALL2015cc<-clonecorrect(ALL2015, strata = ~Location)
ALL2014cc<-clonecorrect(ALL2014, strata = ~Location)
ALL2013cc<-clonecorrect(ALL2013, strata = ~Location)
ALL2016cc<-clonecorrect(ALL2016, strata = ~Location)
ALL2017cc<-clonecorrect(ALL2017, strata = ~Location)



# LOOKING AT YEARLY NETWORKS USING SAMPLE SUBSETTED FROM THE TREE:
#---------------------------------------------------------------------

msn2013<- bruvo.msn(SR2013cc, replen = c(sareplen), add = TRUE, loss = TRUE, showplot = FALSE, include.ties = TRUE)
ggraph.df.2013<-tidygraph::as_tbl_graph(msn2013$graph)
ggraph.df.2013<-mutate(ggraph.df.2013, 
                       Location= FullNames$Location[match(name, FullNames$IDs)],
                       Res= FullNames$Res[match(name, FullNames$IDs)],
                       Year= FullNames$Year[match(name, FullNames$IDs)])

# 2014
msn2014<- bruvo.msn(SR2014cc, replen = c(sareplen), add = TRUE, loss = TRUE, showplot = FALSE, include.ties = TRUE)
ggraph.df.2014<-tidygraph::as_tbl_graph(msn2014$graph)
ggraph.df.2014<-mutate(ggraph.df.2014, 
                       Location= FullNames$Location[match(name, FullNames$IDs)],
                       Res= FullNames$Res[match(name, FullNames$IDs)],
                       Year= FullNames$Year[match(name, FullNames$IDs)])

# 2015
msn2015<- bruvo.msn(SR2015cc, replen = c(sareplen), add = TRUE, loss = TRUE, showplot = FALSE, include.ties = TRUE)
ggraph.df.2015<-tidygraph::as_tbl_graph(msn2015$graph)
ggraph.df.2015<-mutate(ggraph.df.2015, 
                       Location= FullNames$Location[match(name, FullNames$IDs)],
                       Res= FullNames$Res[match(name, FullNames$IDs)],
                       Year= FullNames$Year[match(name, FullNames$IDs)])

#-----------------------------------------------
# Plotting graphs with just the resistance data
#-----------------------------------------------

net2013<-ggraph(ggraph.df.2013, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))
net2013

net2014<-ggraph(ggraph.df.2014, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))
net2014

net2015<-ggraph(ggraph.df.2015, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Res),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(nrow = 2), alpha = guide_legend(nrow = 3))
net2015

#-------------------------------------
# Extract graph attributes for nodes
#-------------------------------------
# 2013
closenessvec<- closeness(msn2013$graph)
degs<-igraph::degree(msn2013$graph, v= V(msn2013$graph))
graphatts2013<- data.frame(closeness=closenessvec, 
                           IDs=names(closenessvec),
                           year=rep("2013", times=length(closenessvec)))
graphatts2013<-graphatts2013[match(names(degs), graphatts2013$IDs),] 
graphatts2013$degrees <-degs

# 2014
closenessvec<- closeness(msn2014$graph)
degs<-igraph::degree(msn2014$graph, v= V(msn2014$graph))
graphatts2014<- data.frame(closeness=closenessvec, 
                           IDs=names(closenessvec),
                           year=rep("2014", times=length(closenessvec)))
graphatts2014<-graphatts2014[match(names(degs), graphatts2014$IDs),] 
graphatts2014$degrees <-degs

# 2015
closenessvec<- closeness(msn2015$graph)
degs<-igraph::degree(msn2015$graph, v= V(msn2015$graph))
graphatts2015<- data.frame(closeness=closenessvec, 
                           IDs=names(closenessvec),
                           year=rep("2015", times=length(closenessvec)))
graphatts2015<-graphatts2015[match(names(degs), graphatts2015$IDs),] 
graphatts2015$degrees <-degs


#-----------------------------------------------
# Extracting adjacency data 
#------------------------------------------------
#merge dfs
allatts<- rbind(graphatts2013,graphatts2014, graphatts2015) %>% 
  merge(., (select(FullNames, IDs, Location, Res)))

# Identify the suceptible IDS
set.seed(26)
SRCladecc<- clonecorrect(SRClade, strata =~Year/Lat)
setPop(SRCladecc) <- ~Res
susIDS <- popsub(SRCladecc, "SS")
setPop(susIDS) <- ~ID
IDvec<- susIDS$pop %>% gsub("ID", "", .)  %>% as.numeric(.)

# Capture the resistant IDS adjacent
get.vertex.attribute(msn2013$graph)
select2013 = (V(msn2013$graph))[names(V(msn2013$graph)) %in% IDvec]
select2013 #775 & 74
adj= adjacent_vertices(msn2013$graph,select2013,mode="all")
N=unlist(adj)
adj2013<- names(V(msn2013$graph)[N])
Res2SusTable_2013<- data.frame(ResIndex=N)
Res2SusTable_2013$SusName <- c("775", "775", "775", "74")
Res2SusTable_2013$ResName <- c("71", "676", "679", "67")  

select2014 = (V(msn2014$graph))[names(V(msn2014$graph)) %in% IDvec]
select2014
adj = adjacent_vertices(msn2014$graph,select2014,mode="all")
N = unlist(adj)
adj2014<- names(V(msn2014$graph)[N])
Res2SusTable_2014<- data.frame(ResIndex=N)
Res2SusTable_2014$SusName <- c("1076", "1076", "1028", "1028")
Res2SusTable_2014$ResName <- c("223", "205", "192", "1001")  

select2015 = (V(msn2015$graph))[names(V(msn2015$graph)) %in% IDvec]
select2015
adj = adjacent_vertices(msn2015$graph,select2015,mode="all")
N = unlist(adj)
adj2015<- names(V(msn2015$graph)[N])
Res2SusTable_2015<- data.frame(ResIndex=N)
Res2SusTable_2015$SusName <- c("379", "379", "868", "868", "931", "386")
Res2SusTable_2015$ResName <- c("338", "835", "833", "386", "834", "868")  

alladj<- c(adj2013,adj2014, adj2015)
TL_SR<-filter(FullNames, IDs %in% alladj & Res=="SR")
TL_SR$IDs

# Add this to the Fullnames data
FullNames$NetGroup <- rep("SS", times=nrow(FullNames))
FullNames$NetGroup[FullNames$Res=="SR"] <- "SR"
FullNames$NetGroup[FullNames$IDs %in% TL_SR$IDs] <- "TL-ADJ-SR"

# Add these ids to the main attributes dataframe
allatts$NetStat<- FullNames$NetGroup[match(allatts$IDs, FullNames$IDs)]

length(which(allatts$year==2013)) #32
length(which(allatts$year==2014)) #74
length(which(allatts$year==2015)) #78
nsamps<- data.frame(year=c(2013,2014,2015),
                    n= c(32,74,78))

allatts<- merge(allatts, nsamps) # add sample numbers

# add in the data that maps res to sus points 
temp <- rbind (Res2SusTable_2013, Res2SusTable_2014, Res2SusTable_2015)
temp <- temp %>% rename(., IDs = ResName)

allatts<-left_join(allatts, temp)

#-----------------------------------------------
# Re-plotting graphs with the adjacency data
#-----------------------------------------------
ggraph.df.2013<-mutate(ggraph.df.2013, 
                       Group= FullNames$NetGroup[match(name, FullNames$IDs)])

net2013<-ggraph(ggraph.df.2013, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Group),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20))  + 
  scale_colour_manual(values=c("#F8766D",  "#00BFC4", "#7CAE00"))
net2013


ggraph.df.2014<-mutate(ggraph.df.2014, 
                       Group= FullNames$NetGroup[match(name, FullNames$IDs)])

net2014<-ggraph(ggraph.df.2014, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Group),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) + 
  scale_colour_manual(values=c("#F8766D",  "#00BFC4", "#7CAE00"))
net2014


ggraph.df.2015<-mutate(ggraph.df.2015, 
                       Group= FullNames$NetGroup[match(name, FullNames$IDs)])

net2015<-ggraph(ggraph.df.2015, layout = "fr") + 
  geom_edge_link(aes(alpha=weight), edge_width=3) + 
  scale_edge_width(range = c(0.05, 1))  + 
  geom_node_point(aes(colour=Group),size=6) + 
  theme(legend.key.size = unit(0.5, 'cm')) + 
  theme(legend.position = "bottom",
        legend.title=element_text(size=30),
        legend.text=element_text(size=20)) + 
  scale_colour_manual(values=c("#F8766D",  "#00BFC4", "#7CAE00"))
net2015


#------------------------------
# merging the network plots 
#------------------------------
nets<- ggarrange(net2013,
                 net2014, 
                 net2015, 
                 nrow=1, 
                 ncol=3,
                 common.legend = TRUE,
                 legend="bottom",
                 labels= c("2013", "2014", "2015"),
                 font.label = list(size = 30))
nets


#-----------------------------------------
# Making network attribute histograms
#-----------------------------------------
dpoints <- allatts %>% filter(., NetStat=="TL-ADJ-SR") %>% select(., degrees, year, NetStat, SusName) 
# calculate an average PER SUS point (account for non-independence)
sumdpoints <- dpoints %>% group_by(SusName, year) %>% summarise(., degrees=mean(degrees))
sumdpoints

degreesHist <- 
  ggplot(allatts, aes(x=degrees, fill=NetStat)) +
  geom_histogram(binwidth=1)+
  geom_vline(data = sumdpoints, aes(xintercept=degrees), colour="#7CAE00", size=2,lty=2) +
  #geom_point(data=points, aes(x=closeness, y=0),stroke= 1, pch=4) + 
  ylab("Frequency")+
  scale_fill_manual(values=c("#F8766D",  "#00BFC4", "#7CAE00")) +
  plot_theme + facet_wrap(~year, scales="free") + 
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=30),
        strip.text.x = element_text(size = 25)) + 
  xlab("Degree centrality")
degreesHist

cpoints <- allatts %>% filter(., NetStat=="TL-ADJ-SR") %>% select(., closeness, year, NetStat, SusName) 
# calculate an average PER SUS point (account for non-independence)
sumcpoints <- cpoints %>% group_by(SusName, year) %>% summarise(., closeness=mean(closeness))
sumcpoints

closenessHist <- 
  ggplot(allatts, aes(x=closeness, fill=NetStat)) +
  geom_histogram(binwidth=0.01)+
  geom_vline(data = cpoints, aes(xintercept=closeness), colour="#7CAE00", size=2,lty=2) +
  #geom_point(data=points, aes(x=closeness, y=0),stroke= 1, pch=4) + 
  ylab("Frequency")+
  scale_fill_manual(values=c("#F8766D",  "#00BFC4", "#7CAE00")) +
  plot_theme + facet_wrap(~year, scales="free") + 
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=30),
        strip.text.x = element_text(size = 25)) + 
  xlab("Closeness centrality")
closenessHist

#================
# FINAL DIAGRAM
#================

RESgraphanalysis<- ggarrange(nets,
                             degreesHist, 
                             closenessHist,
                             nrow=3,
                             common.legend = F,
                             align="hv",
                             heights= c(2,1,1),
                             vjust=1.1, hjust=0)
RESgraphanalysis

png(file = "./Graphs/networkanalysisYEARS_HISTS.png", width = 1200, height = 1000, units = 'px')
print(RESgraphanalysis)
dev.off()


