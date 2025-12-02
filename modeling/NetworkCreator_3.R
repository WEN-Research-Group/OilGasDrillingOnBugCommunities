rm(list = ls())

setwd("C:/Users/rpolivie/OneDrive - Syracuse University/Tao DATA") 

getwd()
list.files()

install.packages(c("netassoc","readxl","dplyr","lubridate","sp","FNN","sf"))
library(netassoc)
library(readxl)
library(dplyr)
library(lubridate)
library(sp)
library(FNN)
library(sf)
library(ggplot2)

bug_data <- read.csv("NetworkAnalysisData", header = TRUE, stringsAsFactors = FALSE)
glimpse(bug_data)

#only considering taxa sampled from wadeable (small) streams, sampled in the spring, with less then 20% development and no aml present in the catchment
bug_data_filtered <- bug_data %>%
  filter(has_aml == 0, season == "spring", Small ==1, DEV<=20)

bug_data_filtered <- bug_data_filtered %>%
  filter(!is.na(Count))

#no na values 
colSums(is.na(bug_data_filtered))


#make networks for huc8
set.seed(99)
sampleMin <- 15
samplingNum <- 10
#for every huc8 with at least 15 samples, randomly sample 10 to make the network
bugdata_clean <- subset(bug_data_filtered, NStations >= sampleMin)
bugdata_clean$U_bool <- as.factor(bugdata_clean$U_bool)
bugdata_clean$C_bool <- as.factor(bugdata_clean$C_bool)

PANetworkAnalysis <- data.frame()

#run 100 times for each qualifying huc8
for(j in 1:100){
  for(WS in unique(bugdata_clean$HUC8)){
    for(U in unique(bugdata_clean$U_bool)){
      for(C in unique(bugdata_clean$C_bool)){
        
        bugdatalocal <- subset(bugdata_clean, HUC8 == WS & U_bool == U & C_bool == C)
        
        if(length(unique(bugdatalocal$STATIONID)) >= sampleMin){
          bugdatalocal <- bugdatalocal[bugdatalocal$STATIONID %in% sample(unique(bugdatalocal$STATIONID),samplingNum),]
          metadata <- unique(bugdatalocal[, c("STATIONID", "IBI", "DEV", "DEP_ASSESS_ID_TV")])
          
          U <- U
          C <- C
          MeanIBI <- mean(metadata$IBI, na.rm = TRUE)
          SDIBI <- sd(metadata$IBI, na.rm = TRUE)
          MeanDEV <- mean(metadata$DEV, na.rm = TRUE)
          SDDEV <- sd(metadata$DEV, na.rm = TRUE)
          MeanPTV <- mean(metadata$DEP_ASSESS_ID_TV, na.rm = TRUE)
          SDPTV <- sd(metadata$DEP_ASSESS_ID_TV, na.rm = TRUE)
          
          
          uniqueTaxa <- unique(bugdatalocal$Taxa.ID)  # All unique taxa in the subset
          station_ids <- unique(bugdatalocal$STATIONID)  # All unique station IDs in the subset
          
          PresenceAbsence <- matrix(0,
                                    nrow = length(uniqueTaxa),
                                    ncol = length(station_ids),
                                    dimnames = list(uniqueTaxa, station_ids))
          
          for (ID in station_ids) {
            
            sampleDF <- subset(bugdatalocal, STATIONID == ID)[, c("Taxa.ID", "Count")]
            
            valid_taxa <- sampleDF$Taxa.ID %in% rownames(PresenceAbsence)
            
            PresenceAbsence[sampleDF$Taxa.ID[valid_taxa], as.character(ID)] <- ifelse(sampleDF$Count[valid_taxa] > 0, 1, 0)}
          
          m_obs <- as.data.frame(PresenceAbsence)
          m_obs[m_obs > 0] <- 1
          
          n <- make_netassoc_network(
            m_obs, vegan::permatfull(m_obs, fixedmar = "both", mtype = "prab", times = 100)$perm[[1]],
            method = "partial_correlation", args = list(method = "shrinkage"), p.method = 'fdr',
            numnulls = 1000, plot = TRUE, alpha = 1e-4, verbose = FALSE)
          
          networkgraph <- as.undirected(n$network_all,mode="collapse") #Generate full graph from co-occurrence patterns.
          networkgraph_pos <- as.undirected(n$network_pos,mode="collapse")
          networkgraph_neg <- as.undirected(n$network_neg, mode="collapse")
          
          if(gsize(networkgraph) > 0 & gsize(networkgraph_pos) > 0 & gsize(networkgraph_neg) > 0){
            networkgraph <- as.undirected(n$network_pos,mode="collapse") #Generate graph from co-occurrence patterns.
            
            networkmatrix <- as.matrix(get.adjacency(networkgraph,attr='weight')) #Generate matrix from co-occurrence graph.
            S_pos <- mean(networkmatrix[networkmatrix!=0]) #Mean edge weight as classified by standardized effect size values
            M_pos <- modularity(cluster_edge_betweenness(networkgraph, weights=NULL,directed=FALSE)) #Modularity
            N_pos <- gorder(networkgraph) #Number of nodes
            C_pos <- edge_density(networkgraph) #Connectance
            taxa_in_network <- V(networkgraph)$name  # This gets all node names
            taxa_in_subset <- paste(taxa_in_network, collapse = ", ")  # Create comma-separated list
            
            row <- t(as.data.frame(c(WS, U, C, MeanIBI, SDIBI, MeanDEV, SDDEV, MeanPTV, SDPTV, S_pos, M_pos, N_pos, C_pos, taxa_in_subset)))
            
            PANetworkAnalysis <- rbind(PANetworkAnalysis, row)
          }
        }
      }
    }
  }
}

colnames(PANetworkAnalysis) <-c("HUC8","U_bool","C_bool","MeanIBI","SDIBI","MeanDEV","SDDEV","MeanPTV","SDPTV","S_pos","M_pos","N_pos", "C_pos",'Taxa_List')
num_cols <- c("MeanIBI", "SDIBI", "MeanDEV", "SDDEV", "MeanPTV", "SDPTV", "S_pos", "M_pos", "N_pos", "C_pos")
PANetworkAnalysis[num_cols] <- lapply(PANetworkAnalysis[num_cols], as.numeric)
rownames(PANetworkAnalysis) <- 1:nrow(PANetworkAnalysis)
write.table(PANetworkAnalysis,"BMI_Network_CooccurrenceAnalysis.txt",quote=FALSE,sep="\t",row.names = TRUE)
