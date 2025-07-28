#Preparing the environment

library(Biostrings)
library(pwalign)
library(igraph)
library(bio3d)
library(ggraph)
library(ggplot2)
library(DSSP)
library(readxl)
library(tidyverse)
library(stringr)

setwd("C:/Users/Jelena/Desktop/zato/dssp/build/Release")




qnx <- read.pdb("1qnx.pdb") #reading the pdb

#Removing all atoms inside the protein:

qnx <- trim(qnx, "protein") #to remove NA and HOH, stopped working for some reason

qnx_dssp <- dssp(qnx, exefile = "mkdssp.exe")  #doesnt work
qnx_dssp <- dssp.pdb(qnx, exefile = "mkdssp.exe")  # also not working

dssp_data <- read_excel("dssp_table.xlsx") #reading the data like this for now
acc <- dssp_data$ACC  #extracting only solvent exposure data

C_alpha <- qnx$atom[qnx$calpha,] #selecting only C alpha atoms

df <- cbind(C_alpha, acc)  #combining the data together
summary(acc)  # could use mean or median for the treshold
#a more precise way is to use Absolute surface area of each residue
#then discard all residues with <25% relative area in contact with the solvent

#making the absolute surface area table for each residue
asa_table <- as.data.frame(AMINO_ACID_CODE)
asa_table$asa <- c(129,274,195,193,167,225,223,104, 224,197,201,236,224,240,159,155,172,285,263,174,0,0,0,0,0,0)
asa_table <- asa_table[1:20,]
asa_table <- asa_table %>% rename(resid = AMINO_ACID_CODE) %>% mutate(resid = toupper(resid))

df <- merge(df, asa_table, by = "resid")
df$rsa <- 100*df$acc/df$asa  #calculating relative surface area
df_outside <- df[df$rsa > 5.0, ] #only selecting residues with big solvent exposure

row.names(df_outside) <- str_c(aa321(df_outside$resid), df_outside$resno, sep="") #changing the labels of atoms


dist_matrix <- as.matrix(dist(df_outside[,c("x", "y", "z")])) #calculating pairwise distance with dist()
dist_matrix[dist_matrix >= 6] <- 0
dist_matrix[dist_matrix > 0] <- 1 # only those close enough are important for the graph

graph <- graph_from_adjacency_matrix(dist_matrix, mode = "undirected", weighted = FALSE)

plot(graph, 
     vertex.size = 10, 
     vertex.label.cex = 0.5,  
     vertex.label.color = "black",
     vertex.color = "orange",
     edge.color = "black",
     edge.width = 1,
     main = "1qnx_graph")

#q111 not connected since A112 has acc1 and G110 has acc4

###########################################################


pattern = "LKPNKY"

#function for making words
extract_path_letters <- function(paths) {
  sapply(paths, function(path) {
    node_names <- path$name
    first_letters <- substring(node_names, 1, 1)
    paste(first_letters, collapse = "")
  })
}


words <- list()

i <- 1
while (i < (length(graph)+1)){
  sub <- induced_subgraph(graph, c(unlist(adjacent_vertices(graph, V(graph)[i])), V(graph)[i]))
  j <- 1
    while (j < (length(sub)+1)){
    paths <- all_simple_paths(sub, from = V(sub)[j], mode = "all")
    words <- append(words , extract_path_letters(paths))
    j = j+1
  }
  i = i+1
}


df_words <- as.data.frame(words)

l <- 1
while(l < length(df_words)) {
  if(nchar(df_words[,l]) < 4) {
    df_words <- df_words[,-l]
  } 
  else { l <- l + 1}
}


pairwiseAlignment(pattern, "LGPRKQ", substitutionMatrix = "BLOSUM50",   gapOpening = 20,  gapExtension = 10)

matches <- data.frame(match = "a", scored = 0, pos = 0)
k = 1
while (k < (length(df_words)+1)){
  pwa <- pairwiseAlignment(pattern, df_words[,k], substitutionMatrix = "BLOSUM50",   gapOpening = 20,  gapExtension = 10)
  if (score(pwa) > 10) {
    matches <- add_row(matches, match = df_words[,k], scored = score(pwa), pos= 0 )
  }
  k <- k + 1
}

getwd()
matches

best <- matches[matches$scored == 17,]
saveRDS(best, "best.RDS")


best <- readRDS(file = "C:/Users/Jelena/Desktop/zato/dssp/build/Release/best.RDS")
