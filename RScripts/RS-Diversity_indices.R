############################################################################################
#Packages
############################################################################################
rm(list=ls())
install.packages(c('ggplot2',  'vegan', 'reshape2', 'fossil', 'corrplot', 'dplyr', 
                   'philentropy', 'gplots', 'cluster', 'dendextend', 'productplots'))
############################################################################################
##Open the dataset
############################################################################################
data<-read.csv2("/Users/alinevaz/Documents/R/Metagenomica_cacau/OTU_TABLE_sem_cordados-sem_streptophyta.csv", header=T)
meta<-read.csv2("/Users/alinevaz/Documents/R/Metagenomica_cacau/metadata_cacao_ME.csv", header=T)
############################################################################################
#Organizing the datasets
# this script refers only to bacteria data
############################################################################################
bac<- data[data$Kingdom == "Bacteria",]
euc<- data[data$Kingdom == "Eukaryota",]
arc<- data[data$Kingdom == "Archaea",]
euc<- euc[-c(273,254,262,289,38,39,55,70,182,187,195,201,209,263,267, 277, 285,44,132,133,191,270, 231), ]
list<- c("Arthropoda", "Cercomonadida[phylum]", "Bacillariophyta", "Chlorophyta", "Colpodellidae[phylum]", "Thaumatomonadida[phylum]","Dinophyceae","Cryptophyta", "Kinetoplastida","Longamoebia")
############################################################################################
####Bacterial datasets - Figure 1
############################################################################################
bac1<-bac[,-c(1:7)]#removing the taxonomy informations
bac2<-aggregate(bac1, by = list(bac[,7]), FUN= sum)
aggreg_gen<- bac2
rownames(aggreg_gen)<- aggreg_gen[,1]
aggreg_gen_1<- aggreg_gen[,-1]
aggreg_gen_2<-t(aggreg_gen_1)
aggreg_gen_2[,1:ncol(aggreg_gen_2)] <- aggreg_gen_2[,1:ncol(aggreg_gen_2)]/ rowSums(aggreg_gen_2[,1:ncol(aggreg_gen_2)])##normalizing the dataset
rownames(aggreg_gen_2)<- c("1.for0h", "8.mix0h", "2.for24h" , "9.mix24h", "3.for48h", "10.mix48h", "4.for72h", "11.mix72h", 
                           "5.for96h", "12.mix96h", "6.for120h", "13.mix120h","7.for144h", "14.mix144h")
#organizing
aggreg_gen_2<- aggreg_gen_2[c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),]
rowSums(aggreg_gen_2)
data_norm<-aggreg_gen_2

#Selecting the ten most prevalent species in each sample
lst_names<- rownames(data_norm)
mostPrev<- setNames(replicate(14,list(0,1),simplify=FALSE),paste0(lst_names[1:14]))
for (i in 1:14){
  mostPrev[i][[1]] <- data.frame(sort(data_norm[i,], decreasing = TRUE) [1:10])
  mostPrev[i][[2]] <- rownames(mostPrev[i])
}
mostPrev<- lapply(mostPrev, setNames, c('Abundance'))#changing the column name in all dataframes in the list

########## creating a dataset with most prevalent species in all datasets
require(plyr)
for(i in 1:length(mostPrev)){
  colnames(mostPrev[[i]]) <- paste0(names(mostPrev)[i])
  mostPrev[[i]]$ROWNAMES  <- rownames(mostPrev[[i]])
}
mostPrevMerge <- join_all(mostPrev, by="ROWNAMES", type="full" )
rownames(mostPrevMerge) <- mostPrevMerge$ROWNAMES; mostPrevMerge$ROWNAMES <- NULL
mostPrevMerge[is.na(mostPrevMerge)] <- 0 #changing the NA for zero

########## 
require(reshape2)
melt_dom0<-melt(t(mostPrevMerge))

library(RColorBrewer)
n <- nrow(mostPrevMerge)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

require(ggplot2)
pdf(file = "/Users/alinevaz/Ecology/Figures/Figure1.pdf" ,  width = 15, height =8)
par(xpd=TRUE)
par(mfrow=c(1,1))
par(mar=c(10,10,5,5))
ggplot(melt_dom0,aes(x=Var1,y=value,fill=factor(Var2)))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(name="OTUs", values = col_vector)+
  xlab("OTUs")+ylab("Relative abundance")+
  theme(legend.text = element_text(colour="black", size=10))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits=c("1.for0h", "2.for24h" , "3.for48h", "4.for72h", "5.for96h", "6.for120h","7.for144h",
                            "8.mix0h", "9.mix24h","10.mix48h", "11.mix72h", "12.mix96h","13.mix120h", "14.mix144h"))+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#############################################################################################
#Diversity indices  - Figure 2
require(vegan)
############################################################################################
div_bac<-data.frame(data[,-c(1:7)])
div_bac2<- div_bac[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]
colnames(div_bac2)<- c("for0h","for24h","for48h","for72h","hib96h","for120h","for144h","mix0h","mix24h","mix48h","mix72h","mix96h","mix120h","mix144h")
div_bac2<-t(div_bac2)
colnames(div_bac2)<-c(1:ncol(div_bac2))
table<- matrix(0, nrow = 14, ncol = 3)
rownames(table)<-rownames(div_bac2)
colnames(table)<-c("richness", "shannon","evenness")

for(i in 1: nrow(table)){
  table[i,1]<-specnumber(which(div_bac2[i,]>0))#richness - vegan package
  table[i,2]<-round(diversity(which(div_bac2[i,]>0)),2)#Shannon
  table[i,3]<-round((1/(1-(diversity(div_bac2[i,], index="simpson"))))/(specnumber(div_bac2[i,])),5)
}

table2<- data.frame(cbind(c("0", "24", "48", "72", "96", "120", "144", "0", "24", "48", "72", "96", "120", "144"),
                          c(rep("for", 7), rep("mix", 7)), table))
colnames(table2) <- c("Time", "Group", "richness", "shannon", "evenness")

pdf(file = "/Users/alinevaz/Ecology/Figures/Figure3a.pdf" ,  width = 15, height =8)
ggplot(data=table2, aes(x=factor(Time, levels = unique(Time)), y=richness, group=Group)) +
  geom_line(aes(color = Group))+
  geom_point(aes(color = Group))+
  xlab("Time")+ylab("OTUs")+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(colour="black", size=10))+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

pdf(file = "/Users/alinevaz/Ecology/Figures/Figure3b.pdf" ,  width = 15, height =8)
ggplot(data=table2, aes(x=factor(Time, levels = unique(Time)), y=shannon, group=Group)) +
  geom_line(aes(color = Group))+
  geom_point(aes(color = Group))+
  xlab("Time")+ylab("OTUs")+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(colour="black", size=10))+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

pdf(file = "/Users/alinevaz/Ecology/Figures/Figure3c.pdf" ,  width = 15, height =8)
ggplot(data=table2, aes(x=factor(Time, levels = unique(Time)), y=evenness, group=Group)) +
  geom_line(aes(color = Group))+
  geom_point(aes(color = Group))+
  xlab("Time")+ylab("OTUs")+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.text = element_text(colour="black", size=10))+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
