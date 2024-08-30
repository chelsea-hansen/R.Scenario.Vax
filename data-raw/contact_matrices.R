rm(list=ls())
library(RColorBrewer)
library(tidyverse)


###POLYMOD
polymod = read.csv("data-raw/POLYMOD_75.csv")

library(RColorBrewer)
n.cols=100
nice.cols <-  colorRampPalette(brewer.pal(9, "YlOrRd"))(n.cols)
heatmap(as.matrix(polymod)/sum(diag(as.matrix(polymod))), Rowv=NA, Colv=NA, scale='none', col=nice.cols)

## create a symmetric contact matrix
symmetricpoly <- 0.5 * (polymod + t(polymod))

## expand the contact under 1 year old
expendpoly_1 <- matrix(data = symmetricpoly[1,1],nrow = 5,ncol=5)
expendpoly_2 <- matrix(data = as.numeric(rep(symmetricpoly[1,],5)),byrow = T,nrow = 5,ncol=75)
expendpoly <- cbind(expendpoly_1,expendpoly_2)
expendpoly <- rbind(expendpoly,cbind(t(expendpoly_2),as.matrix(symmetricpoly)))

expendpoly = symmetricpoly
names(expendpoly)=c(rep("<2m",2),rep("2-3m",2),rep("4-5m",2),rep("6-7m",2),rep("8-9m",2),rep("10-11m",2),rep("1Y",12),rep("2-4Y",36),rep("5-9Y",1),
                    rep("10-19Y",2),rep("20-39Y",4),
                    rep("40-64Y",5),rep("65-74Y",2),rep("75+Y",1))


# aggregate ages into 13 age groups; average the number of contacts within each age group
ave_poly_x = aggregate(expendpoly, by = list(colnames(expendpoly)), FUN='sum' )
ave_poly = aggregate( t(ave_poly_x[,-1]), by = list(colnames(expendpoly)), FUN='mean')
rownames(ave_poly) <- ave_poly[,1]
colnames(ave_poly)[2:length(ave_poly)] <- rownames(ave_poly)

ageorder <- c("Group.1","<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65-74Y","75+Y")
contactpoly <- select(ave_poly,ageorder)
contactpoly <- contactpoly %>% mutate(Group.1=  factor(Group.1, levels = c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65-74Y","75+Y"))) %>%
  arrange(Group.1)
contactpoly <- as.matrix(contactpoly[,-1])
## create a symmetric contact matrix
contactpoly[lower.tri(contactpoly, diag = FALSE)] <- 0
contactpoly <- contactpoly + t(contactpoly)
diag(contactpoly) <- 0.5*diag(contactpoly)
rownames(contactpoly) <-  c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-64Y","65-74Y","75+Y")

# visualize contact matrix
heatmappoly <- as.matrix(contactpoly)
heatmap(heatmappoly/sum(diag(heatmappoly)), Rowv=NA, Colv=NA, scale='none', col=nice.cols)

# scale it to represent per capita contact per day
heatmappoly <- heatmappoly*rowSums(expendpoly)[1]/rowSums(heatmappoly)[1]
heatmappoly
write.csv(heatmappoly, "data-raw/polymod_rsv_net_agegroups.csv")

