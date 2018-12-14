rm(list=ls())

library(geoR)
class(elevation)
plot(elevation, lowess = T)
