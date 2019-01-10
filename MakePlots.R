setwd("/Users/dehaas/Documents/densitydependence")

data1 <- read.csv("data_type1.csv")
data2 <- read.csv("data_type2.csv")
library(plotly)


p1 <- plot_ly(z = as.matrix(data1), type = "heatmap",
             colorscale = 'Viridis')
p2 <- plot_ly(z = as.matrix(data2), type = "heatmap",
             colorscale = 'Viridis')

p1
p2

