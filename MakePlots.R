setwd("/Users/dehaas/Documents/densitydependence")

data1 <- read.csv("data_type1.csv")
data2 <- read.csv("data_type2.csv")
library(plotly)

mathematicadata <- read.csv("./Data/test.csv")


p1 <- plot_ly(z = as.matrix(data1), type = "heatmap",
             colorscale = 'Viridis')
p2 <- plot_ly(z = as.matrix(data2), type = "heatmap",
             colorscale = 'Viridis')

p1
p2


library(ggplot2)

ggplot(data = mathematicadata, aes(x = X0)) + geom_point(aes(y = X1.)) + geom_point(aes(y = X1..1))

    colnames(mathematicadata)
