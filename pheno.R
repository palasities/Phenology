##install packages
library(reshape)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(imager)
library(png)
library(gridExtra)
library(gtable)
library(logisticPCA)
library(dplyr)
library(ca)
library(factoextra)
library(FactoMineR)
library(gridExtra)
library(dendextend)
library(NbClust)
library(pvclust)
library(flexclust)
library(scrime)
library(bayesbio)
library(qvcalc)
library(stringdist)
library(vegan)
library(cluster)
library(purrr)
library(clustree)
library(ggraph)
library(igraph)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
install.packages("ape")
library(ape)
library(magrittr)
library(imager)
library(png)
library(magick)
library(ggplot2)
library(cowplot)


#Para extraer los vectores únicos es necsario que df esté sin coordenadas, sin fenología a visu ni clase sociológica.
df_fen=df_p[,-c(1,2,3,15,16)]

#valores únicos dentro del vector:
fenotipos=unique(df_fen, incomparables = FALSE, fromLast = FALSE,
                 nmax = NA)
#es necesario pasarlo a formato largo, primero se hace un nuevo dataframe:
fenotipos_ID=as.data.frame(cbind(rownames(fenotipos), fenotipos))
#columna uno se llame ID
names(fenotipos_ID)[1]="ID"

#se pasa a formato largo:
linear_fenot_melt=melt(fenotipos_ID,id.vars="ID")
head(linear_fenot_melt,2)

#todos juntos, para que el smooth se vea creo una nueva variable al linear_fenot_melt común a todas
linear_fenot_melt$test="smooth"


##Figura paper: sin colores 
figure1=  ggplot(data = linear_fenot_melt, aes(x = variable, y = value))+
  geom_smooth(method="loess", aes(fill=test, group=test, linetype=test), size=1.2, level=0.99, color="gray20")+
  geom_path(aes(group = ID), linetype=3, size=0.5) +
  #geom_point() +
  geom_jitter(width = 0.45, color="grey19")+
  #geom_line( linetype=4) +
  labs(x = "DOY", y= "phenophase") +
  theme_bw() + theme(legend.position = "none")+
  guides(color=guide_legend(order=1))

plot(figure1)