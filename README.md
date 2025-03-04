# Phenological monitoring through high-resolution multispectral imaging data analysis

Este repositorio describe el análisis de datos fenológicos (observadores de campo), así como datos multiespectrales procesados en índices de vegetación (vuelos de drones) para estimar las diferentes cepas fenológicas de un monte bajo de Quercus pyrenaica mediante análisis de clusterización.

Para información específica del rodal de estudio, el periodo de estudio, la frecuencia de muestreo de los observadores, y la resolución temporal y espacial de los vuelos de los drones consulte el artículo "Phenological monitoring through high-resolution multispectral imaging as a management tool to characterize clonal structure in oak coppices" (Forest Management, doi: XXXXX).

A modo de resumen, se utilizarán varios paquetes de R para el análisis íntegro. Por una parte se analizarán datos categóricos (observadores); y de forma paralela se analizarán datos numéricos (drones). Gran parte de los resultados parciales de este código se procesaron mediante ArcGis (no se incluyen los análisis, pero se describen en el estudio mencionado anteriormente).

### Análisis de datos para los observadores de campo

Primero, se cargan los paquetes necesarios en RStudio:

```r
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
# library(EBImage)
library(magrittr)
library(imager)
library(png)
library(magick)
library(ggplot2)
library(cowplot)
```
Seguidamente, se carga el conjunto de datos "observer_1.csv". Para el resto de observadores es el mismo código.

```r
df_obs1=read.csv(".csv", header = TRUE, sep=";")

```
Quitamos las coordenadas para poder calcular los vectores únicos, que representan los diferentes fenotipos del dataframe. Con esa información, se hace un ggplot para observar las tendencias de los fenotipos en el tiempo de muestreo. Además se añade un "smooth" para observar la tendencia general de los pies estudiados (en negro):

```r

#Para extraer los vectores únicos es necsario que df esté sin coordenadas, sin fenología a visu ni clase sociológica (columnas 1,2,3,15 y 16).
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

#todos juntos, para que el smooth se vea creo una nueva variable al linear_fenot_melt común a todas. 
linear_fenot_melt$smooth="smooth"

##Figura 1 del paper
figure1=  ggplot(data = linear_fenot_melt, aes(x = variable, y = value))+
  geom_smooth(method="loess", aes(fill=smooth, group=smooth, linetype=smooth), size=1.2, level=0.99, color="gray20")+
  geom_path(aes(group = ID), linetype=3, size=0.5) +
  #geom_point() +
  geom_jitter(width = 0.45, color="grey19")+
  #geom_line( linetype=4) +
  labs(x = "DOY", y= "phenophase") +
  theme_bw() + theme(legend.position = "none")+
  guides(color=guide_legend(order=1))
plot(figure1)

```
Análisis clúster:

```r
#análisis de correspondencias múltiples para la base de datos solo con observaciones fenológicas.
set.seed(12)
AC_p=MCA(df_fen)

##varianza explicada por las coordenadas
AC_p$eig

#un nuevo dataframe con las 5 dimensiones por individuo.
MCA_ind=data.frame(AC_p$ind$coord)

##número óptimo de clusters en función de la información fenológica
gap_stat_JP =clusGap(MCA_ind, FUN = kmeans, nstart = 121, K.max = 35, B = 2000)

## plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat_JP)

##como sale 21 clusters de la función clusGap, metemos en la función kmeans 21
k1 = kmeans(MCA_ind, centers = 21, nstart = nrow(MCA_ind), iter.max = 121)

##preparamos la base de datos para ArcGis:
df_k=cbind(df_p[,c(2,3)], k1$cluster) 
names(df_k)=c("X","Y","k")
