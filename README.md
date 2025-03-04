# Phenological monitoring through high-resolution multispectral imaging data analysis

Este repositorio describe el análisis de datos fenológicos (observadores de campo), así como datos multiespectrales procesados en índices de vegetación (vuelos de drones) para estimar las diferentes cepas fenológicas de un monte bajo de Quercus pyrenaica mediante análisis de clusterización.

Para información específica del rodal de estudio, el periodo de estudio, la frecuencia de muestreo de los observadores, y la resolución temporal y espacial de los vuelos de los drones consulte el artículo "Phenological monitoring through high-resolution multispectral imaging as a management tool to characterize clonal structure in oak coppices" (Forest Management, doi: XXXXX).

A modo de resumen, se utilizarán varios paquetes de R para el análisis íntegro. Por una parte se analizarán datos categóricos (observadores); y de forma paralela se analizarán datos numéricos (drones). Gran parte de los resultados parciales de este código se procesaron mediante ArcGis (no se incluyen los análisis, pero se describen en el estudio mencionado anteriormente).

### Análisis de datos para los observadores de campo

Primero, se cargan los paquetes necesarios en RStudio:

```r
# Lista de paquetes necesarios
packages <- c(
  "reshape", "tidyverse", "RColorBrewer", "ggplot2", "patchwork", "imager",
  "png", "gridExtra", "gtable", "logisticPCA", "dplyr", "ca", "factoextra",
  "FactoMineR", "gridExtra", "dendextend", "NbClust", "pvclust", "flexclust",
  "scrime", "bayesbio", "qvcalc", "stringdist", "vegan", "cluster", "purrr",
  "clustree", "ggraph", "igraph", "ape", "magrittr", "magick", "cowplot"
)

# Instalar los paquetes que no estén ya instalados
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Aplicar la instalación a todos los paquetes
lapply(packages, install_if_missing)

# Cargar los paquetes en la sesión de R
invisible(lapply(packages, library, character.only = TRUE))

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
kmax_obs1=nrow(unique(MCA_ind))
gap_stat_obs1 =clusGap(MCA_ind, FUN = kmeans, nstart = 121, K.max = kmax_obs1, B = 2000)

## plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat_obs1)

##como sale 21 clusters de la función clusGap, metemos en la función kmeans 21
k1 = kmeans(MCA_ind, centers = 21, nstart = nrow(MCA_ind), iter.max = 121)

##preparamos la base de datos para ArcGis:
df_k_obs1=cbind(df_p[,c(2,3)], k1$cluster) 
names(df_k_obs1)=c("X","Y","k")
```

Dendrograma

```r

DENDRO=diana(df_k)
pltree(DENDRO, cex = 0.8, hang = -1, main = "", labels=rownames(DENDRO))
rect.hclust(DENDRO, k=21, border=2:10)

```

### Análisis de datos para los índices de vegetación

```r
set.seed(12)
df_ndvi=read.csv("", header = TRUE, sep=";")
df_ndvi=df_ndvi[c(1:121),c(1:11)]
rownames(df_ndvi)=df_ndvi$ID

#le quuito la columna ID
df_ndvi=df_ndvi[,-c(1)]

#análisis PCA sin coordenadas
PCA_ndvi=PCA(df_ndvi[,-c(1,2)])
##incluyo coordenadas
PCA_ndvi_withcoords=cbind(PCA_ndvi$ind$coord, ndvi$X, ndvi$Y)
##varianza explicada por las variables
PCA_ndvi$eig

kmax_NDVI=unique(nrow(PCA_ndvi_withcoords))

##número óptimo de clusters
gap_stat_ndvi=clusGap(PCA_ndvi_withcoords,
FUN = kmeans,
nstart=121,
K.max = kmax_NDVI,
B = 2000)
fviz_gap_stat(gap_stat_ndvi)

#estadísticamente salen 14 clusters vamos a ver cómo los agrupa por dron (sin coordenada, por eso no uso joint:

kmeans_drones=kmeans(PCA_ndvi_withcoords,centers=14,iter.max = 121, nstart = 121 )
kmeans_drones$cluster
```
## Figura errores

```

df_errores=read.csv("errores_obs_split.csv", header = TRUE, sep=";")
rownames(df_errores)=df_errores$ID
errores_melt=melt(df_errores,id.vars="ID")

#cambiar orden del eje X
errores_melt$ID <- factor(errores_melt$ID , levels=c("Observer1","Observer2","Observer3", "NDVI", "NDRE", "GRVI", "RVI"))

##plotear

ggplot(errores_melt, aes(ID, value, fill=ID))+
geom_bar(stat='identity',color="black", alpha=0.6)+
facet_wrap(~variable,  scales = 'free')+
theme(axis.text.x =element_text(angle = 90, hjust = 1))+
theme(axis.text.y =element_text(angle = 90, hjust = 1))+
scale_fill_grey(start = 0.6, end = 0.05)+
theme(axis.text.x = element_blank())+
theme(axis.ticks.x = element_blank())+
theme(axis.text =element_text(size=10.5))+
theme(axis.text = element_text(size=10.5))+
theme(legend.position = "none")+
labs(x = "Obs", y= "Value(%)")+
theme(axis.title.x=element_blank())+ 
theme(axis.title.y = element_text(face="italic", vjust=1.5, colour="black", size=rel(1.1)))

```

##Figura Medida de gestión:

```
df_stool=read.csv("Clonal_stool.csv", header = TRUE, sep=";")
df_distance=read.csv("split_Gestion/Distance.csv", header = TRUE, sep=";")
df_intra=read.csv("Intra_clonal.csv", header = TRUE, sep=";")

#stool
df_stool$ID <- factor(df_stool$ID , levels=c("Observers", "NDVI", "NDRE", "GRVI", "RVI"))

smooth1=24

g1=ggplot(df_stool, aes(x=ID, y=Clonal.stool, fill=ID))+
geom_bar(stat='identity', position = position_dodge(),color="black", alpha=0.6)+
scale_fill_grey(start = 0.6, end = 0.05)+
geom_errorbar(aes(ymin = Clonal.stool -i.c, ymax = Clonal.stool+ i.c), alpha=0.5)+
geom_hline(aes(yintercept=smooth1),  alpha=0.5)+
labs(y= "Clonal density \n (clones 0.25 ha)")+
theme(legend.position = "none")+
theme(axis.title.x=element_blank())+
 theme(axis.text.x = element_blank())+
labs(tag = "A")##€sto es lo que quita el eje X

#distance:
smooth2=4.92005
#intervalo conf
smooth_err2=1.290454703

df_distance$ID <- factor(df_distance$ID , levels=c("Observers", "NDVI", "NDRE", "GRVI", "RVI"))

g2=ggplot(df_distance, aes(x=ID, y=Max..Distance..m., fill=ID))+
geom_bar(stat='identity', position = position_dodge(),color="black", alpha=0.6)+
scale_fill_grey(start = 0.6, end = 0.05)+
geom_errorbar(aes(ymin = Max..Distance..m. -i.c, ymax = Max..Distance..m.+ se), alpha=0.5)+
geom_hline(aes(yintercept=smooth2),  alpha=0.5)+
geom_hline(aes(yintercept=smooth2+smooth_err2), linetype="dashed", alpha=0.3)+
geom_hline(aes(yintercept=smooth2-smooth_err2), linetype="dashed",  alpha=0.3)+
labs(y= "Max. distance \n (m)")+
theme(legend.position = "none")+
theme(axis.title.x=element_blank())+
theme(axis.text.x = element_blank())+
labs(tag = "B")##€sto es lo que quita el eje X

g2

#intra

smooth3=5.409090909

#intervalo de confianza
smooth_err3=1.634536992

df_intra$ID <- factor(df_intra$ID , levels=c("Observers", "NDVI", "NDRE", "GRVI", "RVI"))

g3=ggplot(df_intra, aes(x=ID, y=Intra.clonal.density, fill=ID))+
geom_bar(stat='identity', position = position_dodge(),color="black", alpha=0.6)+
scale_fill_grey(start = 0.6, end = 0.05)+
geom_errorbar(aes(ymin = Intra.clonal.density -i.c, ymax = Intra.clonal.density+ i.c), alpha=0.5)+
geom_hline(aes(yintercept=smooth3),  alpha=0.5)+
geom_hline(aes(yintercept=smooth3+smooth_err3), linetype="dashed", alpha=0.3)+
geom_hline(aes(yintercept=smooth3-smooth_err3), linetype="dashed",  alpha=0.3)+
labs(y= expression( "Intra-clonal density \n  (ramet genet)"))+                         ##(ramet genet"^-1*")")) PARA QUE SALGA COMO EXPRESION
theme(legend.position = "none")+
theme(axis.title.x=element_blank())+
labs(tag = "C")##€sto es lo que quita el eje X


plot_grid(g1, g2, g3, ncol = 1,
                   align="hv",
                   axis="l",
                   rel_widths = 1,
                   rel_heights = 0.5,
                   scale=c(1,1,1,1,1,1))

plot_grid(g1, g2, g3, ncol = 1,
          vjust=12)

```

## Evolución de los vuelos GRVI

```

grvi1=read.csv("1.csv", header = TRUE, sep=";")
grvi2=read.csv("2.csv", header = TRUE, sep=";")
grvi3=read.csv("3.csv", header = TRUE, sep=";")
grvi4=read.csv("4.csv", header = TRUE, sep=";")
grvi5=read.csv("5.csv", header = TRUE, sep=";")
grvi6=read.csv("6.csv", header = TRUE, sep=";")
grvi7=read.csv("7.csv", header = TRUE, sep=";")

grvi=rbind(grvi1,grvi2,grvi3,grvi4,grvi5,grvi6)
grvi$Flight=as.factor(grvi$group)


colores_histo=c("#CD2626","#698B69","#9A32CD","#8B8B00", "#FF7256","#53868B")

#histo lineas, Densidad.

histo_line=ggplot(grvi, aes(GRID_CODE, color=Flight, group=Flight))+
  geom_density()+
  theme_bw()+
  labs(x = "GRVI", y="Density")+
  xlim(1,10)+
  theme(legend.position='top', 
        legend.direction='horizontal')+
  guides(color=guide_legend(nrow=1)) #clave aquí ponerle color, porque si le pones fill como al de arriba no te lo hace porque entiende que es                                       un área

  histo_line
  
  
  
  ##mapas:
  
p1_clip <- ggdraw() + draw_image("grvi_CLIPS/1.jpg", scale = 0.9)
p2_clip <- ggdraw() + draw_image("grvi_CLIPS/2.jpg", scale = 0.9)
p3_clip <- ggdraw() + draw_image("grvi_CLIPS/3.jpg", scale = 0.9)
p4_clip <- ggdraw() + draw_image("grvi_CLIPS/4.jpg", scale = 0.9)
p5_clip <- ggdraw() + draw_image("grvi_CLIPS/5.jpg", scale = 0.9)
p6_clip <- ggdraw() + draw_image("grvi_CLIPS/6.jpg", scale = 0.9)
p7 <- ggdraw() + draw_image("Main.jpg", scale = 0.9)

p_unidad_clip=plot_grid(p1_clip, p2_clip, p3_clip, p4_clip, p5_clip, p6_clip,
                   labels=c("1","2","3","4","5","6"),
                   nrow=2,
                   ncol=3,
                   align="hv",
                   axis="l",
                   rel_widths = 1,
                   rel_heights = 1,
                   scale=c(1,1,1,1,1,1),
                   label_size = 16,
                   hjust=-0.5,
                   vjust=2.5)
                   #byrow = TRUE)

  ptodo_clip=plot_grid(p7, p_unidad_clip,
                rel_widths = c(1,1))


  #combinado: mapas e histograma
  
  plot_grid(ptodo_clip, histo_line, ncol = 1)
  
  ggsave("GRVI_horizontal_histoLine.jpg", units = "in", width = 10, height = 7, bg="white", device = "jpg", dpi = 700)

```
