# Phenological monitoring through high-resolution multispectral imaging data analysis

This repository describes the analysis of phenological data collected by field observers, as well as multispectral data processed into vegetation indices (drone flights) to estimate the different phenological strains of a Quercus pyrenaica coppice through clustering analysis.

For specific information about the study please refer to the article "Linking leaf phenology and clonal structure in oak coppices through multispectral imaging" (Forest Management, doi: XXXXX).

Data Analysis for Field Observers

First, the required R packages are loaded into RStudio:

R version 4.3.0

```r
# List of Required Packages
packages <- c(
  "reshape", "tidyverse", "RColorBrewer", "ggplot2", "patchwork", "imager",
  "png", "gridExtra", "gtable", "logisticPCA", "dplyr", "ca", "factoextra",
  "FactoMineR", "gridExtra", "dendextend", "NbClust", "pvclust", "flexclust",
  "scrime", "bayesbio", "qvcalc", "stringdist", "vegan", "cluster", "purrr",
  "clustree", "ggraph", "igraph", "ape", "magrittr", "magick", "cowplot", "sf"
)

# Install Packages If Not Already Installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply Installation to All Packages
lapply(packages, install_if_missing)

# Load Packages into the R Session
invisible(lapply(packages, library, character.only = TRUE))
```

Next, the dataset "obs1.csv" is loaded. The same code applies to the rest of the observers.

```r
df_obs1=read.csv("obs1.csv", header = TRUE, sep=";")

```
We remove the coordinates to calculate the unique vectors, which represent the different phenotypes in the dataframe. Using this information, a ggplot is created to visualize the trends of the phenotypes over the sampling period. Additionally, a smooth layer is added to observe the overall trend of the studied trees -

```r

#To extract the unique vectors, the dataframe (df) must exclude ID and coordinates
df_fen_obs1=df_obs1[,-c(1,2,3)]

#unique values
fenotipos=unique(df_fen_obs1, incomparables = FALSE, fromLast = FALSE,
        nmax = NA)

fenotipos_ID=as.data.frame(cbind(rownames(fenotipos), fenotipos))

#rename 1st column
names(fenotipos_ID)[1]="ID"

#long format:
linear_fenot_melt=melt(fenotipos_ID,id.vars="ID")
head(linear_fenot_melt,2)

#common variable smooth. 
linear_fenot_melt$smooth="smooth"

##Figure 1
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
Cluster analysis:

## Phenological data. See methods in the article for more details.

```r
#Multiple Correspondence Analysis (MCA) for the database containing only phenological observations.
set.seed(12)

MCA_obs1=MCA(df_fen_obs1) ##no X,Y coords

##new dataframe
MCA_obs1_coords=data.frame(MCA_obs1$ind$coord) ##no X, Y coords (spatial coords)

#define optimal number of clusters
gap_stat_obs1 =clusGap(MCA_obs1_coords, FUN = kmeans, nstart = 121, K.max = 35, B = 10) ##it's going to take long time with B=2000.

## plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat_obs1)

##MCA coords + add spatial coordinates.
for_ID_K_obs1=cbind(MCA_obs1$ind$coord,df_obs1[,c(2,3)])

##OBS1. id-k. kmeans with centers = fviz_gap_stat(gap_stat_obs1)
ID_K_obs1=kmeans(for_ID_K_obs1,centers=21,iter.max = 121, nstart = 121 )

#df for ArcGIS or custom map (next):
db_obs1_k=cbind(obs1[,c(2,3)],ID_K_obs1$cluster)
names(db_obs1_k)=c("X","Y","k")

##load shapefile:
shapefile_path <- "clones_parcela.shp"
aislamiento_cepas=st_read(shapefile_path)

##ID-K map obs1
### from point to sf object
puntos_sf_obs1 <- st_as_sf(db_obs1_k, coords = c("X", "Y"), crs = st_crs(aislamiento_cepas))

#Graph
ggplot() +
  geom_sf(data = aislamiento_cepas, fill = NA, color = "gray40") + 
  geom_sf(data = puntos_sf_obs1, aes(color = as.factor(k)), size = 3.5) +  # k como factor categórico
  geom_text(data = db_obs1_k, aes(x = X, y = Y, label = rownames(db_obs1_k)), 
            size = 4, color = "black") +
  scale_color_viridis_d(option = "turbo", begin = 0, end = 1) +       # Paleta llamativa
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "Cluster k")

```

Dendrograma. Diana function  

```r

db_obs1_k$k=as.factor(db_obs1_k$k)
rownames(db_obs1_k)=df_obs1$ID
DENDRO=diana(db_obs1_k)
pltree(DENDRO, cex = 0.8, hang = -1, main = "Dendrogram diana function", labels=rownames(DENDRO))
rect.hclust(DENDRO, k=21, border=2:10)

```

## Vegetation indexes data analysis. See methods in the article for more details.

```r
set.seed(12)
ndvi=read.csv("NDVI.csv", header = TRUE, sep=";")
ndvi=ndvi[c(1:121),c(1:10)]
rownames(ndvi)=ndvi$ID

ndvi=ndvi[,-c(1)]

#PCA analysis
PCA_ndvi=PCA(ndvi[,-c(1,2)])
PCA_ndvi_withcoords=cbind(PCA_ndvi$ind$coord, ndvi$X, ndvi$Y)


gap_stat_ndvi=clusGap(PCA_ndvi_withcoords,
FUN = kmeans,
nstart=121,
K.max = 35,
B = 10)  ##it's going to take long time with B=2000.
fviz_gap_stat(gap_stat_ndvi)

##id-k for NDVI
ID_K_ndvi=kmeans(PCA_ndvi_withcoords,centers=14,iter.max = 121, nstart = 121 )

##same as obs1

```
## Errors figure
Search de db
```r

df_errores=read.csv("errores_obs_split.csv", header = TRUE, sep=";")
rownames(df_errores)=df_errores$ID
errores_melt=melt(df_errores,id.vars="ID")

#change X-axis order
errores_melt$ID <- factor(errores_melt$ID , levels=c("Observer1","Observer2","Observer3", "NDVI", "NDRE", "GRVI", "RVI"))

##plot
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

## Management tool figure

```r
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
labs(y= expression( "Intra-clonal density \n  (ramet genet)"))+                         ##(ramet genet"^-1*")")) as expresion
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

## GRVI values time-evolution

```r

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

##density lines

histo_line=ggplot(grvi, aes(GRID_CODE, color=Flight, group=Flight))+
  geom_density()+
  theme_bw()+
  labs(x = "GRVI", y="Density")+
  xlim(1,10)+
  theme(legend.position='top', 
        legend.direction='horizontal')+
  guides(color=guide_legend(nrow=1)) #clave aquí ponerle color, porque si le pones fill como al de arriba no te lo hace porque entiende que es                                       un área

  histo_line
  
  
##maps:
  
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

  #maps and histograms combined
  
  plot_grid(ptodo_clip, histo_line, ncol = 1)
  

```
