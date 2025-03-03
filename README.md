# Pheno
Phenological drones 

## 📜 Código en R para el análisis
```r
library(DESeq2)
library(dplyr)

metadata_local$Time <- factor(metadata_local$Time)

dds_local_time <- DESeqDataSetFromMatrix(
  countData = counts_local,
  colData = metadata_local,
  design = ~ Time
```
y ahora sí que puedes hacerlo


