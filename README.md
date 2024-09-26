# Helianthemum marifolium R

## Installing an loading Pakages:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("ggplot2") # Graphics
install.packages("tidyverse")
install.packages('htmltools')
library(htmltools)
install.packages('devtools')
library(devtools)
devtools::install_github("cstubben/trinotateR")
install.packages('patchwork')
install.packages("stringr")
install.packages("FactoMineR")
install.packages("factoextra")
```

``` r
library(ggplot2)
library(edgeR)
library(xfun)
library(tidyverse)
library(FactoMineR)
library(patchwork)
library(trinotateR)
library(stringr)
library(factoextra)
```

## Importing quants:

``` r
counts = read.delim("counts.txt",row.names=1)
```

``` r
group = factor(c(rep("andal", 4), rep("marif", 4), rep("molle", 4), rep("origa", 4)))

replicate= factor(c(rep(c("01","03","06","08","12","15","17","18","22","25","27","29","31","33","38","39"))))
#replicate= factor(c(rep(c("1","2","3","4"), 4)))
names = paste(group, replicate, sep = "")
sample=names
targent= data.frame(group, replicate, sample)
rownames(targent)=names
colnames(counts)= names
```

## DGELIST:

``` r
y=DGEList(counts = counts,
samples = targent$sample,
group = targent$group)
```

### Filtering:

``` r
keep = filterByExpr(y)
y = y[keep, ,keep.lib.sizes=FALSE]
```

### Normalization:

``` r
y = normLibSizes(y,method = "TMM")
```

``` r
design = model.matrix(~0+targent$group)
colnames(design) = levels(group)
y = estimateDisp(y,design)
```

## PCA:

``` r
pca.gene.expression = t(cpm(y))
res.pca = PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )


label = factor(c(rep("andalusicum", 4), rep("marifolium", 4), rep("molle", 4), rep("origanifolium", 4)))

fviz_pca_ind(res.pca, col.ind = group,
pointsize=2, pointshape=21, fill="black",
repel = TRUE,
addEllipses = TRUE,ellipse.type = "confidence",
ellipse.level= 0.95,
legend.title="Conditions",
title="",
show_legend=TRUE, show_guide=TRUE,habillage = label)
```

## Cluster Dendogram:

``` r
res.hcpc = HCPC(res.pca, graph=FALSE,nb.clust = 4)
fviz_dend(res.hcpc,k=4,
cex = 0.75,
palette = "jco",
rect = TRUE, rect_fill = TRUE,
rect_border = "jco",
type="rectangle",
labels_track_height = 1400)
```

## Getting PCA and HCPC multi-panel graphics:

``` r
g1 = fviz_pca_ind(res.pca, col.ind = group,
pointsize=2, pointshape=21, fill="black",
repel = TRUE,
addEllipses = TRUE,ellipse.type = "confidence",
ellipse.level= 0.95,
legend.title="Conditions",
title="PCA",
show_legend=TRUE,show_guide=TRUE,
show_legend=TRUE, show_guide=TRUE,habillage = label)
g2=fviz_dend(res.hcpc,k=4,
cex = 0.75,
palette = "jco",
rect = TRUE, rect_fill = TRUE,
rect_border = "jco",
type="rectangle",
labels_track_height = 1400)



g1 + g2 
```

## Getting GO terms multi-panel graphics:

``` r
x = read_trinotate("trinotate.xls")
summary_trinotate(x)

conteo_cellular_component = str_count(x$gene_ontology_BLASTX, 'cellular_component')
numero_total_cellular_component = sum(na.omit(conteo_cellular_component))

conteo_molecular_function = str_count(x$gene_ontology_BLASTX, 'molecular_function')
numero_total_molecular_function = sum(na.omit(conteo_molecular_function))

conteo_biological_process = str_count(x$gene_ontology_BLASTX, 'biological_process')
numero_total_biological_process = sum(na.omit(conteo_biological_process))

blastx = data.frame(
  GO_term = c("Cellular component", "Molecular function", "Biological process"),
  valor = c(81463, 95621, 104140)
)

ggplot(blastx, aes(x = GO_term, y = valor, fill = GO_term)) +
  geom_bar(stat = "identity") +
  labs(title = "Occurrence of GO terms", x = "GO Term", y = "Value") +
  theme_minimal()
```

``` r
GO = as.factor(na.omit(x$gene_ontology_BLASTX))
GO1= GO[1]
elemento =factor(unlist(strsplit(as.character(GO), "`")))

fac_celular_components <- elemento[grep("cellular_component", elemento)]
caracteres_celular_components <- sub(".*\\^.*\\^(.*)", "\\1", fac_celular_components)
conteos_celular_components <- table(caracteres_celular_components)  
top_10_conteos_celular_components = as.data.frame(sort(conteos_celular_components, decreasing = TRUE)[1:10])

g3=ggplot(top_10_conteos_celular_components, aes(x = caracteres_celular_components, y = Freq)) + 
  geom_col(fill = "skyblue") +  
  coord_flip() +  
  labs(title = "Cellular components", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    plot.title = element_text(size = 16)    
  )
```

``` r
fac_molecular_function <- elemento[grep("molecular_function", elemento)]
caracteres_molecular_function <- sub(".*\\^.*\\^(.*)", "\\1", fac_molecular_function)
conteos_molecular_function <- table(caracteres_molecular_function)  
top_10_conteos_molecular_function = as.data.frame(sort(conteos_molecular_function, decreasing = TRUE)[1:10])

g4=ggplot(top_10_conteos_molecular_function, aes(x = caracteres_molecular_function, y = Freq)) + 
  geom_col(fill = "pink") +  
  coord_flip() +  
  labs(title = "Molecular function", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    plot.title = element_text(size = 16)    
  )
```

``` r
fac_biological_process <- elemento[grep("biological_process", elemento)]
caracteres_biological_process <- sub(".*\\^.*\\^(.*)", "\\1", fac_biological_process)
conteos_biological_process <- table(caracteres_biological_process)  
top_10_conteos_biological_process = as.data.frame(sort(conteos_biological_process, decreasing = TRUE)[1:10])


g5=ggplot(top_10_conteos_biological_process, aes(x = caracteres_biological_process, y = Freq)) + 
  geom_col(fill = "purple") +  
  coord_flip() +  
  labs(title = "Biological process", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    plot.title = element_text(size = 16)    
  )
```

``` r
g3 / g4 /  g5
```
