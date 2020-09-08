
list.of.packages <- c("vcfR",
                      "poppr",
                      "ape",
                      "RColorBrewer",
                      "dplyr",
                      "tidyverse",
                      "adegenet",
                      "igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org/")
sapply(list.of.packages, require, character.only = TRUE)

#Abrir el archivo VCF usando vcfR y verificar que esten las muestras y SNPs cargados.
jaguar.VCF <- read.vcfR("data/populations.snps.vcf")


#Cargar la información sobre la población en un archivo .txt.Este incluye ID y nombre de la población.
#Usar la función read.table():
pop.data <- read.table("data/pop_map_25smp.txt", sep = "\t", header = FALSE, col.names = c("AccesID", "population")) 

pop_data <- readr::read_tsv("data/pop_map_25smp.txt", col_names =  FALSE) %>% 
  dplyr::rename(AccesID = X1) %>% 
  dplyr::rename(population = X2)
populations <- distinct(pop_data, population) %>% pull()
jaguar_vcf <- jaguar.VCF@gt %>% tibble::as_tibble() 

jaguar_vcf_gather <- jaguar_vcf %>% 
  tidyr::gather(AccesID, Info, -FORMAT) %>% 
  na.exclude() 

#Verificar que todas las muestras en el VCF y la tabla de datos de población estén incluidas.
val <- dplyr::distinct(jaguar_vcf_gather, AccesID) %>% 
    dplyr::anti_join(pop.data_tibble) %>%
  nrow 
ifelse(val != 0, "some populations are missing!", "All populations in VFC ok" )




#Convertir el objeto vcfR en un objeto genlight, usando la función: vcfR2genlight:
gl.jaguar <- vcfR::vcfR2genlight(jaguar.VCF)

#Especificar la ploidía del organismo
ploidy(gl.jaguar) <- 2

#Agregar la información sobre las poblaciones al objeto genlight
pop(gl.jaguar) <- pop_data$population

#Crear una matriz de distancia genética a partir de objetos genlight, con la función bitwise.dist().

geneticDistance <- poppr::bitwise.dist(gl.jaguar,
                              percent = TRUE,
                              mat = FALSE,
                              missing_match = TRUE,
                              scale_missing = FALSE,
                              euclidean = FALSE,
                              differences_only = FALSE,
                              threads = 3L
)



png("data/heatmap_distance.png", 640, 640)
heatmap(as.matrix(geneticDistance))
legend(x="bottomright", legend=c("min", "ave", "max"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
title("Mapa de calor entre las distancias bit a bit de todas las muestras")
dev.off()
geneticDistanceTibble <- geneticDistance %>% 
  broom::tidy() %>% 
  rename(population_1 = item1, population_2 = item2)
library(igraph)

png("data/boxplot_distancia.png", 480, 480)

geneticDistanceTibble %>% 
  ggplot(aes(population_1, distance)) + geom_boxplot() + 
  coord_flip() + 
  xlab("Población 1") + 
  ylab("Distancia genética (proporción de loci diferentes") + 
  ggtitle("Distancia genética de una población con respecto a las demás")
dev.off()
#Construir un árbol de distancia genética para representar la relación genética de las muestras. 

tree <- poppr::aboot(gl.jaguar, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = FALSE, cutoff = 50, quiet = T)
cols <- brewer.pal(n = length(populations), name = "Dark2")
?aboot
library(ggtree)
#Representación gráfica del árbol 

png("data/arbol_jaguar.png", 480, 480)
plot(tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.jaguar)], type = "phylogram")
legend('topleft',
       legend = populations,
       fill = cols, border = FALSE, bty = "n", cex = 1)
axis(side = 1)
title(xlab = "Distancia genética (proporción de loci diferentes)")
dev.off()

#Construcción de redes de expasión mínima, a través del objeto genlight objeto y una matriz de distancia.
jaguar.dist <- poppr::bitwise.dist(gl.jaguar)
jaguar.msn <- poppr::poppr.msn(gl.jaguar, jaguar.dist, showplot = TRUE, include.ties = T)
?poppr::poppr.msn

#Representación gráfica de la red.
node.size <- rep(2, times = nInd(gl.jaguar))
names(node.size) <- indNames(gl.jaguar)
vertex.attributes(jaguar.msn$graph)$size <- node.size
#set.seed(1)

png("data/grafo_red_jaguar.png", 480, 480)
plot_poppr_msn(gl.jaguar, jaguar.msn ,
               palette =  cols, gadj = 70)

dev.off()

?plot_poppr_msn
#Realizar un PCA con el objeto genlight usando la función  glPCA.
jaguar.pca <- adegenet::glPca(gl.jaguar, nf = 3)

png("data/PCA_var.png", 480, 480)

tibble(var = 100 * jaguar.pca$eig/sum(jaguar.pca$eig) ) %>% 
  mutate(PCA = row_number()) %>% 
  ggplot() + 
  geom_bar(aes(PCA, rev(var), fill = rev(var)), stat = "identity") + 
  coord_flip() + 
  ggtitle("Porcentaje de varianza explicado por PCA") + 
  xlab("Eigenvalues") + 
  ylab("Porcentaje de varianza explicado")

dev.off()



png("data/PCA_var_acumulada.png", 480, 480)

tibble(var = 100 * jaguar.pca$eig/sum(jaguar.pca$eig) ) %>%
  mutate(var_acummulated = cumsum(var)) %>% 
  mutate(PCA = row_number()) %>% 
  ggplot() + 
  geom_bar(aes(PCA, (var_acummulated), fill = (var_acummulated)), stat = "identity") + 
  #coord_flip() + 
  ggtitle("Porcentaje de varianza explicado por PCA") + 
  xlab("Eigenvalues") + 
  ylab("Varianza acumulada") + 
  scale_fill_continuous(name = "Varianza acumulada")

dev.off()
barplot( 100 * jaguar.pca$eig/sum(jaguar.pca$eig),
         col = heat.colors(50), 
         main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#Ver los resultados de PCA
jaguar.pca.scores <- as_tibble(jaguar.pca$scores) %>% 
  mutate(pop = pop(gl.jaguar)) 
library(ggplot)

p <- jaguar.pca.scores %>% 
  ggplot(aes(x=PC1, y=PC2, col = pop)) +
  geom_point(size=2) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = cols) +
  #geom_hline(yintercept = 0) +
  #geom_vline(xintercept = 0) +
  theme_classic() + 
  ggtitle("PCA de poblaciones de jaguares")

png("data/PCA_plot.png", 480, 480)
p
dev.off()
#Análisis discriminante de componentes principales (DAPC), para maximizar la varianza entre poblaciones.

pnw.dapc <- dapc(gl.jaguar, n.pca = 3, n.da = 2)


#Para confirmar que el DAPC es similar al PCA, trazar los datos en un diagrama de dispersión.

png("data/DAPC.png", 480, 480)
ade4::scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "topright", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, posi.da = "bottomleft" )
title("Análisis discriminante de componentes principales")
dev.off()
#Visualizar las asignaciones posteriores de cada muestra a una población, a traves del diagrama compuesto de barras apiladas.

png("data/compoplot.png", 480, 480)
adegenet::compoplot(pnw.dapc, col = cols, posi = NA, show.lab = TRUE, )
legend("topright",
       inset=c(-0.04,0),
       legend = populations,
       fill = cols, border = FALSE, bty = "n", cex = 1)
title("Probabilidad de pertenencia de una población")

dev.off()

#Separar las muestras por población.
#Observación de probabilidad de pertenencia para una población: nombre de la muestra, la población original y la población asignada


dapc_results <- pnw.dapc$posterior %>% 
  as_tibble() %>% 
  mutate(pop = pop(gl.jaguar)) %>% 
  mutate(Sample = rownames(pnw.dapc$posterior)) %>% 
  tidyr::gather(Assigned_Pop, Posterior_membership_probability, -pop, -Sample) %>% 
  rename(Original_Pop = pop) 
  

#Rpresentación gráfica del marco de estructura poblacional
png("data/estructura_poblacional.png", 480, 480)
dapc_results %>% 
  ggplot(aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cols, name = "Población Asignada") +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + 
  ylab("Probabilidad de pertenencia (posterior)") +
  xlab("Muestra") + 
  #coord_flip() +
  ggtitle("Marco de estructura poblacional")
dev.off()

