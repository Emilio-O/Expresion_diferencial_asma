# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO
# Cargar la etiqueta de los datos que voy a analizar de GEO con getGEO()
gset <- getGEO("GSE43696", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Se hace una columna de nombres
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Con esto vamos a establecer los grupos
# Los grupos que voy a tener son "Control","Moderate Asthma" y "Severe Asthma"
gsms <- paste0("00000000000000000000111111111111111111111111111111",
               "11111111111111111111222222222222222222222222222222",
               "22222222")
sml <- strsplit(gsms, split="")[[1]]
length(sml) # Se tienen en total 108 muestras


# transformacion log2 de los datos
# de modo que estos puedan ser comparables
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# Asignamos muestras a los grpos
gs <- factor(sml)
groups <- make.names(c("Control","Moderate Asthma","Severe Asthma")) # tres grupos
levels(gs) <- groups
gset$group <- gs # Los grupos vam a ser considerdos como niveles
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)
design # Tres columnas extras con los grupos
# 1 indicaria la pertenencia de la muetra al grupo
# 0 indicaria que no pertenece

fit <- lmFit(gset, design)  # Ajustamos a un modelo lineal


cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design) # Establecemos contrastes 1, 0, -1
fit2 <- contrasts.fit(fit, cont.matrix) # Calculamos los coeficientes
head(fit2)


# Aqui es donde hacemos una table de los genes ordenados segun su p-valor
# El orden va del gen cuya expresion diferencial es mas significativa hasta el que es menos significatovo
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250) # Obtener los top 250 genes por su p-valor

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t") # Aqui es donde esta la tabla


p_v <- subset(tT, P.Value < 0.05)
p_v_LogF <- subset(p_v, F > 2)


p_v_LogF

#Esto es algo curioso...
# Esto es post examen
# Y no hice esto
# Hice algo de arriba creyendo que ya era innecesario hacer el subset y etc...

write.csv(p_v_LogF, file = "Ejercicio_5.csv")


# =========================
# Visualizar los resultados
# =========================


# Se hace un histograma de los P-valores detodos los genes
# Se asume que la mayoria de los genes no son expresados diferencialmente
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
# Y como podemos ver se cumple lo antes dicho

# Resumen de resultados como sobre o sub expresados, tambien no expresados
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)


# Para visualizar lo anterior se hace un diagrama de venn
vennDiagram(dT, circle.col=palette())
# En las partes solapadas son los genes significativos entre cada contraste


# Q-Q plot
t.good <- which(!is.na(fit2$F)) # Quitar nas
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
# Hace un plot de los cuantiles de las muestras contra los cuantiles teoricos de una distribuion de T de student
# Este plot ayuda a cerciornos de la calidad de los resultados
# Como los puntos siguen la linea podemos decir que siguen su distribucion teorica predicha




# volcano plot (log P-value vs log fold change)
colnames(fit2) # lista de los nombre de las columnas de los contrastes
ct <- 2        # Elegir el contraste de interes
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
# Visualizar los genes que se sobreexpresan o subexpresan de manera significativa


# MD plot (log fold change vs mean log expression)
# Visualizar aquellos genes que son significativamente expresados
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# boxplot
dev.new(width=3+ncol(gset)/6, height=5) # Esto po alguna razon no funciona en mi computadora
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666")) # Se eligen los colores
par(mar=c(7,4,2,1))
title <- paste ("GSE43696", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord]) # Con esto hacemos el boxplot
legend("topleft", groups, fill=palette(), bty="n")
# La distribucion de los valores de las muestras, el color es de acuerdo al grupo
# Observar la distribucion es util para determinar si las muestras son utiles para un
# analisis de expresion diferencial
# Como comentario mio: esto quiza deberia de estar al principio antes de hacer cualquier cosa...
dev.off()

# Distribucion del valor de la expresion
par(mar=c(4,4,2,1))
title <- paste ("GSE43696", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright") 
# El color segun el grupo y es un complemento del boxplot de arriba para revisar la normalizacion de los datos



# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminar renglones con NAs
ex <- ex[!duplicated(ex), ]  # remover duplicados
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
# Visualizar como las muestras estan relacionadas entre ellas...

legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools") 
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend
plotSA(fit2, main="Mean variance trend, GSE43696") # Revisar la relacion del promedio de la varianza
# después de un fit en el modelo lineal
# Cada puntito es un gen
# La linea roja es el promedio de la varianza
# La linea azul es la varianza constante

